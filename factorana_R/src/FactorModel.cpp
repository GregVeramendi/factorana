#include "FactorModel.h"
#include "distributions.h"
#include "gauss_hermite.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <set>

FactorModel::FactorModel(int n_obs, int n_var, int n_fac, int n_typ,
                         int n_mix, bool correlated, int n_quad)
    : nobs(n_obs), nvar(n_var), nfac(n_fac), ntyp(n_typ),
      nmix(n_mix), fac_corr(correlated), nquad_points(n_quad),
      use_weights(false), use_adaptive(false), adapt_threshold(0.3)
{
    data.resize(nobs * nvar);

    // Initialize with factor parameters
    // For single mixture: nfac factor variances
    nparam = nfac;  // Start with factor variances

    // Add correlation parameters if correlated factors
    // For 2-factor correlated model: 1 correlation parameter
    if (fac_corr && nfac == 2) {
        nparam += 1;  // Add factor_corr_1_2
    }

    // Compute type model parameters
    // Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
    // For ntyp > 1: (ntyp - 1) * nfac loading parameters
    if (ntyp > 1) {
        ntyp_param = (ntyp - 1) * nfac;
        type_param_start = nparam;  // Type params come after factor params
        nparam += ntyp_param;
    } else {
        ntyp_param = 0;
        type_param_start = -1;
    }

    nparam_free = nparam;  // All factor parameters are free by default
}

void FactorModel::AddModel(std::shared_ptr<Model> model, int nparams)
{
    models.push_back(model);
    param_model_start.push_back(nparam);
    param_model_count.push_back(nparams);
    nparam += nparams;
}

void FactorModel::SetData(const std::vector<double>& dat)
{
    if (dat.size() != nobs * nvar) {
        throw std::runtime_error("Data size mismatch");
    }
    data = dat;
}

void FactorModel::SetData(const Eigen::MatrixXd& dat)
{
    if (dat.rows() != nobs || dat.cols() != nvar) {
        throw std::runtime_error("Data dimension mismatch");
    }
    data.resize(nobs * nvar);
    for (int i = 0; i < nobs; i++) {
        for (int j = 0; j < nvar; j++) {
            data[i * nvar + j] = dat(i, j);
        }
    }
}

void FactorModel::SetQuadrature(const std::vector<double>& nodes,
                                const std::vector<double>& weights)
{
    quad_nodes = nodes;
    quad_weights = weights;
}

void FactorModel::SetObservationWeights(const std::vector<double>& weights)
{
    if (weights.size() != nobs) {
        throw std::runtime_error("Observation weights size mismatch: expected " +
                                 std::to_string(nobs) + " but got " +
                                 std::to_string(weights.size()));
    }
    obs_weights = weights;
    use_weights = true;
}

void FactorModel::SetAdaptiveQuadrature(const std::vector<std::vector<double>>& factor_scores,
                                        const std::vector<std::vector<double>>& factor_ses,
                                        const std::vector<double>& factor_vars,
                                        double threshold,
                                        int max_quad)
{
    // Validate inputs
    if (factor_scores.size() != nobs) {
        throw std::runtime_error("factor_scores size mismatch: expected " +
                                 std::to_string(nobs) + " observations");
    }
    if (factor_ses.size() != nobs) {
        throw std::runtime_error("factor_ses size mismatch: expected " +
                                 std::to_string(nobs) + " observations");
    }
    if (factor_vars.size() != nfac) {
        throw std::runtime_error("factor_vars size mismatch: expected " +
                                 std::to_string(nfac) + " factors");
    }

    // Store settings
    adapt_threshold = threshold;
    adapt_factor_var = factor_vars;

    // Resize storage
    obs_nquad.resize(nobs);
    obs_fac_center.resize(nobs);
    obs_weights.resize(nobs, 1.0);

    // Collect all unique nquad values needed
    std::set<int> needed_nquad;

    // Compute per-observation quadrature settings
    for (int iobs = 0; iobs < nobs; iobs++) {
        if (factor_scores[iobs].size() != nfac || factor_ses[iobs].size() != nfac) {
            throw std::runtime_error("Factor score/SE dimensions mismatch for observation " +
                                     std::to_string(iobs));
        }

        obs_nquad[iobs].resize(nfac);
        obs_fac_center[iobs].resize(nfac);

        for (int ifac = 0; ifac < nfac; ifac++) {
            double f_score = factor_scores[iobs][ifac];
            double f_se = factor_ses[iobs][ifac];
            double f_var = factor_vars[ifac];

            // Store factor score center
            obs_fac_center[iobs][ifac] = f_score;

            // Determine number of quadrature points based on SE
            // Formula from legacy code: nquad = 1 + 2 * floor(se / var / threshold)
            // When SE is small, use 1 point (just the factor score)
            // When SE is large, use more points
            double ratio = f_se / f_var / threshold;
            int nq = 1 + 2 * static_cast<int>(std::floor(ratio));

            // Clamp to reasonable range
            if (nq < 1) nq = 1;
            if (nq > max_quad) nq = max_quad;

            // If SE is very large (> sqrt(factor_var)), use default quadrature centered at 0
            // This handles cases where factor score estimation failed
            if (f_se > std::sqrt(f_var)) {
                nq = max_quad;
                obs_fac_center[iobs][ifac] = 0.0;  // Center at prior mean
            }

            obs_nquad[iobs][ifac] = nq;
            needed_nquad.insert(nq);
        }

        // Note: Importance sampling weights are computed per integration point in CalcLkhd,
        // not per observation, because they depend on the quadrature node values.
        obs_weights[iobs] = 1.0;  // Default weight
    }

    // Pre-compute GH quadrature nodes/weights for all needed sizes
    adapt_nodes.clear();
    adapt_weights.clear();

    double sqrt2 = std::sqrt(2.0);
    double sqrt_pi = std::sqrt(M_PI);

    for (int nq : needed_nquad) {
        std::vector<double> nodes, weights;
        calcgausshermitequadrature(nq, nodes, weights);

        // Scale for standard normal: x_scaled = sqrt(2) * x, w_scaled = w / sqrt(pi)
        for (int i = 0; i < nq; i++) {
            nodes[i] *= sqrt2;
            weights[i] /= sqrt_pi;
        }

        adapt_nodes[nq] = nodes;
        adapt_weights[nq] = weights;
    }

    use_weights = true;
    use_adaptive = true;
}

void FactorModel::DisableAdaptiveQuadrature()
{
    use_adaptive = false;
    use_weights = false;  // Also reset weights since they may have been set by adaptive mode
    obs_nquad.clear();
    obs_fac_center.clear();
    obs_weights.clear();
    adapt_factor_var.clear();
    adapt_nodes.clear();
    adapt_weights.clear();
}

void FactorModel::SetParameterConstraints(const std::vector<bool>& fixed)
{
    if (fixed.size() != nparam) {
        throw std::runtime_error("Parameter constraint size mismatch");
    }
    param_fixed = fixed;
    nparam_free = 0;
    freeparlist.clear();
    freeparlist.reserve(nparam);
    for (int i = 0; i < nparam; i++) {
        if (!fixed[i]) {
            freeparlist.push_back(i);
            nparam_free++;
        }
    }

    // Initialize parameter vector with default values
    param.resize(nparam, 0.0);

    // Initialize factor variances to 1.0 by default
    for (int ifac = 0; ifac < nfac; ifac++) {
        param[ifac] = 1.0;
    }
}

void FactorModel::SetParameterConstraints(const std::vector<bool>& fixed,
                                          const std::vector<double>& fixed_values)
{
    if (fixed.size() != nparam) {
        throw std::runtime_error("Parameter constraint size mismatch");
    }
    if (fixed_values.size() != nparam) {
        throw std::runtime_error("Fixed values size mismatch");
    }

    param_fixed = fixed;
    nparam_free = 0;
    freeparlist.clear();
    freeparlist.reserve(nparam);
    for (int i = 0; i < nparam; i++) {
        if (!fixed[i]) {
            freeparlist.push_back(i);
            nparam_free++;
        }
    }

    // Initialize parameter vector with provided values
    // For free parameters, these are defaults (will be overwritten by optimizer)
    // For fixed parameters, these values are used permanently
    param = fixed_values;
}

void FactorModel::MapFreeToFull(const std::vector<double>& free_params)
{
    if (free_params.size() != nparam_free) {
        throw std::runtime_error("Free parameter size mismatch");
    }

    int ifree = 0;
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) {
            param[i] = free_params[ifree++];
        }
        // Fixed parameters keep their current values in param vector
    }
}

void FactorModel::ExtractFreeGradient(const std::vector<double>& full_grad,
                                      std::vector<double>& free_grad)
{
    free_grad.resize(nparam_free);
    int ifree = 0;
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) {
            free_grad[ifree++] = full_grad[i];
        }
    }
}

void FactorModel::ExtractFreeHessian(const std::vector<double>& full_hess,
                                     std::vector<double>& free_hess)
{
    // Extract submatrix corresponding to free parameters
    // Input: full_hess is upper triangle of nparam x nparam matrix
    // Output: free_hess is upper triangle of nparam_free x nparam_free matrix

    free_hess.clear();
    free_hess.reserve(nparam_free * (nparam_free + 1) / 2);

    std::vector<int> free_indices;
    for (int i = 0; i < nparam; i++) {
        if (!param_fixed[i]) free_indices.push_back(i);
    }

    for (int i = 0; i < nparam_free; i++) {
        for (int j = i; j < nparam_free; j++) {
            int full_i = free_indices[i];
            int full_j = free_indices[j];
            int full_idx = full_i * nparam + full_j;
            free_hess.push_back(full_hess[full_idx]);
        }
    }
}

double FactorModel::CalcLogLikelihood(const std::vector<double>& free_params)
{
    double logLkhd;
    std::vector<double> gradL, hessL;
    CalcLkhd(free_params, logLkhd, gradL, hessL, 1);  // flag=1: likelihood only
    return logLkhd;
}

void FactorModel::CalcLkhd(const std::vector<double>& free_params,
                           double& logLkhd,
                           std::vector<double>& gradL,
                           std::vector<double>& hessL,
                           int iflag)
{
    // ===== STEP 1: Map free parameters to full parameter vector =====
    MapFreeToFull(free_params);

    // ===== STEP 2: Extract factor parameters from param vector =====
    // Parameter organization: [factor_vars (nfac) | factor_corr (if correlated) | model_params (rest)]
    // For single mixture (nmix=1):
    //   - param[0..nfac-1] = factor variances (sigma^2)
    //   - param[nfac] = factor correlation (if fac_corr && nfac == 2)
    std::vector<double> factor_var(nfac);
    for (int ifac = 0; ifac < nfac; ifac++) {
        factor_var[ifac] = std::fabs(param[ifac]);  // Enforce positivity
    }

    // Factor correlation (for 2-factor correlated models)
    double rho = 0.0;
    double sqrt_1_minus_rho2 = 1.0;
    if (fac_corr && nfac == 2) {
        rho = param[nfac];  // Correlation parameter at index 2
        // Clamp to valid range to prevent numerical issues
        if (rho > 0.999) rho = 0.999;
        if (rho < -0.999) rho = -0.999;
        sqrt_1_minus_rho2 = std::sqrt(1.0 - rho * rho);
    }

    // Factor means (0 for single mixture)
    std::vector<double> factor_mean(nfac, 0.0);

    // ===== STEP 3: Initialize outputs =====
    logLkhd = 0.0;
    std::vector<double> full_gradL(nparam, 0.0);
    std::vector<double> full_hessL;
    if (iflag == 3) {
        full_hessL.resize(nparam * nparam, 0.0);
    }

    // ===== STEP 4: Compute default number of integration points =====
    // (may be overridden per-observation in adaptive mode)
    int nint_points_default = 1;
    for (int i = 0; i < nfac; i++) {
        nint_points_default *= nquad_points;
    }

    // Temporary storage
    // OPTIMIZATION: Pre-allocate modEval and modHess to maximum needed size
    // This avoids repeated resize operations inside the hot loop
    int max_model_params = 0;
    for (size_t imod = 0; imod < models.size(); imod++) {
        if (param_model_count[imod] > max_model_params) {
            max_model_params = param_model_count[imod];
        }
    }
    int max_ngrad = 1 + nfac + max_model_params;  // Maximum gradient size for any model
    int max_hess_dim = nfac + max_model_params;   // Maximum Hessian dimension
    std::vector<double> modEval(max_ngrad, 0.0);
    std::vector<double> modHess;
    if (iflag == 3) {
        modHess.resize(max_hess_dim * max_hess_dim, 0.0);
    }
    std::vector<double> totalgrad(nparam, 0.0);
    std::vector<double> totalhess;
    if (iflag == 3) totalhess.resize(nparam * nparam, 0.0);

    // Constant for quadrature transformation
    const double sqrt_2 = std::sqrt(2.0);

    std::vector<double> gradilk(nparam, 0.0);  // Gradient at integration point
    std::vector<double> hessilk;
    if (iflag == 3) hessilk.resize(nparam * nparam, 0.0);

    // Pre-allocate working vectors outside the loops for performance
    // These are reused across observations and integration points
    std::vector<double> fac_val(nfac);
    std::vector<double> type_weighted_grad(nparam, 0.0);
    std::vector<double> type_weighted_hess;
    if (iflag == 3) type_weighted_hess.resize(nparam * nparam, 0.0);
    std::vector<double> grad_this_type(nparam, 0.0);
    std::vector<double> hess_this_type;
    if (iflag == 3) hess_this_type.resize(nparam * nparam, 0.0);

    // OPTIMIZATION: Pre-allocate vectors for chain rule factors
    // These store sqrt(factor_var), quadrature nodes, and df/d(sigma^2) per factor
    // Computed once per integration point, reused throughout type/model loops
    std::vector<double> sigma_fac(nfac);
    std::vector<double> x_node_fac(nfac);
    std::vector<double> df_dsigma2(nfac);

    // ===== STEP 5: Loop over observations =====
    for (int iobs = 0; iobs < nobs; iobs++) {
        int iobs_offset = iobs * nvar;

        double totalprob = 0.0;
        if (iflag >= 2) std::fill(totalgrad.begin(), totalgrad.end(), 0.0);
        // OPTIMIZATION: Only zero free parameter entries (like legacy code)
        if (iflag == 3) {
            for (int fi = 0; fi < nparam_free; fi++) {
                int i = freeparlist[fi];
                for (int fj = fi; fj < nparam_free; fj++) {
                    int j = freeparlist[fj];
                    totalhess[i * nparam + j] = 0.0;
                }
            }
        }

        // ===== STEP 6: Multi-dimensional integration over factors =====
        // Determine per-observation integration settings
        int nint_points;
        std::vector<int> obs_nq(nfac);  // Number of quad points per factor
        std::vector<double> fac_center(nfac);  // Center for factor integration

        if (use_adaptive) {
            // Adaptive mode: use per-observation quadrature settings
            nint_points = 1;
            for (int ifac = 0; ifac < nfac; ifac++) {
                obs_nq[ifac] = obs_nquad[iobs][ifac];
                fac_center[ifac] = obs_fac_center[iobs][ifac];
                nint_points *= obs_nq[ifac];
            }
        } else {
            // Standard mode: use global quadrature settings
            nint_points = nint_points_default;
            for (int ifac = 0; ifac < nfac; ifac++) {
                obs_nq[ifac] = nquad_points;
                fac_center[ifac] = factor_mean[ifac];
            }
        }

        // Use "odometer" method: facint tracks indices for each dimension
        std::vector<int> facint(nfac, 0);

        for (int intpt = 0; intpt < nint_points; intpt++) {
            // Reset gradient/Hessian for this integration point
            if (iflag >= 2) std::fill(gradilk.begin(), gradilk.end(), 0.0);
            // OPTIMIZATION: Only zero free parameter entries (like legacy code)
            if (iflag == 3) {
                for (int fi = 0; fi < nparam_free; fi++) {
                    int i = freeparlist[fi];
                    for (int fj = fi; fj < nparam_free; fj++) {
                        int j = freeparlist[fj];
                        hessilk[i * nparam + j] = 0.0;
                    }
                }
            }

            // OPTIMIZATION: Pre-compute sigma_fac and x_node_fac once per integration point
            // These are reused for fac_val computation and gradient/Hessian chain rules
            for (int ifac = 0; ifac < nfac; ifac++) {
                sigma_fac[ifac] = std::sqrt(factor_var[ifac]);
                if (use_adaptive) {
                    int nq = obs_nq[ifac];
                    x_node_fac[ifac] = adapt_nodes.at(nq)[facint[ifac]];
                } else {
                    x_node_fac[ifac] = quad_nodes[facint[ifac]];
                }
            }

            // Get integration weight (product over dimensions)
            double probilk = 1.0;
            for (int ifac = 0; ifac < nfac; ifac++) {
                if (use_adaptive) {
                    int nq = obs_nq[ifac];
                    probilk *= adapt_weights.at(nq)[facint[ifac]];
                } else {
                    probilk *= quad_weights[facint[ifac]];
                }
            }

            // Compute factor values at this integration point (fac_val pre-allocated)
            // OPTIMIZATION: Use pre-computed sigma_fac and x_node_fac instead of recalculating
            if (fac_corr && nfac == 2 && !use_adaptive) {
                // Correlated 2-factor model: use Cholesky transformation
                // (Not supported in adaptive mode - requires independent factors)
                // f1 = σ1 * z1
                // f2 = ρ*σ2*z1 + σ2*√(1-ρ²)*z2
                // where z1, z2 are independent GH quadrature nodes
                fac_val[0] = sigma_fac[0] * x_node_fac[0] + factor_mean[0];
                fac_val[1] = rho * sigma_fac[1] * x_node_fac[0] + sigma_fac[1] * sqrt_1_minus_rho2 * x_node_fac[1] + factor_mean[1];
            } else {
                // Independent factors: standard transformation
                // In adaptive mode: center at factor score instead of 0
                for (int ifac = 0; ifac < nfac; ifac++) {
                    // Factor transformation for GH quadrature
                    // f = sigma * x + center (where center is factor_score in adaptive mode)
                    fac_val[ifac] = sigma_fac[ifac] * x_node_fac[ifac] + fac_center[ifac];
                }
            }

            // Apply importance sampling correction for adaptive quadrature
            // When centering at factor score instead of 0, we need to correct for
            // the different normalizing constant of the prior density.
            // Correction: exp(x²/2) × exp(-f²/(2σ²)) for each factor
            // where x is the GH node and f is the factor value at this point.
            if (use_adaptive) {
                for (int ifac = 0; ifac < nfac; ifac++) {
                    int nq = obs_nq[ifac];
                    double x_node = adapt_nodes.at(nq)[facint[ifac]];
                    double f = fac_val[ifac];
                    double var = factor_var[ifac];

                    // Importance sampling weight correction
                    // This accounts for centering at fac_center instead of prior mean (0)
                    double is_correction = std::exp(x_node * x_node / 2.0 - f * f / (2.0 * var));
                    probilk *= is_correction;
                }
            }

            // ===== STEP 7: Pre-compute chain rule factors for gradients/Hessians =====
            // sigma_fac and x_node_fac are already pre-computed at the start of the integration point loop
            // Here we just compute the chain rule factor for variance derivatives
            for (int ifac = 0; ifac < nfac; ifac++) {
                // Chain rule: df/d(sigma^2) = x_node / (2 * sigma)
                df_dsigma2[ifac] = x_node_fac[ifac] / (2.0 * sigma_fac[ifac]);
            }

            // ===== STEP 8: Evaluate all models with type mixture =====
            // For ntyp > 1: L = Σ_t π_t(f) × L_t(y|f) where L_t = Π_m L_m(y_m|f, intercept_t)
            // For ntyp == 1: Standard case (no type mixture)

            // Compute type probabilities (returns [1.0] for ntyp == 1)
            std::vector<double> type_probs = ComputeTypeProbabilities(fac_val);

            // Accumulators for type-weighted contributions (pre-allocated, zero here)
            double type_weighted_prob = 0.0;
            if (iflag >= 2) std::fill(type_weighted_grad.begin(), type_weighted_grad.end(), 0.0);
            // OPTIMIZATION: Only zero free parameter entries (like legacy code)
            if (iflag == 3) {
                for (int fi = 0; fi < nparam_free; fi++) {
                    int i = freeparlist[fi];
                    for (int fj = fi; fj < nparam_free; fj++) {
                        int j = freeparlist[fj];
                        type_weighted_hess[i * nparam + j] = 0.0;
                    }
                }
            }

            // Loop over types
            for (int ityp = 0; ityp < ntyp; ityp++) {
                double type_prob = type_probs[ityp];

                // Likelihood product for this type
                double prob_this_type = 1.0;

                // Gradient accumulators for this type (pre-allocated, zero here)
                if (iflag >= 2) std::fill(grad_this_type.begin(), grad_this_type.end(), 0.0);
                // OPTIMIZATION: Only zero free parameter entries (like legacy code)
                if (iflag == 3) {
                    for (int fi = 0; fi < nparam_free; fi++) {
                        int i = freeparlist[fi];
                        for (int fj = fi; fj < nparam_free; fj++) {
                            int j = freeparlist[fj];
                            hess_this_type[i * nparam + j] = 0.0;
                        }
                    }
                }

                // Evaluate all models for this type
                for (size_t imod = 0; imod < models.size(); imod++) {
                    int firstpar = param_model_start[imod];

                    // Get type-specific intercept for types > 1
                    // Type 1 (ityp=0) is reference with intercept = 0
                    double type_intercept = 0.0;
                    if (ityp > 0 && ntyp > 1) {
                        // Type intercepts are stored after base model parameters
                        // ityp-1 because type 1 is reference (no intercept)
                        int intercept_offset = GetTypeInterceptIndex(ityp - 1, imod);
                        type_intercept = param[firstpar + intercept_offset];
                    }

                    // If all parameters are fixed for this model, only compute likelihood (flag=1)
                    int model_flag = models[imod]->GetAllParamsFixed() ? 1 : iflag;

                    models[imod]->Eval(iobs_offset, data, param, firstpar, fac_val,
                                      modEval, modHess, model_flag, type_intercept);

                    // Multiply likelihood for this type
                    prob_this_type *= modEval[0];

                    // ===== STEP 9: Accumulate gradients for this type =====
                    if (iflag >= 2) {
                        // Gradients w.r.t. factor variances and correlation
                        if (fac_corr && nfac == 2) {
                            // Correlated 2-factor model: use Cholesky derivatives
                            // Use pre-computed sigma_fac and x_node_fac
                            double dL_df1 = modEval[1];
                            double dL_df2 = modEval[2];

                            grad_this_type[0] += dL_df1 * x_node_fac[0] / (2.0 * sigma_fac[0]);
                            double df2_dsigma2sq = (rho * x_node_fac[0] + sqrt_1_minus_rho2 * x_node_fac[1]) / (2.0 * sigma_fac[1]);
                            grad_this_type[1] += dL_df2 * df2_dsigma2sq;
                            double df2_drho = sigma_fac[1] * (x_node_fac[0] - rho * x_node_fac[1] / sqrt_1_minus_rho2);
                            grad_this_type[2] += dL_df2 * df2_drho;
                        } else {
                            // Independent factors: use pre-computed df_dsigma2
                            for (int ifac = 0; ifac < nfac; ifac++) {
                                grad_this_type[ifac] += modEval[ifac + 1] * df_dsigma2[ifac];
                            }
                        }

                        // Gradients w.r.t. model parameters (base params only, not type intercepts)
                        // param_model_count includes type intercepts, but modEval only has base params
                        int n_type_intercepts = (ntyp > 1) ? (ntyp - 1) : 0;
                        int base_param_count = param_model_count[imod] - n_type_intercepts;
                        for (int iparam = 0; iparam < base_param_count; iparam++) {
                            grad_this_type[firstpar + iparam] += modEval[1 + nfac + iparam];
                        }

                        // Gradient w.r.t. this type's intercept (only for types > 1)
                        // The gradient equals dL/d(linear_predictor) * 1
                        // For models with an intercept covariate (=1), this equals the base intercept gradient
                        if (ntyp > 1 && ityp > 0) {
                            int type_intercept_idx = firstpar + base_param_count + (ityp - 1);
                            // modEval[1 + nfac + 0] is gradient w.r.t. first covariate coefficient
                            // When first covariate is intercept (=1), this equals dL/d(linear_predictor)
                            grad_this_type[type_intercept_idx] += modEval[1 + nfac + 0];
                        }
                    }

                    // ===== STEP 9: Accumulate Hessians for this type =====
                    if (iflag == 3 && modHess.size() > 0) {
                        // Note: modHess from Model::Eval() only includes base parameters,
                        // not type-specific intercepts
                        int n_type_intercepts_hess = (ntyp > 1) ? (ntyp - 1) : 0;
                        int base_param_count_hess = param_model_count[imod] - n_type_intercepts_hess;
                        int nDimModHess = nfac + base_param_count_hess;

                        if (fac_corr && nfac == 2) {
                            // ===== Correlated 2-factor Hessian =====
                            // Use pre-computed sigma_fac and x_node_fac (z1, z2)
                            double df1_dsigma1sq = x_node_fac[0] / (2.0 * sigma_fac[0]);
                            double df2_dsigma2sq = (rho * x_node_fac[0] + sqrt_1_minus_rho2 * x_node_fac[1]) / (2.0 * sigma_fac[1]);
                            double df2_drho = sigma_fac[1] * (x_node_fac[0] - rho * x_node_fac[1] / sqrt_1_minus_rho2);
                            double d2f2_drho2 = -sigma_fac[1] * x_node_fac[1] / (sqrt_1_minus_rho2 * sqrt_1_minus_rho2 * sqrt_1_minus_rho2);
                            double d2f2_drho_dsigma2sq = (x_node_fac[0] - rho * x_node_fac[1] / sqrt_1_minus_rho2) / (2.0 * sigma_fac[1]);

                            double d2L_df1df1 = modHess[0 * nDimModHess + 0];
                            double d2L_df1df2 = modHess[0 * nDimModHess + 1];
                            double d2L_df2df2 = modHess[1 * nDimModHess + 1];
                            double dL_df1 = modEval[1];
                            double dL_df2 = modEval[2];

                            double sigma1_cubed = sigma_fac[0] * sigma_fac[0] * sigma_fac[0];
                            double sigma2_cubed = sigma_fac[1] * sigma_fac[1] * sigma_fac[1];
                            double d2f1_dsigma1sq_sq = -x_node_fac[0] / (4.0 * sigma1_cubed);
                            hess_this_type[0 * nparam + 0] += d2L_df1df1 * df1_dsigma1sq * df1_dsigma1sq
                                                              + dL_df1 * d2f1_dsigma1sq_sq;
                            hess_this_type[0 * nparam + 1] += d2L_df1df2 * df1_dsigma1sq * df2_dsigma2sq;
                            hess_this_type[0 * nparam + 2] += d2L_df1df2 * df1_dsigma1sq * df2_drho;

                            double d2f2_dsigma2sq_sq = -(rho * x_node_fac[0] + sqrt_1_minus_rho2 * x_node_fac[1]) / (4.0 * sigma2_cubed);
                            hess_this_type[1 * nparam + 1] += d2L_df2df2 * df2_dsigma2sq * df2_dsigma2sq
                                                              + dL_df2 * d2f2_dsigma2sq_sq;
                            hess_this_type[1 * nparam + 2] += d2L_df2df2 * df2_dsigma2sq * df2_drho
                                                              + dL_df2 * d2f2_drho_dsigma2sq;
                            hess_this_type[2 * nparam + 2] += d2L_df2df2 * df2_drho * df2_drho
                                                              + dL_df2 * d2f2_drho2;

                            // Use base_param_count_hess for loop bounds since modHess
                            // only contains base parameters, not type-specific intercepts
                            for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                int param_idx = firstpar + iparam;
                                double d2L_df1_dparam = modHess[0 * nDimModHess + (nfac + iparam)];
                                double d2L_df2_dparam = modHess[1 * nDimModHess + (nfac + iparam)];
                                hess_this_type[0 * nparam + param_idx] += d2L_df1_dparam * df1_dsigma1sq;
                                hess_this_type[1 * nparam + param_idx] += d2L_df2_dparam * df2_dsigma2sq;
                                hess_this_type[2 * nparam + param_idx] += d2L_df2_dparam * df2_drho;
                            }

                            for (int i = 0; i < base_param_count_hess; i++) {
                                for (int j = i; j < base_param_count_hess; j++) {
                                    int modhess_idx = (nfac + i) * nDimModHess + (nfac + j);
                                    int full_i = firstpar + i;
                                    int full_j = firstpar + j;
                                    hess_this_type[full_i * nparam + full_j] += modHess[modhess_idx];
                                }
                            }
                        } else {
                            // ===== Independent factors Hessian =====
                            // Use pre-computed df_dsigma2 (chain rule factor) for efficiency
                            for (int i = 0; i < nDimModHess; i++) {
                                for (int j = i; j < nDimModHess; j++) {
                                    int modhess_idx = i * nDimModHess + j;
                                    int full_i = (i < nfac) ? i : (firstpar + i - nfac);
                                    int full_j = (j < nfac) ? j : (firstpar + j - nfac);
                                    int full_idx = full_i * nparam + full_j;

                                    // Use pre-computed df_dsigma2[k] = x_node[k] / (2 * sigma[k])
                                    double chain_factor = 1.0;
                                    if (i < nfac) {
                                        chain_factor *= df_dsigma2[i];
                                    }
                                    if (j < nfac) {
                                        chain_factor *= df_dsigma2[j];
                                    }

                                    hess_this_type[full_idx] += modHess[modhess_idx] * chain_factor;

                                    // Second derivative term for factor variance diagonal
                                    if (i < nfac && i == j) {
                                        double sigma_cubed = sigma_fac[i] * sigma_fac[i] * sigma_fac[i];
                                        double second_deriv_factor = -x_node_fac[i] / (4.0 * sigma_cubed);
                                        hess_this_type[full_idx] += modEval[i + 1] * second_deriv_factor;
                                    }
                                }
                            }
                        }

                        // ===== STEP 9b: Add Hessian terms for type-specific intercepts =====
                        // Type-specific intercepts have the same derivative structure as base intercept.
                        // The Hessian contribution mirrors that of the base intercept (index 0 in model params).
                        if (ntyp > 1 && ityp > 0) {
                            int type_intercept_idx = firstpar + base_param_count_hess + (ityp - 1);

                            // Diagonal term: d²L/d(type_intercept)² = d²L/d(intercept)²
                            // In modHess, intercept is at position (nfac, nfac) -> index nfac * nDimModHess + nfac
                            double hess_intercept_diag = modHess[nfac * nDimModHess + nfac];
                            hess_this_type[type_intercept_idx * nparam + type_intercept_idx] += hess_intercept_diag;

                            // Cross-terms with factor variances: d²L/d(type_intercept)d(factor_var_k)
                            for (int k = 0; k < nfac; k++) {
                                // modHess[k * nDimModHess + nfac] = d²L/df_k d(intercept)
                                double hess_fk_intercept = modHess[k * nDimModHess + nfac];

                                // Use pre-computed df_dsigma2[k] = x_node / (2*sigma)
                                // Hessian cross-term: d²L/d(factor_var_k)d(type_intercept)
                                int full_idx = k * nparam + type_intercept_idx;
                                hess_this_type[full_idx] += hess_fk_intercept * df_dsigma2[k];
                            }

                            // Cross-terms with base model parameters
                            for (int iparam = 0; iparam < base_param_count_hess; iparam++) {
                                int param_idx = firstpar + iparam;
                                // modHess[(nfac + iparam) * nDimModHess + nfac] = d²L/d(param_iparam)d(intercept)
                                // Note: Hessian is symmetric, upper triangular stored
                                int mh_row = nfac;  // intercept row
                                int mh_col = nfac + iparam;  // param column
                                if (mh_row > mh_col) std::swap(mh_row, mh_col);
                                double hess_param_intercept = modHess[mh_row * nDimModHess + mh_col];

                                // Add to cross-term (type_intercept, param)
                                int full_i = std::min(type_intercept_idx, param_idx);
                                int full_j = std::max(type_intercept_idx, param_idx);
                                hess_this_type[full_i * nparam + full_j] += hess_param_intercept;
                            }
                        }
                    }
                } // End of model loop for this type

                // ===== Weight by type probability and accumulate =====
                type_weighted_prob += type_prob * prob_this_type;

                if (iflag >= 2) {
                    // Gradient of type-weighted likelihood:
                    // d(π_t * L_t)/dθ = (dπ_t/dθ) * L_t + π_t * (dL_t/dθ)
                    // For model parameters: dπ_t/dθ = 0, so just π_t * (dL_t/dθ)
                    for (int ipar = 0; ipar < nparam; ipar++) {
                        type_weighted_grad[ipar] += type_prob * prob_this_type * grad_this_type[ipar];
                    }

                    // Gradient w.r.t. type model loadings (for types > 1)
                    // π_t = exp(η_t) / (1 + Σ_s exp(η_s)) where η_t = Σ_k λ_{t,k} * f_k
                    // dπ_ityp/dλ_{s,k} = π_ityp * (δ_{ityp,s} - π_s) * f_k
                    if (ntyp > 1) {
                        for (int t = 0; t < ntyp - 1; t++) {  // Types 2, 3, ... (t+1 in 1-indexed)
                            for (int k = 0; k < nfac; k++) {
                                int param_idx = GetTypeLoadingIndex(t, k);
                                // ityp is 0-based (0 = type 1, 1 = type 2, etc.)
                                // t is 0-based index for types with loadings (0 = type 2, etc.)
                                // dπ_ityp/dη_{t+1} = π_ityp * (δ_{ityp,t+1} - π_{t+1})
                                double dpi_dlambda = type_probs[ityp] *
                                    ((ityp == t + 1 ? 1.0 : 0.0) - type_probs[t + 1]) * fac_val[k];
                                type_weighted_grad[param_idx] += dpi_dlambda * prob_this_type;
                            }
                        }

                        // Additional gradient for factor variance through type probabilities
                        // dπ_t/dσ²_k = dπ_t/df_k * df_k/dσ²_k
                        // where df_k/dσ²_k = x_q / (2 * σ_k) since f = σ * x and dσ/dσ² = 1/(2σ)
                        // and dπ_t/df_k = Σ_s π_t * (δ_{t,s} - π_s) * λ_{s,k}
                        for (int k = 0; k < nfac; k++) {
                            // Compute dπ_t/df_k = Σ_s π_t * (δ_{t,s} - π_s) * λ_{s,k}
                            double dpi_df_k = 0.0;
                            for (int s = 0; s < ntyp - 1; s++) {
                                int lambda_idx = GetTypeLoadingIndex(s, k);
                                double lambda_sk = param[lambda_idx];
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                dpi_df_k += type_probs[ityp] * (delta_ts - type_probs[s + 1]) * lambda_sk;
                            }

                            // Add (dπ_t/dσ²_k) * L_t to factor variance gradient
                            // Use pre-computed df_dsigma2[k]
                            type_weighted_grad[k] += dpi_df_k * df_dsigma2[k] * prob_this_type;
                        }
                    }
                }

                if (iflag == 3) {
                    // ===== STEP 9c: Hessian for model parameters (original formula) =====
                    // d²L_t/dθdφ = L_t * [d²(log L_t)/dθdφ + d(log L_t)/dθ * d(log L_t)/dφ]
                    // Weighted by π_t: contribution to d²L_mix/dθdφ is π_t * d²L_t/dθdφ
                    // OPTIMIZATION: Only compute entries for free parameters (like legacy code)
                    for (int fi = 0; fi < nparam_free; fi++) {
                        int i = freeparlist[fi];
                        for (int fj = fi; fj < nparam_free; fj++) {
                            int j = freeparlist[fj];
                            int idx = i * nparam + j;
                            // Accumulate Hessian weighted by type probability
                            type_weighted_hess[idx] += type_prob * prob_this_type *
                                (hess_this_type[idx] + grad_this_type[i] * grad_this_type[j]);
                        }
                    }

                    // ===== STEP 9d: Additional Hessian terms for type model loadings =====
                    // The mixture L_mix = Σ_t π_t * L_t depends on type loadings through π_t
                    // d²L_mix/dθdφ = Σ_t [d²π_t/dθdφ * L_t + dπ_t/dθ * dL_t/dφ + dπ_t/dφ * dL_t/dθ + π_t * d²L_t/dθdφ]
                    // The last term is already handled above.
                    // For type loading λ_{s,k}: dπ_t/dλ = π_t * (δ_{t,s} - π_s) * f_k
                    // For model parameters: dπ_t/dθ = 0
                    if (ntyp > 1) {
                        // 1. Second derivative terms: d²π_t/dλdλ * L_t
                        // d²π_t/dη_s dη_r = π_t * (δ_{t,s} - π_s) * (δ_{t,r} - π_r) - π_t * π_s * (δ_{s,r} - π_r)
                        for (int s = 0; s < ntyp - 1; s++) {
                            for (int r = s; r < ntyp - 1; r++) {  // r >= s for upper triangle
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                double delta_tr = (ityp == r + 1) ? 1.0 : 0.0;
                                double delta_sr = (s == r) ? 1.0 : 0.0;

                                double d2pi_deta = type_probs[ityp] * (delta_ts - type_probs[s + 1])
                                                                    * (delta_tr - type_probs[r + 1])
                                                 - type_probs[ityp] * type_probs[s + 1]
                                                                    * (delta_sr - type_probs[r + 1]);

                                for (int k = 0; k < nfac; k++) {
                                    for (int l = (s == r ? k : 0); l < nfac; l++) {
                                        int param_idx_sk = GetTypeLoadingIndex(s, k);
                                        int param_idx_rl = GetTypeLoadingIndex(r, l);

                                        int pi = std::min(param_idx_sk, param_idx_rl);
                                        int pj = std::max(param_idx_sk, param_idx_rl);
                                        int idx = pi * nparam + pj;

                                        double d2pi_dlambda = d2pi_deta * fac_val[k] * fac_val[l];
                                        type_weighted_hess[idx] += d2pi_dlambda * prob_this_type;
                                    }
                                }
                            }
                        }

                        // 2. Cross terms: dπ_t/dλ * dL_t/dθ (for λ × model_param)
                        // Note: dL_t/dθ = L_t * d(log L_t)/dθ = prob_this_type * grad_this_type
                        for (int s = 0; s < ntyp - 1; s++) {
                            for (int k = 0; k < nfac; k++) {
                                int type_loading_idx = GetTypeLoadingIndex(s, k);
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                double dpi_dlambda = type_probs[ityp] * (delta_ts - type_probs[s + 1]) * fac_val[k];

                                // Cross with model parameters (not type loadings)
                                for (int ipar = 0; ipar < nparam; ipar++) {
                                    // Skip if ipar is a type loading (handled in section 1)
                                    if (ntyp_param > 0 && ipar >= type_param_start &&
                                        ipar < type_param_start + ntyp_param) {
                                        continue;
                                    }

                                    int pi = std::min(type_loading_idx, ipar);
                                    int pj = std::max(type_loading_idx, ipar);
                                    int idx = pi * nparam + pj;

                                    // dπ_t/dλ * dL_t/dθ = dπ_t/dλ * L_t * d(log L_t)/dθ
                                    type_weighted_hess[idx] += dpi_dlambda * prob_this_type * grad_this_type[ipar];
                                }
                            }
                        }

                        // 3. Hessian terms for factor variance through type probabilities
                        // dπ_t/dσ²_k = dπ_t/df_k * df_k/dσ²_k
                        // where df/dσ² = x_q / (2σ) since f = σ * x and dσ/dσ² = 1/(2σ)
                        // Use pre-computed df_dsigma2[k]
                        for (int k = 0; k < nfac; k++) {
                            // Compute dπ_t/df_k = Σ_s π_t * (δ_{t,s} - π_s) * λ_{s,k}
                            double dpi_df_k = 0.0;
                            for (int s = 0; s < ntyp - 1; s++) {
                                int lambda_idx = GetTypeLoadingIndex(s, k);
                                double lambda_sk = param[lambda_idx];
                                double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                dpi_df_k += type_probs[ityp] * (delta_ts - type_probs[s + 1]) * lambda_sk;
                            }
                            double dpi_dsigma2_k = dpi_df_k * df_dsigma2[k];

                            // 3a. Cross terms: dπ_t/dσ²_k * dL_t/dθ (for σ² × model_param)
                            for (int ipar = 0; ipar < nparam; ipar++) {
                                // Skip type loadings (handled separately below)
                                if (ntyp_param > 0 && ipar >= type_param_start &&
                                    ipar < type_param_start + ntyp_param) {
                                    continue;
                                }
                                // Skip factor variances (handled in 3b)
                                if (ipar < nfac) continue;

                                int pi = std::min(k, ipar);
                                int pj = std::max(k, ipar);
                                int idx = pi * nparam + pj;

                                // dπ_t/dσ² * dL_t/dθ
                                type_weighted_hess[idx] += dpi_dsigma2_k * prob_this_type * grad_this_type[ipar];
                            }

                            // 3b. Cross terms: dπ_t/dσ²_k * dL_t/dσ²_l (for σ²_k × σ²_l)
                            // Plus second derivative d²π_t/dσ²_k dσ²_l * L_t
                            // Use pre-computed df_dsigma2[l]
                            for (int l = k; l < nfac; l++) {
                                // Compute dπ_t/df_l
                                double dpi_df_l = 0.0;
                                for (int s = 0; s < ntyp - 1; s++) {
                                    int lambda_idx = GetTypeLoadingIndex(s, l);
                                    double lambda_sl = param[lambda_idx];
                                    double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                    dpi_df_l += type_probs[ityp] * (delta_ts - type_probs[s + 1]) * lambda_sl;
                                }
                                double dpi_dsigma2_l = dpi_df_l * df_dsigma2[l];

                                int idx = k * nparam + l;

                                // Cross term: dπ_t/dσ²_k * dL_t/dσ²_l + dπ_t/dσ²_l * dL_t/dσ²_k
                                // Note: dL_t/dσ² = L_t * d(log L_t)/dσ² = prob_this_type * grad_this_type[factor_idx]
                                if (k == l) {
                                    // Diagonal: 2 * dπ_t/dσ² * dL_t/dσ² (but we only add once due to symmetry in sum)
                                    type_weighted_hess[idx] += dpi_dsigma2_k * prob_this_type * grad_this_type[l];
                                } else {
                                    // Off-diagonal: add both cross terms
                                    type_weighted_hess[idx] += dpi_dsigma2_k * prob_this_type * grad_this_type[l];
                                    type_weighted_hess[idx] += dpi_dsigma2_l * prob_this_type * grad_this_type[k];
                                }

                                // Second derivative: d²π_t/dσ²_k dσ²_l * L_t
                                // d²π_t/df_k df_l = Σ_s,r (d²π_t/dη_s dη_r) * λ_{s,k} * λ_{r,l}
                                double d2pi_df_k_df_l = 0.0;
                                for (int s = 0; s < ntyp - 1; s++) {
                                    for (int r = 0; r < ntyp - 1; r++) {
                                        int lambda_sk_idx = GetTypeLoadingIndex(s, k);
                                        int lambda_rl_idx = GetTypeLoadingIndex(r, l);
                                        double lambda_sk = param[lambda_sk_idx];
                                        double lambda_rl = param[lambda_rl_idx];

                                        double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                        double delta_tr = (ityp == r + 1) ? 1.0 : 0.0;
                                        double delta_sr = (s == r) ? 1.0 : 0.0;

                                        double d2pi_deta_s_r = type_probs[ityp] * (delta_ts - type_probs[s + 1])
                                                                                * (delta_tr - type_probs[r + 1])
                                                             - type_probs[ityp] * type_probs[s + 1]
                                                                                * (delta_sr - type_probs[r + 1]);
                                        d2pi_df_k_df_l += d2pi_deta_s_r * lambda_sk * lambda_rl;
                                    }
                                }
                                double d2pi_dsigma2_k_l = d2pi_df_k_df_l * df_dsigma2[k] * df_dsigma2[l];
                                type_weighted_hess[idx] += d2pi_dsigma2_k_l * prob_this_type;
                            }

                            // 3c. Cross terms: dπ_t/dσ²_k * dL_t/dλ_{s,l} + d²π_t/dσ²_k dλ_{s,l} * L_t
                            // (for σ²_k × λ_{s,l})
                            for (int s = 0; s < ntyp - 1; s++) {
                                for (int l = 0; l < nfac; l++) {
                                    int type_loading_idx = GetTypeLoadingIndex(s, l);

                                    // Note: dL_t/dλ = 0 (type-specific likelihoods don't depend on type loadings)
                                    // So we only need d²π_t/dσ²_k dλ_{s,l} * L_t

                                    // d²π_t/df_k dλ_{s,l} = d²π_t/df_k dη_s * f_l + dπ_t/dη_s * δ_{k,l} * df/dσ² * dσ²/dσ²
                                    //                     = (d²π_t/dη dη) * λ_{?,k} * f_l + ... (complex)
                                    // Simpler: d(dπ_t/dλ_{s,l})/dσ²_k = d(π_t * (δ_{t,s} - π_s) * f_l)/dσ²_k
                                    //        = (dπ_t/dσ²_k) * (δ_{t,s} - π_s) * f_l
                                    //          + π_t * (-dπ_s/dσ²_k) * f_l
                                    //          + π_t * (δ_{t,s} - π_s) * df_l/dσ²_k

                                    double delta_ts = (ityp == s + 1) ? 1.0 : 0.0;
                                    double delta_kl = (k == l) ? 1.0 : 0.0;

                                    // df_l/dσ²_k = 0 if k != l, else x_node_k / (2 * σ_k)
                                    double df_l_dsigma2_k = delta_kl * df_dsigma2[k];

                                    // dπ_s/dσ²_k (for the type s+1)
                                    double dpi_s_df_k = 0.0;
                                    for (int r = 0; r < ntyp - 1; r++) {
                                        int lambda_rk_idx = GetTypeLoadingIndex(r, k);
                                        double lambda_rk = param[lambda_rk_idx];
                                        double delta_s1_r1 = (s == r) ? 1.0 : 0.0;
                                        dpi_s_df_k += type_probs[s + 1] * (delta_s1_r1 - type_probs[r + 1]) * lambda_rk;
                                    }
                                    double dpi_s_dsigma2_k = dpi_s_df_k * df_dsigma2[k];

                                    // d²(π_t*(δ_{t,s}-π_s)*f_l)/dσ²_k
                                    double d2_term = dpi_dsigma2_k * (delta_ts - type_probs[s + 1]) * fac_val[l]
                                                   - type_probs[ityp] * dpi_s_dsigma2_k * fac_val[l]
                                                   + type_probs[ityp] * (delta_ts - type_probs[s + 1]) * df_l_dsigma2_k;

                                    int pi = std::min(k, type_loading_idx);
                                    int pj = std::max(k, type_loading_idx);
                                    int idx = pi * nparam + pj;

                                    type_weighted_hess[idx] += d2_term * prob_this_type;
                                }
                            }
                        }
                    }
                }
            } // End of type loop

            // Use type-weighted values for the integration point
            // probilk already contains w_q (quadrature weight), multiply by mixture likelihood
            probilk *= type_weighted_prob;

            // gradilk needs to be d(log L_mixture)/dθ = (1/L_mixture) * dL_mixture/dθ
            // type_weighted_grad = dL_mixture/dθ = Σ_t π_t * dL_t/dθ
            // type_weighted_prob = L_mixture = Σ_t π_t * L_t
            if (iflag >= 2 && type_weighted_prob > 1e-100) {
                for (int ipar = 0; ipar < nparam; ipar++) {
                    gradilk[ipar] = type_weighted_grad[ipar] / type_weighted_prob;
                }
            }
            if (iflag == 3 && type_weighted_prob > 1e-100) {
                // hessilk needs d²(log L_mixture)/dθ² which involves complex chain rule
                // For now, convert to log scale: H_log = (1/L)*H - (1/L²)*g*g^T
                // OPTIMIZATION: Only compute entries for free parameters
                double L = type_weighted_prob;
                for (int fi = 0; fi < nparam_free; fi++) {
                    int i = freeparlist[fi];
                    for (int fj = fi; fj < nparam_free; fj++) {
                        int j = freeparlist[fj];
                        int idx = i * nparam + j;
                        hessilk[idx] = type_weighted_hess[idx] / L
                            - (type_weighted_grad[i] * type_weighted_grad[j]) / (L * L);
                    }
                }
            }

            // ===== STEP 10: Accumulate weighted contributions =====
            totalprob += probilk;

            if (iflag >= 2) {
                for (int ipar = 0; ipar < nparam; ipar++) {
                    totalgrad[ipar] += probilk * gradilk[ipar];
                }
            }

            if (iflag == 3) {
                // OPTIMIZATION: Only compute entries for free parameters (like legacy code)
                for (int fi = 0; fi < nparam_free; fi++) {
                    int i = freeparlist[fi];
                    for (int fj = fi; fj < nparam_free; fj++) {
                        int j = freeparlist[fj];
                        int idx = i * nparam + j;
                        totalhess[idx] += probilk * (hessilk[idx] + gradilk[i] * gradilk[j]);
                    }
                }
            }

            // ===== STEP 11: Increment integration point indices (odometer style) =====
            if (intpt < nint_points - 1) {
                for (int ifac = 0; ifac < nfac; ifac++) {
                    facint[ifac]++;
                    if (facint[ifac] < obs_nq[ifac]) break;
                    facint[ifac] = 0;
                }
            }
        }

        // ===== STEP 12: Compute log-likelihood for this observation =====
        // Apply observation weight (default 1.0 if no weights set)
        double obs_weight = use_weights ? obs_weights[iobs] : 1.0;

        if (totalprob > 1e-100) {
            logLkhd += obs_weight * std::log(totalprob);

            // Gradient: d(log L)/dθ = (1/L) * dL/dθ (weighted by observation)
            if (iflag >= 2) {
                for (int ipar = 0; ipar < nparam; ipar++) {
                    full_gradL[ipar] += obs_weight * totalgrad[ipar] / totalprob;
                }
            }

            // Hessian: d²(log L)/dθ² = (1/L)*d²L/dθ² - (1/L²)*(dL/dθ)² (weighted)
            // OPTIMIZATION: Only compute entries for free parameters (like legacy code)
            if (iflag == 3) {
                for (int fi = 0; fi < nparam_free; fi++) {
                    int i = freeparlist[fi];
                    for (int fj = fi; fj < nparam_free; fj++) {
                        int j = freeparlist[fj];
                        int idx = i * nparam + j;
                        double term1 = totalhess[idx] / totalprob;
                        double term2 = (totalgrad[i] * totalgrad[j]) / (totalprob * totalprob);
                        full_hessL[idx] += obs_weight * (term1 - term2);
                    }
                }
            }
        } else {
            // Numerical underflow
            logLkhd += obs_weight * (-1e10);
        }
    }

    // ===== STEP 13: Extract free parameter gradients/Hessians =====
    if (iflag >= 2) {
        ExtractFreeGradient(full_gradL, gradL);
    }
    if (iflag == 3) {
        ExtractFreeHessian(full_hessL, hessL);
    }
}

int FactorModel::GetFactorVarianceIndex(int imix, int ifac)
{
    // Factor variances come first in parameter vector
    return imix * nfac + ifac;
}

int FactorModel::GetFactorMeanIndex(int imix, int ifac)
{
    // Factor means come after variances (if nmix > 1)
    if (nmix == 1) return -1;  // No means for single mixture
    return nmix * nfac + imix * nfac + ifac;
}

int FactorModel::GetMixtureWeightIndex(int imix)
{
    // Mixture weights come after factor parameters
    if (nmix == 1) return -1;  // No weights for single mixture
    return nmix * nfac * 2 + imix;
}

void FactorModel::SetModelParameters(const std::vector<double>& params)
{
    if (params.size() != nparam) {
        throw std::runtime_error("Parameter vector size mismatch in SetModelParameters");
    }
    param = params;
}

void FactorModel::CalcLkhdSingleObs(int iobs,
                                    const std::vector<double>& factor_values,
                                    const std::vector<double>& model_params,
                                    double& logLkhd,
                                    std::vector<double>& gradL,
                                    std::vector<double>& hessL,
                                    int iflag)
{
    // ===== Factor score estimation mode =====
    // Evaluates likelihood for a single observation at given factor values.
    // No quadrature - factor values are the parameters we're optimizing.
    // Includes factor prior: L = p(y|f,θ) * φ(f|0,σ²)

    if (factor_values.size() != nfac) {
        throw std::runtime_error("Factor values size mismatch");
    }
    if (model_params.size() != nparam) {
        throw std::runtime_error("Model parameters size mismatch");
    }

    // Set the model parameters
    param = model_params;

    // Extract factor variances from param vector
    std::vector<double> factor_var(nfac);
    for (int ifac = 0; ifac < nfac; ifac++) {
        factor_var[ifac] = std::fabs(param[ifac]);  // Enforce positivity
    }

    // Initialize outputs
    logLkhd = 0.0;
    if (iflag >= 2) {
        gradL.resize(nfac, 0.0);
        std::fill(gradL.begin(), gradL.end(), 0.0);
    }
    if (iflag == 3) {
        // Upper triangle storage for nfac x nfac Hessian
        hessL.resize(nfac * (nfac + 1) / 2, 0.0);
        std::fill(hessL.begin(), hessL.end(), 0.0);
    }

    int iobs_offset = iobs * nvar;

    // ===== Part 1: Observation likelihood p(y|f,θ) =====
    // We compute log p(y|f) = sum_m log p_m(y_m|f) directly for numerical stability
    double logObs = 0.0;  // Sum of log-likelihoods
    std::vector<double> gradLogObs(nfac, 0.0);  // d(log L_obs)/df = sum_m d(log L_m)/df
    std::vector<double> hessLogObs;
    if (iflag == 3) {
        hessLogObs.resize(nfac * nfac, 0.0);
    }

    // Temporary storage for model evaluation
    std::vector<double> modEval;
    std::vector<double> modHess;

    // Evaluate all models at the given factor values
    for (size_t imod = 0; imod < models.size(); imod++) {
        int firstpar = param_model_start[imod];

        // For factor score estimation, we only need derivatives w.r.t. factor values
        // Pass iflag to get gradients/Hessians w.r.t. factors
        models[imod]->Eval(iobs_offset, data, param, firstpar, factor_values,
                          modEval, modHess, iflag);

        // Get model likelihood
        double Lm = modEval[0];
        if (Lm < 1e-100) {
            Lm = 1e-100;  // Prevent log(0)
        }
        logObs += std::log(Lm);

        // Accumulate gradients: d(log L_obs)/df = sum_m (1/L_m) * dL_m/df
        if (iflag >= 2) {
            for (int ifac = 0; ifac < nfac; ifac++) {
                // modEval[ifac + 1] is dL_m/df_ifac from this model
                gradLogObs[ifac] += modEval[ifac + 1] / Lm;
            }
        }

        // Accumulate Hessians: d²(log L_m)/df_i df_j = (1/L_m)*d²L_m/df² - (1/L_m²)*(dL_m/df_i)*(dL_m/df_j)
        if (iflag == 3 && modHess.size() > 0) {
            int nDimModHess = nfac + param_model_count[imod];
            // Extract only the factor-factor part of the Hessian
            for (int i = 0; i < nfac; i++) {
                for (int j = i; j < nfac; j++) {
                    int modhess_idx = i * nDimModHess + j;
                    int hess_idx = i * nfac + j;

                    double dLi = modEval[i + 1];  // dL_m/df_i
                    double dLj = modEval[j + 1];  // dL_m/df_j
                    double d2L = modHess[modhess_idx];  // d²L_m/df_i df_j

                    // d²(log L_m)/df_i df_j = d²L_m/df_i df_j / L_m - (dL_m/df_i * dL_m/df_j) / L_m²
                    hessLogObs[hess_idx] += d2L / Lm - (dLi * dLj) / (Lm * Lm);
                }
            }
        }
    }

    // ===== Part 2: Factor prior φ(f|0,σ²) =====
    // log φ(f|0,σ²) = -0.5 * Σ [f²/σ² + log(2π*σ²)]
    // = -0.5 * Σ [f²/σ² + log(2π) + log(σ²)]
    double logPrior = 0.0;
    std::vector<double> gradPrior(nfac, 0.0);
    std::vector<double> hessPrior(nfac, 0.0);  // Diagonal only for uncorrelated factors

    for (int ifac = 0; ifac < nfac; ifac++) {
        double f = factor_values[ifac];
        double sigma2 = factor_var[ifac];

        // Log-prior contribution
        logPrior -= 0.5 * (f * f / sigma2 + std::log(2.0 * M_PI * sigma2));

        // Gradient: d(logPrior)/df = -f/σ²
        if (iflag >= 2) {
            gradPrior[ifac] = -f / sigma2;
        }

        // Hessian: d²(logPrior)/df² = -1/σ²
        if (iflag == 3) {
            hessPrior[ifac] = -1.0 / sigma2;
        }
    }

    // ===== Part 3: Combine log-likelihood and log-prior =====
    // log(L_total) = log(L_obs) + log(L_prior)
    logLkhd = logObs + logPrior;

    // Gradient: d(log L_total)/df = d(log L_obs)/df + d(log L_prior)/df
    if (iflag >= 2) {
        for (int ifac = 0; ifac < nfac; ifac++) {
            gradL[ifac] = gradLogObs[ifac] + gradPrior[ifac];
        }
    }

    // Hessian: d²(log L_total)/df² = d²(log L_obs)/df² + d²(log L_prior)/df²
    if (iflag == 3) {
        int hess_out_idx = 0;
        for (int i = 0; i < nfac; i++) {
            for (int j = i; j < nfac; j++) {
                int hess_idx = i * nfac + j;
                double obs_hess = hessLogObs[hess_idx];

                // Add prior Hessian (diagonal only for uncorrelated factors)
                double prior_hess = (i == j) ? hessPrior[i] : 0.0;

                hessL[hess_out_idx++] = obs_hess + prior_hess;
            }
        }
    }
}

int FactorModel::GetTypeLoadingIndex(int ityp, int ifac)
{
    // Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // ifac: 0-based factor index
    // Returns: index in parameter vector for lambda_{ityp+2, ifac+1}
    if (ntyp <= 1 || type_param_start < 0) {
        return -1;  // No type parameters
    }
    return type_param_start + ityp * nfac + ifac;
}

int FactorModel::GetTypeInterceptIndex(int ityp, int model_idx)
{
    // Type-specific intercepts are stored in each model's parameter block
    // They come after the model's base parameters
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // model_idx: 0-based model index
    //
    // Note: The actual implementation depends on how model parameters are organized
    // This returns the index within the model's parameter block
    // The caller needs to add param_model_start[model_idx] to get the full index
    if (ntyp <= 1) {
        return -1;  // No type-specific intercepts
    }
    // Type intercepts come after base model parameters
    // Each model has (ntyp - 1) type-specific intercepts
    // param_model_count[model_idx] includes these, so we need to find the offset
    // Base parameters for model = param_model_count[model_idx] - (ntyp - 1)
    int base_params = param_model_count[model_idx] - (ntyp - 1);
    return base_params + ityp;
}

std::vector<double> FactorModel::ComputeTypeProbabilities(const std::vector<double>& fac)
{
    // Compute type probabilities using multinomial logit
    // Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
    //
    // P(type=1) = 1 / (1 + sum_{t=2}^T exp(sum_k lambda_t_k * f_k))
    // P(type=t) = exp(sum_k lambda_t_k * f_k) / (1 + sum_{t=2}^T exp(sum_k lambda_t_k * f_k))

    std::vector<double> probs(ntyp);

    if (ntyp == 1) {
        probs[0] = 1.0;
        return probs;
    }

    // Compute log-odds for types 2, 3, ... relative to type 1
    std::vector<double> log_odds(ntyp - 1);
    double max_log_odds = 0.0;  // For numerical stability

    for (int t = 0; t < ntyp - 1; t++) {
        log_odds[t] = 0.0;
        for (int k = 0; k < nfac; k++) {
            int param_idx = GetTypeLoadingIndex(t, k);
            log_odds[t] += param[param_idx] * fac[k];
        }
        if (log_odds[t] > max_log_odds) {
            max_log_odds = log_odds[t];
        }
    }

    // Compute denominator with log-sum-exp trick for numerical stability
    double sum_exp = std::exp(-max_log_odds);  // This is exp(0 - max) for type 1
    for (int t = 0; t < ntyp - 1; t++) {
        sum_exp += std::exp(log_odds[t] - max_log_odds);
    }

    // Compute probabilities
    probs[0] = std::exp(-max_log_odds) / sum_exp;  // Type 1 (reference)
    for (int t = 0; t < ntyp - 1; t++) {
        probs[t + 1] = std::exp(log_odds[t] - max_log_odds) / sum_exp;
    }

    return probs;
}
