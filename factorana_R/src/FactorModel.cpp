#include "FactorModel.h"
#include "distributions.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <iostream>

FactorModel::FactorModel(int n_obs, int n_var, int n_fac, int n_typ,
                         int n_mix, bool correlated, int n_quad)
    : nobs(n_obs), nvar(n_var), nfac(n_fac), ntyp(n_typ),
      nmix(n_mix), fac_corr(correlated), nquad_points(n_quad)
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

void FactorModel::SetParameterConstraints(const std::vector<bool>& fixed)
{
    if (fixed.size() != nparam) {
        throw std::runtime_error("Parameter constraint size mismatch");
    }
    param_fixed = fixed;
    nparam_free = 0;
    for (bool f : fixed) {
        if (!f) nparam_free++;
    }

    // Initialize parameter vector with default values
    param.resize(nparam, 0.0);

    // Initialize factor variances to 1.0 by default
    for (int ifac = 0; ifac < nfac; ifac++) {
        param[ifac] = 1.0;
    }
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

    // ===== STEP 4: Compute total number of integration points =====
    int nint_points = 1;
    for (int i = 0; i < nfac; i++) {
        nint_points *= nquad_points;
    }

    // Temporary storage
    std::vector<double> modEval;
    std::vector<double> modHess;
    std::vector<double> totalgrad(nparam, 0.0);
    std::vector<double> totalhess;
    if (iflag == 3) totalhess.resize(nparam * nparam, 0.0);

    std::vector<double> gradilk(nparam, 0.0);  // Gradient at integration point
    std::vector<double> hessilk;
    if (iflag == 3) hessilk.resize(nparam * nparam, 0.0);

    // ===== STEP 5: Loop over observations =====
    for (int iobs = 0; iobs < nobs; iobs++) {
        int iobs_offset = iobs * nvar;

        double totalprob = 0.0;
        if (iflag >= 2) std::fill(totalgrad.begin(), totalgrad.end(), 0.0);
        if (iflag == 3) std::fill(totalhess.begin(), totalhess.end(), 0.0);

        // ===== STEP 6: Multi-dimensional integration over factors =====
        // Use "odometer" method: facint tracks indices for each dimension
        std::vector<int> facint(nfac, 0);

        for (int intpt = 0; intpt < nint_points; intpt++) {
            // Reset gradient/Hessian for this integration point
            if (iflag >= 2) std::fill(gradilk.begin(), gradilk.end(), 0.0);
            if (iflag == 3) std::fill(hessilk.begin(), hessilk.end(), 0.0);

            // Get integration weight (product over dimensions)
            double probilk = 1.0;
            for (int ifac = 0; ifac < nfac; ifac++) {
                probilk *= quad_weights[facint[ifac]];
            }

            // Compute factor values at this integration point
            std::vector<double> fac_val(nfac);

            if (fac_corr && nfac == 2) {
                // Correlated 2-factor model: use Cholesky transformation
                // f1 = σ1 * z1
                // f2 = ρ*σ2*z1 + σ2*√(1-ρ²)*z2
                // where z1, z2 are independent GH quadrature nodes
                double sigma1 = std::sqrt(factor_var[0]);
                double sigma2 = std::sqrt(factor_var[1]);
                double z1 = quad_nodes[facint[0]];
                double z2 = quad_nodes[facint[1]];

                fac_val[0] = sigma1 * z1 + factor_mean[0];
                fac_val[1] = rho * sigma2 * z1 + sigma2 * sqrt_1_minus_rho2 * z2 + factor_mean[1];
            } else {
                // Independent factors: standard transformation
                for (int ifac = 0; ifac < nfac; ifac++) {
                    double sigma = std::sqrt(factor_var[ifac]);
                    double x_node = quad_nodes[facint[ifac]];
                    // Factor transformation for GH quadrature
                    // Empirically, f = sigma * x gives correct parameter recovery
                    // (without sqrt(2) factor that theory suggests)
                    fac_val[ifac] = sigma * x_node + factor_mean[ifac];
                }
            }

            // ===== STEP 7: Evaluate all models =====
            for (size_t imod = 0; imod < models.size(); imod++) {
                int firstpar = param_model_start[imod];

                // If all parameters are fixed for this model, only compute likelihood (flag=1)
                // This saves computation time in multi-stage estimation
                int model_flag = models[imod]->GetAllParamsFixed() ? 1 : iflag;

                models[imod]->Eval(iobs_offset, data, param, firstpar, fac_val,
                                  modEval, modHess, model_flag);

                // Multiply likelihood
                probilk *= modEval[0];

                // ===== STEP 8: Accumulate gradients =====
                if (iflag >= 2) {
                    // Gradients w.r.t. factor variances and correlation
                    if (fac_corr && nfac == 2) {
                        // Correlated 2-factor model: use Cholesky derivatives
                        // f1 = σ1 * z1
                        // f2 = ρ*σ2*z1 + σ2*√(1-ρ²)*z2
                        double sigma1 = std::sqrt(factor_var[0]);
                        double sigma2 = std::sqrt(factor_var[1]);
                        double z1 = quad_nodes[facint[0]];
                        double z2 = quad_nodes[facint[1]];

                        // dL/df1 and dL/df2 from models
                        double dL_df1 = modEval[1];
                        double dL_df2 = modEval[2];

                        // Gradient w.r.t. σ1² (factor variance 1)
                        // df1/dσ1² = z1/(2σ1), df2/dσ1² = 0
                        gradilk[0] += dL_df1 * z1 / (2.0 * sigma1);

                        // Gradient w.r.t. σ2² (factor variance 2)
                        // df2/dσ2² = (ρ*z1 + √(1-ρ²)*z2) / (2σ2)
                        double df2_dsigma2sq = (rho * z1 + sqrt_1_minus_rho2 * z2) / (2.0 * sigma2);
                        gradilk[1] += dL_df2 * df2_dsigma2sq;

                        // Gradient w.r.t. ρ (correlation parameter at index 2)
                        // df2/dρ = σ2*z1 + σ2 * d(√(1-ρ²))/dρ * z2
                        //        = σ2*z1 - σ2*ρ*z2/√(1-ρ²)
                        //        = σ2 * (z1 - ρ*z2/sqrt_1_minus_rho2)
                        double df2_drho = sigma2 * (z1 - rho * z2 / sqrt_1_minus_rho2);
                        gradilk[2] += dL_df2 * df2_drho;
                    } else {
                        // Independent factors: standard chain rule
                        for (int ifac = 0; ifac < nfac; ifac++) {
                            // Chain rule: d/d(sigma^2) = d/df * df/d(sigma^2)
                            // where f = sigma * x + mean
                            // df/d(sigma^2) = df/dσ · dσ/d(sigma^2) = x · 1/(2σ) = x/(2σ)
                            double sigma = std::sqrt(factor_var[ifac]);
                            double x_node = quad_nodes[facint[ifac]];
                            gradilk[ifac] += modEval[ifac + 1] * x_node / (2.0 * sigma);
                        }
                    }

                    // Gradients w.r.t. model parameters
                    for (int iparam = 0; iparam < param_model_count[imod]; iparam++) {
                        gradilk[firstpar + iparam] += modEval[1 + nfac + iparam];
                    }
                }

                // ===== STEP 9: Accumulate Hessians =====
                if (iflag == 3 && modHess.size() > 0) {
                    int nDimModHess = nfac + param_model_count[imod];

                    // Model Hessian contribution (with chain rule for factor variances)
                    for (int i = 0; i < nDimModHess; i++) {
                        for (int j = i; j < nDimModHess; j++) {
                            int modhess_idx = i * nDimModHess + j;
                            int full_i = (i < nfac) ? i : (firstpar + i - nfac);
                            int full_j = (j < nfac) ? j : (firstpar + j - nfac);
                            int full_idx = full_i * nparam + full_j;

                            // Chain rule: d²L/dσᵢ²dσⱼ² = d²L/dfᵢdfⱼ * (dfᵢ/dσᵢ²) * (dfⱼ/dσⱼ²)
                            // where f = sigma * x + mean
                            // df/d(sigma^2) = x/(2σ)
                            double chain_factor = 1.0;

                            if (i < nfac) {
                                // i is a factor variance parameter
                                double sigma_i = std::sqrt(factor_var[i]);
                                double x_node_i = quad_nodes[facint[i]];
                                chain_factor *= x_node_i / (2.0 * sigma_i);
                            }

                            if (j < nfac) {
                                // j is a factor variance parameter
                                double sigma_j = std::sqrt(factor_var[j]);
                                double x_node_j = quad_nodes[facint[j]];
                                chain_factor *= x_node_j / (2.0 * sigma_j);
                            }

                            hessilk[full_idx] += modHess[modhess_idx] * chain_factor;

                            // For diagonal factor variance elements (i == j), add second derivative term
                            // Full chain rule: d²L/d(σ²)² = d²L/df² * (df/d(σ²))² + dL/df * d²f/d(σ²)²
                            // The first term is handled above. The second term is:
                            //   dL/df * d²f/d(σ²)² where d²f/d(σ²)² = -x/(4σ³)
                            // Note: f = σ*x = √(σ²)*x, so:
                            //   df/d(σ²) = x/(2σ)
                            //   d²f/d(σ²)² = d[x/(2σ)]/d(σ²) = -x/(4σ³)
                            if (i < nfac && i == j) {
                                double sigma = std::sqrt(factor_var[i]);
                                double x_node = quad_nodes[facint[i]];
                                double second_deriv_factor = -x_node / (4.0 * sigma * sigma * sigma);
                                // modEval[i + 1] is dL/df_i from this model
                                hessilk[full_idx] += modEval[i + 1] * second_deriv_factor;
                            }
                        }
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
                for (int i = 0; i < nparam; i++) {
                    for (int j = i; j < nparam; j++) {
                        int idx = i * nparam + j;
                        // Chain rule: d²L/dθ² = prob * (d²f/dθ² + (df/dθ)²)
                        totalhess[idx] += probilk * (hessilk[idx] + gradilk[i] * gradilk[j]);
                    }
                }
            }

            // ===== STEP 11: Increment integration point indices (odometer style) =====
            if (intpt < nint_points - 1) {
                for (int ifac = 0; ifac < nfac; ifac++) {
                    facint[ifac]++;
                    if (facint[ifac] < nquad_points) break;
                    facint[ifac] = 0;
                }
            }
        }

        // ===== STEP 12: Compute log-likelihood for this observation =====
        if (totalprob > 1e-100) {
            logLkhd += std::log(totalprob);

            // Gradient: d(log L)/dθ = (1/L) * dL/dθ
            if (iflag >= 2) {
                for (int ipar = 0; ipar < nparam; ipar++) {
                    full_gradL[ipar] += totalgrad[ipar] / totalprob;
                }
            }

            // Hessian: d²(log L)/dθ² = (1/L)*d²L/dθ² - (1/L²)*(dL/dθ)²
            if (iflag == 3) {
                for (int i = 0; i < nparam; i++) {
                    for (int j = i; j < nparam; j++) {
                        int idx = i * nparam + j;
                        double term1 = totalhess[idx] / totalprob;
                        double term2 = (totalgrad[i] * totalgrad[j]) / (totalprob * totalprob);
                        full_hessL[idx] += term1 - term2;
                    }
                }
            }
        } else {
            // Numerical underflow
            logLkhd += -1e10;
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
