#include "Model.h"
#include "distributions.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

Model::Model(ModelType type, int outcome, int missing,
             const std::vector<int>& regs, int nfac, int ntyp,
             const std::vector<double>& fnorm,
             int nchoice, int nrank, bool params_fixed,
             FactorSpec fspec)
    : modtype(type), outcome_idx(outcome), missing_idx(missing),
      regressors(regs), numfac(nfac), numtyp(ntyp),
      facnorm(fnorm), numchoice(nchoice), numrank(nrank),
      ignore(false), all_params_fixed(params_fixed),
      factor_spec(fspec)
{
    nregressors = regressors.size();

    // Compute number of quadratic and interaction loadings based on factor_spec
    if (factor_spec == FactorSpec::QUADRATIC || factor_spec == FactorSpec::FULL) {
        n_quadratic_loadings = numfac;
    } else {
        n_quadratic_loadings = 0;
    }

    if ((factor_spec == FactorSpec::INTERACTIONS || factor_spec == FactorSpec::FULL) && numfac >= 2) {
        n_interaction_loadings = numfac * (numfac - 1) / 2;
    } else {
        n_interaction_loadings = 0;
    }
}

void Model::Eval(int iobs_offset, const std::vector<double>& data,
                 const std::vector<double>& param, int firstpar,
                 const std::vector<double>& fac,
                 std::vector<double>& modEval,
                 std::vector<double>& hess,
                 int flag,
                 double type_intercept)
{
    // Always clear vectors first to ensure no stale data persists
    // This is critical for consistent behavior across different flag values
    modEval.clear();
    hess.clear();

    // Count free factor loadings FIRST (needed for gradient vector sizing)
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    // Determine size of gradient vector
    if (flag >= 2) {
        // Layout: 1 (likelihood) + numfac (d/dtheta for factor variances) +
        //         nregressors (d/dbeta) + ifreefac (d/dalpha for free linear loadings) +
        //         n_quadratic_loadings (d/dalpha_quad) + n_interaction_loadings (d/dalpha_inter) +
        //         model-specific parameters
        int ngrad = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

        // Model-specific parameters
        if (modtype == ModelType::LINEAR) ngrad += 1;  // sigma
        if (modtype == ModelType::LOGIT && numchoice > 2) {
            ngrad = 1 + numfac + (numchoice-1)*(nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings);
        }
        if (modtype == ModelType::OPROBIT) ngrad += (numchoice - 1);  // thresholds

        modEval.resize(ngrad, 0.0);  // Initialize all to 0.0
    } else {
        modEval.resize(1);
    }

    modEval[0] = 1.0;

    // Check missing indicator
    if (missing_idx > -1) {
        if (data[iobs_offset + missing_idx] == 0) return;
    }
    if (ignore) return;

    // Build linear predictor(s)
    int numlogitchoice = 2;
    std::vector<double> expres(1, 0.0);

    if (modtype == ModelType::LOGIT && numchoice > 2) {
        expres.resize(numchoice - 1, 0.0);
        numlogitchoice = numchoice;
    }

    for (int ichoice = 0; ichoice < numlogitchoice - 1; ichoice++) {
        int ifree_local = 0;
        int nparamchoice = nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

        // Add regressor terms
        for (int i = 0; i < nregressors; i++) {
            expres[ichoice] += param[i + firstpar + ichoice*nparamchoice] *
                              data[iobs_offset + regressors[i]];
        }

        // Add factor and type loadings
        // IMPORTANT: Loop bound must respect actual vector sizes to avoid out-of-bounds access
        // When facnorm is provided, use its size; otherwise use numfac
        // The numtyp term is legacy and not currently used in the R interface
        int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
        for (int i = 0; i < fac_loop_bound; i++) {
            if (facnorm.size() == 0) {
                // All loadings are free
                expres[ichoice] += param[ifree_local + firstpar + ichoice*nparamchoice + nregressors] * fac[i];
                ifree_local++;
            } else {
                // Check if this loading is normalized
                if (facnorm[i] > -9998) {
                    // Fixed loading
                    expres[ichoice] += facnorm[i] * fac[i];
                } else {
                    // Free loading
                    expres[ichoice] += param[ifree_local + firstpar + ichoice*nparamchoice + nregressors] * fac[i];
                    ifree_local++;
                }
            }
        }

        // Add quadratic factor terms: sum(lambda_quad_k * f_k^2)
        if (n_quadratic_loadings > 0) {
            int quad_start = nregressors + ifreefac;
            for (int k = 0; k < numfac; k++) {
                double fk_sq = fac[k] * fac[k];
                expres[ichoice] += param[firstpar + ichoice*nparamchoice + quad_start + k] * fk_sq;
            }
        }

        // Add interaction factor terms: sum(lambda_inter_jk * f_j * f_k) for j < k
        if (n_interaction_loadings > 0) {
            int inter_start = nregressors + ifreefac + n_quadratic_loadings;
            int inter_idx = 0;
            for (int j = 0; j < numfac - 1; j++) {
                for (int k = j + 1; k < numfac; k++) {
                    double fj_fk = fac[j] * fac[k];
                    expres[ichoice] += param[firstpar + ichoice*nparamchoice + inter_start + inter_idx] * fj_fk;
                    inter_idx++;
                }
            }
        }
    }

    // Add type-specific intercept to linear predictor(s)
    // For models with multiple types, this shifts the index function
    for (int ichoice = 0; ichoice < numlogitchoice - 1; ichoice++) {
        expres[ichoice] += type_intercept;
    }

    // Dispatch to model-specific evaluation
    if (modtype == ModelType::LINEAR) {
        double Z = 0.0;
        if (outcome_idx > -1) Z = data[iobs_offset + outcome_idx] - expres[0];
        double sigma = std::fabs(param[firstpar + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings]);
        EvalLinear(Z, sigma, fac, param, firstpar, modEval, hess, flag, data, iobs_offset);
    }
    else if (modtype == ModelType::PROBIT) {
        double obsSign = 1.0;
        if (int(data[iobs_offset + outcome_idx]) == 0) obsSign = -1.0;
        EvalProbit(expres[0], obsSign, fac, param, firstpar, modEval, hess, flag, data, iobs_offset);
    }
    else if (modtype == ModelType::LOGIT) {
        double outcome_val = data[iobs_offset + outcome_idx];
        EvalLogit(expres, outcome_val, fac, param, firstpar, modEval, hess, flag, data, iobs_offset);
    }
    else if (modtype == ModelType::OPROBIT) {
        int outcome_val = int(data[iobs_offset + outcome_idx]);
        EvalOprobit(expres[0], outcome_val, fac, param, firstpar, modEval, hess, flag, data, iobs_offset);
    }
}

void Model::EvalLinear(double Z, double sigma, const std::vector<double>& fac,
                       const std::vector<double>& param, int firstpar,
                       std::vector<double>& modEval, std::vector<double>& hess,
                       int flag, const std::vector<double>& data, int iobs_offset)
{
    // Compute likelihood
    modEval[0] = normal_pdf(Z, 0.0, sigma);

    if (flag < 2) return;

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    int npar = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings + 1; // +1 for sigma

    // Initialize Hessian if needed
    if (flag == 3) {
        hess.resize(npar * npar);
        for (int i = 0; i < npar; i++) {
            for (int j = i; j < npar; j++) {
                hess[i*npar + j] = -1.0 / (sigma*sigma);
            }
        }
    }

    // Compute gradients
    double sigma2 = sigma * sigma;
    ifreefac = 0;

    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
    for (int i = 0; i < fac_loop_bound; i++) {
        if (facnorm.size() == 0) {
            // All free
            if (i < numfac) {
                // Gradient w.r.t. factor variance
                modEval[i+1] = Z * param[ifreefac + firstpar + nregressors] / sigma2;
            }
            // Gradient w.r.t. factor loading
            modEval[1 + numfac + nregressors + ifreefac] = Z * fac[i] / sigma2;

            if (flag == 3) {
                if (i < numfac) {
                    for (int j = i; j < npar; j++)
                        hess[i*npar + j] *= param[ifreefac + firstpar + nregressors];
                    for (int j = 0; j <= i; j++)
                        hess[j*npar + i] *= param[ifreefac + firstpar + nregressors];
                }
                int index = numfac + nregressors + ifreefac;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= fac[i];
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= fac[i];
                hess[i*npar + index] += Z / sigma2;
            }
            ifreefac++;
        } else {
            // Mixed fixed/free
            if (facnorm[i] > -9998.0) {
                // Fixed loading
                if (i < numfac) {
                    modEval[i+1] = Z * facnorm[i] / sigma2;
                    if (flag == 3) {
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= facnorm[i];
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= facnorm[i];
                    }
                }
            } else {
                // Free loading
                if (i < numfac) {
                    modEval[i+1] = Z * param[ifreefac + firstpar + nregressors] / sigma2;
                }
                modEval[1 + numfac + nregressors + ifreefac] = Z * fac[i] / sigma2;

                if (flag == 3) {
                    if (i < numfac) {
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= param[ifreefac + firstpar + nregressors];
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= param[ifreefac + firstpar + nregressors];
                    }
                    int index = numfac + nregressors + ifreefac;
                    for (int j = index; j < npar; j++)
                        hess[index*npar + j] *= fac[i];
                    for (int j = 0; j <= index; j++)
                        hess[j*npar + index] *= fac[i];
                    hess[i*npar + index] += Z / sigma2;
                }
                ifreefac++;
            }
        }
    }

    // Gradients w.r.t. regression coefficients
    for (int ireg = 0; ireg < nregressors; ireg++) {
        modEval[ireg + numfac + 1] = Z * data[iobs_offset + regressors[ireg]] / sigma2;
        if (flag == 3) {
            int index = numfac + ireg;
            for (int j = index; j < npar; j++)
                hess[index*npar + j] *= data[iobs_offset + regressors[ireg]];
            for (int j = 0; j <= index; j++)
                hess[j*npar + index] *= data[iobs_offset + regressors[ireg]];
        }
    }

    // Gradients w.r.t. quadratic factor loadings (lambda_quad_k)
    // d(logL)/d(lambda_quad_k) = Z * f_k^2 / sigma^2
    if (n_quadratic_loadings > 0) {
        int quad_grad_start = 1 + numfac + nregressors + ifreefac;
        for (int k = 0; k < numfac; k++) {
            double fk_sq = fac[k] * fac[k];
            modEval[quad_grad_start + k] = Z * fk_sq / sigma2;

            if (flag == 3) {
                int index = numfac + nregressors + ifreefac + k;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= fk_sq;
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= fk_sq;

                // Cross-derivative d^2(logL)/(d(theta_k) d(lambda_quad_k)) = 2 * f_k * lambda_quad_k / sigma^2
                // The factor variance gradient uses chain rule with lambda_quad * 2 * f
                // Adding to existing Hessian structure (d(theta_k)/d(lambda_quad_k) cross term)
                double lambda_quad_k = param[firstpar + nregressors + ifreefac + k];
                hess[k*npar + index] += 2.0 * fac[k] * Z / sigma2;
            }
        }
    }

    // Gradients w.r.t. interaction factor loadings (lambda_inter_jk) for j < k
    // d(logL)/d(lambda_inter_jk) = Z * f_j * f_k / sigma^2
    if (n_interaction_loadings > 0) {
        int inter_grad_start = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings;
        int inter_idx = 0;
        for (int j = 0; j < numfac - 1; j++) {
            for (int k = j + 1; k < numfac; k++) {
                double fj_fk = fac[j] * fac[k];
                modEval[inter_grad_start + inter_idx] = Z * fj_fk / sigma2;

                if (flag == 3) {
                    int index = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                    for (int jj = index; jj < npar; jj++)
                        hess[index*npar + jj] *= fj_fk;
                    for (int jj = 0; jj <= index; jj++)
                        hess[jj*npar + index] *= fj_fk;

                    // Cross-derivatives d^2(logL)/(d(theta) d(lambda_inter))
                    double lambda_inter = param[firstpar + nregressors + ifreefac + n_quadratic_loadings + inter_idx];
                    hess[j*npar + index] += fac[k] * Z / sigma2;  // d(theta_j)/d(lambda_inter_jk)
                    hess[k*npar + index] += fac[j] * Z / sigma2;  // d(theta_k)/d(lambda_inter_jk)
                }
                inter_idx++;
            }
        }
    }

    // Gradient w.r.t. sigma
    int sigma_grad_idx = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;
    modEval[sigma_grad_idx] = (Z*Z / sigma - sigma) / sigma2;

    // Hessian for sigma
    if (flag == 3) {
        int sigma_index = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

        // Multiply sigma row and column by 2*Z/sigma for cross-derivatives
        for (int j = sigma_index; j < npar; j++) {
            hess[sigma_index * npar + j] *= 2.0 * Z / sigma;
        }
        for (int j = 0; j < sigma_index; j++) {
            hess[j * npar + sigma_index] *= 2.0 * Z / sigma;
        }

        // Separately handle diagonal element with correct formula
        hess[sigma_index * npar + sigma_index] = -3.0 * Z * Z / (sigma*sigma*sigma*sigma) + 1.0 / (sigma*sigma);
    }
}

void Model::EvalProbit(double expres, double obsSign, const std::vector<double>& fac,
                       const std::vector<double>& param, int firstpar,
                       std::vector<double>& modEval, std::vector<double>& hess,
                       int flag, const std::vector<double>& data, int iobs_offset)
{
    // Compute likelihood
    modEval[0] = normal_cdf(obsSign * expres);

    if (flag < 2) return;

    double Z = obsSign * expres;
    double pdf = normal_pdf(obsSign * expres);
    double cdf = modEval[0];
    if (obsSign * expres < -35.0) cdf = 1.0e-50;

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    int npar = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

    if (flag == 3) {
        hess.resize(npar * npar, 0.0);
        for (int i = 0; i < npar; i++) {
            for (int j = i; j < npar; j++) {
                hess[i*npar + j] = 1.0;
            }
        }
    }

    // Common term for probit gradients: d(logL)/d(xb) = phi/Phi * obsSign
    double dlogL_dxb = pdf * obsSign / cdf;

    // Compute gradients
    ifreefac = 0;
    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
    for (int i = 0; i < fac_loop_bound; i++) {
        if (facnorm.size() == 0) {
            if (i < numfac) {
                modEval[i+1] = pdf * obsSign * param[ifreefac + firstpar + nregressors] / cdf;
            }
            modEval[1 + numfac + nregressors + ifreefac] = pdf * (obsSign * fac[i]) / cdf;

            if (flag == 3) {
                double lambda_val = pdf / cdf;  // λ = φ/Φ
                if (i < numfac) {
                    double dZ_dtheta = obsSign * param[ifreefac + firstpar + nregressors];
                    double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                    for (int j = i; j < npar; j++)
                        hess[i*npar + j] *= row_mult;
                    for (int j = 0; j <= i; j++)
                        hess[j*npar + i] *= dZ_dtheta;
                }
                int index = numfac + nregressors + ifreefac;
                double dZ_dtheta = obsSign * fac[i];
                double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= row_mult;
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= dZ_dtheta;
            }
            ifreefac++;
        } else {
            if (facnorm[i] > -9998.0) {
                if (i < numfac) {
                    modEval[i+1] = pdf * obsSign * facnorm[i] / cdf;
                    if (flag == 3) {
                        double lambda_val = pdf / cdf;  // λ = φ/Φ
                        double dZ_dtheta = obsSign * facnorm[i];
                        double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= row_mult;
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= dZ_dtheta;
                    }
                }
            } else {
                if (i < numfac) {
                    modEval[i+1] = pdf * obsSign * param[ifreefac + firstpar + nregressors] / cdf;
                }
                modEval[1 + numfac + nregressors + ifreefac] = pdf * (obsSign * fac[i]) / cdf;

                if (flag == 3) {
                    double lambda_val = pdf / cdf;  // λ = φ/Φ
                    if (i < numfac) {
                        double dZ_dtheta = obsSign * param[ifreefac + firstpar + nregressors];
                        double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                        for (int j = i; j < npar; j++)
                            hess[i*npar + j] *= row_mult;
                        for (int j = 0; j <= i; j++)
                            hess[j*npar + i] *= dZ_dtheta;
                    }
                    int index = numfac + nregressors + ifreefac;
                    double dZ_dtheta_loading = obsSign * fac[i];
                    double row_mult_loading = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta_loading;
                    for (int j = index; j < npar; j++)
                        hess[index*npar + j] *= row_mult_loading;
                    for (int j = 0; j <= index; j++)
                        hess[j*npar + index] *= dZ_dtheta_loading;
                }
                ifreefac++;
            }
        }
    }

    // Gradients w.r.t. regression coefficients
    for (int ireg = 0; ireg < nregressors; ireg++) {
        modEval[ireg + numfac + 1] = pdf * (obsSign * data[iobs_offset + regressors[ireg]]) / cdf;
        if (flag == 3) {
            double lambda_val = pdf / cdf;  // λ = φ/Φ
            int index = numfac + ireg;
            double dZ_dtheta = obsSign * data[iobs_offset + regressors[ireg]];
            double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
            for (int j = index; j < npar; j++)
                hess[index*npar + j] *= row_mult;
            for (int j = 0; j <= index; j++)
                hess[j*npar + index] *= dZ_dtheta;
        }
    }

    // Gradients w.r.t. quadratic factor loadings (lambda_quad_k)
    // d(logL)/d(lambda_quad_k) = dlogL_dxb * f_k^2
    if (n_quadratic_loadings > 0) {
        int quad_grad_start = 1 + numfac + nregressors + ifreefac;
        for (int k = 0; k < numfac; k++) {
            double fk_sq = fac[k] * fac[k];
            modEval[quad_grad_start + k] = dlogL_dxb * fk_sq;

            if (flag == 3) {
                double lambda_val = pdf / cdf;
                int index = numfac + nregressors + ifreefac + k;
                double dZ_dtheta = obsSign * fk_sq;
                double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                for (int j = index; j < npar; j++)
                    hess[index*npar + j] *= row_mult;
                for (int j = 0; j <= index; j++)
                    hess[j*npar + index] *= dZ_dtheta;

                // Cross-derivative d(theta_k)/d(lambda_quad_k) = 2 * f_k * lambda_quad_k
                hess[k*npar + index] += 2.0 * fac[k] * lambda_val * obsSign;
            }
        }
    }

    // Gradients w.r.t. interaction factor loadings (lambda_inter_jk) for j < k
    // d(logL)/d(lambda_inter_jk) = dlogL_dxb * f_j * f_k
    if (n_interaction_loadings > 0) {
        int inter_grad_start = 1 + numfac + nregressors + ifreefac + n_quadratic_loadings;
        int inter_idx = 0;
        for (int j = 0; j < numfac - 1; j++) {
            for (int k = j + 1; k < numfac; k++) {
                double fj_fk = fac[j] * fac[k];
                modEval[inter_grad_start + inter_idx] = dlogL_dxb * fj_fk;

                if (flag == 3) {
                    double lambda_val = pdf / cdf;
                    int index = numfac + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                    double dZ_dtheta = obsSign * fj_fk;
                    double row_mult = (-Z * lambda_val - lambda_val * lambda_val) * dZ_dtheta;
                    for (int jj = index; jj < npar; jj++)
                        hess[index*npar + jj] *= row_mult;
                    for (int jj = 0; jj <= index; jj++)
                        hess[jj*npar + index] *= dZ_dtheta;

                    // Cross-derivatives d(theta)/d(lambda_inter)
                    hess[j*npar + index] += fac[k] * lambda_val * obsSign;
                    hess[k*npar + index] += fac[j] * lambda_val * obsSign;
                }
                inter_idx++;
            }
        }
    }

    // Add cross-derivative terms for Hessian
    // d²(log L)/df_k dλ_k = (-Z*λ - λ²) * λ_k * f_k + λ * s
    // where λ = φ/Φ (Mills ratio), s = obsSign
    // The first term comes from the multiplicative Hessian structure,
    // but we need to add the second term λ * s = (pdf/cdf) * obsSign
    if (flag == 3) {
        double lambda_val = pdf / cdf;  // Mills ratio
        ifreefac = 0;
        for (int i = 0; i < numfac; i++) {
            if (facnorm.size() == 0 || facnorm[i] <= -9998) {
                int index = numfac + nregressors + ifreefac;
                hess[i*npar + index] += lambda_val * obsSign;
                ifreefac++;
            }
        }
    }
}

void Model::EvalLogit(const std::vector<double>& expres, double outcome,
                      const std::vector<double>& fac,
                      const std::vector<double>& param, int firstpar,
                      std::vector<double>& modEval, std::vector<double>& hess,
                      int flag, const std::vector<double>& data, int iobs_offset)
{
    // Multinomial logit with K choices (numchoice)
    // Choice 0 is reference category with Z_0 = 0
    // For choices 1 to K-1, we have separate parameters

    // Observed choice (0-indexed: 0, 1, ..., numchoice-1)
    int obsCat = int(outcome) - 1;  // Convert from 1-indexed to 0-indexed

    if (obsCat < 0 || obsCat >= numchoice) {
        std::cerr << "ERROR: Invalid multinomial choice " << int(outcome)
                  << " (must be 1 to " << numchoice << ")" << std::endl;
        modEval[0] = 1e-100;
        return;
    }

    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    // Number of parameters per choice
    int nparamchoice = nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

    // Compute logit denominator: 1 + sum_{k=1}^{K-1} exp(Z_k)
    double logitdenom = 1.0;
    for (int icat = 1; icat < numchoice; icat++) {
        logitdenom += std::exp(expres[icat - 1]);
    }

    // Likelihood: P(Y = obsCat)
    double dens = 1.0 / logitdenom;
    if (obsCat > 0) {
        dens *= std::exp(expres[obsCat - 1]);
    }

    modEval[0] = dens;

    // Numerical stability
    if (modEval[0] < 1e-100) modEval[0] = 1e-100;

    if (flag < 2) return;

    // ===== GRADIENT CALCULATION =====

    // Softmax probabilities for all choices
    std::vector<double> pdf(numchoice);
    pdf[0] = 1.0 / logitdenom;
    for (int icat = 1; icat < numchoice; icat++) {
        pdf[icat] = std::exp(expres[icat - 1]) / logitdenom;
    }

    int npar = numfac + (numchoice - 1) * nparamchoice;

    // Intermediate gradient storage for Hessian (if needed)
    std::vector<double> logitgrad;
    if (flag == 3) {
        logitgrad.resize(numchoice * npar, 0.0);
    }

    // Gradient for factor variance parameters (theta)
    for (int ifac = 0; ifac < numfac; ifac++) {
        // Term from observed choice
        if (obsCat > 0) {
            int param_idx = firstpar + (obsCat - 1) * nparamchoice + nregressors;
            if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                int loading_idx = param_idx + (facnorm.size() > 0 ?
                    std::count_if(facnorm.begin(), facnorm.begin() + ifac,
                                  [](double x) { return x <= -9998; }) : ifac);
                modEval[1 + ifac] += param[loading_idx];

                if (flag == 3) {
                    for (int jcat = 1; jcat < numchoice; jcat++) {
                        int jparam_idx = firstpar + (jcat - 1) * nparamchoice + nregressors;
                        int jloading_idx = jparam_idx + (facnorm.size() > 0 ?
                            std::count_if(facnorm.begin(), facnorm.begin() + ifac,
                                          [](double x) { return x <= -9998; }) : ifac);
                        logitgrad[jcat * npar + ifac] += param[jloading_idx];
                    }
                }
            }
        }

        // Sum over all non-reference choices
        for (int icat = 1; icat < numchoice; icat++) {
            int param_idx = firstpar + (icat - 1) * nparamchoice + nregressors;
            if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                int loading_idx = param_idx + (facnorm.size() > 0 ?
                    std::count_if(facnorm.begin(), facnorm.begin() + ifac,
                                  [](double x) { return x <= -9998; }) : ifac);
                modEval[1 + ifac] += -pdf[icat] * param[loading_idx];

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + ifac] += -pdf[icat] * param[loading_idx];
                    }
                }
            }
        }
    }

    // Gradient for factor loadings and regression coefficients
    // Following legacy TModel.cc pattern exactly

    // Factor loadings - obsCat term and logitgrad
    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;
    int ifree = 0;
    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
            double fval = fac[ifac];

            // obsCat term for gradient
            if (obsCat > 0) {
                int base_idx = numfac + (obsCat - 1) * nparamchoice;
                modEval[1 + base_idx + nregressors + ifree] += fval;
            }

            // logitgrad update (unconditional, not nested in obsCat check!)
            if (flag == 3) {
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                    logitgrad[jcat * npar + jbase_idx + nregressors + ifree] += fval;
                }
            }

            // All categories term
            for (int icat = 1; icat < numchoice; icat++) {
                int base_idx = numfac + (icat - 1) * nparamchoice;
                modEval[1 + base_idx + nregressors + ifree] += -pdf[icat] * fval;

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + base_idx + nregressors + ifree] += -pdf[icat] * fval;
                    }
                }
            }

            ifree++;
        }
    }

    // Regression coefficients - obsCat term and logitgrad
    for (int ireg = 0; ireg < nregressors; ireg++) {
        double xval = data[iobs_offset + regressors[ireg]];

        // obsCat term for gradient
        if (obsCat > 0) {
            int base_idx = numfac + (obsCat - 1) * nparamchoice;
            modEval[1 + base_idx + ireg] += xval;
        }

        // logitgrad update (unconditional, not nested in obsCat check!)
        if (flag == 3) {
            for (int jcat = 1; jcat < numchoice; jcat++) {
                int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                logitgrad[jcat * npar + jbase_idx + ireg] += xval;
            }
        }

        // All categories term
        for (int icat = 1; icat < numchoice; icat++) {
            int base_idx = numfac + (icat - 1) * nparamchoice;
            modEval[1 + base_idx + ireg] += -pdf[icat] * xval;

            if (flag == 3) {
                for (int jcat = 0; jcat < numchoice; jcat++) {
                    logitgrad[jcat * npar + base_idx + ireg] += -pdf[icat] * xval;
                }
            }
        }
    }

    // Quadratic factor loadings - obsCat term and logitgrad
    if (n_quadratic_loadings > 0) {
        for (int k = 0; k < numfac; k++) {
            double fk_sq = fac[k] * fac[k];
            int quad_offset = nregressors + ifreefac + k;

            // obsCat term for gradient
            if (obsCat > 0) {
                int base_idx = numfac + (obsCat - 1) * nparamchoice;
                modEval[1 + base_idx + quad_offset] += fk_sq;
            }

            // logitgrad update
            if (flag == 3) {
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                    logitgrad[jcat * npar + jbase_idx + quad_offset] += fk_sq;
                }
            }

            // All categories term
            for (int icat = 1; icat < numchoice; icat++) {
                int base_idx = numfac + (icat - 1) * nparamchoice;
                modEval[1 + base_idx + quad_offset] += -pdf[icat] * fk_sq;

                if (flag == 3) {
                    for (int jcat = 0; jcat < numchoice; jcat++) {
                        logitgrad[jcat * npar + base_idx + quad_offset] += -pdf[icat] * fk_sq;
                    }
                }
            }
        }
    }

    // Interaction factor loadings - obsCat term and logitgrad
    if (n_interaction_loadings > 0) {
        int inter_idx = 0;
        for (int j = 0; j < numfac - 1; j++) {
            for (int k = j + 1; k < numfac; k++) {
                double fj_fk = fac[j] * fac[k];
                int inter_offset = nregressors + ifreefac + n_quadratic_loadings + inter_idx;

                // obsCat term for gradient
                if (obsCat > 0) {
                    int base_idx = numfac + (obsCat - 1) * nparamchoice;
                    modEval[1 + base_idx + inter_offset] += fj_fk;
                }

                // logitgrad update
                if (flag == 3) {
                    for (int jcat = 1; jcat < numchoice; jcat++) {
                        int jbase_idx = numfac + (jcat - 1) * nparamchoice;
                        logitgrad[jcat * npar + jbase_idx + inter_offset] += fj_fk;
                    }
                }

                // All categories term
                for (int icat = 1; icat < numchoice; icat++) {
                    int base_idx = numfac + (icat - 1) * nparamchoice;
                    modEval[1 + base_idx + inter_offset] += -pdf[icat] * fj_fk;

                    if (flag == 3) {
                        for (int jcat_inner = 0; jcat_inner < numchoice; jcat_inner++) {
                            logitgrad[jcat_inner * npar + base_idx + inter_offset] += -pdf[icat] * fj_fk;
                        }
                    }
                }
                inter_idx++;
            }
        }
    }

    // ===== HESSIAN CALCULATION =====
    if (flag == 3) {
        hess.resize(npar * npar, 0.0);

        // Second-order derivative terms: dZ/dtheta dalpha
        if (obsCat > 0) {
            int ifree = 0;
            for (int ifac = 0; ifac < numfac; ifac++) {
                if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                    int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifree;
                    hess[ifac * npar + index] += 1.0;
                    ifree++;
                }
            }
        }

        // Second term of second-order derivatives
        for (int icat = 1; icat < numchoice; icat++) {
            int ifree = 0;
            for (int ifac = 0; ifac < numfac; ifac++) {
                if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                    int index = numfac + (icat - 1) * nparamchoice + nregressors + ifree;
                    hess[ifac * npar + index] += -pdf[icat];
                    ifree++;
                }
            }
        }

        // Second-order derivative terms for d(theta)/d(lambda_quad): d^2Z/dtheta_k dlambda_quad_k = 2*f_k
        if (n_quadratic_loadings > 0) {
            if (obsCat > 0) {
                for (int k = 0; k < numfac; k++) {
                    int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifreefac + k;
                    hess[k * npar + index] += 2.0 * fac[k];
                }
            }
            for (int icat = 1; icat < numchoice; icat++) {
                for (int k = 0; k < numfac; k++) {
                    int index = numfac + (icat - 1) * nparamchoice + nregressors + ifreefac + k;
                    hess[k * npar + index] += -pdf[icat] * 2.0 * fac[k];
                }
            }
        }

        // Second-order derivative terms for d(theta)/d(lambda_inter): d^2Z/dtheta_j dlambda_inter_jk = f_k (and similar for k)
        if (n_interaction_loadings > 0) {
            if (obsCat > 0) {
                int inter_idx = 0;
                for (int j = 0; j < numfac - 1; j++) {
                    for (int k = j + 1; k < numfac; k++) {
                        int index = numfac + (obsCat - 1) * nparamchoice + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        hess[j * npar + index] += fac[k];  // d(theta_j)/d(lambda_inter_jk)
                        hess[k * npar + index] += fac[j];  // d(theta_k)/d(lambda_inter_jk)
                        inter_idx++;
                    }
                }
            }
            for (int icat = 1; icat < numchoice; icat++) {
                int inter_idx = 0;
                for (int j = 0; j < numfac - 1; j++) {
                    for (int k = j + 1; k < numfac; k++) {
                        int index = numfac + (icat - 1) * nparamchoice + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                        hess[j * npar + index] += -pdf[icat] * fac[k];
                        hess[k * npar + index] += -pdf[icat] * fac[j];
                        inter_idx++;
                    }
                }
            }
        }

        // First-order derivative Hessian terms
        // dtheta d...
        for (int ifac = 0; ifac < numfac; ifac++) {
            for (int icat = 1; icat < numchoice; icat++) {
                int param_idx = firstpar + (icat - 1) * nparamchoice + nregressors;
                double loading_val = 0.0;

                if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                    int loading_idx = param_idx + (facnorm.size() > 0 ?
                        std::count_if(facnorm.begin(), facnorm.begin() + ifac,
                                      [](double x) { return x <= -9998; }) : ifac);
                    loading_val = param[loading_idx];
                }

                // dtheta dtheta
                for (int jfac = ifac; jfac < numfac; jfac++) {
                    hess[ifac * npar + jfac] += -pdf[icat] * logitgrad[icat * npar + jfac] * loading_val;
                }

                // dtheta dalpha, dtheta dbeta, dtheta dquad, dtheta dinter for each category
                for (int jcat = 1; jcat < numchoice; jcat++) {
                    int jbase_idx = numfac + (jcat - 1) * nparamchoice;

                    // dtheta dbeta
                    for (int jreg = 0; jreg < nregressors; jreg++) {
                        int index = jbase_idx + jreg;
                        hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                    }

                    // dtheta dalpha
                    int jfree = 0;
                    for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                        if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                            int index = jbase_idx + nregressors + jfree;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                            jfree++;
                        }
                    }

                    // dtheta dquad
                    if (n_quadratic_loadings > 0) {
                        for (int k = 0; k < numfac; k++) {
                            int index = jbase_idx + nregressors + ifreefac + k;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                        }
                    }

                    // dtheta dinter
                    if (n_interaction_loadings > 0) {
                        for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                            int index = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                            hess[ifac * npar + index] += -pdf[icat] * logitgrad[icat * npar + index] * loading_val;
                        }
                    }
                }
            }
        }

        // Loop over row categories for remaining Hessian terms
        for (int icat = 1; icat < numchoice; icat++) {
            int ibase_idx = numfac + (icat - 1) * nparamchoice;

            for (int jcat = icat; jcat < numchoice; jcat++) {
                int jbase_idx = numfac + (jcat - 1) * nparamchoice;

                // dbeta dbeta
                for (int ireg = 0; ireg < nregressors; ireg++) {
                    double xval_i = data[iobs_offset + regressors[ireg]];
                    for (int jreg = 0; jreg < nregressors; jreg++) {
                        if ((jcat > icat) || (jreg >= ireg)) {
                            int index1 = ibase_idx + ireg;
                            int index2 = jbase_idx + jreg;
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                        }
                    }
                }

                // dbeta dalpha
                for (int ireg = 0; ireg < nregressors; ireg++) {
                    double xval_i = data[iobs_offset + regressors[ireg]];
                    int jfree = 0;
                    for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                        if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                            int index1 = ibase_idx + ireg;
                            int index2 = jbase_idx + nregressors + jfree;
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                            jfree++;
                        }
                    }
                }

                // dalpha dbeta (only for jcat > icat)
                if (jcat > icat) {
                    int ifree = 0;
                    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                            for (int jreg = 0; jreg < nregressors; jreg++) {
                                int index1 = ibase_idx + nregressors + ifree;
                                int index2 = jbase_idx + jreg;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                            ifree++;
                        }
                    }
                }

                // dalpha dalpha
                int ifree = 0;
                for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                    if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                        int jfree = 0;
                        for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                            if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                                if ((jcat > icat) || (jfac >= ifac)) {
                                    int index1 = ibase_idx + nregressors + ifree;
                                    int index2 = jbase_idx + nregressors + jfree;
                                    hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                                }
                                jfree++;
                            }
                        }
                        ifree++;
                    }
                }

                // dbeta dquad
                if (n_quadratic_loadings > 0) {
                    for (int ireg = 0; ireg < nregressors; ireg++) {
                        double xval_i = data[iobs_offset + regressors[ireg]];
                        for (int k = 0; k < numfac; k++) {
                            int index1 = ibase_idx + ireg;
                            int index2 = jbase_idx + nregressors + ifreefac + k;
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                        }
                    }
                }

                // dbeta dinter
                if (n_interaction_loadings > 0) {
                    for (int ireg = 0; ireg < nregressors; ireg++) {
                        double xval_i = data[iobs_offset + regressors[ireg]];
                        for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                            int index1 = ibase_idx + ireg;
                            int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                            hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * xval_i;
                        }
                    }
                }

                // dalpha dquad
                if (n_quadratic_loadings > 0) {
                    ifree = 0;
                    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                            for (int k = 0; k < numfac; k++) {
                                int index1 = ibase_idx + nregressors + ifree;
                                int index2 = jbase_idx + nregressors + ifreefac + k;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                            ifree++;
                        }
                    }
                }

                // dalpha dinter
                if (n_interaction_loadings > 0) {
                    ifree = 0;
                    for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                        if (facnorm.size() == 0 || facnorm[ifac] <= -9998) {
                            for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                                int index1 = ibase_idx + nregressors + ifree;
                                int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fac[ifac];
                            }
                            ifree++;
                        }
                    }
                }

                // dquad dbeta (for jcat > icat), dquad dalpha, dquad dquad, dquad dinter
                if (n_quadratic_loadings > 0) {
                    for (int ik = 0; ik < numfac; ik++) {
                        double fk_sq_i = fac[ik] * fac[ik];
                        int index1 = ibase_idx + nregressors + ifreefac + ik;

                        // dquad dbeta (for jcat > icat)
                        if (jcat > icat) {
                            for (int jreg = 0; jreg < nregressors; jreg++) {
                                int index2 = jbase_idx + jreg;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                            }
                        }

                        // dquad dalpha (for jcat > icat)
                        if (jcat > icat) {
                            int jfree = 0;
                            for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                                if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                                    int index2 = jbase_idx + nregressors + jfree;
                                    hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                                    jfree++;
                                }
                            }
                        }

                        // dquad dquad
                        for (int jk = 0; jk < numfac; jk++) {
                            if ((jcat > icat) || (jk >= ik)) {
                                int index2 = jbase_idx + nregressors + ifreefac + jk;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                            }
                        }

                        // dquad dinter
                        if (n_interaction_loadings > 0) {
                            for (int inter_idx = 0; inter_idx < n_interaction_loadings; inter_idx++) {
                                int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + inter_idx;
                                hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fk_sq_i;
                            }
                        }
                    }
                }

                // dinter dbeta, dinter dalpha, dinter dquad, dinter dinter
                if (n_interaction_loadings > 0) {
                    int i_inter_idx = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            double fj_fk_i = fac[ij] * fac[ik];
                            int index1 = ibase_idx + nregressors + ifreefac + n_quadratic_loadings + i_inter_idx;

                            // dinter dbeta (for jcat > icat)
                            if (jcat > icat) {
                                for (int jreg = 0; jreg < nregressors; jreg++) {
                                    int index2 = jbase_idx + jreg;
                                    hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                }
                            }

                            // dinter dalpha (for jcat > icat)
                            if (jcat > icat) {
                                int jfree = 0;
                                for (int jfac = 0; jfac < fac_loop_bound; jfac++) {
                                    if (facnorm.size() == 0 || facnorm[jfac] <= -9998) {
                                        int index2 = jbase_idx + nregressors + jfree;
                                        hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                        jfree++;
                                    }
                                }
                            }

                            // dinter dquad (for jcat > icat)
                            if (jcat > icat && n_quadratic_loadings > 0) {
                                for (int jk = 0; jk < numfac; jk++) {
                                    int index2 = jbase_idx + nregressors + ifreefac + jk;
                                    hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                }
                            }

                            // dinter dinter
                            int j_inter_idx = 0;
                            for (int jj = 0; jj < numfac - 1; jj++) {
                                for (int jk = jj + 1; jk < numfac; jk++) {
                                    if ((jcat > icat) || (j_inter_idx >= i_inter_idx)) {
                                        int index2 = jbase_idx + nregressors + ifreefac + n_quadratic_loadings + j_inter_idx;
                                        hess[index1 * npar + index2] += -pdf[icat] * logitgrad[icat * npar + index2] * fj_fk_i;
                                    }
                                    j_inter_idx++;
                                }
                            }
                            i_inter_idx++;
                        }
                    }
                }
            }
        }
    }
}

void Model::EvalOprobit(double expres, int outcome_value,
                        const std::vector<double>& fac,
                        const std::vector<double>& param, int firstpar,
                        std::vector<double>& modEval, std::vector<double>& hess,
                        int flag, const std::vector<double>& data, int iobs_offset)
{
    // Count free factor loadings
    int ifreefac = 0;
    for (size_t i = 0; i < facnorm.size(); i++) {
        if (facnorm[i] <= -9998) ifreefac++;
    }
    if (facnorm.size() == 0) ifreefac = numfac + numtyp*(outcome_idx != -2);

    // Use facnorm.size() as loop bound when provided, to avoid out-of-bounds access
    int fac_loop_bound = (facnorm.size() > 0) ? (int)facnorm.size() : numfac;

    // Observed category (1, 2, ..., numchoice)
    int obsCat = outcome_value;

    // Threshold parameters come after betas, alphas, and second-order loadings
    int thresh_idx = firstpar + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings;

    // Build thresholds by accumulating absolute values
    // threshold[0] = lower bound for this category
    // threshold[1] = upper bound for this category
    // Use minimum increment of 0.01 for numerical stability
    const double MIN_THRESH_INCREMENT = 0.01;
    double threshold[2];
    threshold[0] = param[thresh_idx];
    threshold[1] = param[thresh_idx];

    for (int icat = 2; icat <= obsCat; icat++) {
        double incr = std::fabs(param[thresh_idx + icat - 1]);
        if (incr < MIN_THRESH_INCREMENT) incr = MIN_THRESH_INCREMENT;
        if (icat < obsCat) threshold[0] += incr;
        if (icat < numchoice) threshold[1] += incr;
    }

    // Compute CDFs at thresholds
    double CDF[2] = {0.0, 1.0};
    if (obsCat > 1) {
        CDF[0] = normal_cdf(threshold[0] - expres);
    }
    if (obsCat < numchoice) {
        CDF[1] = normal_cdf(threshold[1] - expres);
    }

    // Likelihood: Prob(Y = obsCat) = CDF[upper] - CDF[lower]
    double rawDiffCDF = CDF[1] - CDF[0];

    // Floor the probability to avoid numerical underflow
    // This is critical for multi-factor models where probabilities are multiplied
    const double MIN_PROB = 1.0e-50;
    double diffCDF = rawDiffCDF;
    if (diffCDF < MIN_PROB) diffCDF = MIN_PROB;

    // Return floored probability as the likelihood
    modEval[0] = diffCDF;

    // Fix numerical problems for gradient computation (will divide by CDF later)
    if ((obsCat > 1) && (CDF[0] < MIN_PROB)) CDF[0] = MIN_PROB;
    if (CDF[1] < MIN_PROB) CDF[1] = MIN_PROB;

    if (flag < 2) return;

    // ===== GRADIENT CALCULATION =====

    // Compute PDFs and Z-values at both thresholds
    // Add safeguards for extreme Z-values
    const double MAX_Z = 35.0;  // Beyond this, PDF is essentially 0
    const double MIN_PDF = 1.0e-50;
    double Z[2] = {-9999.0, -9999.0};
    double PDF[2] = {0.0, 0.0};

    if (obsCat > 1) {
        Z[0] = threshold[0] - expres;
        // Bound Z to avoid extreme values in Hessian computation
        if (Z[0] > MAX_Z) Z[0] = MAX_Z;
        if (Z[0] < -MAX_Z) Z[0] = -MAX_Z;
        PDF[0] = normal_pdf(Z[0]);
        if (PDF[0] < MIN_PDF) PDF[0] = MIN_PDF;
    }
    if (obsCat < numchoice) {
        Z[1] = threshold[1] - expres;
        // Bound Z to avoid extreme values in Hessian computation
        if (Z[1] > MAX_Z) Z[1] = MAX_Z;
        if (Z[1] < -MAX_Z) Z[1] = -MAX_Z;
        PDF[1] = normal_pdf(Z[1]);
        if (PDF[1] < MIN_PDF) PDF[1] = MIN_PDF;
    }

    int npar = numfac + nregressors + ifreefac + n_quadratic_loadings + n_interaction_loadings + (numchoice - 1);  // includes thresholds

    if (flag == 3) {
        hess.resize(npar * npar, 0.0);
        for (int i = 0; i < npar; i++) {
            for (int j = i; j < npar; j++) {
                hess[i*npar + j] = 0.0;
            }
        }
    }

    // Loop over two terms: lower threshold (iterm=0) and upper threshold (iterm=1)
    for (int iterm = 0; iterm < 2; iterm++) {
        // Skip if they are end categories
        if (((obsCat > 1) && (iterm == 0)) || ((obsCat < numchoice) && (iterm == 1))) {

            std::vector<double> tmpgrad(npar, 0.0);

            // Gradients w.r.t. factor variance parameters and loadings
            int ifree = 0;
            for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                if (facnorm.size() == 0 || facnorm[ifac] <= -9998.0) {
                    // Free loading
                    if (ifac < numfac) {
                        // d/dtheta (factor variance)
                        tmpgrad[ifac] = (-1.0 * param[ifree + firstpar + nregressors]) * PDF[iterm] / diffCDF;
                    }
                    // d/dalpha (factor loading)
                    tmpgrad[numfac + nregressors + ifree] = (-1.0 * fac[ifac]) * PDF[iterm] / diffCDF;
                    ifree++;
                } else {
                    // Fixed loading
                    if (ifac < numfac) {
                        tmpgrad[ifac] = (-1.0 * facnorm[ifac]) * PDF[iterm] / diffCDF;
                    }
                }
            }

            // Gradients w.r.t. regression coefficients
            for (int ireg = 0; ireg < nregressors; ireg++) {
                tmpgrad[ireg + numfac] = (-1.0 * data[iobs_offset + regressors[ireg]]) * PDF[iterm] / diffCDF;
            }

            // Gradients w.r.t. quadratic factor loadings
            if (n_quadratic_loadings > 0) {
                int quad_grad_start = numfac + nregressors + ifreefac;
                for (int k = 0; k < numfac; k++) {
                    double fk_sq = fac[k] * fac[k];
                    tmpgrad[quad_grad_start + k] = (-1.0 * fk_sq) * PDF[iterm] / diffCDF;
                }
            }

            // Gradients w.r.t. interaction factor loadings
            if (n_interaction_loadings > 0) {
                int inter_grad_start = numfac + nregressors + ifreefac + n_quadratic_loadings;
                int inter_idx = 0;
                for (int j = 0; j < numfac - 1; j++) {
                    for (int k = j + 1; k < numfac; k++) {
                        double fj_fk = fac[j] * fac[k];
                        tmpgrad[inter_grad_start + inter_idx] = (-1.0 * fj_fk) * PDF[iterm] / diffCDF;
                        inter_idx++;
                    }
                }
            }

            // Gradients w.r.t. threshold parameters
            int thres_offset = npar - (numchoice - 1);
            int maxthresloop = obsCat;
            if (iterm == 0) maxthresloop--;

            for (int ithres = 0; ithres < maxthresloop; ithres++) {
                tmpgrad[thres_offset + ithres] = PDF[iterm] / diffCDF;
            }

            // Add to total gradient with appropriate sign
            int obsSign = (iterm == 0) ? -1 : 1;
            for (int i = 0; i < npar; i++) {
                modEval[i + 1] += obsSign * tmpgrad[i];
            }
        }
    }

    // ===== HESSIAN CALCULATION =====
    // Using the same factorized approach as the legacy TModel.cc code
    // Initialize tmphess to 1.0 and multiply by factors for each parameter
    if (flag == 3) {
        for (int iterm = 0; iterm < 2; iterm++) {
            // Skip if end categories
            if (((obsCat > 1) && (iterm == 0)) || ((obsCat < numchoice) && (iterm == 1))) {

                std::vector<double> tmphess(npar * npar, 1.0);
                int obsSign = (iterm == 0) ? -1 : 1;

                // Factor-specific Hessian terms (theta and alpha)
                int ifree = 0;
                for (int ifac = 0; ifac < fac_loop_bound; ifac++) {
                    // No normalizations (all free) or check if this one is free
                    if (facnorm.size() == 0 || facnorm[ifac] <= -9998.0) {
                        if (ifac < numfac) {
                            // lambda^L(theta) - lambda^Prob(theta) (row)
                            for (int j = ifac; j < npar; j++) {
                                tmphess[ifac*npar + j] *= -Z[iterm] * (-1.0 * param[ifree + firstpar + nregressors]) - modEval[1 + ifac];
                            }
                            // dZ/dtheta_i (col)
                            for (int j = 0; j <= ifac; j++) {
                                tmphess[j*npar + ifac] *= -1.0 * param[ifree + firstpar + nregressors];
                            }
                        }

                        // alpha_i index
                        int index = numfac + nregressors + ifree;
                        // lambda^L(alpha) - lambda^Prob(alpha) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * fac[ifac]) - modEval[1 + numfac + nregressors + ifree];
                        }
                        // dZ/dalpha (col)
                        for (int j = 0; j <= index; j++) {
                            tmphess[j*npar + index] *= -1.0 * fac[ifac];
                        }

                        ifree++;
                    } else {
                        // Fixed loading (facnorm[ifac] > -9998)
                        if (ifac < numfac) {
                            // lambda^L(theta) - lambda^Prob(theta) (row)
                            for (int j = ifac; j < npar; j++) {
                                tmphess[ifac*npar + j] *= -Z[iterm] * (-1.0 * facnorm[ifac]) - modEval[1 + ifac];
                            }
                            // dZ/dtheta_i (col)
                            for (int j = 0; j <= ifac; j++) {
                                tmphess[j*npar + ifac] *= -1.0 * facnorm[ifac];
                            }
                        }
                    }
                }

                // Hessian for regression coefficients (X's)
                for (int ireg = 0; ireg < nregressors; ireg++) {
                    int index = numfac + ireg;
                    // lambda^L(X) - lambda^Prob(X) (row)
                    for (int j = index; j < npar; j++) {
                        tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * data[iobs_offset + regressors[ireg]]) - modEval[1 + ireg + numfac];
                    }
                    // dZ/dX_i (col)
                    for (int j = 0; j <= index; j++) {
                        tmphess[j*npar + index] *= -1.0 * data[iobs_offset + regressors[ireg]];
                    }
                }

                // Hessian for quadratic factor loadings
                if (n_quadratic_loadings > 0) {
                    int quad_offset = numfac + nregressors + ifreefac;
                    for (int k = 0; k < numfac; k++) {
                        double fk_sq = fac[k] * fac[k];
                        int index = quad_offset + k;
                        // lambda^L(quad) - lambda^Prob(quad) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * fk_sq) - modEval[1 + quad_offset + k];
                        }
                        // dZ/dquad_k (col)
                        for (int j = 0; j <= index; j++) {
                            tmphess[j*npar + index] *= -1.0 * fk_sq;
                        }
                    }
                }

                // Hessian for interaction factor loadings
                if (n_interaction_loadings > 0) {
                    int inter_offset = numfac + nregressors + ifreefac + n_quadratic_loadings;
                    int inter_idx = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            double fj_fk = fac[ij] * fac[ik];
                            int index = inter_offset + inter_idx;
                            // lambda^L(inter) - lambda^Prob(inter) (row)
                            for (int j = index; j < npar; j++) {
                                tmphess[index*npar + j] *= -Z[iterm] * (-1.0 * fj_fk) - modEval[1 + inter_offset + inter_idx];
                            }
                            // dZ/dinter_jk (col)
                            for (int j = 0; j <= index; j++) {
                                tmphess[j*npar + index] *= -1.0 * fj_fk;
                            }
                            inter_idx++;
                        }
                    }
                }

                // Hessian for thresholds
                int thres_offset = npar - (numchoice - 1);
                int maxthresloop = obsCat;
                if (iterm == 0) maxthresloop--;

                for (int ithres = 0; ithres < numchoice - 1; ithres++) {
                    int index = thres_offset + ithres;
                    if (ithres < maxthresloop) {
                        // lambda^L(chi) - lambda^Prob(chi) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -Z[iterm] - modEval[1 + thres_offset + ithres];
                        }
                    } else {
                        // lambda^L(chi) - lambda^Prob(chi) (row)
                        for (int j = index; j < npar; j++) {
                            tmphess[index*npar + j] *= -modEval[1 + thres_offset + ithres];
                        }
                        // dZ/dchi_i (col) - set to zero
                        for (int j = 0; j <= index; j++) {
                            tmphess[j*npar + index] = 0.0;
                        }
                    }
                }

                // Add cross-derivative term dZ/dtheta dalpha = -1
                ifree = 0;
                for (int i = 0; i < numfac; i++) {
                    if (facnorm.size() == 0 || facnorm[i] <= -9998.0) {
                        int index = numfac + nregressors + ifree;
                        tmphess[i*npar + index] += -1.0;
                        ifree++;
                    }
                }

                // Add cross-derivative term dZ/dtheta_k dlambda_quad_k = -2*f_k
                if (n_quadratic_loadings > 0) {
                    int quad_offset = numfac + nregressors + ifreefac;
                    for (int k = 0; k < numfac; k++) {
                        int index = quad_offset + k;
                        tmphess[k*npar + index] += -2.0 * fac[k];
                    }
                }

                // Add cross-derivative terms dZ/dtheta dlambda_inter
                if (n_interaction_loadings > 0) {
                    int inter_offset = numfac + nregressors + ifreefac + n_quadratic_loadings;
                    int inter_idx = 0;
                    for (int ij = 0; ij < numfac - 1; ij++) {
                        for (int ik = ij + 1; ik < numfac; ik++) {
                            int index = inter_offset + inter_idx;
                            tmphess[ij*npar + index] += -fac[ik];  // d(theta_j)/d(lambda_inter_jk)
                            tmphess[ik*npar + index] += -fac[ij];  // d(theta_k)/d(lambda_inter_jk)
                            inter_idx++;
                        }
                    }
                }

                // Add tmp hessian to totals, scaled by obsSign * PDF / diffCDF
                for (int i = 0; i < npar; i++) {
                    for (int j = i; j < npar; j++) {
                        hess[i*npar + j] += obsSign * tmphess[i*npar + j] * PDF[iterm] / diffCDF;
                    }
                }
            }
        }
    }
}
