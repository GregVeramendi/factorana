#ifndef FACTORMODEL_H
#define FACTORMODEL_H

#include "Model.h"
#include <vector>
#include <memory>
#include <map>
#include <Rcpp.h>
#include <RcppEigen.h>

// Factor structure types
enum class FactorStructure {
    INDEPENDENT,   // Factors are independent
    CORRELATION,   // Correlated factors via Cholesky (2 factors only)
    SE_LINEAR,     // Structural equation: f_k = alpha + alpha_1*f_1 + ... + epsilon
    SE_QUADRATIC   // SE with quadratic terms: f_k = alpha + alpha_1*f_1 + alpha_q1*f_1^2 + ... + epsilon
};

class FactorModel {
private:
    // Data
    int nobs;                          // Number of observations
    int nvar;                          // Number of variables per observation
    std::vector<double> data;          // Flattened data (nobs * nvar)

    // Model components
    std::vector<std::shared_ptr<Model>> models;  // Vector of model pointers

    // Parameters
    std::vector<double> param;         // Full parameter vector
    std::vector<bool> param_fixed;     // Which parameters are fixed
    std::vector<int> freeparlist;      // Indices of free parameters (for Hessian optimization)
    std::vector<int> gradparlist;      // Indices of params needing gradients (free + tied)
    int nparam;                        // Total number of parameters
    int nparam_free;                   // Number of free parameters

    // Factor structure
    int nfac;                          // Number of factors
    int ntyp;                          // Number of types
    int nmix;                          // Number of mixture components (default 1)
    bool fac_corr;                     // Whether factors are correlated (legacy, for backward compat)
    FactorStructure factor_structure;  // Factor dependency structure

    // Structural equation (SE) model parameters
    // For SE_LINEAR: f_k = se_intercept + se_linear[j] * f_j + epsilon
    // Integration is over input factors and residual epsilon
    int n_input_factors;               // Number of input factors (nfac - 1 for SE models)
    int n_outcome_factors;             // Number of outcome factors (1 for SE_LINEAR)
    int se_param_start;                // Starting index for SE parameters in param vector
    int nse_param;                     // Total SE parameters

    // Type model parameters (for n_types > 1)
    // Type probability model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
    // Total: (n_types - 1) * n_factors parameters
    int ntyp_param;                    // Number of type model loading parameters
    int type_param_start;              // Starting index for type model parameters

    // Quadrature
    int nquad_points;                  // Number of quadrature points per dimension
    std::vector<double> quad_nodes;    // GH quadrature nodes (scaled for different npoints)
    std::vector<double> quad_weights;  // GH quadrature weights

    // Observation weights (for weighted likelihood)
    std::vector<double> obs_weights;   // Weight per observation (default: all 1.0)
    bool use_weights;                  // Whether observation weights are set

    // Adaptive integration (for second-stage estimation with factor scores)
    bool use_adaptive;                 // Whether adaptive integration is enabled
    std::vector<std::vector<int>> obs_nquad;     // Per-obs, per-factor quadrature point counts
    std::vector<std::vector<double>> obs_fac_center;  // Per-obs factor score centers [iobs][ifac]
    std::vector<double> adapt_factor_var;       // Factor variances for adaptive mode
    double adapt_threshold;                     // Threshold for determining n_quad per obs
    std::map<int, std::vector<double>> adapt_nodes;   // GH nodes for different nquad values
    std::map<int, std::vector<double>> adapt_weights; // GH weights for different nquad values

    // Parameter organization
    std::vector<int> param_model_start;  // Starting index for each model's parameters
    std::vector<int> param_model_count;  // Number of parameters per model

    // Equality constraints: maps tied param index -> primary (free) param index
    // Value of -1 means parameter is not tied to any other
    std::vector<int> equality_mapping;

public:
    // Constructor (legacy, for backward compatibility)
    FactorModel(int n_obs, int n_var, int n_fac, int n_typ = 0,
                int n_mix = 1, bool correlated = false, int n_quad = 8);

    // Constructor with factor_structure support
    FactorModel(int n_obs, int n_var, int n_fac, int n_typ,
                int n_mix, FactorStructure fac_struct, int n_quad);

    // Add a model component
    void AddModel(std::shared_ptr<Model> model, int nparams);

    // Set data (flattened vector)
    void SetData(const std::vector<double>& dat);
    void SetData(const Eigen::MatrixXd& dat);  // Alternative: matrix form

    // Set quadrature nodes and weights (computed in R)
    void SetQuadrature(const std::vector<double>& nodes,
                      const std::vector<double>& weights);

    // Set observation weights (for weighted likelihood / importance sampling)
    void SetObservationWeights(const std::vector<double>& weights);

    // Set up adaptive integration based on factor scores and standard errors
    // factor_scores: matrix [nobs x nfac] of factor score estimates
    // factor_ses: matrix [nobs x nfac] of standard errors for factor scores
    // factor_vars: vector [nfac] of factor variances from previous stage
    // threshold: threshold for determining quadrature points (smaller = more points)
    // max_quad: maximum number of quadrature points per factor
    void SetAdaptiveQuadrature(const std::vector<std::vector<double>>& factor_scores,
                               const std::vector<std::vector<double>>& factor_ses,
                               const std::vector<double>& factor_vars,
                               double threshold,
                               int max_quad);

    // Disable adaptive integration (revert to standard quadrature)
    void DisableAdaptiveQuadrature();

    // Set which parameters are fixed
    void SetParameterConstraints(const std::vector<bool>& fixed);

    // Set which parameters are fixed and their values
    void SetParameterConstraints(const std::vector<bool>& fixed,
                                 const std::vector<double>& fixed_values);

    // Set equality constraints: tied parameters are mapped to their primary (free) parameter
    // equality_map[i] = j means parameter i is tied to parameter j (j is the primary)
    // equality_map[i] = -1 means parameter i is not tied to any other
    void SetEqualityConstraints(const std::vector<int>& equality_map);

    // Main likelihood calculation
    // Computes likelihood, gradient, and Hessian for current parameter values
    //
    // Parameters:
    //   free_params - Vector of free parameters only
    //   logLkhd - Output: log-likelihood value
    //   gradL - Output: gradient vector (size nparam_free)
    //   hessL - Output: Hessian upper triangle (size nparam_free*(nparam_free+1)/2)
    //   iflag - 1=likelihood only, 2=+gradient, 3=+Hessian
    void CalcLkhd(const std::vector<double>& free_params,
                  double& logLkhd,
                  std::vector<double>& gradL,
                  std::vector<double>& hessL,
                  int iflag);

    // Convenience wrapper that returns likelihood only
    double CalcLogLikelihood(const std::vector<double>& free_params);

    // Factor score estimation: evaluate likelihood for single observation
    // given factor values (no quadrature integration)
    //
    // Parameters:
    //   iobs - Observation index (0-based)
    //   factor_values - Vector of factor values (size nfac)
    //   model_params - Vector of ALL model parameters (size nparam)
    //   logLkhd - Output: log-likelihood value (includes factor prior)
    //   gradL - Output: gradient w.r.t. factor values (size nfac)
    //   hessL - Output: Hessian w.r.t. factor values (upper triangle, size nfac*(nfac+1)/2)
    //   iflag - 1=likelihood only, 2=+gradient, 3=+Hessian
    void CalcLkhdSingleObs(int iobs,
                           const std::vector<double>& factor_values,
                           const std::vector<double>& model_params,
                           double& logLkhd,
                           std::vector<double>& gradL,
                           std::vector<double>& hessL,
                           int iflag);

    // Set model parameters (used for factor score estimation)
    void SetModelParameters(const std::vector<double>& params);

    // Accessors
    int GetNObs() const { return nobs; }
    int GetNVar() const { return nvar; }
    int GetNFac() const { return nfac; }
    int GetNParam() const { return nparam; }
    int GetNParamFree() const { return nparam_free; }
    const std::vector<bool>& GetParamFixed() const { return param_fixed; }

private:
    // Helper: Map free parameters to full parameter vector
    void MapFreeToFull(const std::vector<double>& free_params);

    // Helper: Extract gradient of free parameters from full gradient
    void ExtractFreeGradient(const std::vector<double>& full_grad,
                            std::vector<double>& free_grad);

    // Helper: Extract Hessian of free parameters from full Hessian
    void ExtractFreeHessian(const std::vector<double>& full_hess,
                           std::vector<double>& free_hess);

    // Helper: Get factor parameter indices
    int GetFactorVarianceIndex(int imix, int ifac);
    int GetFactorMeanIndex(int imix, int ifac);
    int GetMixtureWeightIndex(int imix);

    // Helper: Get type model parameter indices
    // Type model: log(P(type=t)/P(type=1)) = sum_k lambda_t_k * f_k
    // Returns index of lambda_t_k (loading for type t on factor k)
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // ifac: 0-based factor index
    int GetTypeLoadingIndex(int ityp, int ifac);

    // Helper: Get type-specific intercept index for a given model
    // ityp: 0-based type index (0 = type 2 since type 1 is reference)
    // model_idx: 0-based model index
    int GetTypeInterceptIndex(int ityp, int model_idx);

    // Helper: Compute type probability from multinomial logit
    // Fills pre-allocated probs vector with probabilities for each type given factor values
    // OPTIMIZATION: Output vector passed by reference to avoid allocation per call
    void ComputeTypeProbabilities(const std::vector<double>& fac, std::vector<double>& probs);

    // Helper: Get SE parameter indices
    // For SE_LINEAR with k factors:
    // - se_intercept at se_param_start
    // - se_linear[j] at se_param_start + 1 + j (for j = 0..n_input_factors-1)
    // - se_residual_var at se_param_start + 1 + n_input_factors
    //
    // For SE_QUADRATIC: same as above but with quadratic coefficients added:
    // - se_linear[j] at se_param_start + 1 + j
    // - se_quadratic[j] at se_param_start + 1 + n_input_factors + j
    // - se_residual_var at se_param_start + 1 + 2*n_input_factors
    int GetSEInterceptIndex() const { return se_param_start; }
    int GetSELinearIndex(int j) const { return se_param_start + 1 + j; }
    int GetSEQuadraticIndex(int j) const { return se_param_start + 1 + n_input_factors + j; }
    int GetSEResidualVarIndex() const {
        if (factor_structure == FactorStructure::SE_QUADRATIC) {
            return se_param_start + 1 + 2 * n_input_factors;
        } else {
            return se_param_start + 1 + n_input_factors;
        }
    }
};

#endif // FACTORMODEL_H
