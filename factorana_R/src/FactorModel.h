#ifndef FACTORMODEL_H
#define FACTORMODEL_H

#include "Model.h"
#include <vector>
#include <memory>
#include <map>
#include <Rcpp.h>
#include <RcppEigen.h>

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
    int nparam;                        // Total number of parameters
    int nparam_free;                   // Number of free parameters

    // Factor structure
    int nfac;                          // Number of factors
    int ntyp;                          // Number of types
    int nmix;                          // Number of mixture components (default 1)
    bool fac_corr;                     // Whether factors are correlated

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

public:
    // Constructor
    FactorModel(int n_obs, int n_var, int n_fac, int n_typ = 0,
                int n_mix = 1, bool correlated = false, int n_quad = 8);

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
    // Returns vector of probabilities for each type given factor values
    std::vector<double> ComputeTypeProbabilities(const std::vector<double>& fac);
};

#endif // FACTORMODEL_H
