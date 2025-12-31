#include <Rcpp.h>
#include <RcppEigen.h>
#include <iomanip>
#include "FactorModel.h"
#include "Model.h"
#include "gauss_hermite.h"

using namespace Rcpp;

// Expose FactorModel class to R using Rcpp modules
RCPP_MODULE(factorana_module) {
    // Use function pointers to disambiguate overloaded methods
    void (FactorModel::*setData1)(const std::vector<double>&) = &FactorModel::SetData;
    void (FactorModel::*setData2)(const Eigen::MatrixXd&) = &FactorModel::SetData;
    void (FactorModel::*setConstraints1)(const std::vector<bool>&) = &FactorModel::SetParameterConstraints;
    void (FactorModel::*setConstraints2)(const std::vector<bool>&, const std::vector<double>&) = &FactorModel::SetParameterConstraints;

    class_<FactorModel>("FactorModel")
        .constructor<int, int, int, int, int, bool, int>("Create a FactorModel")
        .method("SetDataVector", setData1)
        .method("SetDataMatrix", setData2)
        .method("SetQuadrature", &FactorModel::SetQuadrature)
        .method("SetParameterConstraints", setConstraints1)
        .method("SetParameterConstraintsWithValues", setConstraints2)
        .method("CalcLogLikelihood", &FactorModel::CalcLogLikelihood)
        .method("GetNObs", &FactorModel::GetNObs)
        .method("GetNParam", &FactorModel::GetNParam)
        .method("GetNParamFree", &FactorModel::GetNParamFree)
        ;
}

//' Compute Gauss-Hermite quadrature nodes and weights
//'
//' @param n Number of quadrature points
//' @return List with nodes and weights
//' @export
// [[Rcpp::export]]
List gauss_hermite_quadrature(int n) {
    std::vector<double> nodes, weights;
    calcgausshermitequadrature(n, nodes, weights);

    return List::create(
        Named("nodes") = nodes,
        Named("weights") = weights
    );
}

//' Initialize a FactorModel C++ object from R model system
//'
//' @param model_system R model_system object
//' @param data Data frame or matrix with all variables
//' @param n_quad Number of quadrature points
//' @param init_params Optional initial parameter vector (used to set fixed parameter values)
//' @return External pointer to FactorModel object
//' @export
// [[Rcpp::export]]
SEXP initialize_factor_model_cpp(List model_system, SEXP data, int n_quad = 8,
                                  Nullable<NumericVector> init_params = R_NilValue) {

    // Extract factor model information
    List factor_model = model_system["factor"];
    int n_fac = factor_model["n_factors"];
    int n_types = factor_model["n_types"];
    int n_mixtures = factor_model["n_mixtures"];
    bool correlation = factor_model["correlation"];

    // Convert data to matrix (handle both data.frame and matrix)
    NumericMatrix data_mat;
    if (Rf_isMatrix(data)) {
        data_mat = as<NumericMatrix>(data);
    } else {
        // Assume it's a data frame
        DataFrame df = as<DataFrame>(data);
        data_mat = internal::convert_using_rfunction(df, "as.matrix");
    }

    int n_obs = data_mat.nrow();
    int n_var = data_mat.ncol();

    // Get column names for variable indexing
    CharacterVector col_names = colnames(data_mat);
    if (col_names.size() == 0) {
        Rcpp::stop("Data must have column names for variable indexing");
    }

    // Create FactorModel object
    Rcpp::XPtr<FactorModel> fm(
        new FactorModel(n_obs, n_var, n_fac, n_types, n_mixtures, correlation, n_quad),
        true
    );

    // Set data
    std::vector<double> data_vec(n_obs * n_var);
    for (int i = 0; i < n_obs; i++) {
        for (int j = 0; j < n_var; j++) {
            data_vec[i * n_var + j] = data_mat(i, j);
        }
    }
    fm->SetData(data_vec);

    // Compute and set quadrature
    std::vector<double> nodes, weights;

    // Special case: if n_fac == 0, use single integration point at 0 with weight 1
    if (n_fac == 0) {
        nodes.push_back(0.0);
        weights.push_back(1.0);
    } else {
        calcgausshermitequadrature(n_quad, nodes, weights);

        // Scale for standard normal: x_scaled = sqrt(2) * x, w_scaled = w / sqrt(pi)
        double sqrt2 = std::sqrt(2.0);
        double sqrt_pi = std::sqrt(M_PI);
        for (int i = 0; i < n_quad; i++) {
            nodes[i] *= sqrt2;
            weights[i] /= sqrt_pi;
        }
    }

    fm->SetQuadrature(nodes, weights);

    // Add model components
    List components = model_system["components"];
    for (int i = 0; i < components.size(); i++) {
        List comp = components[i];

        // Extract model component information
        std::string model_type_str = as<std::string>(comp["model_type"]);
        std::string outcome_name = as<std::string>(comp["outcome"]);

        // Handle covariates = NULL case
        std::vector<std::string> covariate_names;
        if (!Rf_isNull(comp["covariates"])) {
            covariate_names = as<std::vector<std::string>>(comp["covariates"]);
        }

        // Check if this is a dynamic model first
        bool is_dynamic = comp.containsElementNamed("is_dynamic") &&
                         as<bool>(comp["is_dynamic"]);

        // Find variable indices in data
        int outcome_idx = -1;
        for (int j = 0; j < col_names.size(); j++) {
            if (std::string(col_names[j]) == outcome_name) {
                outcome_idx = j;
                break;
            }
        }
        // For dynamic models, outcome is always zero, so use special marker -2
        // (distinct from -1 which means "missing indicator not present")
        if (outcome_idx == -1) {
            if (is_dynamic) {
                outcome_idx = -2;  // Special marker: outcome is always zero
            } else {
                Rcpp::stop("Outcome variable '" + outcome_name + "' not found in data");
            }
        }

        // Find covariate indices
        std::vector<int> regressor_idx;
        for (const auto& cov_name : covariate_names) {
            int idx = -1;
            for (int j = 0; j < col_names.size(); j++) {
                if (std::string(col_names[j]) == cov_name) {
                    idx = j;
                    break;
                }
            }
            // For dynamic models, 'intercept' covariate is always 1
            // Use special marker -3 if not found in data
            if (idx == -1) {
                if (is_dynamic && cov_name == "intercept") {
                    idx = -3;  // Special marker: intercept is always 1
                } else {
                    Rcpp::stop("Covariate '" + cov_name + "' not found in data");
                }
            }
            regressor_idx.push_back(idx);
        }

        // Find missing/evaluation indicator if present
        int missing_idx = -1;
        if (comp.containsElementNamed("evaluation_indicator") &&
            !Rf_isNull(comp["evaluation_indicator"])) {
            std::string eval_name = as<std::string>(comp["evaluation_indicator"]);
            for (int j = 0; j < col_names.size(); j++) {
                if (std::string(col_names[j]) == eval_name) {
                    missing_idx = j;
                    break;
                }
            }
        }

        // Determine model type
        ModelType mtype;
        if (model_type_str == "linear") mtype = ModelType::LINEAR;
        else if (model_type_str == "probit") mtype = ModelType::PROBIT;
        else if (model_type_str == "logit") mtype = ModelType::LOGIT;
        else if (model_type_str == "oprobit") mtype = ModelType::OPROBIT;
        else {
            Rcpp::stop("Unknown model type: " + model_type_str);
        }

        // Extract factor normalizations from COMPONENT (not factor model)
        // This allows each component to have its own loading constraints
        std::vector<double> facnorm;
        if (comp.containsElementNamed("loading_normalization")) {
            SEXP norm_sexp = comp["loading_normalization"];
            if (!Rf_isNull(norm_sexp)) {
                NumericVector norm_vec = as<NumericVector>(norm_sexp);
                for (int j = 0; j < norm_vec.size(); j++) {
                    if (NumericVector::is_na(norm_vec[j])) {
                        facnorm.push_back(-9999.0);  // Free parameter
                    } else {
                        facnorm.push_back(norm_vec[j]);  // Fixed value
                    }
                }
            }
        }

        int n_choice = comp.containsElementNamed("num_choices") ?
                      int(comp["num_choices"]) : 2;

        // Check if all parameters are fixed (for multi-stage estimation)
        bool all_params_fixed = false;
        if (comp.containsElementNamed("all_params_fixed")) {
            all_params_fixed = as<bool>(comp["all_params_fixed"]);
        }

        // Extract factor_spec (linear, quadratic, interactions, full)
        FactorSpec fspec = FactorSpec::LINEAR;
        if (comp.containsElementNamed("factor_spec")) {
            std::string fs = as<std::string>(comp["factor_spec"]);
            if (fs == "quadratic") fspec = FactorSpec::QUADRATIC;
            else if (fs == "interactions") fspec = FactorSpec::INTERACTIONS;
            else if (fs == "full") fspec = FactorSpec::FULL;
        }

        // Extract dynamic model outcome_factor index (is_dynamic already set above)
        int outcome_factor_idx = -1;
        if (is_dynamic) {
            // outcome_factor is 1-based in R, convert to 0-based for C++
            outcome_factor_idx = as<int>(comp["outcome_factor"]) - 1;
        }

        // Create Model object
        std::shared_ptr<Model> model = std::make_shared<Model>(
            mtype, outcome_idx, missing_idx, regressor_idx,
            n_fac, n_types, facnorm, n_choice, 1, all_params_fixed, fspec,
            is_dynamic, outcome_factor_idx
        );

        // Calculate number of parameters for this model
        int n_free_loadings = 0;
        for (const auto& norm : facnorm) {
            if (norm <= -9998.0) n_free_loadings++;
        }
        if (facnorm.empty()) n_free_loadings = n_fac;

        // Calculate second-order loading counts from model (computed in constructor)
        int n_quad = model->GetNumQuadraticLoadings();
        int n_inter = model->GetNumInteractionLoadings();

        int n_params;
        if (mtype == ModelType::LOGIT && n_choice > 2) {
            // Multinomial logit: each non-reference choice has its own parameters
            n_params = (n_choice - 1) * (regressor_idx.size() + n_free_loadings + n_quad + n_inter);
            // For multinomial logit with types, each choice gets (n_types - 1) type-specific intercepts
            if (n_types > 1) {
                n_params += (n_choice - 1) * (n_types - 1);
            }
        } else if (mtype == ModelType::OPROBIT) {
            // Ordered probit: shared coefficients + thresholds
            n_params = regressor_idx.size() + n_free_loadings + n_quad + n_inter + (n_choice - 1);
            // Add type-specific intercepts for n_types > 1
            if (n_types > 1) {
                n_params += (n_types - 1);
            }
        } else {
            // Binary models (linear, probit, binary logit)
            n_params = regressor_idx.size() + n_free_loadings + n_quad + n_inter;
            if (mtype == ModelType::LINEAR) n_params += 1;  // sigma
            // Add type-specific intercepts for n_types > 1
            if (n_types > 1) {
                n_params += (n_types - 1);
            }
        }

        fm->AddModel(model, n_params);
    }

    // Build parameter constraints based on fixed_coefficients from each component
    int total_params = fm->GetNParam();
    std::vector<bool> param_fixed_vec(total_params, false);

    // Track parameter position as we go through components
    // Start after factor variance parameters (and correlation/type params)
    int param_offset = n_fac;  // Factor variances

    // Add correlation parameter offset if present
    if (correlation && n_fac == 2) {
        param_offset += 1;  // One correlation parameter for 2-factor model
    }

    // Add type model parameters offset if n_types > 1
    if (n_types > 1) {
        param_offset += (n_types - 1) * n_fac;  // Type loadings
    }

    // Process each component's fixed coefficients
    for (int i = 0; i < components.size(); i++) {
        List comp = components[i];

        std::string model_type_str = as<std::string>(comp["model_type"]);
        std::vector<std::string> covariate_names;
        if (!Rf_isNull(comp["covariates"])) {
            covariate_names = as<std::vector<std::string>>(comp["covariates"]);
        }

        int n_choice = comp.containsElementNamed("num_choices") ?
                      int(comp["num_choices"]) : 2;

        // Count free loadings for this component
        int n_free_loadings = 0;
        if (comp.containsElementNamed("loading_normalization")) {
            SEXP norm_sexp = comp["loading_normalization"];
            if (!Rf_isNull(norm_sexp)) {
                NumericVector norm_vec = as<NumericVector>(norm_sexp);
                for (int j = 0; j < norm_vec.size(); j++) {
                    if (NumericVector::is_na(norm_vec[j])) {
                        n_free_loadings++;
                    }
                }
            }
        } else {
            n_free_loadings = n_fac;
        }

        // Check for fixed_coefficients in this component
        if (comp.containsElementNamed("fixed_coefficients") &&
            !Rf_isNull(comp["fixed_coefficients"])) {
            List fixed_coefs = comp["fixed_coefficients"];

            for (int fc_idx = 0; fc_idx < fixed_coefs.size(); fc_idx++) {
                List fc = fixed_coefs[fc_idx];
                std::string cov_name = as<std::string>(fc["covariate"]);

                // Find covariate position
                int cov_pos = -1;
                for (int k = 0; k < covariate_names.size(); k++) {
                    if (covariate_names[k] == cov_name) {
                        cov_pos = k;
                        break;
                    }
                }

                if (cov_pos >= 0) {
                    // Determine parameter index based on model type and choice
                    int param_idx;

                    if (model_type_str == "logit" && n_choice > 2) {
                        // Multinomial logit: check if choice is specified
                        int choice = 1;  // Default to first non-reference choice
                        if (fc.containsElementNamed("choice") && !Rf_isNull(fc["choice"])) {
                            choice = as<int>(fc["choice"]);
                        }
                        // Each choice has: covariates + loadings
                        int params_per_choice = covariate_names.size() + n_free_loadings;
                        param_idx = param_offset + (choice - 1) * params_per_choice + cov_pos;
                    } else {
                        // Binary/linear/probit/oprobit: coefficients come first
                        param_idx = param_offset + cov_pos;
                    }

                    // Mark as fixed
                    if (param_idx < total_params) {
                        param_fixed_vec[param_idx] = true;
                    }
                }
            }
        }

        // Advance param_offset for next component
        // Calculate n_params for this component (same logic as above)
        int n_params_comp;
        ModelType mtype;
        if (model_type_str == "linear") mtype = ModelType::LINEAR;
        else if (model_type_str == "probit") mtype = ModelType::PROBIT;
        else if (model_type_str == "logit") mtype = ModelType::LOGIT;
        else mtype = ModelType::OPROBIT;

        if (mtype == ModelType::LOGIT && n_choice > 2) {
            n_params_comp = (n_choice - 1) * (covariate_names.size() + n_free_loadings);
            if (n_types > 1) n_params_comp += (n_choice - 1) * (n_types - 1);
        } else if (mtype == ModelType::OPROBIT) {
            n_params_comp = covariate_names.size() + n_free_loadings + (n_choice - 1);
            if (n_types > 1) n_params_comp += (n_types - 1);
        } else {
            n_params_comp = covariate_names.size() + n_free_loadings;
            if (mtype == ModelType::LINEAR) n_params_comp += 1;
            if (n_types > 1) n_params_comp += (n_types - 1);
        }
        param_offset += n_params_comp;
    }

    // Set parameter constraints with optional initial values
    if (init_params.isNotNull()) {
        NumericVector ip(init_params);
        std::vector<double> init_params_vec = as<std::vector<double>>(ip);
        fm->SetParameterConstraints(param_fixed_vec, init_params_vec);
    } else {
        fm->SetParameterConstraints(param_fixed_vec);
    }

    return fm;
}

//' Evaluate log-likelihood for given parameters
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param params Vector of parameters
//' @param compute_gradient Whether to compute gradient (default FALSE)
//' @param compute_hessian Whether to compute Hessian (default FALSE)
//' @return List with:
//'   - logLikelihood: scalar log-likelihood value
//'   - gradient: vector of length n_param_free (if requested)
//'   - hessian: vector of length n_param_free*(n_param_free+1)/2 stored as
//'              upper-triangular in row-major order (if requested).
//'              To expand to full symmetric matrix in R:
//'              \code{idx <- 1; for(i in 1:n) for(j in i:n) { H[i,j] <- H[j,i] <- hess[idx]; idx <- idx + 1 }}
//' @export
// [[Rcpp::export]]
List evaluate_likelihood_cpp(SEXP fm_ptr, NumericVector params,
                             bool compute_gradient = false,
                             bool compute_hessian = false) {

    // Get FactorModel object
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    // Convert parameters to std::vector
    std::vector<double> params_vec = as<std::vector<double>>(params);

    // Determine flag
    int iflag = 1;  // likelihood only
    if (compute_gradient) iflag = 2;
    if (compute_hessian) iflag = 3;

    // Compute likelihood
    double logLkhd;
    std::vector<double> gradL, hessL;
    fm->CalcLkhd(params_vec, logLkhd, gradL, hessL, iflag);

    // Return results
    List result = List::create(Named("logLikelihood") = logLkhd);

    if (compute_gradient) {
        result["gradient"] = gradL;
    }
    if (compute_hessian) {
        result["hessian"] = hessL;
    }

    return result;
}

//' Evaluate log-likelihood only (for optimization)
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param params Vector of parameters
//' @return Log-likelihood value
//' @export
// [[Rcpp::export]]
double evaluate_loglik_only_cpp(SEXP fm_ptr, NumericVector params) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);
    std::vector<double> params_vec = as<std::vector<double>>(params);
    return fm->CalcLogLikelihood(params_vec);
}

//' Get parameter counts from FactorModel
//'
//' @param fm_ptr External pointer to FactorModel object
//' @return List with parameter count information
//' @export
// [[Rcpp::export]]
List get_parameter_info_cpp(SEXP fm_ptr) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    return List::create(
        Named("n_obs") = fm->GetNObs(),
        Named("n_param") = fm->GetNParam(),
        Named("n_param_free") = fm->GetNParamFree()
    );
}

//' Extract free parameters from full parameter vector
//'
//' Given a full parameter vector (including fixed parameters),
//' extract only the free parameters based on the model's fixed parameter mask.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param full_params Full parameter vector (size n_param)
//' @return Vector of free parameters only (size n_param_free)
//' @export
// [[Rcpp::export]]
NumericVector extract_free_params_cpp(SEXP fm_ptr, NumericVector full_params) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    int n_param = fm->GetNParam();
    int n_param_free = fm->GetNParamFree();

    if (full_params.size() != n_param) {
        Rcpp::stop("Full params size (%d) doesn't match expected (%d)",
                   full_params.size(), n_param);
    }

    // Get the fixed parameter mask
    const std::vector<bool>& param_fixed = fm->GetParamFixed();

    // Extract free parameters
    NumericVector free_params(n_param_free);
    int ifree = 0;
    for (int i = 0; i < n_param; i++) {
        if (!param_fixed[i]) {
            free_params[ifree++] = full_params[i];
        }
    }

    return free_params;
}

//' Set observation weights for weighted likelihood estimation
//'
//' Sets per-observation weights for the likelihood calculation. When weights
//' are set, each observation's contribution to the log-likelihood is multiplied
//' by its weight. This is used for importance sampling in adaptive integration.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param weights Numeric vector of observation weights (length = n_obs)
//' @export
// [[Rcpp::export]]
void set_observation_weights_cpp(SEXP fm_ptr, NumericVector weights) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);
    std::vector<double> weights_vec = as<std::vector<double>>(weights);
    fm->SetObservationWeights(weights_vec);
}

//' Set up adaptive quadrature based on factor scores and standard errors
//'
//' Enables adaptive integration where the number of quadrature points varies
//' by observation based on the precision of factor score estimates. When factor
//' scores are well-determined (small SE), fewer integration points are used.
//' Importance sampling weights are computed automatically.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param factor_scores Matrix (n_obs x n_factors) of factor score estimates
//' @param factor_ses Matrix (n_obs x n_factors) of standard errors
//' @param factor_vars Vector (n_factors) of factor variances from previous stage
//' @param threshold Threshold for determining quadrature points (default 0.3)
//' @param max_quad Maximum quadrature points per factor (default 16)
//' @param verbose Whether to print summary of adaptive quadrature setup (default TRUE)
//' @export
// [[Rcpp::export]]
void set_adaptive_quadrature_cpp(SEXP fm_ptr,
                                 NumericMatrix factor_scores,
                                 NumericMatrix factor_ses,
                                 NumericVector factor_vars,
                                 double threshold = 0.3,
                                 int max_quad = 16,
                                 bool verbose = true) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    int n_obs = factor_scores.nrow();
    int n_fac = factor_scores.ncol();

    // Convert to nested vectors
    std::vector<std::vector<double>> scores_vec(n_obs);
    std::vector<std::vector<double>> ses_vec(n_obs);
    for (int i = 0; i < n_obs; i++) {
        scores_vec[i].resize(n_fac);
        ses_vec[i].resize(n_fac);
        for (int j = 0; j < n_fac; j++) {
            scores_vec[i][j] = factor_scores(i, j);
            ses_vec[i][j] = factor_ses(i, j);
        }
    }

    std::vector<double> vars_vec = as<std::vector<double>>(factor_vars);

    fm->SetAdaptiveQuadrature(scores_vec, ses_vec, vars_vec, threshold, max_quad);

    // Print summary if verbose
    if (verbose) {
        Rcpp::Rcout << "\nAdaptive Integration Summary\n";
        Rcpp::Rcout << "----------------------------\n";
        Rcpp::Rcout << "Threshold: " << threshold << ", Max quad points: " << max_quad << "\n\n";

        // Compute per-observation total integration points and collect stats
        std::map<int, int> nquad_counts;  // nquad -> count of observations
        double total_points = 0.0;
        double total_points_standard = std::pow(max_quad, n_fac);

        for (int i = 0; i < n_obs; i++) {
            int obs_total = 1;
            for (int j = 0; j < n_fac; j++) {
                // Recompute nquad using same formula as SetAdaptiveQuadrature
                double f_se = factor_ses(i, j);
                double f_var = factor_vars[j];
                double ratio = f_se / f_var / threshold;
                int nq = 1 + 2 * static_cast<int>(std::floor(ratio));
                if (nq < 1) nq = 1;
                if (nq > max_quad) nq = max_quad;
                if (f_se > std::sqrt(f_var)) nq = max_quad;
                obs_total *= nq;
            }
            nquad_counts[obs_total]++;
            total_points += obs_total;
        }

        double avg_points = total_points / n_obs;
        double reduction = 100.0 * (1.0 - avg_points / total_points_standard);

        // Print distribution table
        Rcpp::Rcout << "Integration points per observation:\n";
        Rcpp::Rcout << "  Points   Observations   Percent\n";
        for (auto& kv : nquad_counts) {
            double pct = 100.0 * kv.second / n_obs;
            Rcpp::Rcout << std::setw(8) << kv.first
                        << std::setw(15) << kv.second
                        << std::setw(10) << std::fixed << std::setprecision(1) << pct << "%\n";
        }

        Rcpp::Rcout << "\nAverage integration points: " << std::fixed << std::setprecision(1)
                    << avg_points << " (vs " << static_cast<int>(total_points_standard)
                    << " standard)\n";
        Rcpp::Rcout << "Computational reduction: " << std::fixed << std::setprecision(1)
                    << reduction << "%\n\n";
    }
}

//' Disable adaptive quadrature
//'
//' Reverts to standard (non-adaptive) quadrature integration.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @export
// [[Rcpp::export]]
void disable_adaptive_quadrature_cpp(SEXP fm_ptr) {
    Rcpp::XPtr<FactorModel> fm(fm_ptr);
    fm->DisableAdaptiveQuadrature();
}

//' Evaluate log-likelihood for a single observation at given factor values
//'
//' Used for factor score estimation. The model parameters are held fixed,
//' and the factor values are treated as the parameters to optimize.
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param iobs Observation index (0-based)
//' @param factor_values Vector of factor values (size n_factors)
//' @param model_params Vector of ALL model parameters (from previous estimation)
//' @param compute_gradient Whether to compute gradient (default FALSE)
//' @param compute_hessian Whether to compute Hessian (default FALSE)
//' @return List with log-likelihood, gradient (if requested), and Hessian (if requested)
//' @export
// [[Rcpp::export]]
List evaluate_factorscore_likelihood_cpp(SEXP fm_ptr,
                                         int iobs,
                                         NumericVector factor_values,
                                         NumericVector model_params,
                                         bool compute_gradient = false,
                                         bool compute_hessian = false) {

    // Get FactorModel object
    Rcpp::XPtr<FactorModel> fm(fm_ptr);

    // Convert to std::vector
    std::vector<double> fac_vec = as<std::vector<double>>(factor_values);
    std::vector<double> param_vec = as<std::vector<double>>(model_params);

    // Determine flag
    int iflag = 1;  // likelihood only
    if (compute_gradient) iflag = 2;
    if (compute_hessian) iflag = 3;

    // Compute likelihood for single observation
    double logLkhd;
    std::vector<double> gradL, hessL;
    fm->CalcLkhdSingleObs(iobs, fac_vec, param_vec, logLkhd, gradL, hessL, iflag);

    // Return results
    List result = List::create(Named("logLikelihood") = logLkhd);

    if (compute_gradient) {
        result["gradient"] = gradL;
    }
    if (compute_hessian) {
        result["hessian"] = hessL;
    }

    return result;
}
