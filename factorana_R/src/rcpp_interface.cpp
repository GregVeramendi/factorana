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

    class_<FactorModel>("FactorModel")
        .constructor<int, int, int, int, int, bool, int>("Create a FactorModel")
        .method("SetDataVector", setData1)
        .method("SetDataMatrix", setData2)
        .method("SetQuadrature", &FactorModel::SetQuadrature)
        .method("SetParameterConstraints", &FactorModel::SetParameterConstraints)
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
//' @return External pointer to FactorModel object
//' @export
// [[Rcpp::export]]
SEXP initialize_factor_model_cpp(List model_system, SEXP data, int n_quad = 8) {

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

        // Find variable indices in data
        int outcome_idx = -1;
        for (int j = 0; j < col_names.size(); j++) {
            if (std::string(col_names[j]) == outcome_name) {
                outcome_idx = j;
                break;
            }
        }
        if (outcome_idx == -1) {
            Rcpp::stop("Outcome variable '" + outcome_name + "' not found in data");
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
            if (idx == -1) {
                Rcpp::stop("Covariate '" + cov_name + "' not found in data");
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

        // Create Model object
        std::shared_ptr<Model> model = std::make_shared<Model>(
            mtype, outcome_idx, missing_idx, regressor_idx,
            n_fac, n_types, facnorm, n_choice, 1, all_params_fixed
        );

        // Calculate number of parameters for this model
        int n_free_loadings = 0;
        for (const auto& norm : facnorm) {
            if (norm <= -9998.0) n_free_loadings++;
        }
        if (facnorm.empty()) n_free_loadings = n_fac;

        int n_params;
        if (mtype == ModelType::LOGIT && n_choice > 2) {
            // Multinomial logit: each non-reference choice has its own parameters
            n_params = (n_choice - 1) * (regressor_idx.size() + n_free_loadings);
        } else if (mtype == ModelType::OPROBIT) {
            // Ordered probit: shared coefficients + thresholds
            n_params = regressor_idx.size() + n_free_loadings + (n_choice - 1);
        } else {
            // Binary models (linear, probit, binary logit)
            n_params = regressor_idx.size() + n_free_loadings;
            if (mtype == ModelType::LINEAR) n_params += 1;  // sigma
        }

        fm->AddModel(model, n_params);
    }

    // Set parameter constraints (all free for now)
    int total_params = fm->GetNParam();
    std::vector<bool> param_fixed(total_params, false);
    fm->SetParameterConstraints(param_fixed);

    return fm;
}

//' Evaluate log-likelihood for given parameters
//'
//' @param fm_ptr External pointer to FactorModel object
//' @param params Vector of parameters
//' @param compute_gradient Whether to compute gradient (default FALSE)
//' @param compute_hessian Whether to compute Hessian (default FALSE)
//' @return List with log-likelihood, gradient (if requested), and Hessian (if requested)
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
        Named("n_param_total") = fm->GetNParam(),
        Named("n_param_free") = fm->GetNParamFree()
    );
}
