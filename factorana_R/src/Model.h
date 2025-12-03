#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <string>
#include <RcppEigen.h>

// Model types
enum class ModelType {
    LINEAR = 1,
    PROBIT = 2,
    LOGIT = 3,
    OPROBIT = 4
};

class Model {
private:
    ModelType modtype;        // Model type: 1=linear, 2=probit, 3=logit, 4=ordered probit
    int outcome_idx;          // Index of outcome variable in data
    int missing_idx;          // Index of missing indicator (-1 if none)
    int nregressors;          // Number of regressors
    std::vector<int> regressors; // Indices of regressor variables
    std::vector<double> facnorm;  // Factor loading normalizations (-9999 = free, other = fixed value)
    int numfac;               // Number of factors
    int numtyp;               // Number of types (for type-specific loadings)
    int numchoice;            // Number of choices (for discrete choice models)
    int numrank;              // Number of rankings (for ranked choice models)
    bool ignore;              // If true, skip this model in likelihood
    bool all_params_fixed;    // If true, skip gradient/Hessian computation (use flag=1)

public:
    // Constructor
    Model(ModelType type, int outcome, int missing,
          const std::vector<int>& regs, int nfac, int ntyp = 0,
          const std::vector<double>& fnorm = std::vector<double>(),
          int nchoice = 2, int nrank = 1, bool params_fixed = false);

    // Main evaluation function
    // Computes likelihood, gradient, and Hessian for a single observation
    //
    // Parameters:
    //   iobs_offset - Offset in data vector for this observation
    //   data - Flattened data vector (all observations, all variables)
    //   param - Model-specific parameters (betas, alphas, sigma/thresholds)
    //   firstpar - Index of first parameter for this model in full parameter vector
    //   fac - Factor values for this observation/integration point
    //   modEval - Output: [likelihood, gradient components]
    //   hess - Output: Hessian matrix (upper triangle)
    //   flag - 1=likelihood only, 2=+gradient, 3=+Hessian
    //   type_intercept - Type-specific intercept to add to linear predictor (default 0)
    void Eval(int iobs_offset, const std::vector<double>& data,
              const std::vector<double>& param, int firstpar,
              const std::vector<double>& fac,
              std::vector<double>& modEval,
              std::vector<double>& hess,
              int flag,
              double type_intercept = 0.0);

    // Accessors
    ModelType GetType() const { return modtype; }
    int GetNumFac() const { return numfac; }
    int GetNumReg() const { return nregressors; }
    int GetNumChoice() const { return numchoice; }
    int GetOutcome() const { return outcome_idx; }
    int GetMissing() const { return missing_idx; }
    bool GetIgnore() const { return ignore; }
    void SetIgnore(bool ig) { ignore = ig; }
    bool GetAllParamsFixed() const { return all_params_fixed; }
    void SetAllParamsFixed(bool fixed) { all_params_fixed = fixed; }

private:
    // Helper functions for each model type
    void EvalLinear(double Z, double sigma, const std::vector<double>& fac,
                    const std::vector<double>& param, int firstpar,
                    std::vector<double>& modEval, std::vector<double>& hess,
                    int flag, const std::vector<double>& data, int iobs_offset);

    void EvalProbit(double expres, double obsSign, const std::vector<double>& fac,
                    const std::vector<double>& param, int firstpar,
                    std::vector<double>& modEval, std::vector<double>& hess,
                    int flag, const std::vector<double>& data, int iobs_offset);

    void EvalLogit(const std::vector<double>& expres, double outcome,
                   const std::vector<double>& fac,
                   const std::vector<double>& param, int firstpar,
                   std::vector<double>& modEval, std::vector<double>& hess,
                   int flag, const std::vector<double>& data, int iobs_offset);

    void EvalOprobit(double expres, int outcome_value,
                     const std::vector<double>& fac,
                     const std::vector<double>& param, int firstpar,
                     std::vector<double>& modEval, std::vector<double>& hess,
                     int flag, const std::vector<double>& data, int iobs_offset);
};

#endif // MODEL_H
