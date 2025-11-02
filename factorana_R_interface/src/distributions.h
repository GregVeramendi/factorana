#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <Rcpp.h>
#include <cmath>

// Normal distribution functions
// Uses R's dnorm and pnorm for accuracy and performance

inline double normal_pdf(double x, double mean = 0.0, double sd = 1.0) {
    return R::dnorm(x, mean, sd, 0);  // 0 = return density, not log-density
}

inline double normal_cdf(double x, double mean = 0.0, double sd = 1.0) {
    return R::pnorm(x, mean, sd, 1, 0);  // 1 = lower tail, 0 = not log
}

inline double normal_log_pdf(double x, double mean = 0.0, double sd = 1.0) {
    return R::dnorm(x, mean, sd, 1);  // 1 = return log-density
}

// Inverse normal CDF (quantile function)
inline double normal_quantile(double p, double mean = 0.0, double sd = 1.0) {
    return R::qnorm(p, mean, sd, 1, 0);  // 1 = lower tail, 0 = not log
}

#endif // DISTRIBUTIONS_H
