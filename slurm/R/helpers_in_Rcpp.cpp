#include <Rcpp.h>
#include <string>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <RcppNumerical.h>
using namespace Numer;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

// Recreation of R's rep when two vectors are provided
//
// @param x NumericVector The values to be repeated.
// @param times IntegerVector The number of times to repeat each value.
//   Vector must be of the same length as input `x`.
//
// @returns NumericVector of length sum(times)
// @keywords internal
// [[Rcpp::plugins(cpp11)]]
NumericVector repVec(const NumericVector& x, const IntegerVector& times) {
    std::size_t n = times.size();
    if (n != 1 && n != x.size())
        stop("Invalid 'times' value");
    std::size_t n_out = std::accumulate(times.begin(), times.end(), 0);
    NumericVector res = no_init(n_out);
    auto begin = res.begin();
    for (std::size_t i = 0, ind = 0; i < n; ind += times[i], ++i) {
        auto start = begin + ind;
        auto end = start + times[i];
        std::fill(start, end, x[i]);
    }
    return res;
}

//
// Weibull ////////


// Given rate and shape values, calculate the corresponding scale value
//
// @param rate double The rate parameter of a Weibull distribution. Must be > 0.
// @param shape double The shape parameter of a Weibull distribution. Must be > 0.
//
// @returns double The scale parameter of the Weibull distribution. Will be > 0.
// @keywords internal
double rate2scale(double rate, double shape) {
    NumericVector temp(1, rate);
    return Rcpp::pow(temp, -1.0 / shape)[0];
}

// Given a parameter list with updated rates, calculate new scales
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The input `theta` with updated elements `$scale_H`
//   and `$scale_P`.
// @keywords internal
List update_scales(List theta) {
    theta["scale_H"] = rate2scale(theta["rate_H"], theta["shape_H"]);
    theta["scale_P"] = rate2scale(theta["rate_P"], theta["shape_P"]);
    return theta;
}

// Calculate the difference between a Weibull CDF evaluated at two locations.
//
// Note: Differences are taken between element i of input `b` and element i of
//   input `a`.
//
// @param a NumericVector A vector of one or more quantiles.
// @param b NumericVector A vector of one or more quantiles. Must be the same
//   length as input `a`.
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
//
// @returns NumericVector The ith element is the difference in the Weibull
//   distribution at the ith quantiles of inputs `a` and `b`.
// @keywords internal
NumericVector pweibull_ab(NumericVector a, 
                          NumericVector b, 
                          double shape, 
                          double scale) {
    
    return pweibull(b, shape, scale, true, false) -
        pweibull(a, shape, scale, true, false);
}

// Random sample of x for which the CDF of a Weibull distribution falls between
//   CDF(a) and CDF(b) -- there has to be a more intuitive definition than this!
//
// Note: Upper and lower limits are taken in pairs from inputs `a` and `b`; 
//   i.e., element i of `a` is the lower bound, element i of `b` is the upper
//   bound.
//
// @param a NumericVector A vector of one or more quantiles.
// @param b NumericVector A vector of one or more quantiles. Must be the same
//   length as input `a`.
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
//
// @returns NumericVector The sampled values.
// @keywords internal
NumericVector rweibull_trunc(NumericVector a, 
                             NumericVector b, 
                             double shape, 
                             double scale) {
    
    NumericVector low = pweibull(a, shape, scale);
    NumericVector high = pweibull(b, shape, scale);
    
    NumericVector u(low.size());
    // populate each element of `u`, u_i, with a random draw from U(low_i, high_i)
    std::transform(low.begin(), 
                   low.end(), 
                   high.begin(),
                   u.begin(), [=](double low, double high){ return runif(1, low, high)[0]; }); 
    
    return qweibull(u, shape, scale);
}

// Need an intuitive title here
//
// @param x NumericVector A vector of one or more quantiles at which the
//   density is evaluated.
// @param a NumericVector A vector of one or more quantiles. Must be the
//   same length as input `x`.
// @param b NumericVector A vector of one or more quantiles. Must be the
//   same length as input `x`.
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
// @param uselog bool If TRUE, probabilities p are returned as log(p).
//
// @returns NumericVector The density of the Weibull 
//   distribution constrained to lie between the paired boundaries. (??)
// @keywords internal
NumericVector dweibull_trunc(NumericVector x, 
                             NumericVector a, 
                             NumericVector b, 
                             double shape, 
                             double scale, 
                             bool uselog) {
    if (uselog) {
        return dweibull(x, shape, scale, true) - log(pweibull_ab(a, b, shape, scale));
    } else {
        return dweibull(x, shape, scale, false) / pweibull_ab(a, b, shape, scale);
    }
}

//
// theta ////////

// Store a new rate_H value and update its corresponding scale_H.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param rate_H double The new rate value to be stored.
//
// @returns List The input `theta` with updated elements `$rate_H` and
//  `$scale_H`.
// @keywords internal
List add_rate_H(List theta, double rate_H) {
    theta["rate_H"] = rate_H;
    return update_scales(theta);
}

// Store a new rate_P value and update its corresponding scale_P.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param rate_P double The new rate value to be stored.
//
// @returns List The input `theta` with updated elements `$rate_P` and
//  `$scale_P`.
// @keywords internal
List add_rate_P(List theta, double rate_P) {
    theta["rate_P"] = rate_P;
    return update_scales(theta);
}

// Store a new beta value.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param beta double The new value to be stored in element `$beta`.
//
// @returns List The input `theta` with updated element `$beta`.
// @keywords internal
List add_beta(List theta, double beta) {
    theta["beta"] = beta;
    return theta;
}

// Store a new psi value.
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param psi double The new value to be stored in element `$psi`.
//
// @returns List The input `theta` with updated element `$psi`.
// @keywords internal
List add_psi(List theta, double psi) {
    theta["psi"] = psi;
    return theta;
}

//
// Gibbs theta ////////

// Update beta value
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param prior List The distribution parameters of the prior.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param n_screen_positive integer The total number of positive screens.
//
// @returns List The input `theta` with updated element `$beta`.
List gibbs_beta(List data_objects, 
                List prior, 
                List tau_HPs,
                List theta, 
                int n_screen_positive) { // beta prior and binomial likelihood conjugacy
    
    IntegerVector lengths;
    NumericVector ages, tau, tmp_tau;
    LogicalVector cmp;
    
    List age_screen;
    List obj;
    
    int n_screen = 0;
    for (int i = 0; i < data_objects.size(); ++i) {
        obj = data_objects[i];
        
        age_screen = obj["ages_screen"];
        ages = age_screen["values"];
        lengths = age_screen["lengths"];
        
        tau = tau_HPs[i];
        
        tmp_tau = repVec(tau, lengths);
        cmp = ages > tmp_tau;
        n_screen += sum(as<IntegerVector>(cmp));
    }
    
    double a_beta = prior["a_beta"];
    double b_beta = prior["b_beta"];
    double a_n = a_beta + n_screen_positive;
    double b_n = b_beta + n_screen - n_screen_positive;
    double beta_new = rbeta(1, a_n, b_n)[0];
    
    return add_beta(theta, beta_new);
}

//
// Likelihood ////////

// Derivative of the log-likelihood w.r.t the sojourn time for all individuals
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The current estimated tau for all
//   individuals of the censor_type under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector The derivative of the log-likelihood wrt the
//   sojourn time for each individual.
// @keywords internal
NumericVector dloglik_sojourn_H_obj(List data_object, 
                                    List theta,
                                    NumericVector tau_HP,
                                    double t0) {
    
    NumericVector censor_time = data_object["censor_time"];
    int n = data_object["n"];
    
    NumericVector result(n);
    
    LogicalVector isInfinite = is_infinite(tau_HP);
    
    NumericVector vec1 = censor_time[isInfinite];
    NumericVector result_infinite = pweibull(vec1 - t0,
                                             theta["shape_H"],
                                                  theta["scale_H"],
                                                       false, true);
    result[isInfinite] = result_infinite;
    
    vec1 = tau_HP[!isInfinite];
    NumericVector result_finite = dweibull(vec1 - t0,
                                           theta["shape_H"],
                                                theta["scale_H"],
                                                     true);
    
    result[!isInfinite] = result_finite;
    return result;
}

// Sum of the derivative of the log-likelihood w.r.t the sojourn time for over
//   all individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param tau_HP List The current estimates for tau broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   tau estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns double The sum of the derivative of the log-likelihood wrt
//   the sojourn time over all individuals.
// @keywords internal
double dloglik_sojourn_H_sum(List data_objects, 
                             List tau_HPs,
                             List theta, double t0) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dloglik_sojourn_H_obj(data_objects[i], theta, tau_HPs[i], t0));        
    }
    return result;
}

// Derivative of the log-likelihood w.r.t the sojourn time for all individuals
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The current estimated tau for the
//   censor type under consideration.
// @param indolent IntegerVector 0/1 indicating if is not/is indolent.
//
// @returns NumericVector The derivative of the log-likelihood wrt the
//   sojourn time for each individual.
// @keywords internal
NumericVector dloglik_sojourn_P_obj(List data_object, 
                                    List theta, 
                                    NumericVector tau_HP,
                                    IntegerVector indolent) {
    
    NumericVector censor_time = data_object["censor_time"];
    int censor_type = data_object["censor_type"];
    int n = data_object["n"];
    
    NumericVector result(n, 0.0);
    
    LogicalVector subset = is_finite(tau_HP) & indolent == 0;
    NumericVector vec1 = censor_time[subset] - tau_HP[subset];
    
    if (censor_type == 3) {
        NumericVector res = dweibull(vec1, theta["shape_P"], theta["scale_P"], true);
        result[subset] = res;
    } else {
        NumericVector res = pweibull(vec1, theta["shape_P"], theta["scale_P"], false, true);
        result[subset] = res;
    }
    
    return result;
}

// Sum of the derivative of the log-likelihood w.r.t the sojourn time over all 
//   individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for 0/1 indolent.
//   The 1st contains the indolent estimates for censor_type = "screen" cases; 
//   the 2nd the data for censor_type = "censored" cases; and the 3rd the data
//   for censor_type = "clinical" cases.
// @param tau_HP List The current estimates for tau broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   tau estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns double The sum of the derivative of the log-likelihood wrt
//   the sojourn time over all individuals.
// @keywords internal
double dloglik_sojourn_P_sum(List data_objects, 
                             List indolents, 
                             List tau_HPs, 
                             List theta) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dloglik_sojourn_P_obj(data_objects[i],
                                            theta, 
                                            tau_HPs[i],
                                                   indolents[i]));
    }
    return result;
}

// Derivative of the log-likelihood w.r.t. ... for a single censor type
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param indolent IntegerVector 0/1 indicating is not/is indolent.
//
// @returns NumericVector The derivative of the log-likelihood wrt 
//   indolence for each individual.
// @keywords internal
NumericVector dloglik_indolent_obj(List theta, IntegerVector indolent) {
    double psi = theta["psi"];
    
    int n = indolent.size();
    
    NumericVector result(n);
    result[indolent == 1] = log(psi);
    result[indolent == 0] = log(1.0 - psi);
    return result;    
}

// Derivative of the log-likelihood w.r.t. ... for all individuals
//
// @param indolents List The current estimates for 0/1 indolent.
//   The 1st contains the indolent estimates for censor_type = "screen" cases; 
//   the 2nd the data for censor_type = "censored" cases; and the 3rd the data
//   for censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The sum of the derivative of the log-likelihood wrt
//   the sojourn time over all individuals in each censor type.
// @keywords internal
List dloglik_indolent_List(List indolents, List theta) {
    List result(indolents.size());
    for (int i = 0; i < indolents.size(); ++i) {
        result[i] = dloglik_indolent_obj(theta, indolents[i]);
    }
    return result;
}

// To Be Named
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The current estimates for tau for
//   the censor type under consideration
//
// @returns NumericVector The derivative of the log-likelihood
//   wrt ... for each individual with the censor type under consideration.
// @keywords internal
NumericVector dloglik_screens_obj(List data_object, List theta, NumericVector tau_HP) {
    
    List ages_screen = data_object["ages_screen"];
    NumericVector age_screen = ages_screen["values"];
    IntegerVector starts = ages_screen["starts"];
    IntegerVector ends = ages_screen["ends"];
    IntegerVector lengths = ages_screen["lengths"];
    
    IntegerVector n_screen_positive = data_object["n_screen_positive"];
    
    NumericVector long_taus = repVec(tau_HP, lengths);
    
    IntegerVector comp(age_screen.size(), 0);
    comp[age_screen > long_taus] = 1;
    IntegerVector tmp = cumsum(comp);
    IntegerVector tmp2 = comp[starts];
    IntegerVector n_screen_neg = tmp[ends] - tmp[starts] + tmp2 - n_screen_positive;
    
    // product of Bernoulli RVs
    double beta = theta["beta"];
    return as<NumericVector>(n_screen_positive) * std::log(beta) + 
        as<NumericVector>(n_screen_neg) * std::log(1.0 - beta);
}

// To Be Named
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The list has 3 elements, each pertaining to a single
//   censor type. Each element is a NumericVector containing the derivative of 
//   the log-likelihood wrt ... for all individuals with the specific
//   censor type. Element 1: screen; 2: censored; 3: clinical.
// @keywords internal
List dloglik_screens_List(List data_objects, List tau_HPs, List theta) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = dloglik_screens_obj(data_objects[i], theta, tau_HPs[i]);
    }
    return result;  
}


// Derivative of the log-likelihood wrt sojourn time and indolence for all
//   individuals
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The estimated tau for all individuals.
// @param indolent IntegerVector =1 if indolent.
//
// @returns NumericVector The sum of the derivative of the log-likelihood 
//   wrt sojourn time and wrt to indolence for each individual.
// @keywords internal
NumericVector dloglik_PI_obj(List data_object, 
                             List theta, 
                             NumericVector tau_HP,
                             IntegerVector indolent) {
    NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, theta,
                                                 tau_HP, indolent);
    NumericVector dlog_I = dloglik_indolent_obj(theta, indolent);
    return dlog_P + dlog_I;
}

// Sum of the derivative of the log-likelihood wrt sojourn time and indolence 
//   over all individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   indolent estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns double The sum of the derivative of the log-likelihood 
//   wrt sojourn time and wrt to indolence over all individuals.
// @keywords internal
double dloglik_PI_sum(List data_objects, List indolent, List tau_HPs, List theta) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dloglik_PI_obj(data_objects[i], theta, tau_HPs[i], indolent[i]));
    }
    return result;
}

// Class to facilitate numerical integration
// 
// @param shapeH, scaleH, shapeP, scaleP Each double The parameters
//   of the Weibull distributions.
// @param U double The upper boundary.
//
// Initialization sets shapeH, scaleH, shapeP, scaleP, and U for each integral.
// operator() defines the integral
class WeibPDF: public Func {
private: 
    double shapeH;
    double scaleH;
    double shapeP;
    double scaleP; 
    double U;
    
public: 
    WeibPDF(double shapeH_, double scaleH_, double shapeP_, double scaleP_, double U_) : 
    shapeH(shapeH_), scaleH(scaleH_), shapeP(shapeP_), scaleP(scaleP_), U(U_){}
    double operator()(const double& x) const {
        NumericVector local_vec(1);
        local_vec[0] = x;
        return dweibull(local_vec, shapeH, scaleH, false)[0] * 
            pweibull(U - local_vec, shapeP, scaleP, false, false)[0];
    } 
};

// Provided the lower and upper boundaries, integrate 
//   f(x; k, lambda) F(U - x; k, lambda)
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param L double The lower limit of integration.
// @param U double The upper limit of integration.
//
// @returns double The integral.
// @keywords internal
double compute_integral(List theta, double L, double U) {
    
    // Initialize instance of integral definition for the current upper bound
    WeibPDF f(theta["shape_H"], theta["scale_H"], theta["shape_P"], theta["scale_P"], U);
    
    // error estimate and error code returned from the integration routine.
    // only the error flag is of interest here.
    double err_est;
    int err_code;
    
    // integrate
    double res = integrate(f, L, U, err_est, err_code);
    
    // if error code is not 0, there was a problem. abort.
    if (err_code > 0) stop("Unable to perform integration");
    
    return res;
}

// Compute integral at each ...
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param t0 double The initial time.
//
// @returns An Rcpp::NumericVector of the integrals at each ...
// @keywords internal
NumericVector compute_cp_log(List theta, NumericVector AFS, double t0) {
    
    double L = 0.0;        // lower bound
    NumericVector U = AFS - t0; // upper bound
    
    NumericVector prob_onset_after = pweibull(U, 
                                              theta["shape_H"], theta["scale_H"], 
                                                                     false, false);
    
    NumericVector integral(U.size());
    for (int i = 0; i < U.size(); ++i) {    
        integral[i] = compute_integral(theta, L, U[i]);
    }
    double psi = theta["psi"];
    
    NumericVector cp = psi + (1.0 - psi) * (prob_onset_after + integral); // conditioning probability
    return log(cp);
}

// Sum of the integrals
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param n_AFS NumericVector ...
// @param t0 double The initial time.
//
// @returns double The sum of the integrals over all AFS.
// @keywords internal
double dloglik_cp(List theta, 
                  NumericVector AFS, 
                  NumericVector n_AFS,
                  double t0) {
    return -sum(compute_cp_log(theta, AFS, t0) * n_AFS);
}

// Sum of the derivative of the log-likelihood wrt psi
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   indolent estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param n_AFS NumericVector ...
// @param t0 double The initial time.
//
// @returns double The sum of the derivative wrt psi.
// @keywords internal
double dloglik_psi(List data_objects, 
                   List indolent,
                   List tau_HPs,
                   List theta, 
                   NumericVector AFS, 
                   NumericVector n_AFS, 
                   double t0) {
    double dlog_cp = dloglik_cp(theta, AFS, n_AFS, t0);
    double dlog_PI = dloglik_PI_sum(data_objects, indolent, tau_HPs, theta);
    return dlog_cp + dlog_PI;
}

// Derivative of the log-likelihood wrt tau for a single censor group
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The estimated tau for all individuals.
// @param indolent IntegerVector =1 if indolent.
// @param t0 double Initial time.
//
// @returns NumericVector The derivative for each case in the censor group
//   under consideration.
// @keywords internal
NumericVector dloglik_tau_obj(List data_object,
                              List theta, 
                              NumericVector tau_HP, 
                              IntegerVector indolent, 
                              double t0) {
    NumericVector dlog_H = dloglik_sojourn_H_obj(data_object, theta, tau_HP, t0);
    NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, theta, tau_HP, indolent);
    NumericVector dlog_S = dloglik_screens_obj(data_object, theta, tau_HP);
    return dlog_H + dlog_P + dlog_S;
}

// Derivative of the log-likelihood wrt tau for all censor groups
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   indolent estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List There are 3 elements. Each element holds the derivative
//   for a single censor group.
List dloglik_tau_List(List data_objects, 
                      List indolents, 
                      List tau_HPs, 
                      List theta, 
                      double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = dloglik_tau_obj(data_objects[i], theta, tau_HPs[i], 
                                    indolents[i], t0);
    }
    return result;
}

// Sum of the derivative of the log-likelihood wrt psi
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param n_AFS NumericVector ...
// @param t0 double The initial time.
//
// @returns double The sum of the derivative wrt rate_H.
// @keywords internal
double dloglik_rate_H(List data_objects, 
                      List tau_HPs, 
                      List theta, 
                      NumericVector AFS, 
                      NumericVector n_AFS, 
                      double t0) {
    double dlog_cp = dloglik_cp(theta, AFS, n_AFS, t0);
    double dlog_H = dloglik_sojourn_H_sum(data_objects, tau_HPs, theta, t0);
    return dlog_cp + dlog_H;
}


// Sum of the derivative of the log-likelihood wrt psi
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   indolent estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param n_AFS NumericVector ...
// @param t0 double The initial time.
//
// @returns double The sum of the derivative wrt rate_P
// @keywords internal
double dloglik_rate_P(List data_objects, 
                      List indolents, 
                      List tau_HPs,
                      List theta, 
                      NumericVector AFS, 
                      NumericVector n_AFS, 
                      double t0) {
    double dlog_cp = dloglik_cp(theta, AFS, n_AFS, t0);
    double dlog_P = dloglik_sojourn_P_sum(data_objects, indolents, tau_HPs, theta);
    return dlog_cp + dlog_P;
}

//
// M-H rate_H ////////

// Sample a new rate value
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilong_rate_H double A small value shift.
//
// @returns double The rate shifted by U(-eps, eps). Will never be negative.
// @keywords internal
double rprop_rate_H(List theta, double epsilon_rate_H) {
    double rate_H = theta["rate_H"];
    rate_H += runif(1, - epsilon_rate_H, epsilon_rate_H)[0];
    return std::fabs(rate_H);
}

// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   indolent estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param prior List The distribution parameters of the prior.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param n_AFS NumericVector ...
// @param epsilon_rate_H A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 2 elements. Theta - the accepted new theta; and
//   accept - TRUE if theta was updated; FALSE otherwise.
// @keywords internal
List MH_rate_H(List data_objects,
               List indolents, 
               List prior, 
               List tau_HPs,
               List theta_cur, 
               NumericVector AFS, 
               NumericVector n_AFS, 
               double epsilon_rate_H, 
               double t0) {
    
    // current value
    NumericVector rate_H_cur(1, theta_cur["rate_H"]); // for evaluating prior density
    NumericVector rate_H_new(1, rprop_rate_H(theta_cur, epsilon_rate_H)); // symmetric proposal
    
    List theta_new = clone(theta_cur);
    theta_new = add_rate_H(theta_new, rate_H_new[0]);  // for dloglik_rate_H()
    
    double dlog_lik_cur = dloglik_rate_H(data_objects, tau_HPs, theta_cur, 
                                         AFS, n_AFS, t0);
    double dlog_lik_new = dloglik_rate_H(data_objects, tau_HPs, theta_new, 
                                         AFS, n_AFS, t0);
    
    double prior_rate = prior["rate_H"];
    double dlog_prior_new = dgamma(rate_H_new, prior["shape_H"], 
                                   1.0 / prior_rate, true)[0];
    double dlog_prior_cur = dgamma(rate_H_cur, prior["shape_H"], 
                                   1.0 / prior_rate, true)[0];
    
    double MH_logratio = (dlog_lik_new + dlog_prior_new) - (dlog_lik_cur + dlog_prior_cur);
    
    if (runif(1)[0] < std::exp(MH_logratio)) {
        return List::create(Named("theta") = theta_new, Named("accept") = true);
    } else {
        return List::create(Named("theta") = theta_cur, Named("accept") = false);
    }
}

//
// M-H rate_P ////////

// Sample a new rate value
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilong_rate_H double A small value shift.
//
// @returns double The rate shifted by U(-eps, eps). Will never be negative.
// @keywords internal
double rprop_rate_P(List theta, double epsilon_rate_P) {
    double rate_P = theta["rate_P"];
    rate_P = runif(1, rate_P - epsilon_rate_P, rate_P + epsilon_rate_P)[0];
    return std::fabs(rate_P);
}

// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements. The 1st contains the
//   indolent estimates for censor_type = "screen" cases; the 2nd the data for
//   censor_type = "censored" cases; and the 3rd the data for
//   censor_type = "clinical" cases.
// @param prior List The distribution parameters of the prior.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector ...
// @param n_AFS NumericVector ...
// @param epsilon_rate_P A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 2 elements. Theta - the accepted new theta; and
//   accept - TRUE if theta was updated; FALSE otherwise.
// @keywords internal
List MH_rate_P(List data_objects,
               List indolents, 
               List prior, 
               List tau_HPs,
               List theta_cur, 
               NumericVector AFS, 
               NumericVector n_AFS, 
               double epsilon_rate_P, 
               double t0) {
    
    NumericVector rate_P_cur(1, theta_cur["rate_P"]);
    NumericVector rate_P_new(1, rprop_rate_P(theta_cur, epsilon_rate_P)); 
    
    List theta_new = clone(theta_cur);
    theta_new = add_rate_P(theta_new, rate_P_new[0]);  // for dloglik_rate_P()
    
    // M-H acceptance ratio
    double dlog_lik_cur = dloglik_rate_P(data_objects, indolents, tau_HPs, 
                                         theta_cur, AFS, n_AFS, t0);
    double dlog_lik_new = dloglik_rate_P(data_objects, indolents, tau_HPs,
                                         theta_new, AFS, n_AFS, t0);
    
    double prior_rate = prior["rate_P"];
    double dlog_prior_new = dgamma(rate_P_new, prior["shape_P"], 
                                   1.0 / prior_rate, true)[0];
    double dlog_prior_cur = dgamma(rate_P_cur, prior["shape_P"], 
                                   1.0 / prior_rate, true)[0];
    
    double MH_logratio = (dlog_lik_new + dlog_prior_new) - (dlog_lik_cur + dlog_prior_cur);
    
    if (runif(1)[0] < std::exp(MH_logratio)) {
        return List::create(Named("theta") = theta_new, Named("accept") = true);
    } else {
        return List::create(Named("theta") = theta_cur, Named("accept") = false);
    }
}


//
// M-H tau_HP ////////

// Compute the probability of tau for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @return List Elements include `values`, a NumericVector of all probabilities,
//   `starts`, an IntegerVector containing the first element of `values pertaining
//   to case i; `ends`, an IntegerVector containing the last element pertaining to
//   individual i; and `lengths` the number of element in `values` pertaining to
//   individual i.
List compute_prob_tau_obj(List data_object, List theta, double t0) {
    
    NumericVector censor_time = data_object["censor_time"];
    
    List endpoints = data_object["endpoints"];
    NumericVector ep = endpoints["values"];
    IntegerVector starts = endpoints["starts"];
    IntegerVector ends = endpoints["ends"];
    IntegerVector lengths = endpoints["lengths"];
    
    int censor_type = data_object["censor_type"];
    
    // The number of intervals is 1 less than the number of boundary points
    IntegerVector pt_lengths = lengths - 1;
    // The final index for each individual is 1 fewer than the lengths (0 indexing)
    IntegerVector pt_ends = cumsum(pt_lengths);
    pt_ends = pt_ends - 1;
    // The starts are obtained using the final index and the lengths
    IntegerVector pt_starts = pt_ends - pt_lengths + 1;
    
    // compute probability of interval based on the Weibull waiting time in H 
    // (screen, censored) or in P (clinical)
    NumericVector prob_interval;
    
    // creating a vector of indices to pull all but the last boundary value
    // from endpoints
    // In R this would be endpoints[-ends]
    IntegerVector intervals(sum(pt_lengths));
    IntegerVector idx;
    int st = 0;
    for (int i = 0; i < starts.size(); ++i) {
        idx = seq(starts[i], ends[i] - 1);
        intervals[seq(st, st + idx.size() - 1)] = idx;
        st += idx.size();
    }
    
    // paired points at which the distribution is evaluated 
    NumericVector vec1 = ep[intervals];
    NumericVector vec2 = ep[intervals + 1];
    
    if (censor_type != 3) {
        prob_interval = pweibull_ab(vec1 - t0, vec2 - t0, 
                                    theta["shape_H"], theta["scale_H"]);
        
    } else {
        // need to repeat censor time for each interval
        NumericVector ct = repVec(censor_time, pt_lengths);
        prob_interval = pweibull_ab(ct - vec2,
                                    ct - vec1,
                                    theta["shape_P"], theta["scale_P"]);
    }
    
    // we need to repeat (1-beta)^{K-1}:0
    IntegerVector powers(sum(pt_lengths));
    st = 0;
    for (int i = 0; i < pt_lengths.size(); ++i) {
        idx = seq(0, pt_lengths[i] - 1);
        powers[seq(st, st + idx.size() -1)] = rev(idx);
        st += idx.size();
    }
    
    double beta = theta["beta"];
    NumericVector tbeta(prob_interval.size(), 1.0 - beta);
    
    std::transform(tbeta.begin(),
                   tbeta.end(),
                   powers.begin(),
                   tbeta.begin(), [=] (double x, int p) {return pow(x, p); });
    
    // won't this just be removed by the normalization?
    if (censor_type == 1) tbeta = tbeta * beta;
    
    // combine vectors
    NumericVector prob_tau = prob_interval * tbeta;
    
    // normalize each individual's probabilities
    NumericVector tmp_norm = cumsum(prob_tau);
    NumericVector tmp_pt = prob_tau[pt_starts];
    NumericVector norm = tmp_norm[pt_ends] - tmp_norm[pt_starts] + tmp_pt;
    NumericVector norm_ext = repVec(norm, pt_lengths);
    
    prob_tau = prob_tau / norm_ext; // normalize probabilities
    
    return List::create(Named("values") = prob_tau,
                        Named("starts") = pt_starts,
                        Named("ends") = pt_ends,
                        Named("lengths") = pt_lengths);
}

// Compute probability of tau for all censoring types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is itself a List
//   containing the probabilities of tau details for the specific
//   censoring type.
// exported to initialize these values are start of MC
// [[Rcpp::export]]
List compute_prob_tau_List(List data_objects, List theta, double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = compute_prob_tau_obj(data_objects[i], theta, t0);
    }
    return result;
}

// Needs descriptive title
//
// @param data_object List The data for a single censoring type.
// @param prob_tau List The probabilities of tau for all cases in the
//   censoring type under consideration.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The current estimated tau for all
//   individuals of the censor_type under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector The new estimated tau for each case in the 
//   censoring type under consideration
NumericVector rprop_tau_HP_obj(List data_object, List prob_tau, List theta, 
                               double t0) {
    
    NumericVector censor_time = data_object["censor_time"];
    
    List endpoints = data_object["endpoints"];
    NumericVector ep = endpoints["values"];
    IntegerVector ep_starts = endpoints["starts"];
    IntegerVector ep_ends = endpoints["ends"];
    IntegerVector ep_lengths = endpoints["lengths"];
    
    int censor_type = data_object["censor_type"];
    
    NumericVector pt = prob_tau["values"];
    IntegerVector pt_starts = prob_tau["starts"];
    IntegerVector pt_ends = prob_tau["ends"];
    IntegerVector pt_lengths = prob_tau["lengths"];
    
    // sample each case's interval
    IntegerVector k_new(censor_time.size());
    for (int i = 0; i < pt_starts.size(); ++i) {
        int K = pt_lengths[i];
        NumericVector vec = pt[seq(pt_starts[i], pt_ends[i])];
        k_new[i] = sample(K, 1, false, vec)[0] - 1;
    }
    
    NumericVector tau_HP_new;
    
    NumericVector vec1 = ep[ep_starts + k_new];
    NumericVector vec2 = ep[ep_starts + k_new + 1];
    
    // sample tau_HP in chosen interval
    if (censor_type == 1) {
        NumericVector sojourn_H_new = rweibull_trunc(vec1 - t0, vec2 - t0,
                                                     theta["shape_H"], 
                                                          theta["scale_H"]);
        tau_HP_new = t0 + sojourn_H_new;
    } else if(censor_type == 2) {
        NumericVector sojourn_H_new = rweibull_trunc(vec1 - t0, vec2 - t0,
                                                     theta["shape_H"], 
                                                          theta["scale_H"]);
        tau_HP_new = t0 + sojourn_H_new;
        tau_HP_new[k_new == (pt_lengths - 1)] = R_PosInf;
    } else {
        NumericVector sojourn_P_new = rweibull_trunc(censor_time - vec2, 
                                                     censor_time - vec1,
                                                     theta["shape_P"], 
                                                          theta["scale_P"]);
        tau_HP_new = censor_time - sojourn_P_new;
    }
    return tau_HP_new;
}


// Needs descriptive title
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param prob_tau List The probabilities of tau broken down according to the
//   censor type. There are 3 elements. The 1st contains the probabilities of
//   tau for censor_type = "screen" cases; the 2nd the probabilities of tau for 
//   censor_type = "censored" cases; and the 3rd the probabilities of tau for 
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
// exported to initialize these values are start of MC
// [[Rcpp::export]]
List rprop_tau_HP_List(List data_objects, List prob_tau, List theta, double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = rprop_tau_HP_obj(data_objects[i], prob_tau[i], theta, t0);
    }
    return result;
}

// Needs descriptive title
//
// @param data_object List The data for a single censoring type.
// @param prob_tau List The probabilities of tau for all cases in the
//   censoring type under consideration.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param tau_HP NumericVector The current estimated tau for all
//   individuals of the censor_type under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector 
// @keywords internal
NumericVector dlog_prop_tau_HP_obj(List data_object,
                                   List prob_tau,
                                   List theta,
                                   NumericVector tau_HP, 
                                   double t0) {
    
    NumericVector censor_time = data_object["censor_time"];
    
    List endpoints = data_object["endpoints"];
    NumericVector ep = endpoints["values"];
    IntegerVector ep_starts = endpoints["starts"];
    IntegerVector ep_ends = endpoints["ends"];
    IntegerVector ep_lengths = endpoints["lengths"];
    
    int censor_type = data_object["censor_type"];
    
    NumericVector pt = prob_tau["values"];
    IntegerVector pt_starts = prob_tau["starts"];
    IntegerVector pt_ends = prob_tau["ends"];
    IntegerVector pt_lengths = prob_tau["lengths"];
    
    // contribution of k_new
    NumericVector tmp_tau = repVec(tau_HP, ep_lengths);
    LogicalVector cmp = ep < tmp_tau;
    IntegerVector cmp_csum = cumsum(as<IntegerVector>(cmp));
    LogicalVector cmp_starts = cmp[ep_starts];
    IntegerVector icmp_starts = as<IntegerVector>(cmp_starts);
    
    IntegerVector k_new = cmp_csum[ep_ends] - cmp_csum[ep_starts] + icmp_starts - 1;
    
    NumericVector dlog_k = log(pt[pt_starts + k_new]);
    
    NumericVector vec1 = ep[ep_starts + k_new];
    NumericVector vec2 = ep[ep_starts + k_new + 1];
    
    // contribution of tau_HP_i
    NumericVector dlog_tau;
    if (censor_type == 1) {
        dlog_tau = dweibull_trunc(tau_HP - t0, 
                                  vec1 - t0, 
                                  vec2 - t0,
                                  theta["shape_H"], theta["scale_H"], true);
    } else if (censor_type == 2) {
        dlog_tau = dweibull_trunc(tau_HP - t0, 
                                  vec1 - t0, 
                                  vec2 - t0,
                                  theta["shape_H"], theta["scale_H"], true);
        dlog_tau[(k_new + 1) == pt_lengths] = 0.0;
    } else {
        dlog_tau = dweibull_trunc(censor_time - tau_HP, 
                                  censor_time - vec2, 
                                  censor_time - vec1,
                                  theta["shape_P"], theta["scale_P"], true);
    }
    return dlog_k + dlog_tau;
}

// Needs descriptive title
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements. The 1st contains the data for censor_type = "screen" 
//   cases; the 2nd the data for censor_type = "censored" cases; and the 3rd 
//   the data for censor_type = "clinical" cases.
// @param prob_tau List The probabilities of tau broken down according to the
//   censor type. There are 3 elements. The 1st contains the probabilities of
//   tau for censor_type = "screen" cases; the 2nd the probabilities of tau for 
//   censor_type = "censored" cases; and the 3rd the probabilities of tau for 
//   censor_type = "clinical" cases.
// @param tau_HPs List The current estimates for tau broken down according to
//   the censor type. There are 3 elements. The 1st contains the tau estimates
//   for censor_type = "screen" cases; the 2nd the tau estimates for
//   censor_type = "censored" cases; and the 3rd the tau estimates for
//   censor_type = "clinical" cases.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
List dlog_prop_tau_HP_List(List data_objects, List prob_taus, List tau_HPs, 
                           List theta, double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = dlog_prop_tau_HP_obj(data_objects[i], prob_taus[i], 
                                         theta, tau_HPs[i], t0);
    }
    return result;
}

//
// indolent ////////

NumericVector compute_prob_indolent_obj(List data_object, List theta,
                                        NumericVector tau_HP) {
    
    IntegerVector ind(tau_HP.size());
    
    std::fill(ind.begin(), ind.end(), 0);
    NumericVector L_0 = dloglik_PI_obj(data_object, theta, tau_HP, ind);
    L_0 = exp(L_0);
    
    std::fill(ind.begin(), ind.end(), 1);
    NumericVector L_1 = dloglik_PI_obj(data_object, theta, tau_HP, ind);
    L_1 = exp(L_1);
    
    return L_1 / (L_0 + L_1);
}

// [[Rcpp::export]]
List compute_prob_indolent_List(List data_objects, List tau_HPs, List theta) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = compute_prob_indolent_obj(data_objects[i], theta, tau_HPs[i]);
    }
    return result;
}

IntegerVector rprop_indolent_obj(List data_object, NumericVector prob_indolent) {
    
    int censor_type = data_object["censor_type"];
    
    int n = prob_indolent.size();
    IntegerVector indolent(n);
    
    if (censor_type == 3) {
        indolent.fill(0);
    } else {
        std::transform(prob_indolent.begin(), 
                       prob_indolent.end(), 
                       indolent.begin(), [=](double p){ return rbinom(1, 1, p)[0]; }); 
    }
    
    return indolent;
    
}

// [[Rcpp::export]]
List rprop_indolent_List(List data_objects, List prob_indolents) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = rprop_indolent_obj(data_objects[i], prob_indolents[i]);
    }
    return result;
}


NumericVector dlog_prop_indolent_obj(List data_object, 
                                     NumericVector prob_indolent,
                                     IntegerVector indolent) {
    
    int censor_type = data_object["censor_type"];
    
    NumericVector result = log(1.0 - prob_indolent);
    if (censor_type == 3) {
        result.fill(0.0);
    } else {
        LogicalVector ind1 = indolent == 1;
        NumericVector prob1 = log(prob_indolent[ind1]);
        result[ind1] = prob1;
    }
    
    return result;
}

double dlog_prop_indolent_sum(List data_objects, List indolents, List prob_indolents) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dlog_prop_indolent_obj(data_objects[i], prob_indolents[i], indolents[i]));
    }
    return result;
}

//
// psi ////////
double rprop_psi(List theta, double epsilon_psi) {
    
    double psi = theta["psi"];
    
    // uniform random walk
    double psi_prop = runif(1, psi - epsilon_psi, psi + epsilon_psi)[0];
    
    // reflection on lower bound 0 and upper bound 1 
    if (psi_prop < 0.0) {
        psi_prop = 0.0 + (0.0 - psi_prop);
    } else if (psi_prop > 1.0) {
        psi_prop = 1.0 - (psi_prop - 1.0);
    } 
    
    return psi_prop;
}

//
// M-H tau ////////
List MH_tau_obj(List data_object, List theta, NumericVector tau_HP, 
                IntegerVector indolent, double t0) {
    
    // propose new latent
    List prob_tau = compute_prob_tau_obj(data_object, theta, t0);
    
    NumericVector tau_HP_new = rprop_tau_HP_obj(data_object, prob_tau, theta, t0);
    
    // M-H acceptance ratio
    NumericVector dlog_prop_cur = dlog_prop_tau_HP_obj(data_object,
                                                       prob_tau,
                                                       theta,
                                                       tau_HP, 
                                                       t0);
    
    NumericVector dlog_prop_new = dlog_prop_tau_HP_obj(data_object,
                                                       prob_tau,
                                                       theta,
                                                       tau_HP_new,
                                                       t0);
    
    NumericVector dlog_lik_cur = dloglik_tau_obj(data_object,
                                                 theta,
                                                 tau_HP, 
                                                 indolent, 
                                                 t0);
    
    NumericVector dlog_lik_new = dloglik_tau_obj(data_object,
                                                 theta,
                                                 tau_HP_new, 
                                                 indolent, 
                                                 t0);
    
    NumericVector MH_logratio = dlog_lik_new - dlog_lik_cur + dlog_prop_cur - dlog_prop_new;
    
    LogicalVector test = runif(MH_logratio.size()) < exp(MH_logratio);
    
    NumericVector accepted_tau = tau_HP;
    accepted_tau[test] = tau_HP_new[test];
    
    return List::create(Named("tau_HP") = accepted_tau,
                        Named("accept") = test);
    
}

List MH_tau_List(List data_objects, List indolents, List tau_HPs, List theta, double t0) {
    
    List result_tau(data_objects.size());
    List result_accept(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        List res = MH_tau_obj(data_objects[i], theta, tau_HPs[i], indolents[i], t0);
        result_tau[i] = res["tau_HP"];
        result_accept[i] = res["accept"];
    }
    
    return List::create(Named("tau_HPs") = result_tau,
                        Named("accept") = result_accept);  
}


//
// M-H (psi, indolent) ////////

List rprop_dlog_indolent_obj(List data_object, List theta, NumericVector tau_HP) {
    
    NumericVector prob_indolent_new = compute_prob_indolent_obj(data_object,
                                                                theta,
                                                                tau_HP);
    
    IntegerVector indolent_new = rprop_indolent_obj(data_object, prob_indolent_new);
    
    double dlog_prop_new = sum(dlog_prop_indolent_obj(data_object, 
                                                      prob_indolent_new, 
                                                      indolent_new));
    
    return List::create(Named("indolent") = indolent_new,
                        Named("dlog_prop") = dlog_prop_new);
}

List rprop_dlog_indolent_List(List data_objects, List tau_HPs, List theta) {
    List result_indolent(data_objects.size());
    double dlog_prop = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        List res = rprop_dlog_indolent_obj(data_objects[i], theta, tau_HPs[i]);
        result_indolent[i] = res["indolent"];
        double dlog = res["dlog_prop"];
        dlog_prop += dlog;
    }
    
    return List::create(Named("indolent_new") = result_indolent,
                        Named("dlog_prop_indolent_new") = dlog_prop);  
}

// [[Rcpp::export]]
List MH_psi_indolent(List data_objects,
                     List indolents, 
                     List prior, 
                     List tau_HPs,
                     List theta_cur, 
                     NumericVector AFS, 
                     NumericVector n_AFS,
                     double epsilon_psi, 
                     double t0) {
    
    // propose psi
    NumericVector psi_new(1);
    psi_new[0] = rprop_psi(theta_cur, epsilon_psi); // symmetric proposal
    List theta_new = clone(theta_cur);
    theta_new = add_psi(theta_new, psi_new[0]);
    
    // prior density
    NumericVector cur_psi(1);
    cur_psi[0] = theta_cur["psi"];
    double a_psi = prior["a_psi"];
    double b_psi = prior["b_psi"];
    double dlog_prior_cur = dbeta(cur_psi, a_psi, b_psi, true)[0];
    double dlog_prior_new = dbeta(psi_new, a_psi, b_psi, true)[0];
    
    // propose indolent and compute log proposal density
    List out = rprop_dlog_indolent_List(data_objects, tau_HPs, theta_new);
    List indolent_new = out["indolent_new"];
    
    // proposal density
    double dlog_prop_indolent_new = out["dlog_prop_indolent_new"];
    List prob_indolent_cur = compute_prob_indolent_List(data_objects, tau_HPs, theta_cur);
    double dlog_prop_indolent_cur = dlog_prop_indolent_sum(data_objects,
                                                           indolents, 
                                                           prob_indolent_cur);
    
    // log likelihood
    double dlog_lik_cur = dloglik_psi(data_objects, indolents, tau_HPs, theta_cur, AFS, n_AFS, t0);
    double dlog_lik_new = dloglik_psi(data_objects, indolent_new, tau_HPs, theta_new, AFS, n_AFS, t0);
    
    // M-H acceptance ratio
    double MH_logratio = (dlog_lik_new + dlog_prior_new - dlog_prop_indolent_new) - 
        (dlog_lik_cur + dlog_prior_cur - dlog_prop_indolent_cur);
    if (runif(1)[0] < std::exp(MH_logratio)) {
        return List::create(Named("indolents") = indolent_new, 
                            Named("theta") = theta_new, 
                            Named("accept") = true);
    } else {
        return List::create(Named("indolents") = indolents, 
                            Named("theta") = theta_cur, 
                            Named("accept") = false);
    }
    
}


//
// MCMC ////////
// [[Rcpp::export]]
List MCMC_cpp(List data_objects, List indolents, List prior, List tau_HPs, List theta,
              double epsilon_rate_H, double epsilon_rate_P, double epsilon_psi,
              double t0, int M, int thin, int M_thin, int n_obs,
              int n_screen_positive_total,
              NumericVector AFS, NumericVector n_AFS) {
    
    NumericVector kRATE_H(M_thin), kRATE_P(M_thin), kBETA(M_thin), kPSI(M_thin);
    LogicalVector kACCEPT_PSI(M_thin), kACCEPT_RATE_H(M_thin), kACCEPT_RATE_P(M_thin);
    NumericMatrix kTAU_HP(M_thin, n_obs);
    IntegerMatrix kINDOLENT(M_thin, n_obs);
    LogicalMatrix kACCEPT_LATENT(M_thin, n_obs);
    int ikept = -1;
    
    for (int m = 0; m < M; ++m) {
        // M-H rate_H
        List out_rate_H = MH_rate_H(data_objects, indolents, prior, tau_HPs, theta, 
                                    AFS, n_AFS, epsilon_rate_H, t0);
        theta = out_rate_H["theta"];
        bool accept_rate_H = out_rate_H["accept"];
        
        // M-H rate_P
        List out_rate_P = MH_rate_P(data_objects, indolents, prior, tau_HPs, theta, 
                                    AFS, n_AFS, epsilon_rate_P, t0);
        theta = out_rate_P["theta"];
        bool accept_rate_P = out_rate_P["accept"];
        
        // Gibbs beta
        theta = gibbs_beta(data_objects, prior, tau_HPs, theta, n_screen_positive_total);
        
        // update (psi, indolent)
        List out_psi_indolent = MH_psi_indolent(data_objects, indolents, prior, tau_HPs, theta,
                                                AFS, n_AFS, epsilon_psi, t0);
        indolents = out_psi_indolent["indolents"];
        theta = out_psi_indolent["theta"];
        bool accept_psi = out_psi_indolent["accept"];
        
        // update tau_HP
        List out_tau  = MH_tau_List(data_objects, indolents, tau_HPs, theta, t0);
        tau_HPs = out_tau["tau_HPs"];
        List accept_latent = out_tau["accept"];
        // save output
        if (m % thin == 0) {
            NumericVector tau_screen = tau_HPs[0];
            IntegerVector indolent_screen = indolents[0];
            LogicalVector accept_latent_screen = accept_latent[0];
            
            NumericVector tau_censored = tau_HPs[1];
            IntegerVector indolent_censored = indolents[1];
            LogicalVector accept_latent_censored = accept_latent[1];
            
            NumericVector tau_clinical = tau_HPs[2];
            IntegerVector indolent_clinical = indolents[2];
            LogicalVector accept_latent_clinical = accept_latent[2];
            
            NumericVector tau_HP(tau_screen.size() + 
                tau_censored.size() + tau_clinical.size());
            IntegerVector indolent(indolent_screen.size() + 
                indolent_censored.size() + indolent_clinical.size());
            LogicalVector accept_latent(accept_latent_screen.size() + 
                accept_latent_censored.size() + accept_latent_clinical.size());
            
            int n_screen = tau_screen.size();
            int n_censored = tau_censored.size();
            int n_clinical = tau_clinical.size();
            
            tau_HP[seq(0, n_screen-1)] = tau_screen;
            tau_HP[seq(n_screen, n_screen + n_censored - 1)] = tau_censored;
            tau_HP[seq(n_screen+ n_censored, n_screen + n_censored + n_clinical - 1)] = tau_clinical;
            
            indolent[seq(0, n_screen-1)] = indolent_screen;
            indolent[seq(n_screen, n_screen + n_censored - 1)] = indolent_censored;
            indolent[seq(n_screen+ n_censored, n_screen + n_censored + n_clinical - 1)] = indolent_clinical;
            
            accept_latent[seq(0, n_screen-1)] = accept_latent_screen;
            accept_latent[seq(n_screen, n_screen + n_censored - 1)] = accept_latent_censored;
            accept_latent[seq(n_screen+ n_censored, n_screen + n_censored + n_clinical - 1)] = accept_latent_clinical;
            
            ikept += 1;
            kRATE_H  [ikept] = theta["rate_H"];
            kRATE_P  [ikept] = theta["rate_P"];
            kPSI     [ikept] = theta["psi"];
            kBETA    [ikept] = theta["beta"];
            kTAU_HP  (ikept, _ ) = tau_HP;
            kINDOLENT(ikept, _ ) = indolent;
            kACCEPT_LATENT(ikept, _ ) = accept_latent;
            kACCEPT_PSI   [ikept] = accept_psi;
            kACCEPT_RATE_H[ikept] = accept_rate_H;
            kACCEPT_RATE_P[ikept] = accept_rate_P;
        }
    } // end runtime
    
    // output
    List kTHETA = List::create(Named("RATE_H") = kRATE_H,
                               Named("RATE_P") = kRATE_P,
                               Named("BETA") = kBETA,
                               Named("PSI") = kPSI);
    List epsilon = List::create(Named("epsilon_rate_H") = epsilon_rate_H,
                                Named("epsilon_rate_P") = epsilon_rate_P,
                                Named("epsilon_psi") = epsilon_psi);
    List kACCEPT = List::create(Named("ACCEPT_RATE_H") = kACCEPT_RATE_H,
                                Named("ACCEPT_RATE_P") = kACCEPT_RATE_P,
                                Named("ACCEPT_LATENT") = kACCEPT_LATENT,
                                Named("ACCEPT_PSI") = kACCEPT_PSI);
    return List::create(Named("THETA") = kTHETA, 
                        Named("TAU_HP") = kTAU_HP, 
                        Named("INDOLENT") = kINDOLENT,
                        Named("ACCEPT") = kACCEPT, 
                        Named("epsilon") = epsilon);
}