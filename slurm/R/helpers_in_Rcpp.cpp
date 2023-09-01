#include <Rcpp.h>
#include <string>
#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Eigen>
#include <RcppNumerical.h>
using namespace Numer;
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

// -----------------------------------------------------------------------------
// There are 3 groups:
//   Grp1 satisfy I(T_i < age_at_tau_hp_i < age_at_tau_pc_i)
//   Grp2 satisfy I(age_at_tau_hp_i < T_i < age_at_tau_pc_i)
//   Grp3 satisfy I(age_at_tau_hp_i < age_at_tau_pc_i == T_i)
// I use this notation in comments.
//   Grp1 come from censored cases; 
//   Grp2 come from screen or censored cases;
//   Grp3 come from clinical cases
// -----------------------------------------------------------------------------

// Recreation of R's rep when two vectors are provided
//
// @param x NumericVector The values to be repeated.
// @param times IntegerVector The number of times to repeat each value.
//   Vector must be of the same length as input `x`.
//
// @returns NumericVector of length sum(times)
//
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
//
// @keywords internal
double rate2scale(double rate, double shape) {
    NumericVector temp(1, rate);
    return Rcpp::pow(temp, -1.0 / shape)[0];
}

// Given a parameter List with updated rates, calculate new scales
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The input `theta` with updated elements `$scale_H`
//   and `$scale_P`.
//
// @keywords internal
List update_scales(List theta) {
    theta["scale_H"] = rate2scale(theta["rate_H"], theta["shape_H"]);
    theta["scale_P"] = rate2scale(theta["rate_P"], theta["shape_P"]);
    return theta;
}

// Probability P(a < X <= b) for a Weibull distribution
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
//   CDF at the ith quantiles of inputs `a` and `b`.
//
// @keywords internal
NumericVector pweibull_ab(NumericVector a, 
                          NumericVector b, 
                          double shape, 
                          double scale) {
    
    return pweibull(b, shape, scale, true, false) -
        pweibull(a, shape, scale, true, false);
}

// A bounded random sample from the Weibull distribution.
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
// @returns NumericVector The generated samples.
//
// @keywords internal
NumericVector rweibull_trunc(NumericVector a, 
                             NumericVector b, 
                             double shape, 
                             double scale) {
    
    NumericVector low = pweibull(a, shape, scale);
    NumericVector high = pweibull(b, shape, scale);
    
    NumericVector u = no_init(low.size());
    // populate each element of `u`, u_i, with a random draw from U(low_i, high_i)
    std::transform(low.begin(), 
                   low.end(), 
                   high.begin(),
                   u.begin(), [=](double low, double high){ return runif(1, low, high)[0]; }); 
    
    return qweibull(u, shape, scale);
}

// Truncated Weibull Probability Density
// STH I was expecting an indicator function for values outside of the range
//
// @param x NumericVector A vector of one or more quantiles at which the
//   density is evaluated.
// @param a NumericVector A vector of one or more quantiles. Must be the
//   same length as input `x`. The lower boundary of the reduced support.
// @param b NumericVector A vector of one or more quantiles. Must be the
//   same length as input `x`. The upper boundary of the reduced support
// @param shape double The shape parameter of the Weibull distribution.
//   Must be > 0.
// @param scale double The scale parameter of the Weibull distribution.
//   Must be > 0.
// @param uselog bool If TRUE, probabilities p are returned as log(p).
//
// @returns NumericVector If uselog = false
//     f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)}
//   if uselog = true
//     ln[ f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)} ]
//
// @keywords internal
NumericVector dweibull_trunc(NumericVector x, 
                             NumericVector a, 
                             NumericVector b, 
                             double shape, 
                             double scale, 
                             bool uselog) {
    if (uselog) {
        // ln[ f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)} ]
        return dweibull(x, shape, scale, true) - log(pweibull_ab(a, b, shape, scale));
    } else {
        // f(x; shape, scale) / {F(a; shape, scale) - F(b; shape, scale)}
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
//
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
//
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
//
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
//
// @keywords internal
List add_psi(List theta, double psi) {
    theta["psi"] = psi;
    return theta;
}

//
// Gibbs theta ////////

// Update beta value based on current healthy -> pre-clinical transition times.
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param n_screen_positive A scalar integer specifying the number of 
//   participants whose healthy -> pre-clinical transitions were detected
//   through a screening, i.e., the total number of participants for whom 
//   censor type is "screen".
//
// @returns List The input `theta` with updated element `$beta`.
//
// @keywords internal
List gibbs_beta(List data_objects, 
                List prior, 
                List age_at_tau_hp_hats,
                List theta, 
                int n_screen_positive) { // beta prior and binomial likelihood conjugacy
    
    IntegerVector lengths;
    NumericVector ages, age_at_tau, tmp_age_at_tau;
    LogicalVector cmp;
    
    List age_screen, obj;
    
    int n_screen = 0;
    
    // count the total number of screens that occurred after the current
    // estimated age at time of healthy -> pre-clinical transition
    for (int i = 0; i < data_objects.size(); ++i) {
        obj = data_objects[i];
        
        // extract vectorized age data for current censor type
        age_screen = obj["ages_screen"];
        ages = age_screen["values"];
        lengths = age_screen["lengths"];
        
        // extract age at time of transition and extend to correlate with 
        // vectorized age structure
        age_at_tau = age_at_tau_hp_hats[i];
        tmp_age_at_tau = repVec(age_at_tau, lengths);
        
        // count the number of screens that occurred after the transition time
        cmp = ages > tmp_age_at_tau;
        n_screen += sum(as<IntegerVector>(cmp));
    }
    
    double a_beta = prior["a_beta"];
    double b_beta = prior["b_beta"];
    // a* = a0 + {# of screens that successfully detected pre-clinical cancers}
    double a_n = a_beta + n_screen_positive;
    // b* = b0 + {# of screens that failed to detect pre-clinical cancers}
    double b_n = b_beta + n_screen - n_screen_positive;
    double beta_new = rbeta(1, a_n, b_n)[0];
    
    return add_beta(theta, beta_new);
}

//
// Likelihood ////////

// The log-likelihood terms involving the sojourn time in the healthy state
//   for all individuals of a specific censor type
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type 
//   under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector  Each element provides
//   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
//
// @keywords internal
NumericVector dloglik_sojourn_H_obj(List data_object, 
                                    List theta,
                                    NumericVector age_at_tau_hp_hat,
                                    double t0) {
    
    // this is a bit of a misnomer. more of an event/censor age as this
    // time is either the age at diagnosis of clinical cancer, age at 
    // screen detection of pre-clinical cancers, or age at censoring
    NumericVector censor_time = data_object["censor_time"];
    int n = data_object["n"];
    
    NumericVector result = no_init(n);
    
    LogicalVector isInfinite = is_infinite(age_at_tau_hp_hat);
    
    // for participants that have an infinite age at time of transition, that is
    // so-called "healthy at censoring time",
    // F_h(T_i - t0, shape, scale); Pr(tau_hp <= T_i - t0)
    // S_H(c_i - t_0; \lambda_H, \alpha_H) = 1 - F_h
    NumericVector vec1 = censor_time[isInfinite];
    NumericVector result_infinite = pweibull(vec1 - t0,
                                             theta["shape_H"],
                                             theta["scale_H"],
                                             false, true);
    result[isInfinite] = result_infinite;
    
    // for participants that have a finite estimated age at time of transition,
    // i.e., those participants in Grp2 and Grp3
    // f_h(age_at_tau_hp_hat_i - t0, shape, scale); Pr(tau_hp == tau_hp_hat)
    vec1 = age_at_tau_hp_hat[!isInfinite];
    NumericVector result_finite = dweibull(vec1 - t0,
                                           theta["shape_H"],
                                           theta["scale_H"],
                                           true);
    
    result[!isInfinite] = result_finite;
    return result;
}

// Sum of the log-likelihood terms involving the sojourn time in the
//   healthy state over all individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns double 
// sum_i = 1, n
//   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
//
// @keywords internal
double dloglik_sojourn_H_sum(List data_objects, 
                             List age_at_tau_hp_hats,
                             List theta, 
                             double t0) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dloglik_sojourn_H_obj(data_objects[i], theta, 
                                            age_at_tau_hp_hats[i], t0));        
    }
    return result;
}

// The log-likelihood terms involving the sojourn time in the progressive
//   pre-clinical state for all individuals of a specific censor type
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type 
//   under consideration.
// @param indolent IntegerVector 0/1 indicating if is not/is indolent.
//
// @returns NumericVector  Each element providing
//   Grp1 {0}
//   Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
//   Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
NumericVector dloglik_sojourn_P_obj(List data_object, 
                                    List theta, 
                                    NumericVector age_at_tau_hp_hat,
                                    IntegerVector indolent) {
    
    // this is a bit of a misnomer. more of an event/censor age as this
    // time is either the age at diagnosis of clinical cancer, age at 
    // screen detection of pre-clinical cancers, or age at censoring
    NumericVector censor_time = data_object["censor_time"];
    
    int censor_type = data_object["censor_type"];
    int n = data_object["n"];
    
    NumericVector result(n, 0.0);
    
    // subset includes participants with progressive tumors from Grp2 and Grp3
    LogicalVector subset = is_finite(age_at_tau_hp_hat) & indolent == 0;
    NumericVector tau_p_hat = censor_time[subset] - age_at_tau_hp_hat[subset];
    
    if (censor_type == 3) {
        // Grp3
        // ln[ f_p(T - age_at_tau_hp_hat)]; ln[ Pr(tau_p == tau_p_hat) ]
        NumericVector res = dweibull(tau_p_hat, theta["shape_P"], 
                                     theta["scale_P"], true);
        result[subset] = res;
    } else {
        // Grp2
        // ln[ F_P(T - age_at_tau_hp_hat) ]; ln[ Pr(tau_p <= tau_p_hat) ]
        NumericVector res = pweibull(tau_p_hat, theta["shape_P"], 
                                     theta["scale_P"], false, true);
        result[subset] = res;
    }
    
    return result;
}

// Sum of the log-likelihood terms involving the sojourn time in the
//   progressive pre-clinical state over all individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for 0/1 indolent.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns double 
// sum_i = 1, n
//   Grp1 {0}
//   Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
//   Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
double dloglik_sojourn_P_sum(List data_objects, 
                             List indolents, 
                             List age_at_tau_hp_hats, 
                             List theta) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dloglik_sojourn_P_obj(data_objects[i],
                                            theta, 
                                            age_at_tau_hp_hats[i],
                                            indolents[i]));
    }
    return result;
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of being indolent for a single censor type
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param indolent IntegerVector 0/1 indicating is not/is indolent.
//
// @returns NumericVector Each element provides
//   ln[psi] I(ind_i == 1) + ln[1 - psi] I(ind_i == 0)
//
// @keywords internal
NumericVector dloglik_indolent_obj(List theta, IntegerVector indolent) {
    double psi = theta["psi"];
    
    int n = indolent.size();
    
    NumericVector result = no_init(n);
    result[indolent == 1] = log(psi);
    result[indolent == 0] = log(1.0 - psi);
    return result;    
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of being indolent for all individuals
//
// @param indolents List The current estimates for 0/1 indolent.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List 
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
//   Each element is a NumericVector providing
//     ln[psi] I(ind_i == 1) + ln[1 - psi] I(ind_i == 0)
//
// @keywords internal
List dloglik_indolent_List(List indolents, List theta) {
    List result(indolents.size());
    for (int i = 0; i < indolents.size(); ++i) {
        result[i] = dloglik_indolent_obj(theta, indolents[i]);
    }
    return result;
}

// The log-likelihood terms involving the parameter of the
//   model for the probability of a screen successfully detecting 
//   pre-clinical status for a single censor type
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type 
//   under consideration.
//
// @returns NumericVector Each element provides 
//   n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//
// @keywords internal
NumericVector dloglik_screens_obj(List data_object, 
                                  List theta, 
                                  NumericVector age_at_tau_hp_hat) {
    
    // extracting vectorized age at screening data
    List ages_screen = data_object["ages_screen"];
    NumericVector age_screen = ages_screen["values"];
    IntegerVector starts = ages_screen["starts"];
    IntegerVector ends = ages_screen["ends"];
    IntegerVector lengths = ages_screen["lengths"];
    
    // =1 if screen detected pre-clinical cancer
    IntegerVector n_screen_positive = data_object["n_screen_positive"];
    
    // current age at tau repeated to align with structure of the vectorized
    // age data
    NumericVector long_taus = repVec(age_at_tau_hp_hat, lengths);
    
    // for each participant, count the number of screens that occurred after
    // tau, but did not successfully identify pre-clinical status
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

// The log-likelihood terms involving the parameter of the
//   model for the probability of a screen successfully detecting 
//   pre-clinical status for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List The list has 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
//   Each element is a NumericVector providing
//     n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//
// @keywords internal
List dloglik_screens_List(List data_objects, List age_at_tau_hp_hats, List theta) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = dloglik_screens_obj(data_objects[i], theta, age_at_tau_hp_hats[i]);
    }
    return result;  
}


// The log-likelihood terms involving the sojourn time in the progressive 
//   pre-clinical state and the indolence parameters for all individuals in
//   of single censor type.
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type 
//   under consideration.
// @param indolent IntegerVector =1 if indolent.
//
// @returns NumericVector Each element providing
//   Grp1 ln[psi]
//   Grp2 {I(ind_i == 1) ln[psi] + 
//         I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
//   Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
NumericVector dloglik_PI_obj(List data_object, 
                             List theta, 
                             NumericVector age_at_tau_hp_hat,
                             IntegerVector indolent) {
    // Grp1 {0} +
    // ln[ F_P(tau_hp_hat_i) ] I(ind_i == 0) Grp2 +
    // ln[ f_P(tau_pc_hat_i) ] I(ind_i == 0) Grp3
    NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, theta,
                                                 age_at_tau_hp_hat, indolent);
    
    // ln[psi] I(ind_i == 1) + ln[1 - psi] I(ind_i == 0)
    NumericVector dlog_I = dloglik_indolent_obj(theta, indolent);
    
    return dlog_P + dlog_I;
}

// Sum of the log-likelihood terms involving the sojourn time in the 
//   progressive pre-clinical state and the indolence parameters over all 
//   individuals
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements.
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns double 
// sum i = 1, n
//   Grp1 ln[psi]
//   Grp2 {I(ind_i == 1) ln[psi] + 
//         I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
//   Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
double dloglik_PI_sum(List data_objects, 
                      List indolent, 
                      List age_at_tau_hp_hats, 
                      List theta) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dloglik_PI_obj(data_objects[i], theta, 
                                     age_at_tau_hp_hats[i], indolent[i]));
    }
    return result;
}

// Class to facilitate numerical integration
// 
// @param shapeH, scaleH, shapeP, scaleP double The parameters of the Weibull
//   distributions.
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
        NumericVector local_vec(1, x);
        // f(x; shape_H, scale_H) F(U - x; shape_P, scale_P)
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

// Compute integral at each age at first screening
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector Ages at first screening
// @param t0 double The initial time.
//
// @returns An Rcpp::NumericVector of the integrals at each ...
// @keywords internal
NumericVector compute_cp_log(List theta, NumericVector AFS, double t0) {
    
    double L = 0.0;        // lower bound
    NumericVector U = AFS - t0; // upper bound
    
    // F(S_0 - t_0; shape, scale)
    NumericVector prob_onset_after = pweibull(U, 
                                              theta["shape_H"], 
                                              theta["scale_H"], 
                                              false, false);
    
    NumericVector integral = no_init(U.size());
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
// @param AFS NumericVector Age at first screening
// @param n_AFS NumericVector Frequency of AFS within the data.
// @param t0 double The initial time.
//
// @returns double The sum of the integrals over all AFS.
//
// @keywords internal
double dloglik_cp(List theta, 
                  NumericVector AFS, 
                  NumericVector n_AFS,
                  double t0) {
    return -sum(compute_cp_log(theta, AFS, t0) * n_AFS);
}

// Sum of the log-likelihood terms involving psi, the probability of
//   indolence
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector Age at first screening
// @param n_AFS NumericVector Frequency of AFS within the data.
// @param t0 double The initial time.
//
// @returns double 
// Unk +
// sum i = 1, n
//   Grp1 ln[psi]
//   Grp2 {I(ind_i == 1) ln[psi] + 
//         I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
//   Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
double dloglik_psi(List data_objects, 
                   List indolent,
                   List age_at_tau_hp_hats,
                   List theta, 
                   NumericVector AFS, 
                   NumericVector n_AFS, 
                   double t0) {
    
    // STH I need to better understand this cp term...
    double dlog_cp = dloglik_cp(theta, AFS, n_AFS, t0);
    
    double dlog_PI = dloglik_PI_sum(data_objects, indolent, 
                                    age_at_tau_hp_hats, theta);
    
    return dlog_cp + dlog_PI;
}

// The log-likelihood terms involving tau for a single censor type
//
// @param data_object List The data for a single censor type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type
//   under consideration.
// @param indolent IntegerVector =1 if indolent.
// @param t0 double Initial time.
//
// @returns NumericVector Each element provides
//   Grp1 ln[ F_h(T_i - t0) ] +
//   Grp2 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]} +
//   Grp3 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]} +
//   n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//
// @keywords internal
NumericVector dloglik_tau_obj(List data_object,
                              List theta, 
                              NumericVector age_at_tau_hp_hat, 
                              IntegerVector indolent, 
                              double t0) {
    
    // Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
    NumericVector dlog_H = dloglik_sojourn_H_obj(data_object, theta, 
                                                 age_at_tau_hp_hat, t0);

    // Grp1 {0}
    // Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
    // Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
    NumericVector dlog_P = dloglik_sojourn_P_obj(data_object, theta, 
                                                 age_at_tau_hp_hat, indolent);
    
    // n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
    NumericVector dlog_S = dloglik_screens_obj(data_object, theta, 
                                               age_at_tau_hp_hat);
    
    return dlog_H + dlog_P + dlog_S;
}

// The log-likelihood terms involving tau for all censor types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
//   Each element is a NumericVector providing
//   Grp1 ln[ F_h(T_i - t0) ] +
//   Grp2 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]} +
//   Grp3 {ln[ f_h(tau_hp_hat_i) ] + I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]} +
//   n_successful_screens_i ln[beta] + n_failed_screens_i ln[1 - beta]
//
// @keywords internal
List dloglik_tau_List(List data_objects, 
                      List indolents, 
                      List age_at_tau_hp_hats, 
                      List theta, 
                      double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = dloglik_tau_obj(data_objects[i], theta, age_at_tau_hp_hats[i], 
                                    indolents[i], t0);
    }
    return result;
}

// Sum of the log-likelihood involving the rate of tau_hp distribution
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector Age at first screening
// @param n_AFS NumericVector Frequency of AFS within the data.
// @param t0 double The initial time.
//
// @returns double
// Unk +
// sum_i = 1, n
//   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
//
// @keywords internal
double dloglik_rate_H(List data_objects, 
                      List age_at_tau_hp_hats, 
                      List theta, 
                      NumericVector AFS, 
                      NumericVector n_AFS, 
                      double t0) {
    
    // STH I need to better understand this cp term...
    double dlog_cp = dloglik_cp(theta, AFS, n_AFS, t0);
    
    double dlog_H = dloglik_sojourn_H_sum(data_objects, age_at_tau_hp_hats, 
                                          theta, t0);
    
    return dlog_cp + dlog_H;
}


// The log-likelihood terms involving the rate of tau_p distribution
// 
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector Age at first screening
// @param n_AFS NumericVector Frequency of AFS within the data.
// @param t0 double The initial time.
//
// @returns double
//   Unk +
//   Grp1 {0}
//   Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
//   Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
//
// @keywords internal
double dloglik_rate_P(List data_objects, 
                      List indolents, 
                      List age_at_tau_hp_hats,
                      List theta, 
                      NumericVector AFS, 
                      NumericVector n_AFS, 
                      double t0) {
    // STH I need to better understand this cp term...
    double dlog_cp = dloglik_cp(theta, AFS, n_AFS, t0);
    
    double dlog_P = dloglik_sojourn_P_sum(data_objects, indolents, age_at_tau_hp_hats, theta);
    
    return dlog_cp + dlog_P;
}

//
// M-H rate_H ////////

// Sample a new rate value for the tau_hp distribution
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_rate_H double A small value shift.
//
// @returns double The rate shifted by U(-eps, eps). Will never be negative.
//
// @keywords internal
double rprop_rate_H(List theta, double epsilon_rate_H) {
    double rate_H = theta["rate_H"];
    rate_H += runif(1, - epsilon_rate_H, epsilon_rate_H)[0];
    return std::fabs(rate_H);
}

// Metropolis-Hastings step for rate of tau_HP distribution
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector A grouped vector providing all unique ages at first
//   screening
// @param n_AFS NumericVector  Frequency of AFS within the data.
// @param epsilon_rate_H A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 2 elements. Theta - the accepted new theta; and
//   accept - TRUE if theta was updated; FALSE otherwise.
//
// @keywords internal
List MH_rate_H(List data_objects,
               List indolents, 
               List prior, 
               List age_at_tau_hp_hats,
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
    
    // Unk +
    // sum_i = 1, n
    //   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
    // calculated with current rate_H value
    double dlog_lik_cur = dloglik_rate_H(data_objects, age_at_tau_hp_hats, theta_cur, 
                                         AFS, n_AFS, t0);
    // Unk +
    // sum_i = 1, n
    //   Grp1 ln[ F_h(T_i - t0) ] + {Grp2 + Grp3} ln[ f_h(tau_hp_hat_i) ]
    // calculated with new rate_H value
    double dlog_lik_new = dloglik_rate_H(data_objects, age_at_tau_hp_hats, theta_new, 
                                         AFS, n_AFS, t0);
    
    double prior_rate = prior["rate_H"];
    
    // f_g(rate_H_new; shape_prior, scale_prior)
    double dlog_prior_new = dgamma(rate_H_new, prior["shape_H"], 
                                   1.0 / prior_rate, true)[0];
    // f_g(rate_H_cur; shape_prior, scale_prior)
    double dlog_prior_cur = dgamma(rate_H_cur, prior["shape_H"], 
                                   1.0 / prior_rate, true)[0];
    
    double MH_logratio = (dlog_lik_new + dlog_prior_new) - 
        (dlog_lik_cur + dlog_prior_cur);
    
    if (runif(1)[0] < std::exp(MH_logratio)) {
        return List::create(Named("theta") = theta_new, Named("accept") = true);
    } else {
        return List::create(Named("theta") = theta_cur, Named("accept") = false);
    }
}

//
// M-H rate_P ////////

// Sample a new rate value for the tau_p distribution
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_rate_P double A small value shift.
//
// @returns double The rate shifted by U(-eps, eps). Will never be negative.
//
// @keywords internal
double rprop_rate_P(List theta, double epsilon_rate_P) {
    double rate_P = theta["rate_P"];
    rate_P = runif(1, rate_P - epsilon_rate_P, rate_P + epsilon_rate_P)[0];
    return std::fabs(rate_P);
}

// Metropolis-Hastings step for rate of tau_P distribution
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements.  
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector A grouped vector providing all unique ages at first
//   screening
// @param n_AFS NumericVector  Frequency of AFS within the data.
// @param epsilon_rate_P A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 2 elements. Theta - the accepted new theta; and
//   accept - TRUE if theta was updated; FALSE otherwise.
//
// @keywords internal
List MH_rate_P(List data_objects,
               List indolents, 
               List prior, 
               List age_at_tau_hp_hats,
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
    
    // Unk +
    // Grp1 {0}
    // Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
    // Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
    // evaluated at current rate_P
    double dlog_lik_cur = dloglik_rate_P(data_objects, indolents, age_at_tau_hp_hats, 
                                         theta_cur, AFS, n_AFS, t0);
    // Unk +
    // Grp1 {0}
    // Grp2 I(ind_i == 0) ln[ F_P(tau_p_hat_i) ]
    // Grp3 I(ind_i == 0) ln[ f_P(tau_p_hat_i) ]
    // evaluated at new rate_P
    double dlog_lik_new = dloglik_rate_P(data_objects, indolents, age_at_tau_hp_hats,
                                         theta_new, AFS, n_AFS, t0);
    
    double prior_rate = prior["rate_P"];
    // f_gamma(rate_P_new; shape_prior, scale_prior)
    double dlog_prior_new = dgamma(rate_P_new, prior["shape_P"], 
                                   1.0 / prior_rate, true)[0];
    // f_gamma(rate_P_current; shape_prior, scale_prior)
    double dlog_prior_cur = dgamma(rate_P_cur, prior["shape_P"], 
                                   1.0 / prior_rate, true)[0];
    
    double MH_logratio = (dlog_lik_new + dlog_prior_new) - 
        (dlog_lik_cur + dlog_prior_cur);
    
    if (runif(1)[0] < std::exp(MH_logratio)) {
        return List::create(Named("theta") = theta_new, Named("accept") = true);
    } else {
        return List::create(Named("theta") = theta_cur, Named("accept") = false);
    }
}


//
// M-H age_at_tau_hp_hat ////////

// Compute the probability of tau falling within the observed intervals for a 
//   single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @return List Elements include `values`, a NumericVector of all probabilities,
//   `starts`, an IntegerVector containing the first element of `values
//   pertaining to case i; `ends`, an IntegerVector containing the last element
//   pertaining to individual i; and `lengths` the number of element in `values`
//   pertaining to individual i.
//
// @keywords internal
List compute_prob_tau_obj(List data_object, List theta, double t0) {
    
    // this is a bit of a misnomer. more of an event/censor age as this
    // time is either the age at diagnosis of clinical cancer, age at 
    // screen detection of pre-clinical cancers, or age at censoring
    NumericVector censor_time = data_object["censor_time"];
    
    // extract vectorized observed interval data
    // note that these interval boundaries are actually "age at time of"
    // the interval is obtained by subtracting the lower age if no clinical
    // diagnosis or by subtracting from the age at clinical diagnosis
    List endpoints = data_object["endpoints"];
    NumericVector ep = endpoints["values"];
    IntegerVector starts = endpoints["starts"];
    IntegerVector ends = endpoints["ends"];
    IntegerVector lengths = endpoints["lengths"];
    
    int censor_type = data_object["censor_type"];
    
    // estimating the probability of each interval
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
        // Pr(interval_k_lower < tau_hp <= interval_k_upper)
        prob_interval = pweibull_ab(vec1 - t0, vec2 - t0, 
                                    theta["shape_H"], theta["scale_H"]);
        
    } else {
        // need to repeat censor time to match interval data structure
        NumericVector ct = repVec(censor_time, pt_lengths);
        // Pr(interval_l_lower < tau_p <= interval_l_upper)
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

// Compute probability of tau for all censor types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is itself a List
//   containing the probabilities of tau details for the specific
//   censoring type.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List compute_prob_tau_List(List data_objects, List theta, double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = compute_prob_tau_obj(data_objects[i], theta, t0);
    }
    return result;
}

// Sample current estimated probability of tau for new tau_hp
//  for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param prob_tau List The probabilities of tau for all cases in the
//   censoring type under consideration.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns NumericVector The new estimated tau for each case in the 
//   censoring type under consideration
//
// @keywords internal
NumericVector rprop_age_at_tau_hp_hat_obj(List data_object, 
                                       List prob_tau, 
                                       List theta, 
                                       double t0) {
    
    // this is a bit of a misnomer. more of an event/censor age as this
    // time is either the age at diagnosis of clinical cancer, age at 
    // screen detection of pre-clinical cancers, or age at censoring
    NumericVector censor_time = data_object["censor_time"];
    
    // extract vectorized interval information from data list
    // note that these are "age at time of" boundaries
    List endpoints = data_object["endpoints"];
    NumericVector ep = endpoints["values"];
    IntegerVector ep_starts = endpoints["starts"];
    IntegerVector ep_ends = endpoints["ends"];
    IntegerVector ep_lengths = endpoints["lengths"];
    
    int censor_type = data_object["censor_type"];
    
    // extract vectorized probabilities of tau
    NumericVector pt = prob_tau["values"];
    IntegerVector pt_starts = prob_tau["starts"];
    IntegerVector pt_ends = prob_tau["ends"];
    IntegerVector pt_lengths = prob_tau["lengths"];
    
    // sample each case's interval
    IntegerVector k_new = no_init(censor_time.size());
    for (int i = 0; i < pt_starts.size(); ++i) {
        int K = pt_lengths[i];
        NumericVector vec = pt[seq(pt_starts[i], pt_ends[i])];
        k_new[i] = sample(K, 1, false, vec)[0] - 1;
    }
    
    NumericVector age_at_tau_hp_hat_new;
    
    NumericVector vec1 = ep[ep_starts + k_new];
    NumericVector vec2 = ep[ep_starts + k_new + 1];
    
    if (censor_type == 1) {
        // sample tau_hp from (interval_lower, interval_upper]
        NumericVector sojourn_H_new = rweibull_trunc(vec1 - t0, vec2 - t0,
                                                     theta["shape_H"], 
                                                     theta["scale_H"]);
        age_at_tau_hp_hat_new = t0 + sojourn_H_new;
    } else if(censor_type == 2) {
        // sample tau_hp from (interval_lower, interval_upper]
        // STH At this stage of the game I'm wondering why we don't sample
        // tau_p for those that are not indolent?
        NumericVector sojourn_H_new = rweibull_trunc(vec1 - t0, vec2 - t0,
                                                     theta["shape_H"], 
                                                     theta["scale_H"]);
        age_at_tau_hp_hat_new = t0 + sojourn_H_new;
        
        // if the estimated age at tau_hp > censor_time set tau = Inf
        age_at_tau_hp_hat_new[k_new == (pt_lengths - 1)] = R_PosInf;
    } else {
        // sample tau_p from (interval_lower, interval_upper]
        NumericVector sojourn_P_new = rweibull_trunc(censor_time - vec2, 
                                                     censor_time - vec1,
                                                     theta["shape_P"], 
                                                     theta["scale_P"]);
        // age at tau_hp is age at clinical diagnosis - tau_p
        age_at_tau_hp_hat_new = censor_time - sojourn_P_new;
    }
    return age_at_tau_hp_hat_new;
}


// Sample current estimated probability of tau for new tau_hp for all censor
//  censor types
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_tau List The probabilities of tau broken down according to the
//   censor type. There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List rprop_age_at_tau_hp_hat_List(List data_objects, 
                                  List prob_tau, 
                                  List theta, 
                                  double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = rprop_age_at_tau_hp_hat_obj(data_objects[i], prob_tau[i], 
                                                theta, t0);
    }
    return result;
}

// Contribution to the log-likelihood from probability of tau
//
// @param data_object List The data for a single censoring type.
// @param prob_tau List The probabilities of tau for all cases in the
//   censoring type under consideration.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical for all participants of the censor type under 
//   consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector 
//
// @keywords internal
NumericVector dlog_prop_age_at_tau_hp_hat_obj(List data_object,
                                              List prob_tau,
                                              List theta,
                                              NumericVector age_at_tau_hp_hat, 
                                              double t0) {
    
    // this is a bit of a misnomer. more of an event/censor age as this
    // time is either the age at diagnosis of clinical cancer, age at 
    // screen detection of pre-clinical cancer, or age at censoring
    NumericVector censor_time = data_object["censor_time"];
    
    // extract vectorized interval boundary information
    // note that this an "age at time of" vector
    List endpoints = data_object["endpoints"];
    NumericVector ep = endpoints["values"];
    IntegerVector ep_starts = endpoints["starts"];
    IntegerVector ep_ends = endpoints["ends"];
    IntegerVector ep_lengths = endpoints["lengths"];
    
    int censor_type = data_object["censor_type"];
    
    // extract vectorized probability of tau information
    NumericVector pt = prob_tau["values"];
    IntegerVector pt_starts = prob_tau["starts"];
    IntegerVector pt_ends = prob_tau["ends"];
    IntegerVector pt_lengths = prob_tau["lengths"];
    
    // replicate age_at_tau_hp_hat_i to correlate with vectorized boundaries
    NumericVector tmp_age_at_tau = repVec(age_at_tau_hp_hat, ep_lengths);
    
    // count the number of screens before transition for each patient
    LogicalVector cmp = ep < tmp_age_at_tau;
    IntegerVector cmp_csum = cumsum(as<IntegerVector>(cmp));
    LogicalVector cmp_starts = cmp[ep_starts];
    IntegerVector icmp_starts = as<IntegerVector>(cmp_starts);
    
    // k_new is the last screen before the current estimated transition time
    IntegerVector k_new = cmp_csum[ep_ends] - cmp_csum[ep_starts] + icmp_starts - 1;
    // extract the corresponding probability of tau
    NumericVector dlog_k = log(pt[pt_starts + k_new]);
    
    NumericVector vec1 = ep[ep_starts + k_new];
    NumericVector vec2 = ep[ep_starts + k_new + 1];
    
    NumericVector dlog_tau;
    if (censor_type == 1) {
        // P(tau_hp_i == tau_hp_hat_i | interval_lower < tau_hp <= interval_upper)
        dlog_tau = dweibull_trunc(age_at_tau_hp_hat - t0, 
                                  vec1 - t0, 
                                  vec2 - t0,
                                  theta["shape_H"], theta["scale_H"], true);
    } else if (censor_type == 2) {
        // P(tau_hp_i == tau_hp_hat_i | interval_lower < tau_hp <= interval_upper)
        // STH again wondering why we don't use the tau_p for indolent = 0 cases
        dlog_tau = dweibull_trunc(age_at_tau_hp_hat - t0, 
                                  vec1 - t0, 
                                  vec2 - t0,
                                  theta["shape_H"], theta["scale_H"], true);
        // if estimated age at tau_hp > censor time, set to 0
        dlog_tau[(k_new + 1) == pt_lengths] = 0.0;
    } else {
        // P(tau_p_i == tau_p_hat_i | interval_lower < tau_p <= interval_upper)
        dlog_tau = dweibull_trunc(censor_time - age_at_tau_hp_hat, 
                                  censor_time - vec2, 
                                  censor_time - vec1,
                                  theta["shape_P"], theta["scale_P"], true);
    }
    return dlog_k + dlog_tau;
}

// The log-likelihood terms involving probability of tau for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_tau List The probabilities of tau broken down according to the
//   censor type. There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped
//   according to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
//
// @keywords internal
List dlog_prop_age_at_tau_hp_hat_List(List data_objects, 
                                      List prob_taus, 
                                      List age_at_tau_hp_hats, 
                                      List theta, double t0) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = dlog_prop_age_at_tau_hp_hat_obj(data_objects[i], prob_taus[i], 
                                                    theta, age_at_tau_hp_hats[i], 
                                                    t0);
    }
    return result;
}

//
// indolent ////////

// Probability of being indolent for each case of a single censor type given 
//   current parameter and transition time estimates
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type
//   under consideration.
// @param t0 A scalar double The initial time.
//
// @returns NumericVector 
//
// @keywords internal
NumericVector compute_prob_indolent_obj(List data_object, List theta,
                                        NumericVector age_at_tau_hp_hat) {
    
    IntegerVector ind(age_at_tau_hp_hat.size(), 0);
    
    // Grp1 ln[psi]
    // Grp2 {I(ind_i == 1) ln[psi] + 
    //       I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
    // Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
    // evaluated when all cases are set as non-indolent
    NumericVector L_0 = dloglik_PI_obj(data_object, theta, age_at_tau_hp_hat, ind);
    L_0 = exp(L_0);
    
    std::fill(ind.begin(), ind.end(), 1);
    // Grp1 ln[psi]
    // Grp2 {I(ind_i == 1) ln[psi] + 
    //       I(ind_i == 0) {ln[1 - psi] + ln[ F_P(tau_p_hat_i) ]}} +
    // Grp3 ln[1 - psi] + ln[ f_P(tau_p_hat_i) ]
    // evaluated when all cases are set as indolent
    NumericVector L_1 = dloglik_PI_obj(data_object, theta, age_at_tau_hp_hat, ind);
    L_1 = exp(L_1);
    
    return L_1 / (L_0 + L_1);
}

// Probability of being indolent for all cases given current parameter and
//  transition time estimates
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   according to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List There are 3 elements. Each element is a NumericVector.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List compute_prob_indolent_List(List data_objects, List age_at_tau_hp_hats, List theta) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = compute_prob_indolent_obj(data_objects[i], theta, age_at_tau_hp_hats[i]);
    }
    return result;
}

// Generate vector of indolence for all cases of a single censor type based on
//   current estimated probability of indolence.
//
// Note: indolence is always 0 (false) for clinical cases
//
// @param data_object List The data for a single censoring type.
// @param prob_indolent NumericVector the probability of each case being
//   indolent.
//
// @returns IntegerVector Indicator of indolence.
//
// @keywords internal
IntegerVector rprop_indolent_obj(List data_object, NumericVector prob_indolent) {
    
    int censor_type = data_object["censor_type"];
    
    int n = prob_indolent.size();
    IntegerVector indolent(n);
    
    if (censor_type == 3) {
        // indolence is always false for those with a clinical diagnosis
        indolent.fill(0);
    } else {
        // indolence is a random sample from binomial using the current
        // estimated probability of indolence
        std::transform(prob_indolent.begin(), 
                       prob_indolent.end(), 
                       indolent.begin(), [=](double p){ return rbinom(1, 1, p)[0]; }); 
    }
    
    return indolent;
}

// Generate vector of indolence for all cases.
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_indolents List The current indolence probability broken down
//   according to the censor type. There are 3 elements. 
//   (1) "screen", (2) "censored", and (3) "clinical".
//
// @returns List There are 3 elements. Each element is an IntegerVector.
//
// exported to initialize these values at start of MC
// [[Rcpp::export]]
List rprop_indolent_List(List data_objects, List prob_indolents) {
    List result(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        result[i] = rprop_indolent_obj(data_objects[i], prob_indolents[i]);
    }
    return result;
}


// Case contributions to the log-likelihood for the probability of indolence
//   for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param prob_indolent NumericVector the probability of each case being
//   indolent.
// @param indolent IntegerVector the current indolence indicator
//
// @returns NumericVector 
//
// {Grp1 + Grp2} {I(ind_i == 1) ln[ Pr(Ind == 1) ] + 
//                I(ind_i == 0) ln[ 1 - Pr(Ind == 1) ]} + Grp3 {0}
//
// @keywords internal
NumericVector dlog_prop_indolent_obj(List data_object, 
                                     NumericVector prob_indolent,
                                     IntegerVector indolent) {
    
    int censor_type = data_object["censor_type"];
    
    NumericVector result = log(1.0 - prob_indolent);
    if (censor_type == 3) {
        // dlog-likelihood is 0 for all clinical cases
        result.fill(0.0);
    } else {
        // dlog-likelihood is log(1-p) for 0 cases; log(p) for 1 cases
        LogicalVector ind1 = indolent == 1;
        NumericVector prob1 = log(prob_indolent[ind1]);
        result[ind1] = prob1;
    }
    
    return result;
}

// Log-likelihood for the probability of indolence for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current indolence indicator broken down
//   according to the censor type. There are 3 elements.
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prob_indolents List The current probability of indolence broken down
//   according to the censor type. There are 3 elements.
//   (1) "screen", (2) "censored", and (3) "clinical".
//
// @returns double
// sum i = 1, n
//   {Grp1 + Grp2} {I(ind_i == 1) ln[ Pr(Ind == 1) ] + 
//                  I(ind_i == 0) ln[ 1 - Pr(Ind == 1) ]} + Grp3 {0}
//
// @keywords internal
double dlog_prop_indolent_sum(List data_objects, 
                              List indolents,
                              List prob_indolents) {
    double result = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        result += sum(dlog_prop_indolent_obj(data_objects[i], prob_indolents[i], 
                                             indolents[i]));
    }
    return result;
}

//
// psi ////////

// Random walk in psi space
//
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param epsilon_psi A double. The maximum step size.
//
// @returns A double. The new psi value
//
// @keywords internal
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

// Metropolis-Hastings step for tau distributions for a single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type 
//   under consideration.
// @param indolent IntegerVector Indicator of indolence.
// @param t0 A scalar double The initial time.
//
// @returns List Element `$age_at_tau_hp_hat` is the possibly updated List of
//   tau values and element `$accepted` is a LogicalVector, where true indicates
//   that the corresponding element of `$age_at_tau_hp_hat` has been updated.
//
// @keywords internal
List MH_tau_obj(List data_object, List theta, NumericVector age_at_tau_hp_hat, 
                IntegerVector indolent, double t0) {
    
    // propose new latent variables age_at_tau_hp_hat/tau_CP
    List prob_tau = compute_prob_tau_obj(data_object, theta, t0);
    
    NumericVector age_at_tau_hp_hat_new = rprop_age_at_tau_hp_hat_obj(data_object, 
                                                                      prob_tau, 
                                                                      theta, t0);
    
    // M-H acceptance ratio
    NumericVector dlog_prop_cur = dlog_prop_age_at_tau_hp_hat_obj(data_object,
                                                       prob_tau,
                                                       theta,
                                                       age_at_tau_hp_hat, 
                                                       t0);
    
    NumericVector dlog_prop_new = dlog_prop_age_at_tau_hp_hat_obj(data_object,
                                                       prob_tau,
                                                       theta,
                                                       age_at_tau_hp_hat_new,
                                                       t0);
    
    NumericVector dlog_lik_cur = dloglik_tau_obj(data_object,
                                                 theta,
                                                 age_at_tau_hp_hat, 
                                                 indolent, 
                                                 t0);
    
    NumericVector dlog_lik_new = dloglik_tau_obj(data_object,
                                                 theta,
                                                 age_at_tau_hp_hat_new, 
                                                 indolent, 
                                                 t0);
    
    NumericVector MH_logratio = (dlog_lik_new - dlog_prop_new) - 
        (dlog_lik_cur - dlog_prop_cur);
    
    LogicalVector test = runif(MH_logratio.size()) < exp(MH_logratio);
    
    NumericVector accepted_tau = age_at_tau_hp_hat;
    accepted_tau[test] = age_at_tau_hp_hat_new[test];
    
    return List::create(Named("age_at_tau_hp_hat") = accepted_tau,
                        Named("accept") = test);
    
}

// Metropolis-Hastings step for rate of tau distribution for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current indolence indicator broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical"
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according 
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param t0 A scalar double The initial time.
//
// @returns List Element `$age_at_tau_hp_hats` is the updated estimates for tau 
//   broken down according to the censor type and element `$accept` is a list 
//   broken down according to the censor type indicating if the corresponding
//   elements of `$age_at_tau_hp_hats` were updated.
//
// @keywords internal
List MH_tau_List(List data_objects, 
                 List indolents, 
                 List age_at_tau_hp_hats, 
                 List theta, 
                 double t0) {
    
    List result_tau(data_objects.size());
    List result_accept(data_objects.size());
    for (int i = 0; i < data_objects.size(); ++i) {
        List res = MH_tau_obj(data_objects[i], theta, age_at_tau_hp_hats[i], 
                              indolents[i], t0);
        result_tau[i] = res["age_at_tau_hp_hat"];
        result_accept[i] = res["accept"];
    }
    
    return List::create(Named("age_at_tau_hp_hats") = result_tau,
                        Named("accept") = result_accept);  
}


//
// M-H (psi, indolent) ////////

// Sample current probability of indolence distribution for each case of a 
//   single censor type
//
// @param data_object List The data for a single censoring type.
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param age_at_tau_hp_hat NumericVector The current estimated age at time of 
//   healthy -> pre-clinical transition for all participants of the censor type
//   under consideration.
//
// @returns List Element `$indolent` is the updated indolence indicator and
//   element `$dlog_prop` is the updated derivative evaluated at the updated
//   indolence indicators.
//
// @keywords internal
List rprop_dlog_indolent_obj(List data_object, List theta, NumericVector age_at_tau_hp_hat) {
    
    NumericVector prob_indolent_new = compute_prob_indolent_obj(data_object,
                                                                theta,
                                                                age_at_tau_hp_hat);
    // Generate new indolence vector
    IntegerVector indolent_new = rprop_indolent_obj(data_object, prob_indolent_new);

    double dlog_prop_new = sum(dlog_prop_indolent_obj(data_object, 
                                                      prob_indolent_new, 
                                                      indolent_new));
    
    return List::create(Named("indolent") = indolent_new,
                        Named("dlog_prop") = dlog_prop_new);
}

// Sample current probability of indolence distribution for all cases
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
//
// @returns List Element `$indolent_new` is a List of the the updated 
//   indolence indicator for each censor type and element `$dlog_prop_indolent_new`
//   is the updated derivative for each censor type.
//
// @keywords internal
List rprop_dlog_indolent_List(List data_objects, List age_at_tau_hp_hats, List theta) {
    List result_indolent(data_objects.size());
    double dlog_prop = 0.0;
    for (int i = 0; i < data_objects.size(); ++i) {
        List res = rprop_dlog_indolent_obj(data_objects[i], theta, age_at_tau_hp_hats[i]);
        result_indolent[i] = res["indolent"];
        double dlog = res["dlog_prop"];
        dlog_prop += dlog;
    }
    
    return List::create(Named("indolent_new") = result_indolent,
                        Named("dlog_prop_indolent_new") = dlog_prop);  
}

// Metropolis-Hastings step for psi distribution
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical"
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta_cur List A named List object containing the parameters of the 
//   distributions. 
// @param AFS NumericVector Age at first screening
// @param n_AFS NumericVector Frequency of AFS within the data.
// @param epsilon_psi A small shift value. Must be > 0.
// @param t0 double The initial time.
//
// @returns List There are 3 elements. `$indolence` the possibly updated indolence
//   indicator broken down by censor type; `$theta` the accepted new theta; and
//   `$accept`, true if theta and indolence wer updated; FALSE otherwise.
//
// @keywords internal
List MH_psi_indolent(List data_objects,
                     List indolents, 
                     List prior, 
                     List age_at_tau_hp_hats,
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
    List out = rprop_dlog_indolent_List(data_objects, age_at_tau_hp_hats, 
                                        theta_new);
    List indolent_new = out["indolent_new"];
    
    // proposal density
    double dlog_prop_indolent_new = out["dlog_prop_indolent_new"];
    List prob_indolent_cur = compute_prob_indolent_List(data_objects, 
                                                        age_at_tau_hp_hats, 
                                                        theta_cur);
    double dlog_prop_indolent_cur = dlog_prop_indolent_sum(data_objects,
                                                           indolents, 
                                                           prob_indolent_cur);
    
    // log likelihood
    double dlog_lik_cur = dloglik_psi(data_objects, indolents, 
                                      age_at_tau_hp_hats, 
                                      theta_cur, AFS, n_AFS, t0);
    double dlog_lik_new = dloglik_psi(data_objects, indolent_new, 
                                      age_at_tau_hp_hats, 
                                      theta_new, AFS, n_AFS, t0);
    
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

// Full Markov Chain Monte Carlo Procedure
//
// @param data_objects List The data broken down according to the censor type.
//   There are 3 elements: (1) "screen", (2) "censored", and (3) "clinical".
// @param indolents List The current estimates for indolent broken down
//   according to the censor type. There are 3 elements:
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param prior List The distribution parameters of the prior.
// @param age_at_tau_hp_hats List The current estimates for the age at time of 
//   healthy -> pre-clinical transition for each participant grouped according
//   to the censor type. There are 3 elements: 
//   (1) "screen", (2) "censored", and (3) "clinical".
// @param theta List A named List object containing the parameters of the 
//   distributions. 
// @param double epsilon_rate_H A positive scalar double specifying the maximum
//   step size in the rate_H space.
// @param double epsilon_rate_P A positive scalar double specifying the maximum
//   step size in the rate_P space.
// @param double epsilon_psi A positive scalar double specifying the maximum
//   step size in the psi space.
// @param t0 double The initial time.
// @param M A positive scalar integer specifying the number of MC samples.
// @param thin A positive scalar integer specifying the number of samples to
//   be skipped between accepted samples.
// @param M_thin A positive scalar integer specifying the total number of
//   samples that will be accepted.
// @param n_obs A positive scalar integer specifying the total number of
//   participants. Note this should be equivalent to the sum of element
//   `$n` in each element of `data_objects`.
// @param n_screen_positive_total A positive scalar integer indicating the
//   number of positive screens in the full data.
// @param AFS NumericVector Age at first screening
// @param n_AFS NumericVector Frequency of AFS within the data.
//
// @returns List 
//
// exported main function
// [[Rcpp::export]]
List MCMC_cpp(List data_objects, List indolents, List prior, 
              List age_at_tau_hp_hats, List theta,
              double epsilon_rate_H, double epsilon_rate_P, double epsilon_psi,
              double t0, int M, int thin, int M_thin, int n_obs,
              int n_screen_positive_total,
              NumericVector AFS, NumericVector n_AFS) {
    
    NumericVector kRATE_H(M_thin), kRATE_P(M_thin), kBETA(M_thin), kPSI(M_thin);
    LogicalVector kACCEPT_PSI(M_thin), kACCEPT_RATE_H(M_thin), kACCEPT_RATE_P(M_thin);
    NumericMatrix kage_at_tau_hp_hat(M_thin, n_obs);
    IntegerMatrix kINDOLENT(M_thin, n_obs);
    LogicalMatrix kACCEPT_LATENT(M_thin, n_obs);
    int ikept = -1;
    
    for (int m = 0; m < M; ++m) {
        // M-H rate_H
        List out_rate_H = MH_rate_H(data_objects, indolents, prior, 
                                    age_at_tau_hp_hats, theta, 
                                    AFS, n_AFS, epsilon_rate_H, t0);
        theta = out_rate_H["theta"];
        bool accept_rate_H = out_rate_H["accept"];
        
        // M-H rate_P
        List out_rate_P = MH_rate_P(data_objects, indolents, prior, age_at_tau_hp_hats, theta, 
                                    AFS, n_AFS, epsilon_rate_P, t0);
        theta = out_rate_P["theta"];
        bool accept_rate_P = out_rate_P["accept"];
        
        // Gibbs beta
        theta = gibbs_beta(data_objects, prior, age_at_tau_hp_hats, theta, n_screen_positive_total);
        
        // update (psi, indolent)
        List out_psi_indolent = MH_psi_indolent(data_objects, indolents, prior, age_at_tau_hp_hats, theta,
                                                AFS, n_AFS, epsilon_psi, t0);
        indolents = out_psi_indolent["indolents"];
        theta = out_psi_indolent["theta"];
        bool accept_psi = out_psi_indolent["accept"];
        
        // update age_at_tau_hp_hat
        List out_tau  = MH_tau_List(data_objects, indolents, age_at_tau_hp_hats, theta, t0);
        age_at_tau_hp_hats = out_tau["age_at_tau_hp_hats"];
        List accept_latent = out_tau["accept"];
        // save output
        if (m % thin == 0) {
            NumericVector tau_screen = age_at_tau_hp_hats[0];
            IntegerVector indolent_screen = indolents[0];
            LogicalVector accept_latent_screen = accept_latent[0];
            
            NumericVector tau_censored = age_at_tau_hp_hats[1];
            IntegerVector indolent_censored = indolents[1];
            LogicalVector accept_latent_censored = accept_latent[1];
            
            NumericVector tau_clinical = age_at_tau_hp_hats[2];
            IntegerVector indolent_clinical = indolents[2];
            LogicalVector accept_latent_clinical = accept_latent[2];
            
            NumericVector age_at_tau_hp_hat(tau_screen.size() + 
                tau_censored.size() + tau_clinical.size());
            IntegerVector indolent(indolent_screen.size() + 
                indolent_censored.size() + indolent_clinical.size());
            LogicalVector accept_latent(accept_latent_screen.size() + 
                accept_latent_censored.size() + accept_latent_clinical.size());
            
            int n_screen = tau_screen.size();
            int n_censored = tau_censored.size();
            int n_clinical = tau_clinical.size();
            
            age_at_tau_hp_hat[seq(0, n_screen-1)] = tau_screen;
            age_at_tau_hp_hat[seq(n_screen, n_screen + n_censored - 1)] = tau_censored;
            age_at_tau_hp_hat[seq(n_screen+ n_censored, n_screen + n_censored + n_clinical - 1)] = tau_clinical;
            
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
            kage_at_tau_hp_hat  (ikept, _ ) = age_at_tau_hp_hat;
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
                        Named("age_at_tau_hp_hat") = kage_at_tau_hp_hat, 
                        Named("INDOLENT") = kINDOLENT,
                        Named("ACCEPT") = kACCEPT, 
                        Named("epsilon") = epsilon);
}