#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec ComputeExpBetaX(arma::vec beta, arma::mat X){
  long n = X.n_rows;
  arma::vec exp_beta_x(n);
  exp_beta_x = exp(X * beta);
  return exp_beta_x;
}

arma::vec ComputeS(arma::vec Time,
                   arma::vec Time_bounds,
                   arma::vec lambda,
                   arma::vec exp_beta_x,
                   arma::mat simpdat){
  long n = Time_bounds.size();
  long m = Time.size();
  arma::vec S(n); S.zeros();
  for(int i = 0; i < n; ++i){
    if(Time_bounds(i) == INFINITY){
      S(i) = INFINITY;
      continue;
    }
    long id_i = simpdat(i,4);
    for(int k = 0; k < m; ++k){
      if(Time(k) <= Time_bounds(i)){
        S(i) += lambda(k);
      }
    }
    S(i) = S(i) * exp_beta_x(id_i-1);
  }
  
  return S;
}

arma::mat Estep(arma::vec exp_beta_x, 
                arma::vec lambda, 
                arma::mat X, 
                arma::mat simpdat,
                arma::vec Time){
  long n = X.n_rows;
  long m = Time.size();
  arma::mat w(n, m); w.zeros();
  arma::vec L = simpdat.col(1);
  arma::vec R = simpdat.col(2);
  arma::vec S_L = ComputeS(Time, L, lambda, exp_beta_x, simpdat);
  arma::vec S_R = ComputeS(Time, R, lambda, exp_beta_x, simpdat);
  // use sorted simple data
  for (long i = 0; i < n; ++i) {
    double L_i = simpdat(i,1);
    double R_i = simpdat(i,2);
    long id_i = simpdat(i,4);
    
    for (long k = 0; k < m; ++k) {
      double t_k = Time(k);
      if (R_i != INFINITY && t_k > L_i && t_k <= R_i) {
        w(i,k) = lambda(k) * exp_beta_x(id_i-1) / (1.0 - exp(S_L(i) - S_R(i)));
      }
    }
  }
  return w;
}

double AbsoluteDifferenceSafe(double prev, double next, double delta) {
  double min = std::abs(prev);
  if(std::abs(next) < std::abs(prev)) { min = std::abs(next); }
  if (min <= delta) {
    return std::abs(prev - next);
  }else{
    return std::abs(prev - next) / min;
  }
}

double VectorAbsoluteDifferenceSafe(arma::vec prev, arma::vec next, double delta) {
  
  double max_difference = 0.0;
  for (int i = 0; i < prev.size(); ++i) {
    double current_difference = AbsoluteDifferenceSafe(prev(i), next(i), delta);
    if (current_difference > max_difference) max_difference = current_difference;
  }
  return max_difference;
}

arma::vec CumulativeLambda(arma::vec lambda){
  arma::vec cumulative_lambda(lambda.size()); cumulative_lambda.zeros();
  for (int i = 0; i < lambda.size(); ++i) {
    if (i != 0) {
      cumulative_lambda(i) = cumulative_lambda(i - 1);
    }
    cumulative_lambda(i) += lambda(lambda.size()-i-1);
  }
  return cumulative_lambda;
}

arma::vec profile_likelihood(arma::vec beta, 
                             arma::vec lambda_initial, 
                             arma::mat X, 
                             arma::mat simpdat, 
                             double eps, 
                             arma::vec Time){
  long m = Time.size();
  long n = X.n_rows;
  arma::vec lambda_hat(m); lambda_hat = lambda_initial;
  arma::vec lambda_old; lambda_old = lambda_initial;
  arma::vec exp_beta_x_new = ComputeExpBetaX(beta, X);
  
  double error = 100;
  int Iter = 0;
  int MaxIter = 1000;
  double delta = 0.01;
  
  while (error > eps && Iter < MaxIter){
    arma::mat w = Estep(exp_beta_x_new, lambda_old, X, simpdat, Time);
    //update lambda
    for (int k = 0; k < m; ++k) {
      double num = 0;
      double den = 0;
      for (int i = 0; i < n; ++i) {
        if (simpdat(i,3)>=Time(k)) {
          int id = simpdat(i,0);
          num += w(i,k);
          den += exp_beta_x_new(id-1);
        }
      }
      lambda_hat(k) = num / den;
    }
    
    arma::vec cumulative_lambda_old = CumulativeLambda(lambda_old);
    arma::vec cumulative_lambda_hat = CumulativeLambda(lambda_hat);
    error = VectorAbsoluteDifferenceSafe(cumulative_lambda_hat, cumulative_lambda_old, delta);
    lambda_old = lambda_hat;
    Iter = Iter + 1;
    //std::cout << "Lambda Iter = " << Iter << ", Error = " << error << "\n" << std::endl;
  }
  
  if(error <= eps){
    std::cout << "Lambda converged in " << Iter << " interations" << ", Error = " << error << std::endl;
  } else{
    std::cout << "Lambda did not converge!\n" << std::endl;
  }
  
  return lambda_hat;
}

double likelihood_i(arma::vec beta, 
                    arma::vec lambda, 
                    arma::mat X, 
                    arma::mat simpdat, 
                    arma::vec Time, 
                    long i){
  double pl;
  long id_i = simpdat(i,0);
  double L_i = simpdat(i,1);
  double R_i = simpdat(i,2);
  arma::rowvec X_i = X.row(id_i-1);
  double exp_beta_Xi = arma::as_scalar(exp(X_i * beta));
  arma::uvec index1 = arma::find(Time <= L_i);
  arma::uvec index2 = arma::find(Time <= R_i);
  
  if(R_i == INFINITY){
    pl = - exp_beta_Xi * arma::sum(lambda(index1));
  } else{
    double sum1 = exp(-exp_beta_Xi * arma::sum(lambda(index1)));
    double sum2 = exp(-exp_beta_Xi * arma::sum(lambda(index2)));
    pl = log(sum1 - sum2);
  }
  return(pl);
}


// [[Rcpp::export]]
Rcpp::List get_W(arma::vec beta_new, arma::mat X_new, arma::vec lambda_old,
                 arma::mat simpdat, double eps, arma::vec Time, double constant){
  Rcpp::List ret;
  long n = X_new.n_rows;
  long p = X_new.n_cols;
  double h = constant/sqrt(n);
  
  arma::mat Dpl(n, p);
  
  //first order
  for(long j = 0; j < p; ++j){
    arma::vec beta_j1 = beta_new;
    arma::vec beta_j2 = beta_new;
    beta_j1(j) = beta_j1(j) + h;
    beta_j2(j) = beta_j2(j) - h;
    arma::vec lambda_j1 = profile_likelihood(beta_j1, lambda_old, X_new, simpdat, eps, Time);
    arma::vec lambda_j2 = profile_likelihood(beta_j2, lambda_old, X_new, simpdat, eps, Time);
    for(long i = 0; i < n; ++i){
      Dpl(i,j) = (likelihood_i(beta_j1, lambda_j1, X_new, simpdat, Time, i) - 
        likelihood_i(beta_j2, lambda_j2, X_new, simpdat, Time, i)) / (2*h);
    }
  }
  
  double W = arma::sum(Dpl.col(p-1));
  arma::mat Info = arma::trans(Dpl) * Dpl;
  
  ret["W"] = W;
  ret["Dpl"] = Dpl;
  ret["Info"] = Info/n;
  
  return ret;
}

