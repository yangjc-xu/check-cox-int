#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

struct MstepResult{
  arma::vec beta;
  arma::vec lambda;
  arma::vec cumulative_lambda;
  arma::vec exp_beta_x;
};

struct EMResult{
  arma::vec beta;
  arma::vec lambda;
  arma::vec cumulative_lambda;
};

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

struct MstepResult Mstep( arma::vec beta, 
                                arma::mat w, 
                                arma::mat X, 
                                arma::mat simpdat,
                                arma::vec Time,
                                arma::vec exp_beta_x){
  long p = beta.size();
  long n = X.n_rows;
  long m = Time.size();
  arma::vec beta_new(p);
  arma::vec lambda_new(m);
  
  double temp0 = 0;
  arma::vec temp1(p); temp1.zeros();
  arma::mat temp2(p, p); temp2.zeros();
  
  arma::vec Score(p); Score.zeros();
  arma::mat InfMat(p, p); InfMat.zeros();
  arma::vec colsum_w(m); colsum_w.zeros();
  
  struct MstepResult parameter;
  
  long index_star = 0;
  for (long k = 0; k < m; ++k) {
    double t_k = Time(k);
    colsum_w(k) = arma::sum(w.col(k));
    arma::vec Score_X(p); Score_X.zeros();
    
    for (long j = index_star; j < n && simpdat(j, 3) >= t_k; ++j) {
      long id_j = simpdat(j,4);
      arma::rowvec X_jk = X.row(id_j-1);
      double temp = exp(arma::as_scalar(X_jk * beta));
      temp0 = temp0 + temp;
      temp1 = temp1 + temp * arma::trans(X_jk);
      temp2 = temp2 + temp * arma::trans(X_jk) * X_jk;
      index_star = index_star + 1;
    }
    
    if(temp0 > 0){
      for (long i = 0 ; i < n; ++i) {
        long id_i = simpdat(i,4);
        Score_X = Score_X + w(i, k) * arma::trans(X.row(id_i-1));
      }
      InfMat = InfMat + colsum_w(k) * (temp2/temp0 - (temp1 / temp0) * arma::trans(temp1 / temp0));
      Score = Score - colsum_w(k) * (temp1 / temp0) + Score_X;
    }
  }
  
  //update beta
  beta_new = beta + arma::solve(InfMat, Score);
  //update lambda
  arma::vec exp_beta_x_new = ComputeExpBetaX(beta_new, X);
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
    lambda_new(k) = num / den;
  }
  
  
  parameter.beta = beta_new;
  parameter.lambda = lambda_new;
  parameter.cumulative_lambda = CumulativeLambda(lambda_new);
  parameter.exp_beta_x = exp_beta_x_new;
  return parameter;
}



struct EMResult EM(arma::mat X, 
                         arma::mat simpdat, 
                         double eps,
                         arma::vec Time){
  //initial values of parameters
  long p = X.n_cols;
  long m = Time.size();
  arma::vec beta_hat(p); beta_hat.zeros();
  arma::vec lambda_hat(m); lambda_hat.ones(); lambda_hat = lambda_hat / m;
  arma::vec cumulative_lambda_hat(m); cumulative_lambda_hat = CumulativeLambda(lambda_hat);
  
  //arma::vec diff(1000); diff.zeros();
  
  //temporary variables
  arma::vec beta_old;
  arma::vec lambda_old;
  arma::vec cumulative_lambda_old;
  struct EMResult estimate;
  struct MstepResult Mstep_temp;
  double error = 100;
  int Iter = 0;
  int MaxIter = 5000;
  double delta = 0.01;
  
  arma::vec exp_beta_x_old = ComputeExpBetaX(beta_hat, X);
  while (error > eps && Iter < MaxIter) {
    beta_old = beta_hat;
    lambda_old = lambda_hat;
    cumulative_lambda_old = cumulative_lambda_hat;
    arma::mat w = Estep(exp_beta_x_old, lambda_hat, X, simpdat, Time);
    Mstep_temp = Mstep(beta_hat, w, X, simpdat, Time, exp_beta_x_old);
    beta_hat = Mstep_temp.beta;
    lambda_hat = Mstep_temp.lambda;
    cumulative_lambda_hat = Mstep_temp.cumulative_lambda;
    exp_beta_x_old = Mstep_temp.exp_beta_x;
    
    double beta_max = VectorAbsoluteDifferenceSafe(beta_hat, beta_old, delta);
    double cumulative_lambda_max = VectorAbsoluteDifferenceSafe(cumulative_lambda_hat, cumulative_lambda_old, delta);
    error = beta_max + cumulative_lambda_max;
    
    ++Iter;
    //if(Iter == 1){ std::cout << "\nAfter " << Iter  <<" iteration, beta=\n" << beta_hat << std::endl;}
    //if(Iter % 50 == 1){
    //  std::cout << "Iter = " << Iter << ", Error = " << error << std::endl;
    //}
  }
  
  if(error <= eps){
    std::cout << "Algorithm converged in " << Iter << " interations" << ", Error = " << error << std::endl;
  } else{
    std::cout << "Algorithm did not converge!\n" << std::endl;
  }
  
  estimate.beta = beta_hat;
  estimate.lambda = lambda_hat;
  estimate.cumulative_lambda = CumulativeLambda(lambda_hat);
  return estimate;
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
  int MaxIter = 5000;
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

double likelihood_beta(arma::vec beta, 
                       arma::vec lambda, 
                       arma::mat X,
                       arma::mat simpdat,
                       arma::vec Time){
  double pl = 0;
  long n = X.n_rows;
  arma::vec exp_beta_x = ComputeExpBetaX(beta, X);
  for(int i = 0; i < n; ++i){
    long id_i = simpdat(i,0);
    double L_i = simpdat(i,1);
    double R_i = simpdat(i,2);
    arma::uvec index1 = arma::find(Time <= L_i);
    arma::uvec index2 = arma::find(Time <= R_i);
    
    if(R_i == INFINITY){
      pl += - exp_beta_x(id_i-1) * arma::sum(lambda(index1));
    } else{
      double sum1 = exp(-exp_beta_x(id_i-1) * arma::sum(lambda(index1)));
      double sum2 = exp(-exp_beta_x(id_i-1) * arma::sum(lambda(index2)));
      pl += log(sum1 - sum2);
    }
  }
  return(pl);
}


// [[Rcpp::export]]
Rcpp::List unireg_indep_EM(arma::mat X, arma::mat simpdat, double eps, arma::vec Time){
  Rcpp::List ret;
  long p = X.n_cols;
  long n = X.n_rows;
  double h = 5/sqrt(n);
  double eps_h = eps;
  struct EMResult res = EM(X, simpdat, eps, Time);
  
  arma::vec pl_orig(n);
  arma::mat Dpl(n,p);
  arma::mat Cov1(p,p);
  arma::mat Info1(p,p);
  
  arma::vec maximizing_lambda = profile_likelihood(res.beta, res.lambda, X, simpdat, eps, Time);
  double PL = likelihood_beta(res.beta, maximizing_lambda, X, simpdat, Time);
  
  //first order
  for(long i = 0; i < n; ++i){
    pl_orig(i) = likelihood_i(res.beta, maximizing_lambda, X, simpdat, Time, i);
  }
  
  for(long j = 0; j < p; ++j){
    arma::vec beta_j = res.beta;
    beta_j(j) = beta_j(j) + h;
    arma::vec lambda_j = profile_likelihood(beta_j, res.lambda, X, simpdat, eps_h, Time);
    for(long i = 0; i < n; ++i){
      Dpl(i,j) = (likelihood_i(beta_j, lambda_j, X, simpdat, Time, i) - pl_orig(i)) / h;
    }
  }
  
  Info1 = arma::trans(Dpl) * Dpl;
  Cov1 = arma::inv(Info1); 
  
  ret["beta"] = res.beta;
  ret["lambda"] = res.lambda;
  ret["cumulative_lambda"] = res.cumulative_lambda;
  ret["eff_score"] = Dpl;
  ret["Cov1"] = Cov1;
  ret["I1"] = Info1/n;
  ret["PL"] = PL;
  
  return ret;
}

