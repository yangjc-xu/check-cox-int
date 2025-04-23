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

arma::mat Estep(arma::vec beta, 
                arma::vec lambda, 
                arma::mat X_dep, 
                arma::mat simpdat,
                arma::vec Time){
  long n = simpdat.n_rows;
  long m = Time.size();
  arma::mat w(n, m); w.zeros();
  // use sorted simple data
  for (long i = 1; i < n+1; ++i) {
    double L_i = simpdat(i-1,1);
    double R_i = simpdat(i-1,2);
    long id_i = simpdat(i-1,4);
    
    for (long k = 1; k < m+1; ++k) {
      double t_k = Time(k-1);
      if (R_i < INFINITY && t_k > L_i && t_k <= R_i) {
        arma::rowvec X_ik = X_dep.row((id_i-1)*m+k-1);
        arma::uvec index_t = arma::find(Time > L_i && Time <= R_i);
        
        double num = lambda(k-1) * exp(arma::as_scalar(X_ik * beta));
        double sum = 0;
        for (long l = 0; l < index_t.size() ; ++l) {
          arma::rowvec X_il = X_dep.row((id_i-1)*m+index_t(l));
          sum = sum + lambda(index_t(l)) * exp(arma::as_scalar(X_il * beta));
        }
        double den = 1 - exp(-sum);
        w(i-1,k-1) = num / den;
      }
    }
  }
  return w;
}

arma::mat Estep_indep(arma::vec exp_beta_x, 
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

arma::vec profile_likelihood(arma::vec beta, arma::vec lambda_initial, arma::mat X_dep, 
                             arma::mat simpdat, double eps, arma::vec Time){
  long m = Time.size();
  long n = simpdat.n_rows;
  arma::vec lambda_hat(m); lambda_hat = lambda_initial;
  arma::vec lambda_old; lambda_old = lambda_initial;
  
  double error = 100;
  int Iter = 0;
  int MaxIter = 1000;
  double delta = 0.01;
  
  while (error > eps && Iter < MaxIter){
    arma::mat w = Estep(beta, lambda_old, X_dep, simpdat, Time);
    arma::rowvec colsum_w = arma::sum(w, 0);
    //update lambda
    
    for (long p = 0; p < m; ++p) {
      double den = 0;
      for (long i = 0; i < n && simpdat(i, 3) >= Time(p); ++i) {
        long id_i_new = simpdat(i,4);
        arma::rowvec X_ip = X_dep.row((id_i_new-1)*m+p);
        den = den + exp(arma::as_scalar(X_ip * beta));
      }
      lambda_hat(p) = colsum_w(p) / den;
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

arma::vec profile_likelihood_indep(arma::vec beta, 
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
    arma::mat w = Estep_indep(exp_beta_x_new, lambda_old, X, simpdat, Time);
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
                    arma::mat X_dep,
                    arma::vec exp_beta_x,
                    arma::mat simpdat, 
                    arma::vec Time, 
                    long i){
  double pl = 0;
  long m = Time.size();
  long p = beta.n_elem;
  long id_i = simpdat(i,4);
  double L_i = simpdat(i,1);
  double R_i = simpdat(i,2);
  arma::uvec index1 = arma::find(Time <= L_i);
  arma::uvec index2 = arma::find(Time <= R_i);
  
  if(R_i == INFINITY){
    
    arma::uvec X_index1(index1.size());
    for(long k = 0; k < index1.size(); ++k){
      X_index1(k) = (id_i-1) * m + index1(k);
    }
    arma::mat X_sub = X_dep.rows(X_index1);
    pl = - exp_beta_x(id_i-1) * arma::sum(exp(X_sub.col(p-1) * beta(p-1)) % lambda(index1));
    
  } else{
    
    arma::uvec X_index1(index1.size());
    arma::uvec X_index2(index2.size());
    for(long k = 0; k < index1.size(); ++k){
      X_index1(k) = (id_i-1) * m + index1(k);
    }
    for(long k = 0; k < index2.size(); ++k){
      X_index2(k) = (id_i-1) * m + index2(k);
    }
    arma::mat X_sub1 = X_dep.rows(X_index1);
    arma::mat X_sub2 = X_dep.rows(X_index2);
    double sum1 = exp_beta_x(id_i-1) * arma::sum(exp(X_sub1.col(p-1) * beta(p-1)) % lambda(index1));
    double sum2 = exp_beta_x(id_i-1) * arma::sum(exp(X_sub2.col(p-1) * beta(p-1)) % lambda(index2));
    
    pl = log(exp(-sum1) - exp(-sum2));
  }
  return(pl);
}

double likelihood_i_indep(arma::vec beta, 
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
Rcpp::List get_W(arma::vec beta_new, arma::mat X_dep_new, arma::mat X_indep, arma::vec lambda_old,
                 arma::mat simpdat, double eps, arma::vec Time, double constant){
  Rcpp::List ret;
  long n = simpdat.n_rows;
  long p = X_dep_new.n_cols;
  double h = constant/sqrt(n);
  
  arma::mat Dpl(n, p);
  
  //first order
  for(long j = 0; j < p; ++j){
    arma::vec beta_j1 = beta_new;
    arma::vec beta_j2 = beta_new;
    beta_j1(j) = beta_j1(j) + h;
    beta_j2(j) = beta_j2(j) - h;
    
    if(j < p-1){
      arma::vec lambda_j1 = profile_likelihood_indep(beta_j1.subvec(0,p-2), lambda_old, X_indep, simpdat, eps, Time);
      arma::vec lambda_j2 = profile_likelihood_indep(beta_j2.subvec(0,p-2), lambda_old, X_indep, simpdat, eps, Time);
      for(long i = 0; i < n; ++i){
        Dpl(i,j) = (likelihood_i_indep(beta_j1.subvec(0,p-2), lambda_j1, X_indep, simpdat, Time, i) - 
          likelihood_i_indep(beta_j2.subvec(0,p-2), lambda_j2, X_indep, simpdat, Time, i)) / (2*h);
      }
      std::cout << "Completed for j = " << j << "\n" << std::endl; 
    } else{
      arma::vec lambda_j1 = profile_likelihood(beta_j1, lambda_old, X_dep_new, simpdat, eps, Time);
      arma::vec lambda_j2 = profile_likelihood(beta_j2, lambda_old, X_dep_new, simpdat, eps, Time);
      arma::vec exp_beta_X1 = ComputeExpBetaX(beta_j1.subvec(0, p-2), X_indep);
      arma::vec exp_beta_X2 = ComputeExpBetaX(beta_j2.subvec(0, p-2), X_indep);
      for(long i = 0; i < n; ++i){
        Dpl(i,j) = (likelihood_i(beta_j1, lambda_j1, X_dep_new, exp_beta_X1, simpdat, Time, i) - 
          likelihood_i(beta_j2, lambda_j2, X_dep_new, exp_beta_X2, simpdat, Time, i)) / (2*h);
      }
      std::cout << "Completed for j = " << j << "\n" << std::endl; 
    }
  }
  
  double W = arma::sum(Dpl.col(p-1));
  arma::mat Info = arma::trans(Dpl) * Dpl;
  
  ret["W"] = W;
  ret["Dpl"] = Dpl;
  ret["Info"] = Info/n;
  
  return ret;
}
