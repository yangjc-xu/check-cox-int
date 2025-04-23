#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

struct MstepResult{
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

struct MstepResult Mstep(arma::vec beta, arma::mat w, arma::mat X_dep, arma::mat simpdat, arma::vec Time){
  long p = beta.size();
  long n = simpdat.n_rows;
  long m = Time.size();
  arma::vec beta_new(p);
  arma::vec lambda_new(m);
  
  arma::vec Score(p); Score.zeros();
  arma::mat InfMat(p, p); InfMat.zeros();
  arma::vec colsum_w(m); colsum_w.zeros();
  
  struct MstepResult parameter;
  
  for (long k = 0; k < m; ++k) {
    double t_k = Time(k);
    colsum_w(k) = arma::sum(w.col(k));
    arma::vec Score_X(p); Score_X.zeros();
    //arma::vec sum_tmp(p); sum_tmp.zeros();
    //arma::vec sum_tmp1(p); sum_tmp1.zeros();
    //arma::vec sum_tmp2(p); sum_tmp2.zeros();
    //arma::vec quotient_k(p); quotient_k.zeros();
    
    double temp0 = 0;
    arma::vec temp1(p); temp1.zeros();
    arma::mat temp2(p, p); temp2.zeros();
    
    for (long j = 0; j < n; ++j) {
      if(simpdat(j, 3) >= t_k){
        long id_j = simpdat(j, 4);
        arma::rowvec X_jk = X_dep.row((id_j-1)*m+k);
        double temp = exp(arma::as_scalar(X_jk * beta));
        temp0 = temp0 + temp;
        temp1 = temp1 + temp * arma::trans(X_jk);
        temp2 = temp2 + temp * arma::trans(X_jk) * X_jk;
      }
    }
    
    if(temp0 > 0){
      for (long i = 0 ; i < n; ++i) {
        if(simpdat(i, 3) >= t_k){
          long id_i = simpdat(i,4);
          Score_X = Score_X + w(i, k) * arma::trans(X_dep.row((id_i-1)*m+k));
        }
      }
      InfMat = InfMat + colsum_w(k) * (temp2/temp0 - (temp1 / temp0) * arma::trans(temp1 / temp0));
      Score = Score - colsum_w(k) * (temp1 / temp0) + Score_X;
      //sum_tmp = - colsum_w(k) * (temp1 / temp0) + Score_X;
      //sum_tmp1 = Score_X;
      //sum_tmp2 = colsum_w(k) * (temp1 / temp0);
      //quotient_k = temp1 / temp0;
      //if(k>444){
      //  std::cout << "k = " << k << std::endl;
      //  std::cout << "sum_tmp = " << sum_tmp(3) << std::endl;
      //  std::cout << "sum_tmp1 = " << sum_tmp1(3) << std::endl;
      //  std::cout << "sum_tmp2 = " << sum_tmp2(3) << std::endl;
      //  std::cout << "colsum_w = " << colsum_w(k) << std::endl;
      //  std::cout << "quotient_k = " << quotient_k(3) << std::endl;
      //  std::cout << "temp0 = " << temp0 << std::endl;
      //  std::cout << "temp1 = " << temp1 << std::endl;
      //}
    }
  }
  
  //update beta
  beta_new = beta + 0.5 * arma::solve(InfMat, Score);
  //std::cout << "\nScore = \n" << Score << std::endl;
  //std::cout << "\nInfMat = \n" << InfMat << std::endl;
  //update lambda
  for (long p = 0; p < m; ++p) {
    double den = 0;
    for (long i = 0; i < n && simpdat(i, 3) >= Time(p); ++i) {
      long id_i_new = simpdat(i,4);
      arma::rowvec X_ip = X_dep.row((id_i_new-1)*m+p);
      den = den + exp(arma::as_scalar(X_ip * beta_new));
    }
    lambda_new(p) = colsum_w(p) / den;
  }
  
  parameter.beta = beta_new;
  parameter.lambda = lambda_new;
  parameter.cumulative_lambda = CumulativeLambda(lambda_new);
  return parameter;
}

struct MstepResult EM(arma::mat X_dep, arma::mat simpdat, double eps, arma::vec Time){
  //initial values of parameters
  long p = X_dep.n_cols;
  long m = Time.size();
  arma::vec beta_hat(p); beta_hat.zeros();
  arma::vec lambda_hat(m); lambda_hat.ones(); lambda_hat = lambda_hat / m;
  arma::vec cumulative_lambda_hat(m); cumulative_lambda_hat = CumulativeLambda(lambda_hat);
  
  //temporary variables
  arma::vec beta_old;
  arma::vec lambda_old;
  arma::vec cumulative_lambda_old;
  struct MstepResult estimate;
  struct MstepResult Mstep_temp;
  double error = 100;
  int Iter = 0;
  int MaxIter = 5000;
  double delta = 0.01;
  
  while (error > eps && Iter < MaxIter) {
    beta_old = beta_hat;
    lambda_old = lambda_hat;
    cumulative_lambda_old = cumulative_lambda_hat;
    arma::mat w = Estep(beta_hat, lambda_hat, X_dep, simpdat, Time);
    Mstep_temp = Mstep(beta_hat, w, X_dep, simpdat, Time);
    beta_hat = Mstep_temp.beta;
    lambda_hat = Mstep_temp.lambda;
    cumulative_lambda_hat = Mstep_temp.cumulative_lambda;
    
    //double beta_max = std::abs(beta_hat-beta_old).max();
    //double lambda_max = std::abs(lambda_hat-lambda_old).max();
    //if(beta_max > lambda_max){
    //  error = beta_max;
    //} else{
    //  error= lambda_max;
    //}
    
    double beta_max = VectorAbsoluteDifferenceSafe(beta_hat, beta_old, delta);
    double cumulative_lambda_max = VectorAbsoluteDifferenceSafe(cumulative_lambda_hat, cumulative_lambda_old, delta);
    error = beta_max + cumulative_lambda_max;
    
    ++Iter;
    //if(Iter == 1){
      //std::cout << "\nAfter first iteration, beta=\n" << beta_hat << std::endl;
      //std::cout << "\nAfter first iteration, sum_w=\n" << arma::sum(arma::sum(w)) << std::endl;
    //}
    if(Iter % 2 == 1){
     std::cout << "Iter = " << Iter << ", Error = " << error << std::endl;
    }
  }
  
  if(error <= eps){
    std::cout << "Algorithm converged in " << Iter << " interations" << ", Error = " << error << "\n" << std::endl;
  } else{
    std::cout << "Algorithm did not converge!\n" << std::endl;
  }
  
  estimate.beta = beta_hat;
  estimate.lambda = lambda_hat;
  estimate.cumulative_lambda = cumulative_lambda_hat;
  return estimate;
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

double likelihood_i(arma::vec beta, 
                    arma::vec lambda, 
                    arma::mat X_dep,
                    arma::mat simpdat, 
                    arma::vec Time, 
                    long i){
  double pl = 0;
  long m = Time.size();
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
    pl = - arma::sum(exp(X_sub * beta) % lambda(index1));
    
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
    double sum1 = arma::sum(exp(X_sub1 * beta) % lambda(index1));
    double sum2 = arma::sum(exp(X_sub2 * beta) % lambda(index2));
    
    pl = log(exp(-sum1) - exp(-sum2));
  }
  return(pl);
}

double likelihood_beta(arma::vec beta, 
                       arma::vec lambda,
                       arma::mat X_dep,
                       arma::mat simpdat,
                       arma::vec Time){
  double pl = 0;
  long n = simpdat.n_rows;
  long m = Time.size();
  arma::vec exp_beta_x = ComputeExpBetaX(beta, X_dep);
  for(int i = 0; i < n; ++i){
    long id_i = simpdat(i,4);
    double L_i = simpdat(i,1);
    double R_i = simpdat(i,2);
    arma::uvec index1 = arma::find(Time <= L_i);
    arma::uvec index2 = arma::find(Time <= R_i);
    
    if(R_i == INFINITY){
      for(int k = 0; k < index1.size(); ++k){
        pl += - lambda(index1(k)) * exp_beta_x((id_i-1)*m + index1(k));
      }
    } else{
      double sum1 = 0;
      double sum2 = 0;
      for(int k = 0; k < index1.size(); ++k){
        sum1 += - lambda(index1(k)) * exp_beta_x((id_i-1)*m + index1(k));
      }
      for(int k = 0; k < index2.size(); ++k){
        sum2 += - lambda(index2(k)) * exp_beta_x((id_i-1)*m + index2(k));
      }
      pl += log(exp(sum1) - exp(sum2));
    }
  }
  return(pl);
}

// [[Rcpp::export]]
Rcpp::List unireg_dep_EM(arma::mat X_dep, arma::mat simpdat, double eps, arma::vec Time){
  Rcpp::List ret;
  long p = X_dep.n_cols;
  long n = simpdat.n_rows;
  double h = 5/sqrt(n);
  double eps_h = eps;
  struct MstepResult res = EM(X_dep, simpdat, eps, Time);
  
  arma::vec pl_orig(n);
  arma::mat Dpl(n,p);
  arma::mat Cov1(p,p);
  arma::mat Info1(p,p);
  
  double pl;
  arma::vec pl1(p);
  arma::mat pl2(p,p); pl2.zeros();
  arma::mat Cov(p,p);
  arma::mat Info(p,p);
  
  arma::vec maximizing_lambda = profile_likelihood(res.beta, res.lambda, X_dep, simpdat, eps_h, Time);
  pl = likelihood_beta(res.beta, maximizing_lambda, X_dep, simpdat, Time);
  
  //first order
  for(long i = 0; i < n; ++i){
    pl_orig(i) = likelihood_i(res.beta, maximizing_lambda, X_dep, simpdat, Time, i);
  }
  
  for(long j = 0; j < p; ++j){
    arma::vec beta_j = res.beta;
    beta_j(j) = beta_j(j) + h;
    arma::vec lambda_j = profile_likelihood(beta_j, res.lambda, X_dep, simpdat, eps_h, Time);
    for(long i = 0; i < n; ++i){
      Dpl(i,j) = (likelihood_i(beta_j, lambda_j, X_dep, simpdat, Time, i) - pl_orig(i)) / h;
    }
  }
  
  Info1 = arma::trans(Dpl) * Dpl;
  Cov1 = arma::inv(Info1);
  
  // second order
  // for(long j = 0; j < p; ++j){
  //   arma::vec beta_j = res.beta;
  //   beta_j(j) = beta_j(j) + h;
  //   arma::vec lambda_j = profile_likelihood(beta_j, res.lambda, X_dep, simpdat, eps_h, Time);
  //   pl1(j) = likelihood_beta(beta_j, lambda_j, X_dep, simpdat, Time);
  // }
  // 
  // for(long i = 0; i < p; ++i){
  //   for(long j = 0; j <= i; ++j){
  //     arma::vec beta_ij = res.beta;
  //     beta_ij(i) = beta_ij(i) + h;
  //     beta_ij(j) = beta_ij(j) + h;
  //     arma::vec lambda_ij = profile_likelihood(beta_ij, res.lambda, X_dep, simpdat, eps_h, Time);
  //     pl2(i,j) = likelihood_beta(beta_ij, lambda_ij, X_dep, simpdat, Time);
  //     Info(i,j) = (pl - pl1(i) - pl1(j) + pl2(i,j)) / (h * h);
  //     Info(j,i) = Info(i,j);
  //   }
  // }
  // 
  // Cov = -arma::inv(Info);
  
  ret["beta"] = res.beta;
  ret["lambda"] = res.lambda;
  ret["cumulative_lambda"] = res.cumulative_lambda;
  ret["eff_score"] = Dpl;
  ret["Cov1"] = Cov1;
  ret["I1"] = Info1/n;
  //ret["Cov2"] = Cov;
  //ret["I2"] = -Info/n;
  ret["PL"] = pl;

  return ret;
}
