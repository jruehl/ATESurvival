// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::vec EBS(const std::vector< arma::mat >& Z,
              const arma::vec &event,
              const arma::vec& time,
              const NumericVector& t,
              const std::vector< arma::vec >& beta,
              const std::vector< int >& index_A,
              const Rcpp::Nullable<Rcpp::NumericVector>& factor_init = R_NilValue,
              const Rcpp::Nullable<Rcpp::NumericVector>& trunc_time_init = R_NilValue,
              const int& cause = 1){
  
  int n = Z[cause-1].n_rows;
  int ncauses = Z.size();
  
  arma::vec factor(n, fill::value(0));
  if(factor_init.isNotNull()){
    factor = Rcpp::as<arma::vec>(factor_init);
  }else{
    // find number of duplicates of each observation
    int iter = 0;
    while(iter < n){
      int m_iter = 1;
      for(int iter2=iter+1; iter2<n; iter2++){
        if(time[iter] == time[iter2]){
          m_iter++;
        }else{
          break;
        }
      }
      factor[iter] = m_iter;
      iter += m_iter;
    }
  }
  arma::uvec nodup = find(factor > 0);
  
  arma::uvec eval_ind = find((event != 0) && (factor > 0) && (time <= max(t)));
  
  // left-truncation
  arma::vec trunc_time = zeros(n);
  if(trunc_time_init.isNotNull()){
    trunc_time = Rcpp::as<arma::vec>(trunc_time_init);
  }
  
  // dN, Y, S0, Z_bar & omega
  std::vector< arma::vec > dN(ncauses);
  arma::vec Y(eval_ind.size());
  std::vector < arma::mat > Z_A1(ncauses);
  std::vector < arma::mat > Z_A0(ncauses);
  std::vector< arma::vec > exp_Zb(ncauses);
  std::vector< arma::vec > exp_Zb_A1(ncauses);
  std::vector< arma::vec > exp_Zb_A0(ncauses); 
  std::vector< arma::vec > S0(ncauses);
  for(int c=0; c<ncauses; c++){
    dN[c] = factor.elem(eval_ind) % ((trunc_time.elem(eval_ind) < time.elem(eval_ind)) && (event.elem(eval_ind) == c+1));
    Z_A1[c] = Z[c];
    Z_A0[c] = Z[c];
    if(index_A[c] > 0){
      Z_A1[c].col(index_A[c]-1) = ones(n);
      Z_A0[c].col(index_A[c]-1) = zeros(n);
    }
    exp_Zb[c] = exp(Z[c].rows(nodup) * beta[c]);
    exp_Zb_A1[c] = exp(Z_A1[c] * beta[c]);
    exp_Zb_A0[c] = exp(Z_A0[c] * beta[c]);
    S0[c].set_size(eval_ind.size());
  }
  for(int j=0; j < int(eval_ind.size()); j++){
    arma::uvec Y_i = find((trunc_time.elem(nodup) < arma::vec(time.elem(eval_ind))[j]) && (arma::vec(time.elem(eval_ind))[j] <= time.elem(nodup)));
    for(int c=0; c<ncauses; c++){
      S0[c][j] = 1/double(n) * sum(exp_Zb[c].elem(Y_i) % arma::mat(factor.elem(nodup)).elem(Y_i));
    }
  }
  
  // estimated ATE
  arma::vec S_A1 = ones(n);
  arma::vec S_A0 = ones(n);
  arma::vec F_A1 = zeros(n);
  arma::vec F_A0 = zeros(n);
  int counter = 0;
  arma::vec ATE(t.size(), fill::value(0.0));
  for(int s=1; s <= int(eval_ind.size()); s++){
    while(t[counter] < time[eval_ind[s-1]]){
      ATE[counter] = sum(F_A1 - F_A0)/double(n);
      counter++;
    }
    if(dN[cause-1][s-1] > 0){
      F_A1 += S_A1 % (dN[cause-1][s-1] / (double(n)*S0[cause-1][s-1]) * exp_Zb_A1[cause-1]);
      F_A0 += S_A0 % (dN[cause-1][s-1] / (double(n)*S0[cause-1][s-1]) * exp_Zb_A0[cause-1]);
      F_A1 = min(ones(n), max(zeros(n), F_A1));
      F_A0 = min(ones(n), max(zeros(n), F_A0));
    }
    if(s == int(eval_ind.size())){
      for(int u=counter; u < t.size(); u++){
        ATE[u] = sum(F_A1 - F_A0)/double(n);
      }
      break;
    }
    if((time[eval_ind[s-1]] <= t[counter]) && (time[eval_ind[s]] > t[counter])){
      ATE[counter] = sum(F_A1 - F_A0)/double(n);
      counter++;
      if(counter == t.size()) break;
    }
    arma::vec dLambda_A1_all(n, fill::value(0.0));
    arma::vec dLambda_A0_all(n, fill::value(0.0));
    for(int c=0; c<ncauses; c++){
      dLambda_A1_all += dN[c][s-1] / (double(n)*S0[c][s-1]) * exp_Zb_A1[c];
      dLambda_A0_all += dN[c][s-1] / (double(n)*S0[c][s-1]) * exp_Zb_A0[c];
    }
    S_A1 %= (1.0 - dLambda_A1_all);
    S_A0 %= (1.0 - dLambda_A0_all);
    S_A1 = min(ones(n), max(zeros(n), S_A1));
    S_A0 = min(ones(n), max(zeros(n), S_A0));
  }
  
  return(ATE);
  
}