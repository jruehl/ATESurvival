// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppClock)]]
#include <RcppClock.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List ATE_IF_WBS(const std::vector< arma::mat >& Z,
                const arma::vec &event, 
                const arma::vec& time, 
                const NumericVector& t, 
                const std::vector< arma::vec >& beta,
                const std::vector< int >& index_A,
                const Rcpp::Nullable<Rcpp::NumericVector>& trunc_time_init = R_NilValue, 
                const int& cause = 1,
                const bool& IF = true,
                const bool& WBS = true,
                const int& bs_iter = 1e4, 
                const double& level = 0.95,
                const bool& arms = false){
  
  Rcpp::Clock clock;
  clock.tick("ATE");
  
  int n = Z[cause-1].n_rows;
  int ncauses = Z.size();
  
  arma::uvec event_ind = find((event != 0));
  arma::uvec eval_ind = find((event != 0) && (time <= max(t)));
  
  // left-truncation
  arma::vec trunc_time = zeros(n);
  if(trunc_time_init.isNotNull()){
    trunc_time = Rcpp::as<arma::vec>(trunc_time_init);
  }
  
  // dN, Y, S0, Z_bar & omega
  std::vector< arma::vec > dN(ncauses);
  arma::vec Y(event_ind.size());
  std::vector < arma::mat > Z_A1(ncauses);
  std::vector < arma::mat > Z_A0(ncauses);
  std::vector< arma::vec > exp_Zb(ncauses);
  std::vector< arma::vec > exp_Zb_A1(ncauses);
  std::vector< arma::vec > exp_Zb_A0(ncauses); 
  std::vector< arma::vec > S0(ncauses);
  std::vector< arma::mat > Z_bar(ncauses);
  std::vector< arma::mat > omega(ncauses);
  std::vector< arma::mat > omega_inv(ncauses);
  for(int c=0; c<ncauses; c++){
    dN[c] = ones(event_ind.size()) % ((trunc_time.elem(event_ind) < time.elem(event_ind)) && (event.elem(event_ind) == c+1));
    Z_A1[c] = Z[c];
    Z_A0[c] = Z[c];
    if(index_A[c] > 0){
      Z_A1[c].col(index_A[c]-1) = ones(n);
      Z_A0[c].col(index_A[c]-1) = zeros(n);
    }
    exp_Zb[c] = exp(Z[c] * beta[c]);
    exp_Zb_A1[c] = exp(Z_A1[c] * beta[c]);
    exp_Zb_A0[c] = exp(Z_A0[c] * beta[c]);
    S0[c].set_size(event_ind.size());
    if(IF || WBS){
      Z_bar[c].set_size(event_ind.size(), Z[c].n_cols);
      omega[c].set_size(Z[c].n_cols, Z[c].n_cols);
      omega[c].fill(0.0);
    }
  }
  for(int j=0; j < int(event_ind.size()); j++){
    arma::uvec Y_i = find((trunc_time < arma::vec(time.elem(event_ind))[j]) && (arma::vec(time.elem(event_ind))[j] <= time));
    for(int c=0; c<ncauses; c++){
      S0[c][j] = 1/double(n) * sum(exp_Zb[c].elem(Y_i));
      if(IF || WBS){
        Z_bar[c].row(j) = 1/double(n) * (exp_Zb[c].elem(Y_i)).t() * Z[c].rows(Y_i) / S0[c][j];
        arma::mat exp_Zb_Z = Z[c].rows(Y_i);
        exp_Zb_Z.each_col() %= exp_Zb[c].elem(Y_i);
        arma::mat S2 = 1/double(n) * (exp_Zb_Z).t() * Z[c].rows(Y_i);
        omega[c] += 1/double(n) * (S2/S0[c][j] - Z_bar[c].row(j).t() * Z_bar[c].row(j)) * dN[c][j];
      }
    }
    Y[j] = Y_i.size();
  }
  if(IF || WBS){
    for(int c=0; c<ncauses; c++){
      omega_inv[c] = inv_sympd(omega[c]);
    }
  }
  
  // dLambda
  std::vector < arma::mat > dLambda_A1(ncauses);
  std::vector < arma::mat > dLambda_A0(ncauses);
  for(int c=0; c<ncauses; c++){
    dLambda_A1[c].set_size(n, eval_ind.size());
    dLambda_A0[c].set_size(n, eval_ind.size());
    for(int s=0; s<int(eval_ind.size()); s++){
      dLambda_A1[c].col(s) = dN[c][s] / (double(n)*S0[c][s]) * exp_Zb_A1[c];
      dLambda_A0[c].col(s) = dN[c][s] / (double(n)*S0[c][s]) * exp_Zb_A0[c];
    }
  }
  
  // S & F
  arma::mat S_A1(n, int(eval_ind.size()+1));
  arma::mat S_A0(n, int(eval_ind.size()+1));
  arma::mat F_A1(n, int(eval_ind.size()+1));
  arma::mat F_A0(n, int(eval_ind.size()+1));
  S_A1.col(0) = ones(n);
  S_A0.col(0) = ones(n);
  F_A1.col(0) = zeros(n);
  F_A0.col(0) = zeros(n);
  for(int s=1; s <= int(eval_ind.size()); s++){
    arma::vec dLambda_A1_all(n, fill::value(0.0));
    arma::vec dLambda_A0_all(n, fill::value(0.0));
    for(int c=0; c<ncauses; c++){
      dLambda_A1_all += dLambda_A1[c].col(s-1);
      dLambda_A0_all += dLambda_A0[c].col(s-1);
    }
    S_A1.col(s) = min(ones(n), max(zeros(n), S_A1.col(s-1) % (1.0 - dLambda_A1_all)));
    S_A0.col(s) = min(ones(n), max(zeros(n), S_A0.col(s-1) % (1.0 - dLambda_A0_all)));
    //S_A1.col(s) = S_A1.col(s-1) % (1.0 - dLambda_A1_all);
    //S_A0.col(s) = S_A0.col(s-1) % (1.0 - dLambda_A0_all);
    arma::vec dF_A1(n, fill::value(0.0));
    arma::vec dF_A0(n, fill::value(0.0));
    if(dN[cause-1][s-1] > 0){
      dF_A1 = S_A1.col(s-1) % dLambda_A1[cause-1].col(s-1);
      dF_A0 = S_A0.col(s-1) % dLambda_A0[cause-1].col(s-1);
    }
    F_A1.col(s) = min(ones(n), max(zeros(n), F_A1.col(s-1) + dF_A1));
    F_A0.col(s) = min(ones(n), max(zeros(n), F_A0.col(s-1) + dF_A0));
    //F_A1.col(s) = F_A1.col(s-1) + dF_A1;
    //F_A0.col(s) = F_A0.col(s-1) + dF_A0;
  }
  
  
  // mean F
  arma::vec F_A1_mean(t.size(), fill::value(0.0));
  arma::vec F_A0_mean(t.size(), fill::value(0.0));
  for(int s=0; s<t.size(); s++){
    if(t[s] >= min(time)){
      arma::uvec eval_ind_t = find((event != 0) && (time <= t[s]));
      F_A1_mean[s] = sum(F_A1.col(eval_ind_t.size()))/double(n);
      F_A0_mean[s] = sum(F_A0.col(eval_ind_t.size()))/double(n);
    }
  }
  
  clock.tock("ATE");
  
  
  clock.tick("IF");
  arma::mat res_CI_IF;
  arma::mat res_CB_IF;
  
  if(IF){
    
    // influence functions
    
    // IF_beta
    std::vector< arma::mat > IF_beta(ncauses);
    for(int c=0; c<ncauses; c++){
      IF_beta[c].set_size(n, Z[c].n_cols);
      IF_beta[c].fill(0.0);
      for(int i=0; i<n; i++){
        arma::uvec i_event_ind = find(i == event_ind);
        if(i_event_ind.size() > 0){
          IF_beta[c].row(i) = (Z[c].row(i) - Z_bar[c].rows(i_event_ind)) * as_scalar(dN[c].elem(i_event_ind)) * 1/double(n)*omega_inv[c].t();
        }
        arma::rowvec IF_beta_2(Z[c].n_cols, fill::value(0.0));
        int j=0;
        while(time[event_ind[j]] <= time[i]){
          IF_beta_2 += (Z[c].row(i) - Z_bar[c].row(j)) * dN[c][j] / (double(n)*S0[c][j]);
          j++;
          if(j == int(event_ind.size())) break;
        }
        IF_beta[c].row(i) -= exp_Zb[c][i] * IF_beta_2 * 1/double(n)*omega_inv[c].t();
      }
    }
    
    // dIF_Lambda_0
    std::vector< arma::mat > dIF_Lambda_0(ncauses);
    for(int c=0; c<ncauses; c++){
      dIF_Lambda_0[c].set_size(n, int(eval_ind.size()+1));
      dIF_Lambda_0[c].fill(0.0);
      for(int s=1; s <= int(eval_ind.size()); s++){
        dIF_Lambda_0[c].col(s) = -IF_beta[c] * Z_bar[c].row(s-1).t() * dN[c][s-1] / (double(n)*S0[c][s-1]);
        if(dN[c][s-1] > 0){
          dIF_Lambda_0[c].col(s) -= (time[eval_ind[s-1]] <= time) % exp_Zb[c] * dN[c][s-1] / pow(double(n)*S0[c][s-1], 2);
          dIF_Lambda_0[c](eval_ind[s-1], s) += dN[c][s-1] / (double(n)*S0[c][s-1]);
        }
      }
    }
    
    // IF_F
    arma::mat IF_F_A1_mean(n, t.size(), fill::value(0.0));
    arma::mat IF_F_A0_mean(n, t.size(), fill::value(0.0));
    arma::vec IF_F_A1_mean_count(n, fill::value(0.0));
    arma::vec IF_F_A0_mean_count(n, fill::value(0.0));
    int counter = 0;
    for(int s=0; s < int(t.size()); s++){
      if(counter == int(eval_ind.size())){
        IF_F_A1_mean.col(s) = IF_F_A1_mean_count;
        IF_F_A0_mean.col(s) = IF_F_A0_mean_count;
        continue;
      }
      while(time[int(eval_ind[counter])] <= t[s]){
        if(dN[cause-1][counter] > 0){
          for(int c=0; c<ncauses; c++){
            IF_F_A1_mean_count -= 1/double(n) * sum(S_A1.col(counter) % exp_Zb_A1[c] % dLambda_A1[cause-1].col(counter)) * sum(dIF_Lambda_0[c].head_cols(counter+1), 1);
            IF_F_A0_mean_count -= 1/double(n) * sum(S_A0.col(counter) % exp_Zb_A0[c] % dLambda_A0[cause-1].col(counter)) * sum(dIF_Lambda_0[c].head_cols(counter+1), 1);
            if(counter > 0){
              arma::mat S_exp_Zb_dLambda_Z_A1 = Z_A1[c];
              arma::mat S_exp_Zb_dLambda_Z_A0 = Z_A0[c];
              S_exp_Zb_dLambda_Z_A1.each_col() %= S_A1.col(counter) % exp_Zb_A1[c] % dLambda_A1[cause-1].col(counter);
              S_exp_Zb_dLambda_Z_A0.each_col() %= S_A0.col(counter) % exp_Zb_A0[c] % dLambda_A0[cause-1].col(counter);
              double Lambda_0 = sum(dN[c].head(counter)/(double(n)*S0[c].head(counter)));
              IF_F_A1_mean_count -= 1/double(n) * IF_beta[c] * sum(S_exp_Zb_dLambda_Z_A1, 0).t() * Lambda_0;
              IF_F_A0_mean_count -= 1/double(n) * IF_beta[c] * sum(S_exp_Zb_dLambda_Z_A0, 0).t() * Lambda_0;
            }
          }
          IF_F_A1_mean_count += 1/double(n) * sum(S_A1.col(counter) % exp_Zb_A1[cause-1]) * dIF_Lambda_0[cause-1].col(counter+1);
          IF_F_A0_mean_count += 1/double(n) * sum(S_A0.col(counter) % exp_Zb_A0[cause-1]) * dIF_Lambda_0[cause-1].col(counter+1);
          arma::mat S_dLambda_A1 = Z_A1[cause-1];
          arma::mat S_dLambda_A0 = Z_A0[cause-1];
          S_dLambda_A1.each_col() %= S_A1.col(counter) % dLambda_A1[cause-1].col(counter);
          S_dLambda_A0.each_col() %= S_A0.col(counter) % dLambda_A0[cause-1].col(counter);
          IF_F_A1_mean_count += 1/double(n) * IF_beta[cause-1] * sum(S_dLambda_A1, 0).t();
          IF_F_A0_mean_count += 1/double(n) * IF_beta[cause-1] * sum(S_dLambda_A0, 0).t();
        }
        counter++;
        if(counter == int(eval_ind.size())) break;
      }
      IF_F_A1_mean.col(s) = IF_F_A1_mean_count;
      IF_F_A0_mean.col(s) = IF_F_A0_mean_count;
    }
    
    // IF_ATE
    arma::mat IF_ATE(n, t.size(), fill::value(0.0));
    for(int s=0; s<t.size(); s++){
      if(t[s] >= min(time)){
        arma::uvec eval_ind_t = find((event != 0) && (time <= t[s]));
        IF_ATE.col(s) = 1/double(n)*F_A1.col(int(eval_ind_t.size())) - 1/double(n)*F_A0.col(int(eval_ind_t.size())) -
          (1/double(n)*F_A1_mean[s] - 1/double(n)*F_A0_mean[s]) * ones(n) + 
          IF_F_A1_mean.col(s) - IF_F_A0_mean.col(s);
      }
    }
    
    // SE_IF
    arma::vec SE_IF(t.size());
    for(int s=0; s<t.size(); s++){
      SE_IF[s] = sqrt(as_scalar(sum(pow(IF_ATE.col(s), 2))));
    }
    
    
    res_CI_IF = join_rows(max(F_A1_mean - F_A0_mean - SE_IF * R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0), -1 * ones(t.size())),
                          min(F_A1_mean - F_A0_mean + SE_IF * R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0), ones(t.size()))).t();
    
    // BS multipliers
    arma::mat G(bs_iter, n, fill::randn);
    arma::mat IF_ATE_by_SE_IF = IF_ATE.cols(find(SE_IF > 0));
    IF_ATE_by_SE_IF.each_row() /= SE_IF.elem(find(SE_IF > 0)).t();
    arma::vec q_IF = max(abs(G * IF_ATE_by_SE_IF), 1);
    res_CB_IF = join_rows(max(F_A1_mean - F_A0_mean - SE_IF * arma::vec(sort(q_IF))[round(bs_iter * level)], -1 * ones(t.size())),
                          min(F_A1_mean - F_A0_mean + SE_IF * arma::vec(sort(q_IF))[round(bs_iter * level)], ones(t.size()))).t();
    
  }
  
  clock.tock("IF");
  
  
  clock.tick("WBS");
  arma::mat res_CI_WBS;
  arma::mat res_CB_WBS;
  
  if(WBS){
    
    // BS multipliers
    arma::mat G_Lin(event_ind.size(), bs_iter, fill::randn);
    arma::mat G_Bey(event_ind.size(), bs_iter);
    arma::mat G_Weird(event_ind.size(), bs_iter);
    for(int j=0; j<bs_iter; j++){
      for(int i=0; i<int(event_ind.size()); i++){
        G_Bey(i,j) = R::rpois(1) - 1;
        G_Weird(i,j) = R::rbinom(Y[i], 1/Y[i]) - 1;
      }
    }
    
    // phi & psi
    arma::mat phi_A1(t.size(), Z[cause-1].n_cols);
    arma::mat phi_A0(t.size(), Z[cause-1].n_cols);
    arma::rowvec phi_A1_count(Z[cause-1].n_cols, fill::value(0.0));
    arma::rowvec phi_A0_count(Z[cause-1].n_cols, fill::value(0.0));
    std::vector < arma::mat > psi_A1(ncauses);
    std::vector < arma::mat > psi_A0(ncauses);
    std::vector < arma::rowvec > psi_A1_count_sum2(ncauses);
    std::vector < arma::rowvec > psi_A0_count_sum2(ncauses);
    for(int c=0; c<ncauses; c++){
      psi_A1[c].set_size(t.size(), Z[c].n_cols);
      psi_A0[c].set_size(t.size(), Z[c].n_cols);
      psi_A1_count_sum2[c].set_size(Z[c].n_cols);
      psi_A1_count_sum2[c].fill(0.0);
      psi_A0_count_sum2[c].set_size(Z[c].n_cols);
      psi_A0_count_sum2[c].fill(0.0);
    }
    int counter = 0;
    for(int s=0; s<t.size(); s++){
      if(counter < int(eval_ind.size())){
        while(time[int(eval_ind[counter])] <= t[s]){
          for(int c=0; c<ncauses; c++){
            if(dN[c][counter] == 1){
              arma::mat Z_minus_Z_bar_A1 = Z_A1[c];
              arma::mat Z_minus_Z_bar_A0 = Z_A0[c];
              Z_minus_Z_bar_A1.each_row() -= Z_bar[c].row(counter);
              Z_minus_Z_bar_A0.each_row() -= Z_bar[c].row(counter);
              if(c == cause-1){
                arma::mat phi_A1_count_j = Z_minus_Z_bar_A1;
                arma::mat phi_A0_count_j = Z_minus_Z_bar_A0;
                phi_A1_count_j.each_col() %= S_A1.col(counter) % dLambda_A1[cause-1].col(counter);
                phi_A0_count_j.each_col() %= S_A0.col(counter) % dLambda_A0[cause-1].col(counter);
                phi_A1_count += sum(phi_A1_count_j, 0);
                phi_A0_count += sum(phi_A0_count_j, 0);
              }
              arma::mat psi_A1_count_sum2_j = Z_minus_Z_bar_A1;
              arma::mat psi_A0_count_sum2_j = Z_minus_Z_bar_A0;
              psi_A1_count_sum2_j.each_col() %= F_A1.col(counter+1) % dLambda_A1[c].col(counter);
              psi_A0_count_sum2_j.each_col() %= F_A0.col(counter+1) % dLambda_A0[c].col(counter);
              psi_A1_count_sum2[c] += sum(psi_A1_count_sum2_j, 0);
              psi_A0_count_sum2[c] += sum(psi_A0_count_sum2_j, 0);
            }
          }
          counter++;
          if(counter == int(eval_ind.size())) break;
        }
      }
      std::vector < arma::rowvec > psi_A1_count_sum1(ncauses);
      std::vector < arma::rowvec > psi_A0_count_sum1(ncauses);
      for(int c=0; c<ncauses; c++){
        psi_A1_count_sum1[c].set_size(Z[c].n_cols);
        psi_A1_count_sum1[c].fill(0.0);
        psi_A0_count_sum1[c].set_size(Z[c].n_cols);
        psi_A0_count_sum1[c].fill(0.0);
      }
      if(t[s] >= min(time)){
        arma::uvec eval_ind_t = find((event != 0) && (time <= t[s]));
        for(int r=0; r<int(eval_ind_t.size()); r++){
          for(int c=0; c<ncauses; c++){
            if(dN[c][r] == 1){
              arma::mat psi_A1_count_sum1_j = Z_A1[c];
              arma::mat psi_A0_count_sum1_j = Z_A0[c];
              psi_A1_count_sum1_j.each_row() -= Z_bar[c].row(r);
              psi_A0_count_sum1_j.each_row() -= Z_bar[c].row(r);
              psi_A1_count_sum1_j.each_col() %= F_A1.col(eval_ind_t.size()) % dLambda_A1[c].col(r);
              psi_A0_count_sum1_j.each_col() %= F_A0.col(eval_ind_t.size()) % dLambda_A0[c].col(r);
              psi_A1_count_sum1[c] += sum(psi_A1_count_sum1_j, 0);
              psi_A0_count_sum1[c] += sum(psi_A0_count_sum1_j, 0);
            }
          }
        }
      }
      phi_A1.row(s) = 1/double(n) * phi_A1_count;
      phi_A0.row(s) = 1/double(n) * phi_A0_count;
      for(int c=0; c<ncauses; c++){
        psi_A1[c].row(s) = 1/double(n) * (psi_A1_count_sum1[c] - psi_A1_count_sum2[c]);
        psi_A0[c].row(s) = 1/double(n) * (psi_A0_count_sum1[c] - psi_A0_count_sum2[c]);
      }
    }
    
    std::vector< arma::mat > sum34(ncauses);
    for(int c=0; c<ncauses; c++){
      arma::mat Z_minus_Z_bar_dN = Z[c].rows(event_ind) - Z_bar[c];
      Z_minus_Z_bar_dN.each_col() %= dN[c];
      sum34[c] = omega_inv[c] * Z_minus_Z_bar_dN.t();
    }
    
    arma::vec sum1_A1(event_ind.size(), fill::value(0.0));
    arma::vec sum1_A0(event_ind.size(), fill::value(0.0));
    arma::mat sum2_A1(event_ind.size(),t.size(), fill::value(0.0));
    arma::mat sum2_A0(event_ind.size(),t.size(), fill::value(0.0));
    sum1_A1.head(eval_ind.size()) += sum(S_A1.head_cols(eval_ind.size()) % dLambda_A1[cause-1].head_cols(eval_ind.size()), 0).t();
    sum1_A0.head(eval_ind.size()) += sum(S_A0.head_cols(eval_ind.size()) % dLambda_A0[cause-1].head_cols(eval_ind.size()), 0).t();
    for(int s=0; s<t.size(); s++){
      if(t[s] >= min(time)){
        arma::uvec eval_ind_t = find((event != 0) && (time <= t[s]));
        arma::mat dLambda_A1_all(n, eval_ind_t.size());
        arma::mat dLambda_A0_all(n, eval_ind_t.size());
        for(int c=0; c<ncauses; c++){
          dLambda_A1_all += dLambda_A1[c].cols(0, eval_ind_t.size()-1);
          dLambda_A0_all += dLambda_A0[c].cols(0, eval_ind_t.size()-1);
        }
        sum2_A1.submat(0,s,eval_ind_t.size()-1,s) += sum((repelem(F_A1.col(eval_ind_t.size()), 1, eval_ind_t.size()) - F_A1.cols(1, eval_ind_t.size())) % dLambda_A1_all, 0).t();
        sum2_A0.submat(0,s,eval_ind_t.size()-1,s) += sum((repelem(F_A0.col(eval_ind_t.size()), 1, eval_ind_t.size()) - F_A0.cols(1, eval_ind_t.size())) % dLambda_A0_all, 0).t();
      }
    }
    
    arma::mat res_Lin_A1(bs_iter, t.size());
    arma::mat res_Lin_A0(bs_iter, t.size());
    arma::mat res_Bey_A1(bs_iter, t.size());
    arma::mat res_Bey_A0(bs_iter, t.size());
    arma::mat res_Weird_A1(bs_iter, t.size());
    arma::mat res_Weird_A0(bs_iter, t.size());
    
    arma::rowvec SD_A1(t.size(), fill::value(0.0));
    arma::rowvec SD_A0(t.size(), fill::value(0.0));
    arma::rowvec SD(t.size(), fill::value(0.0));
    
    for(int s=0; s<t.size(); s++){
      
      arma::vec t_ind = zeros(event_ind.size());
      if(t[s] >= min(time.elem(eval_ind))){
        int index_t = max(find(time.elem(eval_ind) <= t[s]));
        t_ind.head_rows(index_t+1) = 1*ones(index_t+1);
      }
      
      arma::vec sum3_A1 = (phi_A1.row(s)*sum34[cause-1]).t();
      arma::vec sum3_A0 = (phi_A0.row(s)*sum34[cause-1]).t();
      
      arma::vec sum4_A1(event_ind.size(), fill::value(0.0));
      arma::vec sum4_A0(event_ind.size(), fill::value(0.0));
      for(int c=0; c<ncauses; c++){
        sum4_A1 += (psi_A1[c].row(s) * sum34[c]).t();
        sum4_A0 += (psi_A0[c].row(s) * sum34[c]).t();
      }
      
      for(int i=0; i<bs_iter; i++){
        res_Lin_A1(i,s) = 1/sqrt(n) * sum(G_Lin.col(i) % (t_ind%(sum1_A1-sum2_A1.col(s))+sum3_A1-sum4_A1));
        res_Lin_A0(i,s) = 1/sqrt(n) * sum(G_Lin.col(i) % (t_ind%(sum1_A0-sum2_A0.col(s))+sum3_A0-sum4_A0));
        res_Bey_A1(i,s) = 1/sqrt(n) * sum(G_Bey.col(i) % (t_ind%(sum1_A1-sum2_A1.col(s))+sum3_A1-sum4_A1));
        res_Bey_A0(i,s) = 1/sqrt(n) * sum(G_Bey.col(i) % (t_ind%(sum1_A0-sum2_A0.col(s))+sum3_A0-sum4_A0));
        res_Weird_A1(i,s) = 1/sqrt(n) * sum(G_Weird.col(i) % (t_ind%(sum1_A1-sum2_A1.col(s))+sum3_A1-sum4_A1));
        res_Weird_A0(i,s) = 1/sqrt(n) * sum(G_Weird.col(i) % (t_ind%(sum1_A0-sum2_A0.col(s))+sum3_A0-sum4_A0));
      }
      
      SD_A1[s] = 1/sqrt(n) * sqrt(sum(pow(t_ind%(sum1_A1-sum2_A1.col(s))+sum3_A1-sum4_A1, 2)));
      SD_A0[s] = 1/sqrt(n) * sqrt(sum(pow(t_ind%(sum1_A0-sum2_A0.col(s))+sum3_A0-sum4_A0, 2)));
      SD[s] = 1/sqrt(n) * sqrt(sum(pow((t_ind%(sum1_A1-sum2_A1.col(s))+sum3_A1-sum4_A1) - (t_ind%(sum1_A0-sum2_A0.col(s))+sum3_A0-sum4_A0), 2)));
      
    }
    
    arma::vec lower_calc(t.size(), fill::value(0.0));
    arma::vec upper_calc(t.size(), fill::value(0.0));
    arma::vec lower_Lin(t.size(), fill::value(0.0));
    arma::vec upper_Lin(t.size(), fill::value(0.0));
    arma::vec lower_Bey(t.size(), fill::value(0.0));
    arma::vec upper_Bey(t.size(), fill::value(0.0));
    arma::vec lower_Weird(t.size(), fill::value(0.0));
    arma::vec upper_Weird(t.size(), fill::value(0.0));
    
    arma::rowvec SD_Lin(t.size(), fill::value(0.0));
    arma::rowvec SD_Bey(t.size(), fill::value(0.0));
    arma::rowvec SD_Weird(t.size(), fill::value(0.0));
    
    int quant_ind = round(bs_iter * level);
    
    for(int s=0; s<t.size(); s++){
      if(t[s] >= min(time)){
        
        lower_calc[s] = max(F_A1_mean[s] - F_A0_mean[s] - R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0) * SD[s] / sqrt(n), -1.0);
        upper_calc[s] = min(F_A1_mean[s] - F_A0_mean[s] + R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0) * SD[s] / sqrt(n), 1.0);
        double q_Lin = arma::vec(sort(abs(res_Lin_A1.col(s) - res_Lin_A0.col(s))))[quant_ind];
        lower_Lin[s] = max(F_A1_mean[s] - F_A0_mean[s] - q_Lin / sqrt(n), -1.0);
        upper_Lin[s] = min(F_A1_mean[s] - F_A0_mean[s] + q_Lin / sqrt(n), 1.0);
        double q_Bey = arma::vec(sort(abs(res_Bey_A1.col(s) - res_Bey_A0.col(s))))[quant_ind];
        lower_Bey[s] = max(F_A1_mean[s] - F_A0_mean[s] - q_Bey / sqrt(n), -1.0);
        upper_Bey[s] = min(F_A1_mean[s] - F_A0_mean[s] + q_Bey / sqrt(n), 1.0);
        double q_Weird = arma::vec(sort(abs(res_Weird_A1.col(s) - res_Weird_A0.col(s))))[quant_ind];
        lower_Weird[s] = max(F_A1_mean[s] - F_A0_mean[s] - q_Weird / sqrt(n), -1.0);
        upper_Weird[s] = min(F_A1_mean[s] - F_A0_mean[s] + q_Weird / sqrt(n), 1.0);
        
        SD_Lin[s] = stddev(res_Lin_A1.col(s) - res_Lin_A0.col(s));
        SD_Bey[s] = stddev(res_Bey_A1.col(s) - res_Bey_A0.col(s));
        SD_Weird[s] = stddev(res_Weird_A1.col(s) - res_Weird_A0.col(s));
      }
    }
    
    arma::uvec eval_ind_new = find(SD != 0);
    //double q_sqrt = arma::vec(sort(max(abs(res_A1.cols(eval_ind_new) - res_A0.cols(eval_ind_new)), 1)))[quant_ind];
    double q_Lin_calc = arma::vec(sort(max(abs(res_Lin_A1.cols(eval_ind_new) - res_Lin_A0.cols(eval_ind_new)) / repelem(SD.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
    double q_Lin_EP = arma::vec(sort(max(abs(res_Lin_A1.cols(eval_ind_new) - res_Lin_A0.cols(eval_ind_new)) / repelem(SD_Lin.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
    double q_Bey_calc = arma::vec(sort(max(abs(res_Bey_A1.cols(eval_ind_new) - res_Bey_A0.cols(eval_ind_new)) / repelem(SD.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
    double q_Bey_EP = arma::vec(sort(max(abs(res_Bey_A1.cols(eval_ind_new) - res_Bey_A0.cols(eval_ind_new)) / repelem(SD_Bey.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
    double q_Weird_calc = arma::vec(sort(max(abs(res_Weird_A1.cols(eval_ind_new) - res_Weird_A0.cols(eval_ind_new)) / repelem(SD.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
    double q_Weird_EP = arma::vec(sort(max(abs(res_Weird_A1.cols(eval_ind_new) - res_Weird_A0.cols(eval_ind_new)) / repelem(SD_Weird.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
    ////double q_HW = arma::vec(sort(max(abs(res_A1.cols(eval_ind_new) - res_A0.cols(eval_ind_new)) / repelem(1 + SD.elem(eval_ind_new)^2, 1, bs_iter).t(), 1)))[quant_ind];
    
    if(arms){
      
      arma::vec lower_calc_A1(t.size(), fill::value(0.0));
      arma::vec upper_calc_A1(t.size(), fill::value(0.0));
      arma::vec lower_calc_A0(t.size(), fill::value(0.0));
      arma::vec upper_calc_A0(t.size(), fill::value(0.0));
      arma::vec lower_Lin_A1(t.size(), fill::value(0.0));
      arma::vec upper_Lin_A1(t.size(), fill::value(0.0));
      arma::vec lower_Lin_A0(t.size(), fill::value(0.0));
      arma::vec upper_Lin_A0(t.size(), fill::value(0.0));
      arma::vec lower_Bey_A1(t.size(), fill::value(0.0));
      arma::vec upper_Bey_A1(t.size(), fill::value(0.0));
      arma::vec lower_Bey_A0(t.size(), fill::value(0.0));
      arma::vec upper_Bey_A0(t.size(), fill::value(0.0));
      arma::vec lower_Weird_A1(t.size(), fill::value(0.0));
      arma::vec upper_Weird_A1(t.size(), fill::value(0.0));
      arma::vec lower_Weird_A0(t.size(), fill::value(0.0));
      arma::vec upper_Weird_A0(t.size(), fill::value(0.0));
      
      arma::rowvec SD_Lin_A1(t.size(), fill::value(0.0));
      arma::rowvec SD_Lin_A0(t.size(), fill::value(0.0));
      arma::rowvec SD_Bey_A1(t.size(), fill::value(0.0));
      arma::rowvec SD_Bey_A0(t.size(), fill::value(0.0));
      arma::rowvec SD_Weird_A1(t.size(), fill::value(0.0));
      arma::rowvec SD_Weird_A0(t.size(), fill::value(0.0));
      
      for(int s=0; s<t.size(); s++){
        if(t[s] >= min(time)){
          
          lower_calc_A1[s] = max(F_A1_mean[s] - R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0) * SD_A1[s] / sqrt(n), 0.0);
          upper_calc_A1[s] = min(F_A1_mean[s] + R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0) * SD_A1[s] / sqrt(n), 1.0);
          lower_calc_A0[s] = max(F_A0_mean[s] - R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0) * SD_A0[s] / sqrt(n), 0.0);
          upper_calc_A0[s] = min(F_A0_mean[s] + R::qnorm(1-(1-level)/2, 0.0, 1.0, 1, 0) * SD_A0[s] / sqrt(n), 1.0);
          double q_Lin_A1 = arma::vec(sort(abs(res_Lin_A1.col(s))))[quant_ind];
          lower_Lin_A1[s] = max(F_A1_mean[s] - q_Lin_A1 / sqrt(n), 0.0);
          upper_Lin_A1[s] = min(F_A1_mean[s] + q_Lin_A1 / sqrt(n), 1.0);
          double q_Lin_A0 = arma::vec(sort(abs(res_Lin_A0.col(s))))[quant_ind];
          lower_Lin_A0[s] = max(F_A0_mean[s] - q_Lin_A0 / sqrt(n), 0.0);
          upper_Lin_A0[s] = min(F_A0_mean[s] + q_Lin_A0 / sqrt(n), 1.0);
          double q_Bey_A1 = arma::vec(sort(abs(res_Bey_A1.col(s))))[quant_ind];
          lower_Bey_A1[s] = max(F_A1_mean[s] - q_Bey_A1 / sqrt(n), 0.0);
          upper_Bey_A1[s] = min(F_A1_mean[s] + q_Bey_A1 / sqrt(n), 1.0);
          double q_Bey_A0 = arma::vec(sort(abs(res_Bey_A0.col(s))))[quant_ind];
          lower_Bey_A0[s] = max(F_A0_mean[s] - q_Bey_A0 / sqrt(n), 0.0);
          upper_Bey_A0[s] = min(F_A0_mean[s] + q_Bey_A0 / sqrt(n), 1.0);
          double q_Weird_A1 = arma::vec(sort(abs(res_Weird_A1.col(s))))[quant_ind];
          lower_Weird_A1[s] = max(F_A1_mean[s] - q_Weird_A1 / sqrt(n), 0.0);
          upper_Weird_A1[s] = min(F_A1_mean[s] + q_Weird_A1 / sqrt(n), 1.0);
          double q_Weird_A0 = arma::vec(sort(abs(res_Weird_A0.col(s))))[quant_ind];
          lower_Weird_A0[s] = max(F_A0_mean[s] - q_Weird_A0 / sqrt(n), 0.0);
          upper_Weird_A0[s] = min(F_A0_mean[s] + q_Weird_A0 / sqrt(n), 1.0);
          
          SD_Lin_A1[s] = stddev(res_Lin_A1.col(s));
          SD_Lin_A0[s] = stddev(res_Lin_A0.col(s));
          SD_Bey_A1[s] = stddev(res_Bey_A1.col(s));
          SD_Bey_A0[s] = stddev(res_Bey_A0.col(s));
          SD_Weird_A1[s] = stddev(res_Weird_A1.col(s));
          SD_Weird_A0[s] = stddev(res_Weird_A0.col(s));
          
        }
      }
      
      double q_Lin_calc_A1 = arma::vec(sort(max(abs(res_Lin_A1.cols(eval_ind_new)) / repelem(SD_A1.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Lin_calc_A0 = arma::vec(sort(max(abs(res_Lin_A0.cols(eval_ind_new)) / repelem(SD_A0.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Lin_EP_A1 = arma::vec(sort(max(abs(res_Lin_A1.cols(eval_ind_new)) / repelem(SD_Lin_A1.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Lin_EP_A0 = arma::vec(sort(max(abs(res_Lin_A0.cols(eval_ind_new)) / repelem(SD_Lin_A0.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Bey_calc_A1 = arma::vec(sort(max(abs(res_Bey_A1.cols(eval_ind_new)) / repelem(SD_A1.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Bey_calc_A0 = arma::vec(sort(max(abs(res_Bey_A0.cols(eval_ind_new)) / repelem(SD_A0.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Bey_EP_A1 = arma::vec(sort(max(abs(res_Bey_A1.cols(eval_ind_new)) / repelem(SD_Bey_A1.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Bey_EP_A0 = arma::vec(sort(max(abs(res_Bey_A0.cols(eval_ind_new)) / repelem(SD_Bey_A0.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Weird_calc_A1 = arma::vec(sort(max(abs(res_Weird_A1.cols(eval_ind_new)) / repelem(SD_A1.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Weird_calc_A0 = arma::vec(sort(max(abs(res_Weird_A0.cols(eval_ind_new)) / repelem(SD_A0.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Weird_EP_A1 = arma::vec(sort(max(abs(res_Weird_A1.cols(eval_ind_new)) / repelem(SD_Weird_A1.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      double q_Weird_EP_A0 = arma::vec(sort(max(abs(res_Weird_A0.cols(eval_ind_new)) / repelem(SD_Weird_A0.elem(eval_ind_new).t(), bs_iter, 1), 1)))[quant_ind];
      
      res_CI_WBS.set_size(24, t.size());
      res_CI_WBS.row(0) = lower_calc_A1.t();
      res_CI_WBS.row(1) = upper_calc_A1.t();
      res_CI_WBS.row(2) = lower_calc_A0.t();
      res_CI_WBS.row(3) = upper_calc_A0.t();
      res_CI_WBS.row(4) = lower_calc.t();
      res_CI_WBS.row(5) = upper_calc.t();
      res_CI_WBS.row(6) = lower_Lin_A1.t();
      res_CI_WBS.row(7) = upper_Lin_A1.t();
      res_CI_WBS.row(8) = lower_Lin_A0.t();
      res_CI_WBS.row(9) = upper_Lin_A0.t();
      res_CI_WBS.row(10) = lower_Lin.t();
      res_CI_WBS.row(11) = upper_Lin.t();
      res_CI_WBS.row(12) = lower_Bey_A1.t();
      res_CI_WBS.row(13) = upper_Bey_A1.t();
      res_CI_WBS.row(14) = lower_Bey_A0.t();
      res_CI_WBS.row(15) = upper_Bey_A0.t();
      res_CI_WBS.row(16) = lower_Bey.t();
      res_CI_WBS.row(17) = upper_Bey.t();
      res_CI_WBS.row(18) = lower_Weird_A1.t();
      res_CI_WBS.row(19) = upper_Weird_A1.t();
      res_CI_WBS.row(20) = lower_Weird_A0.t();
      res_CI_WBS.row(21) = upper_Weird_A0.t();
      res_CI_WBS.row(22) = lower_Weird.t();
      res_CI_WBS.row(23) = upper_Weird.t();
      
      res_CB_WBS.set_size(36, t.size());
      res_CB_WBS.fill(0.0);
      res_CB_WBS.submat(0,eval_ind_new[0],0,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - q_Lin_calc_A1 * SD_A1.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(1,eval_ind_new[0],1,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) + q_Lin_calc_A1 * SD_A1.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(2,eval_ind_new[0],2,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - q_Lin_EP_A1 * SD_Lin_A1.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(3,eval_ind_new[0],3,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) + q_Lin_EP_A1 * SD_Lin_A1.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(4,eval_ind_new[0],4,t.size()-1) = max(F_A0_mean.elem(eval_ind_new) - q_Lin_calc_A0 * SD_A0.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(5,eval_ind_new[0],5,t.size()-1) = min(F_A0_mean.elem(eval_ind_new) + q_Lin_calc_A0 * SD_A0.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(6,eval_ind_new[0],6,t.size()-1) = max(F_A0_mean.elem(eval_ind_new) - q_Lin_EP_A0 * SD_Lin_A0.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(7,eval_ind_new[0],7,t.size()-1) = min(F_A0_mean.elem(eval_ind_new) + q_Lin_EP_A0 * SD_Lin_A0.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(8,eval_ind_new[0],8,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Lin_calc * SD.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(9,eval_ind_new[0],9,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Lin_calc * SD.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(10,eval_ind_new[0],10,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Lin_EP * SD_Lin.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(11,eval_ind_new[0],11,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Lin_EP * SD_Lin.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      
      res_CB_WBS.submat(12,eval_ind_new[0],12,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - q_Bey_calc_A1 * SD_A1.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(13,eval_ind_new[0],13,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) + q_Bey_calc_A1 * SD_A1.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(14,eval_ind_new[0],14,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - q_Bey_EP_A1 * SD_Bey_A1.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(15,eval_ind_new[0],15,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) + q_Bey_EP_A1 * SD_Bey_A1.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(16,eval_ind_new[0],16,t.size()-1) = max(F_A0_mean.elem(eval_ind_new) - q_Bey_calc_A0 * SD_A0.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(17,eval_ind_new[0],17,t.size()-1) = min(F_A0_mean.elem(eval_ind_new) + q_Bey_calc_A0 * SD_A0.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(18,eval_ind_new[0],18,t.size()-1) = max(F_A0_mean.elem(eval_ind_new) - q_Bey_EP_A0 * SD_Bey_A0.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(19,eval_ind_new[0],19,t.size()-1) = min(F_A0_mean.elem(eval_ind_new) + q_Bey_EP_A0 * SD_Bey_A0.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(20,eval_ind_new[0],20,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Bey_calc * SD.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(21,eval_ind_new[0],21,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Bey_calc * SD.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(22,eval_ind_new[0],22,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Bey_EP * SD_Bey.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(23,eval_ind_new[0],23,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Bey_EP * SD_Bey.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      
      res_CB_WBS.submat(24,eval_ind_new[0],24,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - q_Weird_calc_A1 * SD_A1.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(25,eval_ind_new[0],25,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) + q_Weird_calc_A1 * SD_A1.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(26,eval_ind_new[0],26,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - q_Weird_EP_A1 * SD_Weird_A1.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(27,eval_ind_new[0],27,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) + q_Weird_EP_A1 * SD_Weird_A1.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(28,eval_ind_new[0],28,t.size()-1) = max(F_A0_mean.elem(eval_ind_new) - q_Weird_calc_A0 * SD_A0.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(29,eval_ind_new[0],29,t.size()-1) = min(F_A0_mean.elem(eval_ind_new) + q_Weird_calc_A0 * SD_A0.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(30,eval_ind_new[0],30,t.size()-1) = max(F_A0_mean.elem(eval_ind_new) - q_Weird_EP_A0 * SD_Weird_A0.elem(eval_ind_new) / sqrt(n), zeros(eval_ind_new.size())).t();
      res_CB_WBS.submat(31,eval_ind_new[0],31,t.size()-1) = min(F_A0_mean.elem(eval_ind_new) + q_Weird_EP_A0 * SD_Weird_A0.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(32,eval_ind_new[0],32,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Weird_calc * SD.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(33,eval_ind_new[0],33,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Weird_calc * SD.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(34,eval_ind_new[0],34,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Weird_EP * SD_Weird.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(35,eval_ind_new[0],35,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Weird_EP * SD_Weird.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      
    }else{
      
      res_CI_WBS.set_size(8, t.size());
      res_CI_WBS.row(0) = lower_calc.t();
      res_CI_WBS.row(1) = upper_calc.t();
      res_CI_WBS.row(2) = lower_Lin.t();
      res_CI_WBS.row(3) = upper_Lin.t();
      res_CI_WBS.row(4) = lower_Bey.t();
      res_CI_WBS.row(5) = upper_Bey.t();
      res_CI_WBS.row(6) = lower_Weird.t();
      res_CI_WBS.row(7) = upper_Weird.t();
      
      res_CB_WBS.set_size(12, t.size());
      res_CB_WBS.fill(0.0);
      //res_CB_WBS.submat(0,eval_ind_new[0],0,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_sqrt/sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      //res_CB_WBS.submat(0,eval_ind_new[0],0,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_sqrt/sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(0,eval_ind_new[0],0,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Lin_calc * SD.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(1,eval_ind_new[0],1,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Lin_calc * SD.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(2,eval_ind_new[0],2,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Lin_EP * SD_Lin.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(3,eval_ind_new[0],3,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Lin_EP * SD_Lin.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      
      res_CB_WBS.submat(4,eval_ind_new[0],4,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Bey_calc * SD.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(5,eval_ind_new[0],5,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Bey_calc * SD.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(6,eval_ind_new[0],6,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Bey_EP * SD_Bey.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(7,eval_ind_new[0],7,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Bey_EP * SD_Bey.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      
      res_CB_WBS.submat(8,eval_ind_new[0],8,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Weird_calc * SD.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(9,eval_ind_new[0],9,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Weird_calc * SD.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(10,eval_ind_new[0],10,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_Weird_EP * SD_Weird.elem(eval_ind_new) / sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      res_CB_WBS.submat(11,eval_ind_new[0],11,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_Weird_EP * SD_Weird.elem(eval_ind_new) / sqrt(n), ones(eval_ind_new.size())).t();
      //res_CB_WBS.submat(0,eval_ind_new[0],0,t.size()-1) = max(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) - q_HW*(1 + SD.elem(eval_ind_new)^2)/sqrt(n), -1.0*ones(eval_ind_new.size())).t();
      //res_CB_WBS.submat(0,eval_ind_new[0],0,t.size()-1) = min(F_A1_mean.elem(eval_ind_new) - F_A0_mean.elem(eval_ind_new) + q_HW*(1 + SD.elem(eval_ind_new)^2)/sqrt(n), ones(eval_ind_new.size())).t();
      
    }
    
  }
  
  clock.tock("WBS");
  clock.stop("clock");
  
  List res = List::create(F_A1_mean - F_A0_mean,
                          res_CI_IF = res_CI_IF,
                          res_CB_IF = res_CB_IF,
                          res_CI_WBS = res_CI_WBS,
                          res_CB_WBS = res_CB_WBS);
  
  return(res);
  
  
}