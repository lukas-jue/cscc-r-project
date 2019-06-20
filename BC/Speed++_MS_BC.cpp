// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace arma; // use the Armadillo library for matrix computations
using namespace Rcpp;

//Market Share Computations (simple, pa=1)
//[[Rcpp::export]]
vec probabilities_cpp(mat const& beta, mat const& design_joint){
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  mat Xbeta = zeros(draws,nplay);
  //Xbeta
  Xbeta = beta * trans(design_joint);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    vec Xbeta_temp = zeros(nplay); 
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r)); //first row of matrix is NOT interpreted as vec
    vec Xbeta_temp_stab = zeros(nplay);
    vec max_temp = zeros(nplay);
    max_temp.fill(max(Xbeta_temp));
    Xbeta_temp_stab = Xbeta_temp - max_temp;
    vec numerator = exp(Xbeta_temp_stab);
    vec denominator = zeros(nplay);
    denominator.fill(sum(exp(Xbeta_temp_stab)));
    market_share_draws(r,span()) = trans(numerator/denominator);
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}

//Market Share Computations (simple, pa=1)
//[[Rcpp::export]]
vec probabilities_log_cpp(mat const& beta, mat const& design_joint){
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  mat Xbeta = zeros(draws,nplay);
  //Xbeta
  Xbeta = beta * trans(design_joint);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    vec Xbeta_temp = zeros(nplay); 
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r)); //first row of matrix is NOT interpreted as vec
    vec Xbeta_temp_stab = zeros(nplay);
    vec max_temp = zeros(nplay);
    max_temp.fill(max(Xbeta_temp));
    Xbeta_temp_stab = Xbeta_temp - max_temp;
    vec numerator = Xbeta_temp_stab;
    vec denominator = zeros(nplay);
    denominator.fill(log(sum(exp(Xbeta_temp_stab))));
    market_share_draws(r,span()) = trans(exp((numerator-denominator)));
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}


//Market Share Computations (simple, pa=1) with BC
//FIRST entry in constrained vector corresponds to BC parameter
//pr-th row of design contains prices
//[[Rcpp::export]]
mat probabilities_BC_cpp(mat const& beta, mat const& design_joint, int const& pr){
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  mat Xbeta = zeros(draws,nplay);
  //Xbeta
  mat beta_red = beta(span::all,span(1,nvar-1)); //Compute Xbeta
  Xbeta = beta_red * trans(design_joint);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    vec Xbeta_temp = zeros(nplay); 
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r)); //first row of matrix is NOT interpreted as vec
    vec Xbeta_temp_stab = zeros(nplay);
    vec max_temp = zeros(nplay);
    max_temp.fill(max(Xbeta_temp));
    Xbeta_temp_stab = Xbeta_temp - max_temp;
    vec constraint_ind = ones(nplay);
    for(int n=0; n<nplay; n++){
      if (design_joint(n,pr_cpp) > budget) {
        constraint_ind(n) = 0; //set indicator to zero if price is larger than budget
      }
    }
    vec numerator = exp(Xbeta_temp_stab)%constraint_ind;
    vec denominator = zeros(nplay);
    denominator.fill(sum(numerator));
    market_share_draws(r,span()) = trans(numerator/denominator);
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}

//Market Share Computations (simple, pa=1) with BC
//FIRST entry in constrained vector corresponds to BC parameter
//pr-th row of design contains prices
//[[Rcpp::export]]
mat probabilities_BC_log_cpp(mat const& beta, mat const& design_joint, int const& pr){
  int pr_cpp = pr - 1; //c++ indexing
  int draws = beta.n_rows;
  int nplay = design_joint.n_rows;
  int nvar = beta.n_cols;
  mat Xbeta = zeros(draws,nplay);
  //Xbeta
  mat beta_red = beta(span::all,span(1,nvar-1)); //Compute Xbeta
  Xbeta = beta_red * trans(design_joint);
  mat market_share_draws = zeros(draws,nplay);
  //Compute stabilized choice probs now...
  for(int r = 0; r<draws; r++){
    double budget = beta(r,0); //update budget for r-th draw
    vec constraint_ind = ones(nplay);
    for(int n=0; n<nplay; n++){
      if (design_joint(n,pr_cpp) > budget) {
        constraint_ind(n) = 0; //set indicator to zero if price is larger than budget
      }
    }
    vec Xbeta_temp = zeros(nplay); 
    Xbeta_temp = arma::conv_to< vec >::from(Xbeta.row(r))%constraint_ind; //Do NOT choose maximum xbeta not affordable
    vec Xbeta_temp_stab = zeros(nplay);
    vec max_temp = zeros(nplay);
    max_temp.fill(max(Xbeta_temp));
    Xbeta_temp_stab = Xbeta_temp - max_temp;
    
    vec numerator = log(exp(Xbeta_temp_stab)%constraint_ind);
    vec denominator = zeros(nplay);
    denominator.fill(log(sum(exp(numerator))));
    market_share_draws(r,span()) = trans(exp((numerator - denominator)));
  }
  vec exp_ms = arma::conv_to< vec >::from(mean(market_share_draws)); //computes mean for each column
  return exp_ms;
}


//Market Share Computations
//[[Rcpp::export]]
vec exp_market_share_cpp(mat const& beta, cube const& design_joint){
  int draws = beta.n_rows;
  int pa = design_joint.n_slices;
  int nplay = design_joint.n_rows;
  Cube<double> Xbeta = zeros(draws,nplay,pa);
  //Xbeta of all product possibilities
  for(int p=0; p<pa; p++){
    mat sub_design = design_joint.slice(p);
    Xbeta.slice(p) = beta * trans(sub_design);
  }
  Mat<double> market_share_draws = zeros(draws,pa);
  //Compute stabilized choice probs now...
  for(int p = 0; p<pa; p++){
    for(int r = 0; r<draws; r++){
      vec Xbeta_temp = zeros(nplay);
      Xbeta_temp = trans(Xbeta.slice(p).row(r));
      vec Xbeta_temp_stab = zeros(nplay);
      vec max_temp = zeros(nplay);
      max_temp.fill(max(Xbeta_temp));
      Xbeta_temp_stab = Xbeta_temp - max_temp;
      double numerator = exp(Xbeta_temp_stab(0));
      double denominator = sum(exp(Xbeta_temp_stab));
      market_share_draws(r,p) = numerator/denominator;
    }
  }
  vec exp_ms = zeros(pa);
  for(int p=0; p<pa; p++){
    vec pa_draws = market_share_draws(span(),p);
    exp_ms(p) = mean(pa_draws);
  }
  return exp_ms;
}


