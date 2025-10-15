data{

	int<lower=0> succ_exp; 	// number of successes in exp arm

	int<lower=0> succ_cntl; 	// number of successes in cntl arm
	
	int<lower=0> fail_exp; 	// number of failure in exp arm
	
	int<lower=0> fail_cntl; 	// number of failures in cntl arm


  //beta priors for probabilities
  real<lower=0> a_cntl; 
  real<lower=0> b_cntl;
  real<lower=0> a_exp;
  real<lower=0> b_exp; 

}parameters{
  
  real beta0;  
  real beta1;//fixed effects//fixed effects

}model{

  //priors
  target += beta_lpdf(inv_logit(beta0) | a_cntl, b_cntl);
  target += beta_lpdf( inv_logit(beta0+beta1) | a_exp, b_exp);
  
  //likelihood
  target += bernoulli_lpmf(1 | inv_logit(beta0+beta1))*succ_exp;
  target += bernoulli_lpmf(1 | inv_logit(beta0))*succ_cntl;
  target += bernoulli_lpmf(0 | inv_logit(beta0+beta1))*fail_exp;
  target += bernoulli_lpmf(0 | inv_logit(beta0))*fail_cntl;
}


