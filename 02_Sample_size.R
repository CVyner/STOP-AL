library(rstan)
library(binom)

##local run
local=TRUE
Output_path <- "<PATH_TO_OUTPUT>" #make sure to end with a /
Read_path <- "<PATH_TO_RSTAN_MODEL>" #make sure to end with a /

##read in analysis model
compiled_model <- readRDS(paste0(Read_path,"logistic_priors.rds"))

simulate_trial_simple <- function(N,null_leak_rate,alt_leak_rate,allocation=0.5){
  
  #add intercept
  X0 <- rep(1,length.out=N)
  
  ##number in experimental and control
  exp_N <- ceiling(allocation*N)
  cntl_N <- N-exp_N
  
  ##indicator of treatment
  X1 <- c(rep(1,exp_N),rep(0,cntl_N))
  
  ##probability of complication
  p <- null_leak_rate*(1-X1) + alt_leak_rate*X1
  
  ##simulate outcome
  y <- rbinom(N,1,p)
  
  return(data.frame(y=y,X=cbind(X0,X1)))
}

run_simple <- function(N_2,N_3,alpha_interim,alpha_final,null_leak,alt_leak,interim_check){
  
  if(a_cntl > 1){
    ##allocation so that the two beta distributions have roughly even information at the end of stage_I
    allocation_stage_I <- (N_2+a_cntl+b_cntl-1)/(2*N_2)
  }else{
    allocation_stage_I <- 0.5
  }
  
  stage_I_data <- simulate_trial_simple(
    N=N_2,
    null_leak_rate=null_leak,
    alt_leak_rate=alt_leak,
    allocation=allocation_stage_I)
  
  success_in_exp_interim <- sum(1*(stage_I_data$y==1 & stage_I_data$X.X1==1))
  success_in_cntl_interim <- sum(1*(stage_I_data$y==1 & stage_I_data$X.X1==0))
  failure_in_exp_interim <- sum(1*(stage_I_data$y==0 & stage_I_data$X.X1==1))
  failure_in_cntl_interim <- sum(1*(stage_I_data$y==0 & stage_I_data$X.X1==0))
  
  stan_dat <- 
    list(
      succ_exp=as.integer(success_in_exp_interim),
      succ_cntl=as.integer(success_in_cntl_interim),
      fail_exp=as.integer(failure_in_exp_interim),
      fail_cntl=as.integer(failure_in_cntl_interim),
      a_cntl=a_cntl,
      b_cntl=b_cntl,
      a_exp=a_exp,
      b_exp=b_exp
    )
  
  ##analysis at interim
  stan_out_interim <- extract(sampling(object=compiled_model,data=stan_dat,chains=1,iter=4000,show_messages=FALSE,check_data=FALSE,verbose=FALSE,refresh=0))$beta1
  
  test_stat_interim <- quantile(stan_out_interim,alpha_interim)
  
  s1 <- 0
  ##interim check
  if(test_stat_interim > interim_check ){
    s1 <- 1
  }
  
  ##Regardless of if we proceed, simulate rest of trial
  stage_II_data <- simulate_trial_simple(
    N=N_3,
    null_leak_rate=null_leak,
    alt_leak_rate=alt_leak)
  
  ##stage 2 data
  success_in_exp <- success_in_exp_interim + sum(1*(stage_II_data$y==1 & stage_II_data$X.X1==1))
  success_in_cntl <- success_in_cntl_interim + sum(1*(stage_II_data$y==1 & stage_II_data$X.X1==0))
  failure_in_exp <- failure_in_exp_interim + sum(1*(stage_II_data$y==0 & stage_II_data$X.X1==1))
  failure_in_cntl <- failure_in_cntl_interim + sum(1*(stage_II_data$y==0 & stage_II_data$X.X1==0))
  
  stan_dat <- 
    list(
      succ_exp=as.integer(success_in_exp),
      succ_cntl=as.integer(success_in_cntl),
      fail_exp=as.integer(failure_in_exp),
      fail_cntl=as.integer(failure_in_cntl),
      a_cntl=a_cntl,
      b_cntl=b_cntl,
      a_exp=a_exp,
      b_exp=b_exp
    )
  
  stan_out <- extract(sampling(object=compiled_model,data=stan_dat,chains=1,iter=4000,show_messages=FALSE,check_data=FALSE,verbose=FALSE,refresh=0))$beta1
 
  test_stat <- quantile(stan_out,1-alpha_final)
  
  if(test_stat < 0){
    s2 <- 0
  }else{
    s2 <- 1
  }
  s11 <-0
  s01 <- 0
  s10 <- 0
  s00 <- 0
  if(s1==1 & s2==1){
    ##would have stopped at interim and final
    s11 <- 1
  }else if(s1==0 & s2==1){
    ##would have stopped at final only
    s01 <-1
  }else if(s1==1 & s2==0){
    ##would have stopped at interim, but if continued would not have stopped at final
    s10 <- 1
  }else if (s1==0 & s2==0){
    ##never stop
    s00 <- 1
  }
  
  return(c(s1,s2,s00,s01,s10,s11))
}

##wrapper for BOSSS

sim_trial <- function(n=1000,N_int=200,alpha_interim=0.2,alpha_final=0.06,interim_check=0,cntl_leak=0.12,exp_leak=0.06){
 
  
  N_2 <- N_int
  N_3 <- n - N_2
  
  
  out <- run_simple(
    N_2=N_2,
    N_3=N_3,
    alpha_interim=alpha_interim,
    alpha_final=alpha_final,
    interim_check=interim_check,
    null_leak = cntl_leak,
    alt_leak = exp_leak
  )
  
  return(out)
  
}

dat_display <- function(sim_out){
  
  return(data.frame(
    stop_at_both=sim_out[,6],
    stop_at_final_only=sim_out[,4],
    stop_at_interim_only=sim_out[,5],
    no_stop=sim_out[,3],
    reject_binding=1-sim_out[,6]-sim_out[,4]-sim_out[,5],
    reject_non_binding=1-sim_out[,6]-sim_out[,4],
    sense_check=sim_out[,3]+sim_out[,4]+sim_out[,5]+sim_out[,6]
  ))
  
}






power_out_int <- rep(0,6)
type_I_out_int<- rep(0,6)

if(local){
  n_sims <- 10000
  alpha_final <- 0.05
  seed <- 91212946
  set.seed(91212946)
  
  n <- 750
  n_int <- 300
  alpha_int <- 0.3
  alpha_final <- 0.05
  a_cntl <- 5.788
  b_cntl <- 48.5
  a_exp <- 1/2
  b_exp <- 1/2
  null_rate <- 0.095
}

for(i in 1:n_sims){
    

    ##simulate under alternate hypothesis
    out <- sim_trial(n=n,N_int=n_int,alpha_interim=alpha_int,alpha_final = alpha_final,interim_check = 0,cntl_leak=null_rate,exp_leak=null_rate/2)
    power_out_int[1:6] <- power_out_int[1:6]+1*out[1:6]/n_sims
    
    ##simulate under null hypothesis
    out <- sim_trial(n=n,N_int=n_int,alpha_interim=alpha_int,alpha_final = alpha_final,interim_check = 0,cntl_leak=null_rate,exp_leak=null_rate)
    type_I_out_int[1:6] <- type_I_out_int[1:6]+1*out[1:6]/n_sims
 
  }
  #s1,s2,s00,s01,s10,s11
  curr_out <- data.frame(n_sims=n_sims,
                         n=n,
                         n_int=n_int,
                         null_rate=null_rate,
                         alt_rate=null_rate/2,
                         alpha_final=alpha_final,
                         alpha_int=alpha_int,
                         a_cntl=a_cntl,
                         b_cntl=b_cntl,
                         a_exp=a_exp,
                         b_exp=b_exp,
                         seed=seed,
                         power_s1=power_out_int[1],
                         power_s2=power_out_int[2],
                         power_s00=power_out_int[3],
                         power_s01=power_out_int[4],
                         power_s10=power_out_int[5],
                         power_s11=power_out_int[6],
                         type_I_s1=type_I_out_int[1],
                         type_I_s2=type_I_out_int[2],
                         type_I_s00=type_I_out_int[3],
                         type_I_s01=type_I_out_int[4],
                         type_I_s10=type_I_out_int[5],
                         type_I_s11=type_I_out_int[6]
  )

##power
pow <- curr_out$power_s00
pow_CI <- binom::binom.confint(x=pow*n_sims,n=n_sims,conf.level=0.95,methods=c("exact"))
power <- c(pow_CI$lower,pow,pow_CI$upper)

##type-I
tI <- curr_out$type_I_s00
tI_CI <- binom::binom.confint(x=tI*n_sims,n=n_sims,conf.level=0.95,methods=c("exact"))
type_I <- c(tI_CI$lower,tI,tI_CI$upper)

##non-binding type-I
tI_non_binding <- curr_out$type_I_s00 + curr_out$type_I_s10
tI_non_binding_CI <- binom::binom.confint(x=tI_non_binding*n_sims,n=n_sims,conf.level=0.95,methods=c("exact"))
type_I_non_binding <- c(tI_non_binding_CI$lower,tI_non_binding,tI_non_binding_CI$upper)

##combine for output
operating_char_mat <- rbind(power,type_I,type_I_non_binding)

summary <- data.frame(
  summary=c("power","type_I","non_binding_type_I"),
  lower_CI=operating_char_mat[,1],
  estimate=operating_char_mat[,2],
  upper_CI=operating_char_mat[,3]
)

  
write.csv(curr_out,paste0(Output_path,"full_output.csv"))
write.csv(summary,paste0(Output_path,"summary.csv"))
