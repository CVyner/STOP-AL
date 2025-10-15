library(rstan)

output_path <- "<PATH_TO_STAN_FILE>" #make sure to end with /

compiled_model = stan_model(paste0( output_path,"Logistic_priors.stan"),
                      auto_write =TRUE,model_name="logistic_priors")


