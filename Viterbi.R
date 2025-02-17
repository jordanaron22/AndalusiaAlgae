
library(AndalusiaAlgaePackage)
library(tidyverse)
library(matrixStats)


##### NEW FUNCTION CALLS 
##### NOW RETURNS ALGAE_TOX_DFS SO WE CAN COMPARE TRUE MC STATES
SimData <- function(init_true,tran_true,mu_a_true,k_true,betas_true,missing_perc,data_len,cell_counts,bin_prob){
  states <- c(0:1)
  toxin_states <- c(0:3)
  
  tox_em_true <- Ord2Mat(betas_true[[1]],betas_true[[2]],states,toxin_states)
  
  algae_tox_dfs <- GenerateSimulatedMC(init_true,tran_true,data_len,mu_a_true,k_true,tox_em_true,missing_perc,cell_counts,bin_prob)
  algae_toxin_sim_true_df <-algae_tox_dfs[[1]]
  algae_toxin_sim_obs_df <- algae_tox_dfs[[2]]
  
  algae_data_true <- algae_toxin_sim_true_df$algae
  toxin_data_true <- algae_toxin_sim_true_df$toxin
  markov_chain <- algae_toxin_sim_true_df$mc
  
  algae_data <- algae_toxin_sim_obs_df$algae
  toxin_data <- algae_toxin_sim_obs_df$toxin
  return(list(algae_data,toxin_data,algae_tox_dfs))
}


###SAME AS BEFORE BUT NEED FOR VITERBI
ClassificationAlgae <- function(observed_val,latent_val,mu_a,k_a,cell_counts,bin_prob){
  if (is.na(latent_val)){
    return(1)
  } else if (is.na(observed_val)){
    return(1)
  } else {
    if (latent_val == 1){
      if (cell_counts){
        return(dnbinom(observed_val,mu = mu_a, size = k_a))
      } else {
        return(dbinom(observed_val,1,bin_prob))
      }
    } else {
      if (observed_val == 0){
        return(1)
      } else {
        return(0)
      }
    }

  }
}


##SAME AS BEFORE BUT NEED FOR VITERB
ClassificationToxin <- function(observed_val,last_observed_val,latent_val,tox_em){
  if (is.na(latent_val)){
    return(1)
  } else if (is.na(observed_val)){
    return(1)
  } else if (is.na(last_observed_val)){
    return(1)
  } else {
    return(tox_em[last_observed_val+1,observed_val+1,latent_val+1])
  }
}

########


##Still need to discretize this
# jun_data <- read.csv("juneau_alex_psp (1).csv") %>% mutate(CellCount = as.numeric(Cells.in.3.min.Net.Tow))

#Is algae data in cell/liter, then set to true
#Is algae data binary absence/presence, then set to false
cell_counts <- T

#Lmit of detection for algae counts
#Set to 0 if algae count is binary
# threshold <- 40
threshold <- 0

#Stopping threshold for EM
epsilon <- .01

#starting initial Markov state values
init_true <- c(.75,.25)

#starting Markov transition probabilities
tran_true <- matrix(c(.60,.40,
                      .20,.80), 2,2, byrow = T)

#starting mean paramater for negative binomial
# mu_a_true <- 300
mu_a_true <- 300

#starting size paramater for negative binomial
# k_true <- 1.1
k_true <- 1.1

#starting probability parameter for binomial
# bin_prob_true <- NA
bin_prob_true <- .7

#starting regression coefficients for ordinal logistic regression
#first grouping is for effects of last DST state and current Markov state, respectivly
#second grouping is for intercepts
betas_true <- list(c(1,3),c(3,4,5))

######Generate Simulated Data######
#Percent of data to remove
missing_perc <- 1/2

#Number of days to simulate
sample_size <- 2000

#Creates simulated data and organizes it
simulated_data <- SimData(init_true,tran_true,mu_a_true,k_true,betas_true,missing_perc,sample_size,cell_counts,bin_prob_true)
algae_data <- simulated_data[[1]]
toxin_data <- simulated_data[[2]]
true_mc <- simulated_data[[3]][[1]]$mc
###################################

# algae_data <- jun_data$CellCount
# toxin_data <- jun_data$psp_result


#Runs bivariate autoregressive HMM on algae_data and toxin data
#Uses init_true as starting initial Markov state values
#Uses tran_true as starting Markov transition probabilities
#Uses mu_a_true as starting mean paramater for negative binomial
#Uses k_true as starting size paramater for negative binomial
#Uses betas_true as starting regression coefficients for ordinal logistic regression
#Uses threshold as threshold for algae data
#Uses epsilon as stopping value for EM

est_param <- RunEM(algae_data,toxin_data,
                   init_true,tran_true,mu_a_true,k_true,betas_true,
                   threshold,epsilon,cell_counts,bin_prob_true)

#Output est_param is a list of lists
#est_param[[1]] are the estimated initial Markov state values
#est_param[[2]] are the estimated Markov transition probabilities
#est_param[[3]] is the estimated mean paramater for negative binomial
#est_param[[4]] is the estimated size paramater for negative binomial
#est_param[[5]] is the estimated probability parameter for the binomial
#est_param[[6]][[1]][[1]] is the estimated regression effect from previous DST state on current DST state 
#est_param[[6]][[1]][[2]] is the estimated regression effect from current Markov state on current DST state 
#est_param[[6]][[2]] are the estimated intercept coefficients

########

init <- est_param[[1]]
tran <- est_param[[2]]
mu_a <- est_param[[3]]
k_a <- est_param[[4]]
bin_prob <- est_param[[5]]
mc_coef <- est_param[[6]][[1]]
inter_coef <- est_param[[6]][[2]]

PosteriorDecoding <- function(algae_data, toxin_data, est_param){
  states <- c(0:1)
  toxin_states <- c(0:3)
  
  mc_coef <- est_param[[6]][[1]]
  inter_coef <- est_param[[6]][[2]]
  
  tox_em <- Ord2Mat(mc_coef,inter_coef,states,toxin_states)
  
  forw <- ForwardLinear(algae_data,toxin_data,sample_size,est_param[[1]],est_param[[2]],est_param[[3]],est_param[[4]],tox_em,states,toxin_states,cell_counts = T,est_param[[5]])
  backw <- BackwardLinear(algae_data,toxin_data, 1,est_param[[2]],est_param[[3]],est_param[[4]],tox_em,states,toxin_states,cell_counts = T,est_param[[5]])
  times <- c(1:sample_size)
  
  forbackw_decode <- lapply(times, ProbSt,forw = forw, backw = backw, states = states, toxin_states = toxin_states)
  posterior_decode <- unlist(lapply(forbackw_decode, which.max)) - 1
  # ProbSt(forw,backw, 30, states, toxin_states)
  
  return(posterior_decode)
}

Viterbi <- function(algae_data,toxin_data,est_param){
  init <- est_param[[1]]
  tran <- est_param[[2]]
  mu_a <- est_param[[3]]
  k_a <- est_param[[4]]
  bin_prob <- est_param[[5]]
  mc_coef <- est_param[[6]][[1]]
  inter_coef <- est_param[[6]][[2]]
  
  states <- c(0:1)
  toxin_states <- c(0:3)
  tox_em <- Ord2Mat(mc_coef,inter_coef,states,toxin_states)
  
  viterbi_mat <- matrix(NA,2,length(algae_data))
  
  viterbi_mat[1,1] <- log(init[1]) + 
    log(ClassificationAlgae(algae_data[1],0,mu_a,k_a,T,bin_prob)) + 
    log(ClassificationToxin(toxin_data[1],toxin_data[1],0,tox_em))
  
  viterbi_mat[2,1] <- log(init[2]) + 
    log(ClassificationAlgae(algae_data[1],1,mu_a,k_a,T,bin_prob)) + 
    log(ClassificationToxin(toxin_data[1],toxin_data[1],1,tox_em))
  
  
  for (time in 2:sample_size){
    viterbi_mat[1,time] <- log(ClassificationAlgae(algae_data[time],0,mu_a,k_a,T,bin_prob)) + 
      log(ClassificationToxin(toxin_data[time],toxin_data[time-1],0,tox_em)) + 
      max(viterbi_mat[1,time-1] + log(tran[1,1]),
          viterbi_mat[2,time-1] + log(tran[2,1]))
    
    
    viterbi_mat[2,time] <- log(ClassificationAlgae(algae_data[time],1,mu_a,k_a,T,bin_prob)) + 
      log(ClassificationToxin(toxin_data[time],toxin_data[time-1],1,tox_em)) + 
      max(viterbi_mat[1,time-1] + log(tran[1,2]),
          viterbi_mat[2,time-1] + log(tran[2,2]))
  }
  
  decoded_mc <- apply(viterbi_mat,2,which.max) - 1
  
  return(decoded_mc)
  
}

viterbi_decode <- Viterbi(algae_data,toxin_data,est_param)
posterior_decode <- PosteriorDecoding(algae_data, toxin_data, est_param)
  

sum(decoded_mc != true_mc) / sample_size
sum(posterior_decode != true_mc) / sample_size
