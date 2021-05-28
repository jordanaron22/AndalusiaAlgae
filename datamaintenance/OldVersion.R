###########################
library(matrixStats)
library(MASS)
ClassificationAlgae <- function(observed_val,latent_val,mu_a,k_a){
  if (is.na(latent_val)){
    return(1)
  } else if (is.na(observed_val)){
    return(1)
  } else {
    if (latent_val == 1){
      return(dnbinom(observed_val,mu = mu_a, size = k_a))
    } else {
      if (observed_val == 0){
        return(1)
      } else {
        return(0)
      } 
    }
    
  }
}

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

ForwardLinearHelper <- function(old_alpha,alg_data,tox_data,last_tox,tran,mu_a,k_a,tox_em){
  states <- c(0:1)
  toxin_states <- c(0:3)
  pseudo_alpha_matrix <- numeric(length(states))
  
  for (new_mc in 1:length(states)){
    new_alpha <- numeric(length(states))
    for (last_mc in 1:length(states)){
      new_alpha[last_mc] <- old_alpha[last_mc] + log(tran[last_mc,new_mc]) + log(ClassificationAlgae(alg_data,states[new_mc],mu_a,k_a)) + log(ClassificationToxin(tox_data,last_tox,states[new_mc],tox_em))
    }
    pseudo_alpha_matrix[new_mc] <- logSumExp(new_alpha)
  }
  return(pseudo_alpha_matrix)
}

ForwardLinearHelperMissing <- function(old_alpha_mat,alg_data,tox_data,tran,mu_a,k_a,tox_em){
  states <- c(0:1)
  toxin_states <- c(0:3)
  pseudo_alpha_matrix <- numeric(length(states))
  
  for (new_mc in 1:length(states)){
    new_alpha <- numeric(length(states))
    for (last_mc in 1:length(states)){
      to_sum_toxin_vec <- numeric(length(toxin_states))
      for (last_tox in 1:length(toxin_states)){
        old_alpha <- old_alpha_mat[last_tox,]
        to_sum_toxin_vec[last_tox] <- old_alpha[last_mc] + log(tran[last_mc,new_mc]) + log(ClassificationAlgae(alg_data,states[new_mc],mu_a,k_a)) + log(ClassificationToxin(tox_data,toxin_states[last_tox],states[new_mc],tox_em))
      }
      new_alpha[last_mc] <- logSumExp(to_sum_toxin_vec)
    }
    pseudo_alpha_matrix[new_mc] <- logSumExp(new_alpha)
  }
  return(pseudo_alpha_matrix)
}

ForwardLinear <- function(algae_data, toxin_data,time,init,tran,mu_a,k_a,tox_em){
  states <- c(0:1)
  toxin_states <- c(0:3)
  
  alpha_matrix<-vector("list",time)
  alpha_i <- numeric(length(states))
  
  for (i in 1:time){
    if (is.na(toxin_data[i])){
      for (j in 1:length(toxin_states)){
        alpha_matrix[[i]][[j]] <- alpha_i
      }
    } else {
      alpha_matrix[[i]][[1]] <- alpha_i
    }
  }
  
  
  #ASSUME OBS TOX AT TIME 1 IS 0
  if (!is.na(toxin_data[1])){
    for (j in 1:length(states)){
      alpha_matrix[[1]][[1]][[j]] <- log(init[j]) + log(ClassificationAlgae(algae_data[1],states[j],mu_a,k_a)) + log(ClassificationToxin(toxin_data[1],0,states[j],tox_em))
    }
  } else {
    for (j in 1:length(toxin_states)){
      for (i in 1:length(states)){
        alpha_matrix[[1]][[j]][[i]] <- log(init[i]) + log(ClassificationAlgae(algae_data[1],states[i],mu_a,k_a)) + log(ClassificationToxin(toxin_states[j],0,states[i],tox_em))
      }
    }
  }
  
  if (time > 1){
    for (i in 2:time){
      if (!is.na(toxin_data[i-1])){
        old_alpha <- alpha_matrix[[i-1]][[1]]
        old_alpha_mat <- NA
        if (!is.na(toxin_data[i])){
          alpha_matrix[[i]][[1]] <- ForwardLinearHelper(old_alpha,algae_data[i],toxin_data[i],toxin_data[i-1],tran,mu_a,k_a,tox_em)
        } else {
          for (curr_tox_states in 1:length(toxin_states)){
            alpha_matrix[[i]][[curr_tox_states]] <- ForwardLinearHelper(old_alpha,algae_data[i],toxin_states[curr_tox_states],toxin_data[i-1],tran,mu_a,k_a,tox_em)
          }
        }
        
      } else if (is.na(toxin_data[i-1])){
        old_alpha_mat <- matrix(0,length(toxin_states),length(states))
        for (j in 1:length(alpha_matrix[[i-1]])){
          old_alpha_mat[j,] <- alpha_matrix[[i-1]][[j]]
        }
        old_alpha <- NA
        if (!is.na(toxin_data[i])){
          alpha_matrix[[i]][[1]] <- ForwardLinearHelperMissing(old_alpha_mat,algae_data[i],toxin_data[i],tran,mu_a,k_a,tox_em)
        } else if (is.na(toxin_data[i])){
          for (curr_tox_states in 1:length(toxin_states)){
            alpha_matrix[[i]][[curr_tox_states]] <- ForwardLinearHelperMissing(old_alpha_mat,algae_data[i],toxin_states[curr_tox_states],tran,mu_a,k_a,tox_em)
          }
        }
      }
    }
  }
  
  return(alpha_matrix)
}

BackwardLinear <- function(algae_data,toxin_data,time,tran,mu_a,k_a,tox_em){
  states <- c(0:1)
  toxin_states <- c(0:3)
  
  beta_i <- numeric(length(states))
  beta_matrix <- vector("list",length(algae_data))
  
  
  for (i in 1:(length(algae_data))){
    if (is.na(toxin_data[i])){
      for (j in 1:length(toxin_states)){
        beta_matrix[[i]][[j]] <- beta_i
      } 
    } else {
      beta_matrix[[i]][[1]] <- beta_i
    }
  }
  
  for (j in 1:length(beta_matrix[[length(algae_data)]])){
    for (i in 1:length(states)){
      beta_matrix[[length(algae_data)]][[j]][i] <- log(1)
    }
  }
  
  if (time != length(algae_data)) {
    for (i in (length(algae_data)-1):time){
      if (!is.na(toxin_data[i])){
        if (!is.na(toxin_data[i+1])){
          beta_matrix[[i]][[1]] <- BackwardLinearHelper(beta_matrix[[i+1]],algae_data[i+1],toxin_data[i+1],toxin_data[i],tran,mu_a,k_a,tox_em)
        } else {
          beta_matrix[[i]][[1]] <- BackwardLinearHelperMissing(beta_matrix[[i+1]],algae_data[i+1],toxin_data[i],tran,mu_a,k_a,tox_em)
        }
      } else {
        for (j in 1:length(toxin_states)){
          if (!is.na(toxin_data[i+1])){
            beta_matrix[[i]][[j]] <- BackwardLinearHelper(beta_matrix[[i+1]],algae_data[i+1],toxin_data[i+1],toxin_states[j],tran,mu_a,k_a,tox_em)
          } else {
            #NEED TO CHECK THIS CONDITION
            beta_matrix[[i]][[j]] <- BackwardLinearHelperMissing(beta_matrix[[i+1]],algae_data[i+1],toxin_states[j],tran,mu_a,k_a,tox_em)
          }
        }
      }
    }
  }
  
  
  return(beta_matrix)
}

BackwardLinearHelper <- function(new_beta,alg_data,tox_data,obs_tox,tran,mu_a,k_a,tox_em){
  new_beta <- new_beta[[1]]
  pseudo_beta_matrix <- numeric(length(states))
  states <- c(0:1)
  toxin_states <- c(0:3)
  
  for (j in 1:length(states)){
    log_beta <- numeric(length(states))
    for (j_2 in 1:length(states)){
      log_beta[j_2] <-  new_beta[j_2] + log(tran[j,j_2]) + log(ClassificationAlgae(alg_data,states[j_2],mu_a,k_a)) + log(ClassificationToxin(tox_data,obs_tox,states[j_2],tox_em))
    }
    
    pseudo_beta_matrix[j] <- logSumExp(log_beta)
  }
  return(pseudo_beta_matrix)
}

BackwardLinearHelperMissing <- function(new_beta,alg_data,obs_tox,tran,mu_a,k_a,tox_em){
  pseudo_beta_matrix <- numeric(length(states))
  new_beta_mat <- matrix(0,length(toxin_states),length(states))
  states <- c(0:1)
  toxin_states <- c(0:3)
  
  #tox_data is NA
  
  for (i2 in 1:length(new_beta)){
    new_beta_mat[i2,] <- new_beta[[i2]]
  }
  
  for (j in 1:length(states)){
    log_beta <- numeric(length(states))
    for (j_2 in 1:length(states)){
      log_beta_vec <- numeric(length(toxin_states))
      for (k_1 in 1:length(toxin_states)){
        log_beta_vec[k_1] <- new_beta_mat[k_1,j_2] + log(tran[j,j_2]) + log(ClassificationAlgae(alg_data,states[j_2],mu_a,k_a)) + log(ClassificationToxin(toxin_states[k_1],obs_tox,states[j_2],tox_em))
      }
      log_beta[j_2] <- logSumExp(log_beta_vec)
    }
    
    pseudo_beta_matrix[j] <- logSumExp(log_beta)
  }
  return(pseudo_beta_matrix)
}

CollapseBackw <- function(backw,time){
  states <- c(0:1)
  toxin_states <- c(0:3)
  if (length(backw[[time]]) == 1){
    backw_prob <- backw[[time]][[1]]
  } else {
    backw_prob <- numeric(length(states))
    for (i in 1:length(states)){
      to_logsumexp <- numeric(length(toxin_states))
      for (j in 1:length(toxin_states)){
        to_logsumexp[j] <- backw[[time]][[j]][i]
      }
      backw_prob[i] <- logSumExp(to_logsumexp)
    }
  }
  return(backw_prob)
}

CollapseForw <- function(forw,time){
  states <- c(0:1)
  toxin_states <- c(0:3)
  if (length(forw[[time]]) == 1){
    forw_prob <- forw[[time]][[1]]
  } else {
    forw_prob <- numeric(length(states))
    for (i in 1:length(states)){
      to_logsumexp <- numeric(length(toxin_states))
      for (j in 1:length(toxin_states)){
        to_logsumexp[j] <- forw[[time]][[j]][i]
      }
      forw_prob[i] <- logSumExp(to_logsumexp)
    }
  }
  return(forw_prob)
}

CalcInit <- function(forw,backw){
  states <- c(0:1)
  toxin_states <- c(0:3)

  forw_prob <- CollapseForw(forw,1)
  backw_prob <- CollapseBackw(backw,1)

  denom <- logSumExp(forw_prob + backw_prob)
  return((forw_prob + backw_prob) - denom)
}

CalcTran <- function(forw,backw,tran,algae_data,toxin_data,mu_a,k_a,tox_em){
  states <- c(0:1)
  toxin_states <- c(0:3)
  tran_matrix_current <- array(NA,c(length(states),length(states),length(algae_data)-1,length(toxin_states),length(toxin_states)))

  for (time in 2:length(algae_data)){

    for (initial_state in 1:length(states)){
      for(new_state in 1:length(states)){


        if (!is.na(toxin_data[time-1])){
          if (!is.na(toxin_data[time])){
            forw_prob <- CollapseForw(forw,time-1)
            backw_prob <- CollapseBackw(backw,time)
            tran_matrix_current[initial_state,new_state,time-1,1,1] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
              log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a)) + log(ClassificationToxin(toxin_data[time],toxin_data[time-1],new_state - 1,tox_em))
          } else {
            forw_prob <- CollapseForw(forw,time-1)
            for (curr_tox in 1:length(toxin_states)){
              backw_prob <- backw[[time]][[curr_tox]]
              tran_matrix_current[initial_state,new_state,time-1,1,curr_tox] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
                log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a)) + log(ClassificationToxin(toxin_states[curr_tox],toxin_data[time-1],new_state - 1,tox_em))
            }
          }
        } else {
          if (!is.na(toxin_data[time])){
            backw_prob <- CollapseBackw(backw,time)
            for(last_tox in 1:length(toxin_states)){
              forw_prob <- forw[[time-1]][[last_tox]]
              tran_matrix_current[initial_state,new_state,time-1,last_tox,1] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
                log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a)) + log(ClassificationToxin(toxin_data[time],toxin_states[last_tox],new_state - 1,tox_em))
            }
          } else {
            for (last_tox in 1:length(toxin_states)){
              for(curr_tox in 1:length(toxin_states)){
                forw_prob <- forw[[time-1]][[last_tox]]
                backw_prob <- backw[[time]][[curr_tox]]
                tran_matrix_current[initial_state,new_state,time-1,last_tox,curr_tox] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
                  log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a)) + log(ClassificationToxin(toxin_states[curr_tox],toxin_states[last_tox],new_state - 1,tox_em))
              }
            }
          }
        }
      }
    }
  }
  tran_matrix <- matrix(0,2,2)
  for (i in 1:length(states)){
    for (j in 1:length(states)){
      tran_matrix[i,j] <- logSumExp(tran_matrix_current[i,j,,,],na.rm = T)
    }
  }

  for (i in 1:length(states)){
    row_sum <- logSumExp(tran_matrix[i,])
    for (j in 1:length(states)){
      tran_matrix[i,j] <- exp(tran_matrix[i,j] - row_sum)
    }
  }

  return(tran_matrix)
}

ProbSt <- function(forw,backw,time){
  forw_prob <- CollapseForw(forw,time)
  backw_prob <- CollapseBackw(backw,time)
  prob <- forw_prob + backw_prob
  prob <- exp(prob - logSumExp(prob))
  return(prob)
}

ProbWeights <- function(data, forw,backw){
  weights <- numeric(length(data))
  for (i in 1:length(weights)){
    weights[i] <- ProbSt(forw,backw,i)[2]
  }
  return(weights)
}

LogLikenbinom <- function(x,data,weights) {
  rs <- 0
  for (i in 1:length(data)){
    if (!is.na(data[i])){
      if (!is.na(log(dnbinom(data[i], mu = x[1], size = x[2])))){
        rs <- rs + (weights[i] *  log(dnbinom(data[i], mu = x[1], size = x[2])))
      } else {
        return(-Inf)
      }
    }
  }
  return(-rs)
}

Phi2K <- function(mu,phi){
  return(1 / (((phi*mu) - mu) / mu^2))
}

K2Phi <- function(mu,k_a){
  return((mu + (mu^2/k_a)) / mu)
}

GenerateSimulatedMC <- function(init, tran, toxin_data_len,mu_a,k_a,tox_em){
  sim_markov_chain <- numeric(toxin_data_len)
  sim_toxin_data <- numeric(toxin_data_len)
  sim_algae_data <- numeric(toxin_data_len)

  sim_markov_chain[1] <- which(rmultinom(1,1,init) == 1) - 1
  sim_toxin_data[1] <- which(rmultinom(1,1,tox_em[1,,sim_markov_chain[1]+1]) == 1) - 1
  if (sim_markov_chain[1] == 1){
    sim_algae_data[1] <- rnbinom(1,mu = mu_a, size = k_a)
  }

  for (i in 2:toxin_data_len){
    sim_markov_chain[i] <- which(rmultinom(1,1,tran[sim_markov_chain[i-1]+1,]) == 1) - 1
    sim_toxin_data[i] <- which(rmultinom(1,1,tox_em[sim_toxin_data[i-1]+1,,sim_markov_chain[i]+1]) == 1) - 1
    if (sim_markov_chain[i] == 1){
      sim_algae_data[i] <- rnbinom(1,mu = mu_a, size = k_a)
    }
  }

  simulated_true_df <- data.frame(time = c(1:toxin_data_len),
                                  algae = sim_algae_data,
                                  toxin = sim_toxin_data,
                                  mc = sim_markov_chain)

  algae_miss <- rbinom(toxin_data_len,1,.87)
  toxin_miss <- rbinom(toxin_data_len,1,.85)
  # algae_miss <- rbinom(toxin_data_len,1,1/3)
  # toxin_miss <- rbinom(toxin_data_len,1,1/3)
  for (i in 1:toxin_data_len){
    if (algae_miss[i] == 1){
      sim_algae_data[i] <- NA
    }
    if (toxin_miss[i] == 1){
      sim_toxin_data[i] <- NA
    }
  }

  simulated_observed_df <- data.frame(time = c(1:toxin_data_len),
                                      algae = sim_algae_data,
                                      toxin = sim_toxin_data)
  return(list(simulated_true_df,simulated_observed_df))
}

Ord2Mat <- function(coeffs, zetas){
  states <- c(0:1)
  toxin_states <- c(0:3)
  tox_em <- array(1, dim = c(length(toxin_states), length(toxin_states),length(states)))
  for (curr_tox in 1:(length(toxin_states) - 1)){
    for (last_tox in 1:length(toxin_states)){
      for (curr_mc in 1:length(states)){
        
        inter <- zetas[[curr_tox]]
        
        if (last_tox > 1){
          tox_coef <- coeffs[[1]] * (last_tox-1)
        } else{
          tox_coef <- 0
        }
        
        if (curr_mc == 1){
          mc_coef <- 0
        } else {
          mc_coef <- coeffs[[2]]
        }
        
        prob <- inter - tox_coef - mc_coef
        tox_em[last_tox,curr_tox,curr_mc] <- exp(prob) / (1 + exp(prob))
      }
    }
  }
  
  for (last_tox in 1:length(toxin_states)){
    for (curr_mc in 1:length(states)){
      tox_em[last_tox,4,curr_mc] <- tox_em[last_tox,4,curr_mc] - tox_em[last_tox,3,curr_mc]
      tox_em[last_tox,3,curr_mc] <- tox_em[last_tox,3,curr_mc] - tox_em[last_tox,2,curr_mc]
      tox_em[last_tox,2,curr_mc] <- tox_em[last_tox,2,curr_mc] - tox_em[last_tox,1,curr_mc]
    }
  }
  return(tox_em)
}

CalcBetas <- function(algae_data,toxin_data,forw,backw,tran,mu_a,k,tox_em,denom){
  states <- c(0:1)
  toxin_states <- c(0:3)
  ord_log_df <- data.frame("current_toxin"=integer(),
                           "last_toxin"=integer(),
                           "Markov_chain" = integer(),
                           "weight"=integer(),
                           "time"=integer(),
                           stringsAsFactors=T) 
  
  QuantCalcMissing <- function(curr_time,last_toxin,curr_toxin,curr_algae,last_mc,curr_mc,tran,mu_a,k,tox_em){
    return(forw[[curr_time-1]][[last_toxin+1]][last_mc+1] + log(tran[last_mc+1,curr_mc+1]) +
             log(ClassificationAlgae(curr_algae,curr_mc,mu_a,k)) + log(ClassificationToxin(curr_toxin,last_toxin,curr_mc,tox_em)) +
             backw[[curr_time]][[1]][curr_mc+1])
  }
  
  QuantCalcCompletelyMissing <- function(curr_time,last_toxin,curr_toxin,curr_algae,last_mc,curr_mc,tran,mu_a,k,tox_em){
    return(forw[[curr_time-1]][[last_toxin+1]][last_mc+1] + log(tran[last_mc+1,curr_mc+1]) +
             log(ClassificationAlgae(curr_algae,curr_mc,mu_a,k)) + log(ClassificationToxin(curr_toxin,last_toxin,curr_mc,tox_em)) +
             backw[[curr_time]][[curr_toxin+1]][curr_mc+1])
  }
  
  if(!is.na(toxin_data[1])){
    weights_vec <- forw[[1]][[1]] + backw[[1]][[1]]
    weights_vec <- exp(weights_vec - denom)
    for (curr_mc in 1:length(states)){
      ord_log_df <- rbind(ord_log_df,c(toxin_data[1],0,curr_mc-1,weights_vec[curr_mc],1))
    }
  } else {
    weights_vec <- c(forw[[1]][[1]] + backw[[1]][[1]],
                     forw[[1]][[2]] + backw[[1]][[2]],
                     forw[[1]][[3]] + backw[[1]][[3]],
                     forw[[1]][[4]] + backw[[1]][[4]])
    weights_vec <- exp(weights_vec - denom)
    for (curr_toxin_ind in 1:length(toxin_states)){
      for (curr_mc_ind in 1:length(states)){
        ord_log_df <- rbind(ord_log_df,c(curr_toxin_ind-1,0,curr_mc_ind-1,weights_vec[2*(curr_toxin_ind-1) + curr_mc_ind],1))
      }
    }
  }
  
  
  for (curr_time in 2:length(toxin_data)){
    if (!is.na(toxin_data[curr_time])){
      if (!is.na(toxin_data[curr_time-1])){
        weights_vec <- forw[[curr_time]][[1]] + backw[[curr_time]][[1]]
        weights_vec <- exp(weights_vec - denom)
        for (curr_mc in 1:length(states)){
          ord_log_df <- rbind(ord_log_df,c(toxin_data[curr_time],toxin_data[curr_time-1],curr_mc-1,weights_vec[curr_mc],curr_time))
        }
        
      } else {
        weights_vec <- c()
        for (curr_mc in states){
          for (last_toxin in toxin_states){
            q1 <- QuantCalcMissing(curr_time,last_toxin,toxin_data[curr_time],algae_data[curr_time],0,curr_mc,tran,mu_a,k,tox_em)
            q2 <- QuantCalcMissing(curr_time,last_toxin,toxin_data[curr_time],algae_data[curr_time],1,curr_mc,tran,mu_a,k,tox_em)
            weights_vec <- c(weights_vec,logSumExp(c(q1,q2)))
            # if (is.na(q1)){print("A")}
          }
        }
        
        
        weights_vec <- exp(weights_vec - denom)
        for (curr_mc_ind in 1:length(states)){
          for (last_toxin_ind in 1:length(toxin_states)){
            ord_log_df <- rbind(ord_log_df,c(toxin_data[curr_time],last_toxin_ind-1,curr_mc_ind-1,weights_vec[4 * (curr_mc_ind - 1) + last_toxin_ind],curr_time))
          }
        }
      }
      
    } else if (is.na(toxin_data[curr_time])){
      if (!is.na(toxin_data[curr_time - 1])){
        weights_vec <- c(forw[[curr_time]][[1]] + backw[[curr_time]][[1]],
                         forw[[curr_time]][[2]] + backw[[curr_time]][[2]],
                         forw[[curr_time]][[3]] + backw[[curr_time]][[3]],
                         forw[[curr_time]][[4]] + backw[[curr_time]][[4]])
        weights_vec <- exp(weights_vec - denom)
        for (curr_toxin_ind in 1:length(toxin_states)){
          for (curr_mc_ind in 1:length(states)){
            ord_log_df <- rbind(ord_log_df,c(curr_toxin_ind-1,toxin_data[curr_time-1],curr_mc_ind-1,weights_vec[2*(curr_toxin_ind-1) + curr_mc_ind],curr_time))
          }
        }
      } else {
        weights_vec <- c()
        for (curr_toxin in toxin_states){
          for (curr_mc in states){
            for (last_toxin in toxin_states){
              q1 <- QuantCalcCompletelyMissing(curr_time,last_toxin,curr_toxin,algae_data[curr_time],0,curr_mc,tran,mu_a,k,tox_em)
              q2 <- QuantCalcCompletelyMissing(curr_time,last_toxin,curr_toxin,algae_data[curr_time],1,curr_mc,tran,mu_a,k,tox_em)
              weights_vec <- c(weights_vec,logSumExp(c(q1,q2)))
              # if (is.na(q1)){print("B")}
            }
          }
        }
        
        weights_vec <- exp(weights_vec - denom)
        for (curr_toxin_ind in 1:length(toxin_states)){
          for (curr_mc_ind in 1:length(states)){
            for (last_toxin_ind in 1:length(toxin_states)){
              ord_log_df <- rbind(ord_log_df,c(curr_toxin_ind-1,last_toxin_ind-1,curr_mc_ind-1,weights_vec[8*(curr_toxin_ind-1) + 4 * (curr_mc_ind - 1) + last_toxin_ind],curr_time))
            }
          }
        }
      }
    }
  }
  
  
  colnames(ord_log_df) <- c("CurrentTox","LastTox","MarkovChain","Weights","Time")
  ord_log_df$CurrentTox <- factor(ord_log_df$CurrentTox)
  betas <- polr(CurrentTox ~ LastTox + MarkovChain, data = ord_log_df, weights = Weights, method = "logistic")
  # betas <- polr(CurrentTox ~ MarkovChain, data = ord_log_df, weights = Weights, method = "logistic")
  return(list(betas$coefficients,betas$zeta))
}

Cont2Ord <- function(toxin_data,co1,co2,co3){
  for (i in 1:length(toxin_data)){
    if (!is.na(toxin_data[i])){
      if (toxin_data[i] <= co1){
        toxin_data[i] <- 0
      } else if (toxin_data[i] <= co2 & toxin_data[i] > co1){
        toxin_data[i] <- 1
      } else if (toxin_data[i] <= co3 & toxin_data[i] > co2){
        toxin_data[i] <- 2
      } else {
        toxin_data[i] <- 3
      }
    }
  }
  return(toxin_data)
}

GenerateSimulatedMC <- function(init, tran, toxin_data_len,mu_a,k,tox_em,missing_perc){
  sim_markov_chain <- numeric(toxin_data_len)
  sim_toxin_data <- numeric(toxin_data_len)
  sim_algae_data <- numeric(toxin_data_len)
  
  sim_markov_chain[1] <- which(rmultinom(1,1,init) == 1) - 1
  sim_toxin_data[1] <- which(rmultinom(1,1,tox_em[1,,sim_markov_chain[1]+1]) == 1) - 1
  if (sim_markov_chain[1] == 1){
    sim_algae_data[1] <- rnbinom(1,mu = mu_a, size = k)
  } 
  
  for (i in 2:toxin_data_len){
    sim_markov_chain[i] <- which(rmultinom(1,1,tran[sim_markov_chain[i-1]+1,]) == 1) - 1
    sim_toxin_data[i] <- which(rmultinom(1,1,tox_em[sim_toxin_data[i-1]+1,,sim_markov_chain[i]+1]) == 1) - 1
    if (sim_markov_chain[i] == 1){
      sim_algae_data[i] <- rnbinom(1,mu = mu_a, size = k)
    }
  }
  
  simulated_true_df <- data.frame(time = c(1:toxin_data_len),
                                  algae = sim_algae_data,
                                  toxin = sim_toxin_data,
                                  mc = sim_markov_chain)
  
  algae_miss <- rbinom(toxin_data_len,1,missing_perc)
  toxin_miss <- rbinom(toxin_data_len,1,missing_perc)
  for (i in 1:toxin_data_len){
    if (algae_miss[i] == 1){
      sim_algae_data[i] <- NA
    }
    if (toxin_miss[i] == 1){
      sim_toxin_data[i] <- NA
    }
  }
  
  simulated_observed_df <- data.frame(time = c(1:toxin_data_len),
                                      algae = sim_algae_data,
                                      toxin = sim_toxin_data)
  return(list(simulated_true_df,simulated_observed_df))
}

RunEM <- function(algae_data,toxin_data,init,tran,mu_a,k_a,betas,threshold){
  algae_data <- replace(algae_data, algae_data < threshold,0)
  denom_vec <- numeric()
  like_decrease <- F
  last_denom <- -Inf
  tox_em <- Ord2Mat(betas[[1]],betas[[2]])
  forw <- ForwardLinear(algae_data,toxin_data,length(algae_data),init,tran,mu_a,k_a,tox_em)
  backw <- BackwardLinear(algae_data,toxin_data,1,tran,mu_a,k_a,tox_em)
  denom <- logSumExp(CollapseForw(forw,length(forw)))
  init <- exp(CalcInit(forw,backw))
  tran <- CalcTran(forw,backw,tran, algae_data,toxin_data,mu_a,k_a,tox_em)
  
  betas <- CalcBetas(algae_data,toxin_data,forw,backw,tran,mu_a,k_a,tox_em,denom)
  tox_em <- Ord2Mat(betas[[1]],betas[[2]])
  
  weights <- ProbWeights(algae_data,forw,backw)
  nbin_par <- suppressWarnings(optim(c(mu_a,k_a), LogLikenbinom, data = algae_data, weights = weights)$par)
  mu_a <- nbin_par[1]
  k_a <- nbin_par[2]
  
  forw <- ForwardLinear(algae_data,toxin_data,length(algae_data),init,tran,mu_a,k_a,tox_em)
  backw <- BackwardLinear(algae_data,toxin_data,1,tran,mu_a,k_a,tox_em)
  new_denom <- logSumExp(CollapseForw(forw,length(forw)))
  print(new_denom - denom)
  while ((new_denom - denom) > .001){
    denom <- new_denom
    init <- exp(CalcInit(forw,backw))
    tran <- CalcTran(forw,backw,tran, algae_data,toxin_data,mu_a,k_a,tox_em)
    betas <- CalcBetas(algae_data,toxin_data,forw,backw,tran,mu_a,k_a,tox_em,denom)
    tox_em <- Ord2Mat(betas[[1]],betas[[2]])
    
    weights <- ProbWeights(algae_data,forw,backw)
    nbin_par <- suppressWarnings(optim(c(mu_a,k_a), LogLikenbinom, data = algae_data, weights = weights)$par)
    mu_a <- nbin_par[1]
    k_a <- nbin_par[2]
    
    forw <- ForwardLinear(algae_data,toxin_data,length(algae_data),init,tran,mu_a,k_a,tox_em)
    backw <- BackwardLinear(algae_data,toxin_data,1,tran,mu_a,k_a,tox_em)
    new_denom <- logSumExp(CollapseForw(forw,length(forw)))
    print(new_denom - denom)
    denom_vec <- c(denom_vec,new_denom)
  }
  
  estimated_parameters <- list(init, tran, mu_a, k_a, betas)
  return(estimated_parameters)
}

SimData <- function(init_true,tran_true,mu_a_true,k_true,betas_true,missing_perc, sample_size){
  tox_em_true <- Ord2Mat(betas_true[[1]],betas_true[[2]])
  
  algae_tox_dfs <- GenerateSimulatedMC(init_true,tran_true,sample_size,mu_a_true,k_true,tox_em_true,missing_perc)
  algae_toxin_sim_true_df <-algae_tox_dfs[[1]]
  algae_toxin_sim_obs_df <- algae_tox_dfs[[2]]
  
  algae_data_true <- algae_toxin_sim_true_df$algae
  toxin_data_true <- algae_toxin_sim_true_df$toxin
  markov_chain <- algae_toxin_sim_true_df$mc
  
  algae_data <- algae_toxin_sim_obs_df$algae
  toxin_data <- algae_toxin_sim_obs_df$toxin
  return(list(algae_data,toxin_data))
}

##########Simulated##########

states <- c(0:1)
toxin_states <- c(0:3)

init_true <- c(.75,.25)
tran_true <- matrix(c(.60,.40,
                      .20,.80), 2,2, byrow = T)
mu_a_true <- 64
k_true <- 0.25
betas_true <- list(c(1,3),c(3,4,5))
missing_perc <- .3
sample_size <- 2000

simulated_data <- SimData(init_true,tran_true,mu_a_true,k_true,betas_true,missing_perc,sample_size)
algae_data <- simulated_data[[1]]
toxin_data <- simulated_data[[2]]
est_param <- RunEM(algae_data,toxin_data,init_true,tran_true,mu_a_true,k_true,betas_true,50)



  
