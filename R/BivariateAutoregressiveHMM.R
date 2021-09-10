###########################
#' Calculates probability of observed algae
#'
#' Calculates probability of observed algae based on latent state and negative binomial parameters
#' @param observed_val Observed algae measurement
#' @param latent_val Latent Markov state
#' @param mu_a mean negative binomial parameter
#' @param k_a size negative binomial parameter
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Probability of observed algae
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

#' Calculates probability of observed toxin
#'
#' Calculates probability of observed toxin based on latent state and ordinal regression parameters
#' @param observed_val Observed toxin measurement
#' @param last_observed_val Toxin measurement at last time point
#' @param latent_val Latent Markov state
#' @param tox_em Classification matrix from ordinal regression parameters
#' @importFrom stats dnbinom
#' @return Probability of observed toxin
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

ForwardLinearHelper <- function(old_alpha,alg_data,tox_data,last_tox,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){
  pseudo_alpha_matrix <- numeric(length(states))

  for (new_mc in 1:length(states)){
    new_alpha <- numeric(length(states))
    for (last_mc in 1:length(states)){
      new_alpha[last_mc] <- old_alpha[last_mc] + log(tran[last_mc,new_mc]) + log(ClassificationAlgae(alg_data,states[new_mc],mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(tox_data,last_tox,states[new_mc],tox_em))
    }
    pseudo_alpha_matrix[new_mc] <- logSumExp(new_alpha)
  }
  return(pseudo_alpha_matrix)
}

ForwardLinearHelperMissing <- function(old_alpha_mat,alg_data,tox_data,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){
  pseudo_alpha_matrix <- numeric(length(states))

  for (new_mc in 1:length(states)){
    new_alpha <- numeric(length(states))
    for (last_mc in 1:length(states)){
      to_sum_toxin_vec <- numeric(length(toxin_states))
      for (last_tox in 1:length(toxin_states)){
        old_alpha <- old_alpha_mat[last_tox,]
        to_sum_toxin_vec[last_tox] <- old_alpha[last_mc] + log(tran[last_mc,new_mc]) + log(ClassificationAlgae(alg_data,states[new_mc],mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(tox_data,toxin_states[last_tox],states[new_mc],tox_em))
      }
      new_alpha[last_mc] <- logSumExp(to_sum_toxin_vec)
    }
    pseudo_alpha_matrix[new_mc] <- logSumExp(new_alpha)
  }
  return(pseudo_alpha_matrix)
}

#' Autoregressive Forward Algorithm adapted from Stanculescu et al
#'
#' Calculates forward probabilities accounting for autoregressive dependence in the toxin measurements
#' @param algae_data Algae count data
#' @param toxin_data Toxin data measured in micro grams / OA equivalent
#' @param time Day to calculate forward probability at
#' @param init Initial Markov state probabilities
#' @param tran Transition Markov state probabilities
#' @param mu_a mean negative binomial parameter
#' @param k_a size negative binomial parameter
#' @param tox_em Classification matrix from ordinal regression parameters
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Return forward probabilities from day 1 to day time(variable name)
#' @export
ForwardLinear <- function(algae_data, toxin_data,time,init,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){
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
      alpha_matrix[[1]][[1]][[j]] <- log(init[j]) + log(ClassificationAlgae(algae_data[1],states[j],mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_data[1],0,states[j],tox_em))
    }
  } else {
    for (j in 1:length(toxin_states)){
      for (i in 1:length(states)){
        alpha_matrix[[1]][[j]][[i]] <- log(init[i]) + log(ClassificationAlgae(algae_data[1],states[i],mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_states[j],0,states[i],tox_em))
      }
    }
  }

  if (time > 1){
    for (i in 2:time){
      if (!is.na(toxin_data[i-1])){
        old_alpha <- alpha_matrix[[i-1]][[1]]
        old_alpha_mat <- NA
        if (!is.na(toxin_data[i])){
          alpha_matrix[[i]][[1]] <- ForwardLinearHelper(old_alpha,algae_data[i],toxin_data[i],toxin_data[i-1],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
        } else {
          for (curr_tox_states in 1:length(toxin_states)){
            alpha_matrix[[i]][[curr_tox_states]] <- ForwardLinearHelper(old_alpha,algae_data[i],toxin_states[curr_tox_states],toxin_data[i-1],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
          }
        }

      } else if (is.na(toxin_data[i-1])){
        old_alpha_mat <- matrix(0,length(toxin_states),length(states))
        for (j in 1:length(alpha_matrix[[i-1]])){
          old_alpha_mat[j,] <- alpha_matrix[[i-1]][[j]]
        }
        old_alpha <- NA
        if (!is.na(toxin_data[i])){
          alpha_matrix[[i]][[1]] <- ForwardLinearHelperMissing(old_alpha_mat,algae_data[i],toxin_data[i],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
        } else if (is.na(toxin_data[i])){
          for (curr_tox_states in 1:length(toxin_states)){
            alpha_matrix[[i]][[curr_tox_states]] <- ForwardLinearHelperMissing(old_alpha_mat,algae_data[i],toxin_states[curr_tox_states],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
          }
        }
      }
    }
  }

  return(alpha_matrix)
}

#' Autoregressive Backward Algorithm adapted from Stanculescu et al
#'
#' Calculates backward probabilities accounting for autoregressive dependence in the toxin measurements
#' @param algae_data Algae count data
#' @param toxin_data Toxin data measured in micro grams / OA equivalent
#' @param time Day to calculate forward probability at
#' @param tran Transition Markov state probabilities
#' @param mu_a mean negative binomial parameter
#' @param k_a size negative binomial parameter
#' @param tox_em Classification matrix from ordinal regression parameters
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Return backward probabilities from day time to end of data
#' @export
BackwardLinear <- function(algae_data,toxin_data,time,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){

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
          beta_matrix[[i]][[1]] <- BackwardLinearHelper(beta_matrix[[i+1]],algae_data[i+1],toxin_data[i+1],toxin_data[i],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
        } else {
          beta_matrix[[i]][[1]] <- BackwardLinearHelperMissing(beta_matrix[[i+1]],algae_data[i+1],toxin_data[i],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
        }
      } else {
        for (j in 1:length(toxin_states)){
          if (!is.na(toxin_data[i+1])){
            beta_matrix[[i]][[j]] <- BackwardLinearHelper(beta_matrix[[i+1]],algae_data[i+1],toxin_data[i+1],toxin_states[j],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
          } else {
            beta_matrix[[i]][[j]] <- BackwardLinearHelperMissing(beta_matrix[[i+1]],algae_data[i+1],toxin_states[j],tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
          }
        }
      }
    }
  }


  return(beta_matrix)
}

BackwardLinearHelper <- function(new_beta,alg_data,tox_data,obs_tox,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){
  new_beta <- new_beta[[1]]
  pseudo_beta_matrix <- numeric(length(states))

  for (j in 1:length(states)){
    log_beta <- numeric(length(states))
    for (j_2 in 1:length(states)){
      log_beta[j_2] <-  new_beta[j_2] + log(tran[j,j_2]) + log(ClassificationAlgae(alg_data,states[j_2],mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(tox_data,obs_tox,states[j_2],tox_em))
    }

    pseudo_beta_matrix[j] <- logSumExp(log_beta)
  }
  return(pseudo_beta_matrix)
}

BackwardLinearHelperMissing <- function(new_beta,alg_data,obs_tox,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){
  pseudo_beta_matrix <- numeric(length(states))
  new_beta_mat <- matrix(0,length(toxin_states),length(states))

  #tox_data is NA

  for (i2 in 1:length(new_beta)){
    new_beta_mat[i2,] <- new_beta[[i2]]
  }

  for (j in 1:length(states)){
    log_beta <- numeric(length(states))
    for (j_2 in 1:length(states)){
      log_beta_vec <- numeric(length(toxin_states))
      for (k_1 in 1:length(toxin_states)){
        log_beta_vec[k_1] <- new_beta_mat[k_1,j_2] + log(tran[j,j_2]) + log(ClassificationAlgae(alg_data,states[j_2],mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_states[k_1],obs_tox,states[j_2],tox_em))
      }
      log_beta[j_2] <- logSumExp(log_beta_vec)
    }

    pseudo_beta_matrix[j] <- logSumExp(log_beta)
  }
  return(pseudo_beta_matrix)
}

#' Collapses backward probabilities
#'
#' Collapses backward probabilities over missing data
#' @param backw Backward probabilities to collapse
#' @param time Time to collapse at
#' @import matrixStats
#' @return Returns collapsed backward quantity
#' @export
CollapseBackw <- function(backw,time,states,toxin_states){
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

#' Collapses forward probabilities
#'
#' Collapses forward probabilities over missing data
#' @param forw Forward probabilities to collapse
#' @param time Time to collapse at
#' @return Returns collapsed forward quantity
#' @export
CollapseForw <- function(forw,time,states,toxin_states){
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

#' Calculates latent state initial probabilities
#'
#' @param forw Forward probabilities
#' @param backw Backward probabilities
#' @return Returns initial probabilities
#' @export
CalcInit <- function(forw,backw,states,toxin_states){
  forw_prob <- CollapseForw(forw,1,states,toxin_states)
  backw_prob <- CollapseBackw(backw,1,states,toxin_states)

  denom <- logSumExp(forw_prob + backw_prob)
  return((forw_prob + backw_prob) - denom)
}

#' Calculates latent state transition probabilities
#'
#' @param forw Forward probabilities
#' @param backw Backward probabilities
#' @param tran Old Markov state transition probabilities
#' @param algae_data Algae count data
#' @param toxin_data Toxin data measured in micro grams / OA equivalent
#' @param mu_a mean negative binomial parameter
#' @param k_a size negative binomial parameter
#' @param tox_em Classification matrix from ordinal regression parameters
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Returns initial probabilities
#' @export
CalcTran <- function(forw,backw,tran,algae_data,toxin_data,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob){
  tran_matrix_current <- array(NA,c(length(states),length(states),length(algae_data)-1,length(toxin_states),length(toxin_states)))

  for (time in 2:length(algae_data)){

    for (initial_state in 1:length(states)){
      for(new_state in 1:length(states)){


        if (!is.na(toxin_data[time-1])){
          if (!is.na(toxin_data[time])){
            forw_prob <- CollapseForw(forw,time-1,states,toxin_states)
            backw_prob <- CollapseBackw(backw,time,states,toxin_states)
            tran_matrix_current[initial_state,new_state,time-1,1,1] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
              log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_data[time],toxin_data[time-1],new_state - 1,tox_em))
          } else {
            forw_prob <- CollapseForw(forw,time-1,states,toxin_states)
            for (curr_tox in 1:length(toxin_states)){
              backw_prob <- backw[[time]][[curr_tox]]
              tran_matrix_current[initial_state,new_state,time-1,1,curr_tox] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
                log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_states[curr_tox],toxin_data[time-1],new_state - 1,tox_em))
            }
          }
        } else {
          if (!is.na(toxin_data[time])){
            backw_prob <- CollapseBackw(backw,time,states,toxin_states)
            for(last_tox in 1:length(toxin_states)){
              forw_prob <- forw[[time-1]][[last_tox]]
              tran_matrix_current[initial_state,new_state,time-1,last_tox,1] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
                log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_data[time],toxin_states[last_tox],new_state - 1,tox_em))
            }
          } else {
            for (last_tox in 1:length(toxin_states)){
              for(curr_tox in 1:length(toxin_states)){
                forw_prob <- forw[[time-1]][[last_tox]]
                backw_prob <- backw[[time]][[curr_tox]]
                tran_matrix_current[initial_state,new_state,time-1,last_tox,curr_tox] <- forw_prob[initial_state] + backw_prob[new_state] + log(tran[initial_state,new_state]) +
                  log(ClassificationAlgae(algae_data[time],new_state - 1,mu_a,k_a,cell_counts,bin_prob)) + log(ClassificationToxin(toxin_states[curr_tox],toxin_states[last_tox],new_state - 1,tox_em))
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

#' Calculates latent state probability at specific time
#'
#' @param forw Forward probabilities
#' @param backw Backward probabilities
#' @param time Time to calculate Markov at
#' @return Returns latent state probability at specific time
#' @export
ProbSt <- function(forw,backw,time,states,toxin_states){
  forw_prob <- CollapseForw(forw,time,states,toxin_states)
  backw_prob <- CollapseBackw(backw,time,states,toxin_states)
  prob <- forw_prob + backw_prob
  prob <- exp(prob - logSumExp(prob))
  return(prob)
}

#' Calculates relative weights of the binary Markov states
#'
#' @param data Either algae or toxin data to get the length of time in study
#' @param forw Forward probabilities
#' @param backw Backward probabilities
#' @return Returns relative weights of Markov latent states
#' @export
ProbWeights <- function(data, forw,backw,states,toxin_states){
  weights <- numeric(length(data))
  for (i in 1:length(weights)){
    weights[i] <- ProbSt(forw,backw,i,states,toxin_states)[2]
  }
  return(weights)
}

#' Calculates log likelihood of negative binomial for algae data
#'
#' @param x vector of (mean and size) or probability parameter
#' @param data Algae data to calculate likelihood of
#' @param weights Weights from previous function ProbWeights
#' @param cell_counts T/F for using algae cell counts or binary
#' @return Returns log likelihood for negative binomial of algae data
#' @export
LogLikenbinom <- function(x,data,weights,cell_counts) {
  rs <- 0
  for (i in 1:length(data)){
    if (!is.na(data[i])){
      
      if (cell_counts){
        if (!is.na(log(dnbinom(data[i], mu = x[1], size = x[2])))){
          rs <- rs + (weights[i] *  log(dnbinom(data[i], mu = x[1], size = x[2])))
        } else {
          return(-Inf)
        }
      }
      
      if (!cell_counts){
        if (!is.na(log(dbinom(data[i], size = 1,prob = x[1])))){
          rs <- rs + (weights[i] *  log(dbinom(data[i], size = 1,prob = x[1])))
        } else {
          return(-Inf)
        }
      }
      
    }
  }
  return(-rs)
}

#' Transforms negative binomial phi parameter to size
#'
#' Transforms negative binomial phi parameter to size parameter. There are two different possible parameterizations of negative binomial
#' Size (k) is easier to interpret but R uses phi
#' @param mu mean paramater
#' @param phi other parameterization parameter
#' @return Returns size paramaeter for negative binomial
Phi2K <- function(mu,phi){
  return(1 / (((phi*mu) - mu) / mu^2))
}

#' Transforms negative binomial size parameter to phi
#'
#' Transforms negative binomial size parameter to phi parameter. Opposite of Phi2K
#' @param mu mean paramater
#' @param k_a size parameter
#' @return Returns phi paramaeter for negative binomial
K2Phi <- function(mu,k_a){
  return((mu + (mu^2/k_a)) / mu)
}

#' Transforms ordinal regression parameters to classification matrix
#'
#' @param coeffs Regression coefficents
#' @param zetas Regression intercepts
#' @return Returns classification matrix
#' @export
Ord2Mat <- function(coeffs, zetas,states,toxin_states){
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

#' Calculates ordinal logistic regression parameters
#'
#' @param algae_data Algae count data
#' @param toxin_data Toxin data measured in micro grams / OA equivalent
#' @param forw Forward probabilities
#' @param backw Backward probabilities
#' @param tran Old Markov state transition probabilities
#' @param mu_a mean negative binomial parameter
#' @param k size negative binomial parameter
#' @param tox_em Classification matrix from ordinal regression parameters
#' @param denom Denominator used in expected value calculations, total likelihood of all data
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Returns betas for ordinal logistic regression
#' @import MASS
#' @export
CalcBetas <- function(algae_data,toxin_data,forw,backw,tran,mu_a,k,tox_em,denom,states,toxin_states,cell_counts,bin_prob){
  ord_log_df <- data.frame("current_toxin"=integer(),
                           "last_toxin"=integer(),
                           "Markov_chain" = integer(),
                           "weight"=integer(),
                           "time"=integer(),
                           stringsAsFactors=T)

  QuantCalcMissing <- function(curr_time,last_toxin,curr_toxin,curr_algae,last_mc,curr_mc,tran,mu_a,k,tox_em){
    return(forw[[curr_time-1]][[last_toxin+1]][last_mc+1] + log(tran[last_mc+1,curr_mc+1]) +
             log(ClassificationAlgae(curr_algae,curr_mc,mu_a,k,cell_counts,bin_prob)) + log(ClassificationToxin(curr_toxin,last_toxin,curr_mc,tox_em)) +
             backw[[curr_time]][[1]][curr_mc+1])
  }

  QuantCalcCompletelyMissing <- function(curr_time,last_toxin,curr_toxin,curr_algae,last_mc,curr_mc,tran,mu_a,k,tox_em){
    return(forw[[curr_time-1]][[last_toxin+1]][last_mc+1] + log(tran[last_mc+1,curr_mc+1]) +
             log(ClassificationAlgae(curr_algae,curr_mc,mu_a,k,cell_counts,bin_prob)) + log(ClassificationToxin(curr_toxin,last_toxin,curr_mc,tox_em)) +
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
  betas <- polr(CurrentTox ~ LastTox + MarkovChain, data = ord_log_df, weights = ord_log_df$Weights, method = "logistic")
  return(list(betas$coefficients,betas$zeta))
}

#' Continuous to discrete for toxin data
#'
#' @param toxin_data Toxin data measured in micro grams / OA equivalent
#' @param co1 Cutoff 1 for toxin data
#' @param co2 Cutoff 1 for toxin data
#' @param co3 Cutoff 1 for toxin data
#' @return Returns discretized toxin data
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

#' Generates simulated Markov chain data
#'
#' @param init Initial Markov state probabilities
#' @param tran Transition Markov state probabilities
#' @param toxin_data_len Length of days to simulate
#' @param mu_a mean negative binomial parameter
#' @param k size negative binomial parameter
#' @param tox_em Classification matrix from ordinal regression parameters
#' @param missing_perc Percent of data to remove from simulated data
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Returns simulated data
#' @import stats
#' @export
GenerateSimulatedMC <- function(init, tran, toxin_data_len,mu_a,k,tox_em,missing_perc,cell_counts,bin_prob){
  sim_markov_chain <- numeric(toxin_data_len)
  sim_toxin_data <- numeric(toxin_data_len)
  sim_algae_data <- numeric(toxin_data_len)

  sim_markov_chain[1] <- which(rmultinom(1,1,init) == 1) - 1
  sim_toxin_data[1] <- which(rmultinom(1,1,tox_em[1,,sim_markov_chain[1]+1]) == 1) - 1
  if (sim_markov_chain[1] == 1){
    if (cell_counts){
      sim_algae_data[1] <- rnbinom(1,mu = mu_a, size = k)
    } else{
      sim_algae_data[1] <- rbinom(1,prob = bin_prob, size = 1)
    }
  }

  for (i in 2:toxin_data_len){
    sim_markov_chain[i] <- which(rmultinom(1,1,tran[sim_markov_chain[i-1]+1,]) == 1) - 1
    sim_toxin_data[i] <- which(rmultinom(1,1,tox_em[sim_toxin_data[i-1]+1,,sim_markov_chain[i]+1]) == 1) - 1
    if (sim_markov_chain[i] == 1){
      if (cell_counts){
        sim_algae_data[i] <- rnbinom(1,mu = mu_a, size = k)
      } else{
        sim_algae_data[i] <- rbinom(1,prob = bin_prob, size = 1)
        
      }
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


#' Runs EM
#'
#'Runs EM algorithm and outputs estimated parameters
#' @param algae_data Algae count data
#' @param toxin_data Toxin data measured in micro grams / OA equivalent
#' @param init Initial Markov state probabilities
#' @param tran Old Markov state transition probabilities
#' @param mu_a mean negative binomial parameter
#' @param k_a size negative binomial parameter
#' @param betas Ordinal logistic regression coefficients
#' @param threshold Threshold for algae data. Values lower than threshold are set to 0
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Returns initial probabilities
#' @export
RunEM <- function(algae_data,toxin_data,init,tran,mu_a,k_a,betas,threshold,epsilon, cell_counts,bin_prob){

  states <- c(0:1)
  toxin_states <- c(0:3)
  if (!cell_counts){threshold <- 0}
  algae_data <- replace(algae_data, algae_data < threshold,0)
  denom_vec <- numeric()
  like_decrease <- F
  last_denom <- -Inf
  tox_em <- Ord2Mat(betas[[1]],betas[[2]],states,toxin_states)
  forw <- ForwardLinear(algae_data,toxin_data,length(algae_data),init,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
  backw <- BackwardLinear(algae_data,toxin_data,1,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
  denom <- logSumExp(CollapseForw(forw,length(forw),states,toxin_states))
  init <- exp(CalcInit(forw,backw,states,toxin_states))
  tran <- CalcTran(forw,backw,tran, algae_data,toxin_data,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)

  betas <- CalcBetas(algae_data,toxin_data,forw,backw,tran,mu_a,k_a,tox_em,denom,states,toxin_states,cell_counts,bin_prob)
  tox_em <- Ord2Mat(betas[[1]],betas[[2]],states,toxin_states)

  weights <- ProbWeights(algae_data,forw,backw,states,toxin_states)
  if (cell_counts){
    nbin_par <- suppressWarnings(optim(c(mu_a,k_a), LogLikenbinom, data = algae_data, weights = weights, cell_counts = cell_counts)$par)
    mu_a <- nbin_par[1]
    k_a <- nbin_par[2]
  } else {
    nbin_par <- suppressWarnings(optim(c(bin_prob), LogLikenbinom, data = algae_data, weights = weights, cell_counts = cell_counts)$par)
    bin_prob <- nbin_par[1]
  }
    

  forw <- ForwardLinear(algae_data,toxin_data,length(algae_data),init,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
  backw <- BackwardLinear(algae_data,toxin_data,1,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
  new_denom <- logSumExp(CollapseForw(forw,length(forw),states,toxin_states))
  print(new_denom - denom)
  while ((new_denom - denom) > epsilon){
    denom <- new_denom
    init <- exp(CalcInit(forw,backw,states,toxin_states))
    tran <- CalcTran(forw,backw,tran, algae_data,toxin_data,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
    betas <- CalcBetas(algae_data,toxin_data,forw,backw,tran,mu_a,k_a,tox_em,denom,states,toxin_states,cell_counts,bin_prob)
    tox_em <- Ord2Mat(betas[[1]],betas[[2]],states,toxin_states)

    weights <- ProbWeights(algae_data,forw,backw,states,toxin_states)
    if (cell_counts){
      nbin_par <- suppressWarnings(optim(c(mu_a,k_a), LogLikenbinom, data = algae_data, weights = weights, cell_counts = cell_counts)$par)
      mu_a <- nbin_par[1]
      k_a <- nbin_par[2]
    } else {
      nbin_par <- suppressWarnings(optim(c(bin_prob), LogLikenbinom, data = algae_data, weights = weights, cell_counts = cell_counts)$par)
      bin_prob <- nbin_par[1]
    }

    forw <- ForwardLinear(algae_data,toxin_data,length(algae_data),init,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
    backw <- BackwardLinear(algae_data,toxin_data,1,tran,mu_a,k_a,tox_em,states,toxin_states,cell_counts,bin_prob)
    new_denom <- logSumExp(CollapseForw(forw,length(forw),states,toxin_states))
    print(new_denom - denom)
    denom_vec <- c(denom_vec,new_denom)
  }

  estimated_parameters <- list(init, tran, mu_a, k_a, bin_prob, betas)
  return(estimated_parameters)
}

#' Generates simulated data
#'
#' @param init_true Initial Markov state probabilities for simulated data
#' @param tran_true Transition Markov state probabilities for simulated data
#' @param mu_a_true mean negative binomial parameter for simulated data
#' @param k_true size negative binomial parameter
#' @param betas_true Betas for simulated data
#' @param missing_perc Percent of data to remove from simulated data
#' @param data_len Number of days to simulate
#' @param cell_counts T/F for using algae cell counts or binary
#' @param bin_prob probability parameter for binomial
#' @return Returns simulated data
#' @export
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
  return(list(algae_data,toxin_data))
}
