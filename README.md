# R package `AndalusiaAlgae`

This package is used to model algae toxin and algae counts simultaneously using a bivariate hidden Markov model.
It was developed focusing on western Andalusia. Algae counts (or binary presence/absence) are modeled using a negative binomial (or just binomial) model.
Toxin measurements are discretized into 4 states and then modeled with ordinal logistic regression. 

## Installation

```{r}
install.packages("devtools")
devtools::install_github("jordanaron22/AndalusiaAlgae")
library(AndalusiaAlgaePackage)
```

## Simulated Data

```{r}
#Is algae data in cell/liter, then set to true
#Is algae data binary absence/presence, then set to false
cell_counts <- F

#Lmit of detection for algae counts
#Set to 0 if algae count is binary
threshold <- 0

#Stopping threshold for EM
epsilon <- 1

#starting initial Markov state values
init_true <- c(.75,.25)

#starting Markov transition probabilities
tran_true <- matrix(c(.60,.40,
                      .20,.80), 2,2, byrow = T)

#starting mean paramater for negative binomial
mu_a_true <- 64


#starting size paramater for negative binomial
k_true <- 0.25

#starting probability parameter for binomial
bin_prob_true <- .7

#starting regression coefficients for ordinal logistic regression
#first grouping is for effects of last DST state and current Markov state, respectivly
#second grouping is for intercepts
betas_true <- list(c(1,3),c(3,4,5))

######Generate Simulated Data######
#Percent of data to remove
missing_perc <- .3

#Number of days to simulate
sample_size <- 2000

#Creates simulated data and organizes it
simulated_data <- SimData(init_true,tran_true,mu_a_true,k_true,betas_true,missing_perc,sample_size,cell_counts,bin_prob_true)
algae_data <- simulated_data[[1]]
toxin_data <- simulated_data[[2]]
###################################

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
```
