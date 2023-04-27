library(tidyverse)
library(ggplot2)

# load and clean data 
df <- read_csv("covid_working_data.csv")
df <- df %>%
  mutate(week_start = week(week_start)) %>%
  drop_na()



# function for calculating joint log likelihood
# Yij: vector of case counts
# Xij: design matrix
# nij: vector of population sizes
# Pij: vector of population density
# Eij: vector of percentage over 65
# uij: vector of uij (note that this is obtained from ui by identifying the ui for each observation in the data frame)
# eij: vector of eij
# alpha: alpha coefficient
# beta: beta coefficient 
# gamma: gamma coefficient
# delta: delta coefficient
loglik <- function(Yij, Xij, nij, Pij, Eij, uij, eij, alpha, beta, gamma, delta) {
  lambda_ij <- exp(alpha + Xij %*% beta + gamma * Pij + delta * Eij + uij + eij)
  loglik_ij <- - lambda_ij * nij + Yij * log(lambda_ij) + Yij * log(nij)
  loglik_ij <- loglik_ij[loglik_ij != -Inf] # remove -Inf values
  return(sum(loglik_ij, na.rm = TRUE))
}

# function for calculating joint log prior
logprior <- function(alpha, beta, gamma, delta, sigma_u, sigma_e, ui, eij) {
  return(
    -1/2 * log(sigma_u^2*sigma_e^2) -
      (1/200) * (sum(beta^2) + gamma^2 + delta^2 + sigma_u^2 + sigma_e^2) -
      (1/(2 * sigma_u^2)) * sum(ui) - (1/(2 * sigma_e^2)) * sum(eij) 
  )
}
# function for calculating log posterior
logpost <- function(Yij, Xij, nij, Pij, Eij, ui, uij, eij, alpha, beta, gamma, delta, sigma_u, sigma_e) {
  if (min(sigma_u, sigma_e) <= 0) # make sure sigma_u and sigma_e are positive
    return(-Inf)
  else
    return(
      loglik(Yij, Xij, nij, Pij, Eij, uij, eij, alpha, beta, gamma, delta) + 
        logprior(alpha, beta, gamma, delta, sigma_u, sigma_e, ui, eij)
    )
}


MH <- function(df, Yij, Xij, nij, Pij, Eij, pars, avec) {
  res <- pars
  npars <- length(res)
  res_ui <- res[1:50] 
  res_uij <- df %>%
    left_join(data.frame(state = unique(df$state), ui = res_ui), by = "state") %>%
    pull(ui)
  res_eij <- res[51:2377]
  res_alpha <- res[2378]
  res_beta <-  res[2379:2386]
  res_gamma <- res[2387]
  res_delta <- res[2388]
  res_sigma_u <- res[2389]
  res_sigma_e <- res[2390]
    
  for (i in 1:npars) {
    prop <- res
    prop[i] <- res[i] + 2 * avec[i] * (runif(1) - 0.5)
    prop_ui <- prop[1:50] 
    prop_uij <- df %>%
      left_join(data.frame(state = unique(df$state), ui = prop_ui), by = "state") %>%
      pull(ui)
    prop_eij <- prop[51:2377]
    prop_alpha <- prop[2378]
    prop_beta <-  prop[2379:2386]
    prop_gamma <- prop[2387]
    prop_delta <- prop[2388]
    prop_sigma_u <- prop[2389]
    prop_sigma_e <- prop[2390]
    
    if (
      log(runif(1)) < (
        logpost(Yij, Xij, nij, Pij, Eij, prop_ui, prop_uij, prop_eij, prop_alpha, prop_beta, prop_gamma, prop_delta, prop_sigma_u, prop_sigma_e) - 
        logpost(Yij, Xij, nij, Pij, Eij, res_ui, res_uij, res_eij, res_alpha, res_beta, res_gamma, res_delta, res_sigma_u, res_sigma_e) 
      )
    ) {
      res[i] <- prop[i]
    }
        
  }
  return(res)
}

# set up MH algorithm stuff
set.seed(123)
# number of reps 
nrep <- 1000
# information that doesn't change
Yij <- df$weekly_cases
Xij <- cbind(
  df$week_start,
  df$retail_and_recreation_percent_change_from_baseline, 
  df$parks_percent_change_from_baseline,
  df$transit_stations_percent_change_from_baseline,
  df$government_response_index,
  df$containment_index,
  df$economic_support_index,
  df$stringency_index
)
nij <- df$pop2019
Pij <- df$population_density
Eij <- df$Percentage_over_65
# information that does change (i.e. parameters)
ui_start <- rep(0.25, unique(df$state) %>% length())
uij_start <- df %>%
  left_join(data.frame(state = unique(df$state), ui = ui_start), by = "state") %>%
  pull(ui)
eij_start <- rep(0.5, nrow(df))
alpha_start <- -10
beta_start <- rep(0, ncol(Xij))
gamma_start <- 0
delta_start <- 0
sigma_u_start <- 1
sigma_e_start <- 1
pars_start <- c(ui_start, eij_start, alpha_start, beta_start, gamma_start, delta_start, sigma_u_start, sigma_e_start)
# set "a" for each parameter 
avec <- c(rep(0.75, length(ui_start)), rep(0.90, length(eij_start)), 0.15, rep(0.01, length(beta_start)), 0.01, 0.20, rep(0.98, 2))
# create matrix to hold parameter results
mchain <- matrix(NA, nrow = nrep, ncol = length(pars_start))
mchain[1,] <- pars_start

# run MH algorithm 
pb <- txtProgressBar(min = 0, max = nrep, initial = 0, style = 3)
for (i in 2:nrep) {
  setTxtProgressBar(pb, i)
  mchain[i,] <- MH(df, Yij, Xij, nij, Pij, Eij, mchain[i-1,], avec)
}
# create matrices for the different parameters 
mchain_param <- mchain[,2378:2390] 
colnames(mchain_param) <- c("alpha", paste0("beta", rep(1:8)), "gamma", "delta", "sigma_u", "sigma_e")
mchain_ui <- mchain[,1:50]
mchain_eij <- mchain[,51:2377]
# save(mchain, mchain_param, mchain_ui, mchain_eij, file = "mchain.Rdata")
# load("mchain.Rdata")

# plot convergence results
# line plots
data.frame(rep = 1:nrep, mchain_param) %>%
  ggplot(aes(x = rep, y = beta5)) +
  geom_line() +
  theme_bw()
data.frame(rep = 1:nrep, mchain_ui) %>%
  ggplot(aes(x = rep, y = X1)) +
  geom_line() +
  theme_bw()
data.frame(rep = 1:nrep, mchain_eij) %>%
  ggplot(aes(x = rep, y = X1)) +
  geom_line() +
  theme_bw()
# histograms
data.frame(rep = 1:nrep, mchain_param) %>%
  ggplot(aes(x = delta)) +
  geom_histogram() +
  theme_bw()
data.frame(rep = 1:nrep, mchain_ui) %>%
  ggplot(aes(x = X1)) +
  geom_histogram() +
  theme_bw()
data.frame(rep = 1:nrep, mchain_eij) %>%
  ggplot(aes(x = X1)) +
  geom_histogram() +
  theme_bw()

# function to investigate if choices of "a" are good
numunique <- function(mat){
  for (i in 1:ncol(mat))
    cat(i, "\t", length(unique(mat[,i])), "\n")
}
# look at how well "a" works for the different parameters 
numunique(mchain_param)
numunique(mchain_eij)
numunique(mchain_ui)



