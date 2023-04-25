library(tidyverse)
library(ggplot2)

## prepare data, add log(infection_rate)
data <- read.csv("covid_working_data.csv")
data <- data %>%
  mutate(week_start = as.Date(week_start),
         log_infection_rate = log(infection_rate))

## log likelihood
#### OBSERVED DATA:
## log_infection_rate: vector of length 2,548
## x: 7 covariates related to government interventions and mobility changes, 2,548*7 matrix
## pop_dens: population_density, vector of length 2,548
## eld_pop: Percentage_over_65, vector of length 2,548
#### PARAMETERS:
## alpha: intercept
## beta: coefficients for x, vector of length 7
## gamma: coefficient for population_density
## delta: coefficient for Percentage_over_65
## sigma_u: std. dev. of state-level RE
## sigma_e: std. dev. of residual error
loglik <- function(log_infection_rate, x, pop_dens, eld_pop, alpha, beta, gamma, delta, sigma_u, sigma_e) {
  mean <- alpha + x%*%beta + gamma*pop_dens + delta*eld_pop
  var <- sigma_u^2 + sigma_e^2
  loglik_ij <- -log(2*pi*var)/2 - ((log_infection_rate - mean)^2)/(2*pi*var)
  loglik_ij <- loglik_ij[loglik_ij != -Inf] ## remove -Inf values
  return(sum(loglik_ij, na.rm = T))
}

## log of prior dist.
logprior <- function(alpha, beta, gamma, delta, sigma_u, sigma_e) {
  return((-1/(2*10)) * (alpha^2 + sum(beta^2) + gamma^2 + delta^2 + sigma_u^2 + sigma_e^2))
}

## log of posterior dist.
logpost <- function(log_infection_rate, x, pop_dens, eld_pop, alpha, beta, gamma, delta, sigma_u, sigma_e) {
  if (min(sigma_u, sigma_e) <= 0) ## make sure sigma_u and sigma_e are positive
    return(-Inf)
  else
    return(loglik(log_infection_rate, x, pop_dens, eld_pop, alpha, beta, gamma, delta, sigma_u, sigma_e) + logprior(alpha, beta, gamma, delta, sigma_u, sigma_e))
}

## MH algorithm
## pars: initial guess of parameters, vector of length 12
## avec: window length for each parameter, vector of length 12
MH <- function(log_infection_rate, x, pop_dens, eld_pop, pars, avec) {
  res <- pars
  npars <- length(pars)
  for (i in 1:npars) {
    prop <- res
    ## proposed step of i^th element =  current i^th element + U(-avec[i], avec[i]) 
    prop[i] <- res[i] + 2 * avec[i] * (runif(1) - 0.5)
    if (log(runif(1)) < logpost(log_infection_rate, x, pop_dens, eld_pop, alpha = prop[1], beta = prop[2:8], gamma = prop[9], delta = prop[10], sigma_u = prop[11], sigma_e = prop[12]) - logpost(log_infection_rate, x, pop_dens, eld_pop, alpha = res[1], beta = res[2:8], gamma = res[9], delta = res[10], sigma_u = res[11], sigma_e = res[12]))
        res[i] <- prop[i]
  }
  return(res)
}

## prepare covariates in x
x <- data %>%
  select(-c(state, week_start, weekly_cases, pop2019, LandArea, Percentage_over_65, infection_rate, population_density, log_infection_rate))
x <- as.matrix(x)

## run MH algorithm
set.seed(123)
nrep <- 90000 
npars <- 12
avec <- rep(0.5, npars) 
mchain <- matrix(NA, nrow = nrep, ncol = npars)
mchain[1,] <- rep(1, npars) 
for (i in 2:nrep) {
  mchain[i,] <- MH(data$log_infection_rate, x, data$population_density, data$Percentage_over_65, mchain[i-1,], avec)
}

## save mchain
mchain <- as.data.frame(mchain) 
colnames(mchain) <- c("alpha", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "gamma", "delta", "sigma_u", "sigma_e")
#save(mchain, file = "mchain.Rda")
#load("mchain.Rda")

## plot convergence result
ggplot(mchain, aes(x = 1:nrow(mchain), y = delta)) +
  geom_line() +
  theme_bw()

## histogram
ggplot(mchain, aes(x = delta)) +
  geom_histogram() +
  theme_bw()

