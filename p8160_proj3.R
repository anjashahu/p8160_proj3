library(tidyverse)
library(ggplot2)
library(lme4)

## prepare data, add log(infection_rate)
data <- read.csv("covid_working_data.csv")
data <- data %>%
  mutate(week_start = as.Date(week_start),
         week = week(week_start),
         log_infection_rate = log(infection_rate))

## fit lmer 
fit_lmer <- lmer(
  log_infection_rate ~ retail_and_recreation_percent_change_from_baseline +
    parks_percent_change_from_baseline + transit_stations_percent_change_from_baseline +
    government_response_index + containment_index + economic_support_index +
    stringency_index + week + population_density + Percentage_over_65 + (1 | state),
  data = data %>%
    mutate(log_infection_rate = log(infection_rate)) %>%
    filter(log_infection_rate != -Inf)
)

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
  loglik_ij <- -log(2*pi*var)/2 - ((log_infection_rate - mean)^2)/(2*var)
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
    if (log(runif(1)) < logpost(log_infection_rate, x, pop_dens, eld_pop, alpha = prop[1], beta = prop[2:9], gamma = prop[10], delta = prop[11], sigma_u = prop[12], sigma_e = prop[13]) - logpost(log_infection_rate, x, pop_dens, eld_pop, alpha = res[1], beta = res[2:9], gamma = res[10], delta = res[11], sigma_u = res[12], sigma_e = res[13]))
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
nrep <- 100000 
npars <- 13
avec <- c(rep(0.01, 10), 0.1, 0.01, 0.01)
mchain <- matrix(NA, nrow = nrep, ncol = npars)
mchain[1,] <- c(summary(fit_lmer)$coefficients[,1], 0.8213, 0.7826) 
for (i in 2:nrep) {
  mchain[i,] <- MH(data$log_infection_rate, x, data$population_density, data$Percentage_over_65, mchain[i-1,], avec)
}

## save mchain
mchain <- as.data.frame(mchain) 
colnames(mchain) <- c("alpha", "beta1", "beta2", "beta3", "beta4", "beta5", "beta6", "beta7", "week", "gamma", "delta", "sigma_u", "sigma_e")
#save(mchain, file = "mchain.Rda")
#load("mchain.Rda")
mchain_last <- mchain[50001:nrep,]

## line plots (all)
par(mfrow = c(3,5))
for(i in 1:13){
  plot(x = 1:nrow(mchain), y = mchain[,i], type = "l", lty = 1, main = colnames(mchain)[i])
  lines(x = 1:nrow(mchain), y = mchain[,i], type = "l", lty = 1)
}

## line plots (after burnout)
par(mfrow = c(3,5))
for(i in 1:13){
  plot(x = 1:nrow(mchain_last), y = mchain_last[,i], type = "l", lty = 1, main = colnames(mchain_last)[i])
  lines(x = 1:nrow(mchain_last), y = mchain_last[,i], type = "l", lty = 1)
}

## histogram
ggplot(mchain, aes(x = alpha)) +
  geom_histogram() +
  theme_bw()

## function to investigate if choices of “a” are good
numunique <- function(mat){
  for (i in 1:ncol(mat))
    cat(i, "\t", length(unique(mat[,i])), "\n")
}

## look at how well “a” works for the different parameters
numunique(mchain)





