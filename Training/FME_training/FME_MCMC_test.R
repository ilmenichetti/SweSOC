
library(FME)

mu <- 10
std <- 1
Nfun <- function(p){
 -2*log(dnorm(max(p), mean = mu, sd = std))
}


pri <- function(p) -2*log(dnorm(p, 8, 1))

MCMC3 <- modMCMC (f = Nfun, p = c( 9.5), niter = 2000, jump = 5,
                     updatecov = 10, prior = pri)


MCMC3 <- modMCMC (f = Nfun, p = c(2, 9.5), niter = 2000, jump = 5,
                  updatecov = 10, prior = pri)





NN <- function(p) {
  mu <- c(1,2,3)
  -2*sum(log(dnorm(p, mean = mu, sd = 0.1))) #-2*log(probability)
}
# simple Metropolis-Hastings
MCMC <- modMCMC(f = NN, p =1, niter = 50000,
                outputlength = 1000, jump = 0.5)

summary(MCMC)
plot(MCMC) # check convergence



NN <- function(p) {
  mu <- c(1,2,3)
  -2*sum(log(dnorm(p, mean = mu, sd = 0.1))) #-2*log(probability)
}
# simple Metropolis-Hastings
MCMC <- modMCMC(f = NN, p =1, niter = 50000,
                outputlength = 1000, jump = 0.5)

summary(MCMC)
plot(MCMC) # check convergence







func <- function(pars) {
  mu <- c(1,2,3)
  -2*sum(log(dnorm(pars, mean = mu, sd = 0.1))) #-2*log(probability)
}


func <- function(pars) {
  mu <- c(1,2,3)
  -2*sum(log(dnorm(pars, mean = mu, sd = 0.1))) #-2*log(probability)
}




# simple Metropolis-Hastings
MCMC <- modMCMC(f = func,
                p = 0:2,
                niter = 5000,
                outputlength = 1000, jump = 0.5)


summary(MCMC)
plot(MCMC) # check convergence
pairs(MCMC)

