# MCMC example in R.
# Adapted from: https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/
# RS13/02/2017

# Set the "true" parameters for the model (we want to simulate a linear relationship).
beta1 <- 5
beta0 <- 0
trueSd <- 10
sampleSize <- 31

# Create independent x-values.
x <- (-(sampleSize - 1) / 2) : ((sampleSize - 1) / 2)

# Create dependent values according to beta1 * x + beta0 + N(0,sd)
y <- beta1 * x + beta0 + rnorm(n = sampleSize,
                               mean = 0,
                               sd = trueSd)
plot(x, y, main = "Scatterplot of test data")

# Calculate the likelihood of the data, given the parameters.
likelihood <- function(param) {
  beta1 = param[1]
  beta0 = param[2]
  sd = param[3]
  
  pred = beta1 * x + beta0
  
  # Draw from a normal density, but logarithmically. The result is in terms of log(p).
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
  
  # This way, the product turns into a sum (numerically more stable).
  sumll = sum(singlelikelihoods)
  
  return(sumll)   
}

# Return the likelihood function values.
slopevalues <- function(x) {
  return(likelihood(c(x, beta0, trueSd)))
}

# Plot likelihood function.
slopelikelihoods <- lapply(seq(3, 7, by = .05), slopevalues)
plot(seq(3, 7, by=.05), slopelikelihoods, type="l", 
     xlab = "Values of slope parameter beta1", ylab = "Log likelihood")

# Calculate the prior densities.
prior <- function(param) {
  beta1 = param[1]
  beta0 = param[2]
  sd = param[3]
 
  # Use a flat prior for beta1, sd, and a normal prior for beta0.
  # Again, work in terms of log(p).
  beta1.prior = dunif(beta0, min = 0, max = 10, log = T)
  beta0.prior = dnorm(beta1, sd = 5, log = T)
  sdprior = dunif(sd, min = 0, max = 30, log = T)
  
  return(beta1.prior + beta0.prior + sdprior)
}

# Calculate the posterior density.
posterior <- function(param) {
  return(likelihood(param) + prior(param))
}

# MH Proposal function.
# Use a normal RV as proposal function.
proposalfunction <- function(param) {
  return(rnorm(3, mean = param, sd = c(0.1, 0.5, 0.3)))
}

# MH MCMC algorithm.
run_metropolis_MCMC <- function(startvalue, iterations) {
  # Set up the chain of samples.
  chain = array(dim = c(iterations + 1, 3))
  chain[1, ] = startvalue
  
  for (i in 1:iterations){
    # Draw from the proosal distribution q.
    proposal = proposalfunction(chain[i, ])
    
    # Calculate the acceptance probability alpha for move.
    alpha = exp(posterior(proposal) - posterior(chain[i, ]))
    
    # Update the chain.
    if (runif(1) < alpha){
      chain[i + 1, ] = proposal
    } else {
      chain[i + 1, ] = chain[i, ]
    }
  }
  
  return(chain)
}

# MCMC simulation.
N = 10000
startvalue = c(4, 0, 10)
set.seed(1233)
chain = run_metropolis_MCMC(startvalue, N)

burnIn = 5000
acceptance = 1 - mean(duplicated(chain[-(1:burnIn), ]))

## Plot histogram of posterior densitites (marginal), true value and posterior mean.
par(mfrow = c(2, 3))
hist(chain[-(1:burnIn), 1], nclass = 30, main="Posterior of beta1", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn), 1]))
abline(v = beta1, col = "red" )
hist(chain[-(1:burnIn), 2], nclass = 30, main="Posterior of beta0", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn), 2]))
abline(v = beta0, col = "red" )
hist(chain[-(1:burnIn), 3], nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn), 3]))
abline(v = trueSd, col = "red" )
plot(chain[-(1:burnIn), 1], type = "l", xlab="True value = red line" , main = "Chain values of beta1")
abline(h = beta1, col = "red" )
plot(chain[-(1:burnIn), 2], type = "l", xlab="True value = red line" , main = "Chain values of beta0")
abline(h = beta0, col = "red" )
plot(chain[-(1:burnIn), 3], type = "l", xlab="True value = red line" , main = "Chain values of sd")
abline(h = trueSd, col = "red" )

# Fit a normal regression model for comparison.
summary(lm(y ~ x))
