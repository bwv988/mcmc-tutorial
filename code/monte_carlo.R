# Monte Carlo integration example.
# RS14022017

# Draw N random samples normal.
N <- 10000
mu <- 0
sd <- 1
data <- rnorm(n = N, mean = mu, sd = sd)

# P(X <= p)
p.wanted <- 0.5
prob.in <- data <= p.wanted

# Monte Carlo estimate.
print(sum(prob.in) / N)

# Compare to true value.
print(pnorm(p.wanted, mean = mu, sd = sd))