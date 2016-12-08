# simulate data
set.seed(1)
n <- 500 # obs
x <- runif(n) # covariate

## predictor dependent mixture weights
weights <- cbind(exp(-2 * x), 
                 1-exp(-2 * x))

## mean of each regression
means <- cbind(x, x^4)

## variance of each regression
vars <- cbind(0.01, 0.04)

## mean of function
mean.function <- weights[, 1] * means[, 1] + weights[, 2] * means[, 2]

## variance of function
var.function <- weights[, 1] ^ 2 * vars[, 1] + weights[, 2] ^ 2 * vars[, 2]

## data observed
y <- rnorm(n, mean.function, sqrt(var.function))

## visualize data
library(ggplot2)
theme_set(theme_classic())
data <- data.frame(x=x, y=y, mean=mean.function)

ggplot(data, aes(x, y)) +
    geom_point() +
    geom_line(aes(y=mean), colour="red", size=2) +
    # same limits as figure 4 lower righthand corner
    xlim(0, 1) +
    ylim(-0.5, 1.5) +
    ylab("E(y|x)")

# hyperparameters
nu.1 <- nu.2 <- 0.01
eta.1 <- eta.2 <- 2
nu.0 <- 1
Sigma.0 <- diag(1)
mu.0 <- 0
Sigma.mu <- n * solve(t(x) %*% x)
a.gamma <- -0.05
b.gamma <- 1.05
psi <- 0.05