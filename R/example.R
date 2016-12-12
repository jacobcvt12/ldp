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

# mcmc parameters
burn <- 1000
post <- 1000
thin <- 2

# DP parameters
N <- 50 # number of stick pieces

# initialize parameters from prior
alpha <- rgamma(1, eta.1, eta.2)
tau <- rgamma(1, nu.1, nu.2)
Sigma.beta <- solve(rWishart(1, nu.0, solve(nu.0 * Sigma.0))[, , 1])
mu.beta <- rnorm(1, mu.0, Sigma.mu)
beta.h.star <- rnorm(N, mu.beta, Sigma.beta)
Gamma.h <- runif(N, a.gamma, b.gamma)
V.h <- rbeta(N, 1, alpha)
K <- integer(length(y))

for (i in 1:n) {
    # predictor dependent set indexing the locations belonging
    # to the psi-nbd of x, eta_x^{psi}
    L.x <- which(abs(x[i] - Gamma.h) < psi) 
    pi <- L.x
    # cardinality of L.x
    N.x <- length(L.x)
    
    # allocate memory for p
    p <- numeric(N.x)
    
    p[1] <- V.h[pi[1]]
    
    # weights/stick pieces
    for (l in 2:(N.x -1)) {
        p[l] <- V.h[pi[l]] * prod(1-V.h[pi[1:l]])
    }
    
    p[N.x] <- prod(1-V.h[pi])
    
    # now sample K 
    K[i] <- sample(N.x, 1, replace=TRUE, prob=p)
}

for (s in 1:(burn+post)) {
    # conditional for K_i
    for (i in 1:n) {
        # predictor dependent set indexing the locations belonging
        # to the psi-nbd of x, eta_x^{psi}
        L.x <- which(abs(x[i] - Gamma.h) < psi) 
        pi <- L.x
        # cardinality of L.x
        N.x <- length(L.x)
        
        # allocate memory for p
        p <- numeric(N.x)
        
        p[1] <- V.h[pi[1]]
        
        # weights/stick pieces
        for (l in 2:(N.x -1)) {
            p[l] <- V.h[pi[l]] * prod(1-V.h[pi[1:l]])
        }
        
        p[N.x] <- prod(1-V.h[pi])
        p.prime <- dnorm(y, x * beta.h.star[K], sqrt(tau)) * p
        
        # now sample K 
        K[i] <- sample(N.x, 1, replace=TRUE, prob=p)
    }
    
    # conditional for V_h, h=1, ..., N
    N.x <- sapply(x, function(y) max(which(abs(y-Gamma.h) < psi)))
    for (h in 1:N) {
        a <- sum(K == h & K != N.x)
        b <- sum(K > h)
        
        V.h[h] <- rbeta(1, 1 + a, alpha + b)
    }
    
    # conditional for Gamma_h
    for (h in 1:N) {
        a <- max(max(x - psi), a.gamma)
        Gamma.h[h] <- runif(a, b)
    }
    
    # conditional for beta_h^*, h=1, ..., N
    for (h in 1:N) {
        beta.h.star <- rnorm(1, mu.betah, sqrt(Sigma.betah))
        mu.betah <- Sigma.betah %*% (solve(Sigma.beta) %*% mu.beta + tau * x[K == h] %*% y[K == h])
        Sigma.betah <- solve(solve(Sigma.beta) + tau * x[K == h] %*% t(x[K == h]))
    }
    
    # conditional for mu_beta
    Sigma.mu.hat <- solve(solve(Sigma.mu) + N * solve(Sigma.beta))
    mu.0.hat <- Sigma.mu.hat %*% (solve(Sigma.mu %*% mu.0) + 
                                  solve(Sigma.beta) * sum(beta.h.star))
    mu.beta <- rnorm(1, mu.0.hat, Sigma.mu.hat)
    
    # conditional for Sigma_{beta}^{-1}
    Sigma.beta <- solve(rWishart(1, N + nu.0, 
                                 solve()))
    
    # conditional for alpha
    alpha <- rgamma(1, eta.1 + N, eta.2 - sum(log(1-V.h)))
}