# simulate data
set.seed(1)

n <- 500
x <- runif(n)
phi <- cbind(exp(-2 * x), 
             1-exp(-2 * x))
z <- apply(phi, 1, function(x) sample(0:1, 1, prob=x))
y <- rnorm(n, 
           x * z + (x ^ 4) * (1-z), 
           sqrt(0.01 * z) + sqrt(0.04 * (1-z)))
