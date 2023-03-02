
{
    n <- 1e4
    a <- 3  # gamma shape
    b <- 2 # gamma rate
    shape <- 2 # Weibull shape
    
    x <- seq(0.01, 10, by = 0.01)
    f      <- shape * x^(shape-1) * a / b * (b/(b+x^shape))^(a+1)
    f_prop <- x^(shape-1) / (b+x^shape)^(a+1) # proportional to
    plot(x,f, type = "l")
}

lambda <- rgamma(n, a, b)
hist(lambda, freq = F)

sigma <- rweibull(n, shape, lambda^(-1/shape)) # random Weibull rate
hist(sigma, breaks = 100, freq = F)
lines(x,f)

hist(rweibull(n, shape, (a/b)^(-1/shape)), breaks = 100, freq=F) # keeping the Weibul fixed
lines(x,f)
