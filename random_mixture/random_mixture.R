pdf("plots.pdf")

library(TMB)
compile("random_mixture.cpp")
dyn.load(dynlib("random_mixture"))

############################################## CASE 1: Many data points

nobs <- 100
ngroup <- 2

## Template objects for MakeADFun
## * We fix all parameters except one mu parameter to get identifiabilty
map <- list(mu    = factor( c(1,  NA)),
            logsd = factor( c(NA, NA)) )
parameters <- list(
    u = matrix( c(-1, 1), nobs, ngroup, byrow=TRUE),
    mu = c(-2, -2),
    logsd = c(0, 0)
)
data <- list( N=rep(0, nobs) )

## Simulate data (using simulation code from C++)
obj <- MakeADFun(data, parameters, random="u", map=map)
set.seed(123)
data <- obj$simulate()

## Pass simulated data to new object
obj <- MakeADFun(data, parameters, random="u", map=map)

## Use fixed initialization of inner problem
obj$env$random.start <- expression(par[random])

opt <- nlminb(obj$par + 2, obj$fn, obj$gr)
plot(Vectorize(obj$fn), opt$par - 1, opt$par + 1)
title(paste(nobs,"observations", "mu_true=", parameters$mu[1]))

######################################### CASE 2: A minimal example (one observation)
nobs <- 1
ngroup <- 2
map <- list(mu    = factor( c(1,  NA)),
            logsd = factor( c(NA, NA)) )
parameters <- list(
    u = matrix( c(-1, 1), nobs, ngroup, byrow=TRUE),
    mu = c(-2, -2),
    logsd = c(0, 0)
)
data <- list( N=2.7 )  ## 2*(exp(mu+1)+Q) where Q=1, mu=-2
obj <- MakeADFun(data, parameters, random="u", map=map)
obj$env$random.start <- expression(par[random])
opt <- nlminb(obj$par + 2, obj$fn, obj$gr)
plot(Vectorize(obj$fn), opt$par - 3, opt$par + 1, n=500)
title(paste(nobs,"observations", "mu_true=", parameters$mu[1]))


########################################## Understand what causes the spike:

plotImage <- function(mu){
    grid <- seq(-3,4.5,by=0.05) - 2
    df <- expand.grid(grid,grid)
    df <- cbind(df, mu)
    mat <- t(as.matrix(df))
    val <- apply(mat,2,obj$env$f)
    dim(val) <- rep(length(grid),2)
    image(grid,grid,exp(-val))
    contour(grid,grid,exp(-val),add=TRUE)
    lap <- exp(-obj$fn(mu))
    p <- obj$env$last.par[1:2]
    abline(0,1,col="green",lty="dashed")
    points(p[1],p[2],pch=16)
    prior <- obj$env$parList()$mu
    points(prior[1],prior[2],pch=16,col="blue")
    euler <- sum(exp(-val))*(grid[2]-grid[1])^2
    txt <- paste("euler:",round(euler,6),"laplace:",round(lap,6),"mu1:",mu)
    title(txt)
    print(txt)
}

gr <- seq(-.1,.1,by=0.01) + -2
for(mu in gr){
    plotImage(mu)
}

dev.off()
