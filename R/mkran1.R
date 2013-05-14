## Combine random effects for mixed-effect models
mkran1 <- function(ran1,ran2)
{
    z <- cbind(ran1$z,ran2$z)
    env <- list(sz1=dim(ran1$z)[2],sig1=ran1$sigma,nz1=length(ran1$init),
                sz2=dim(ran2$z)[2],sig2=ran2$sigma,nz2=length(ran2$init))
    fun <- function(zeta,env) {
        idx1 <- 1:env$sz1
        idx2 <- env$sz1+(1:env$sz2)
        sig <- matrix(0,env$sz1+env$sz2,env$sz1+env$sz2)
        sig[idx1,idx1] <- env$sig1$fun(zeta[1:env$nz1],env$sig1$env)
        sig[idx2,idx2] <- env$sig2$fun(zeta[env$nz1+(1:env$nz2)],env$sig2$env)
        sig
    }
    sigma <- list(fun=fun,env=env)
    list(z=z,sigma=sigma,init=c(ran1$init,ran2$init))
}
