## S3 method
predict9 <- function (object,...) UseMethod("predict9")
## Calculate prediction and Bayesian SE from ssanova objects
predict9.gssanova <- function(object,newdata,ci=FALSE,level=.95,nu=NULL,...)
{
    est <- predict(object,newdata,se.fit=ci)
    if (object$family=="binomial") {
        if (!ci) return(plogis(est))
        else {
            z <- est$fit+outer(est$se.fit*qnorm(1-(1-level)/2),c(0,-1,1))
            z <- plogis(z)
            return(list(fit=z[,1],lcl=z[,2],ucl=z[,3]))
        }
    }
    if (object$family%in%c("poisson","Gamma","inverse.gaussian")) {
        if (!ci) return(exp(est))
        else {
            z <- est$fit+outer(est$se.fit*qnorm(1-(1-level)/2),c(0,-1,1))
            z <- exp(z)
            return(list(fit=z[,1],lcl=z[,2],ucl=z[,3]))
        }
    }
    if (object$family=="nbinomial") {
        if (!is.vector(model.response(object$mf))) {
            if (is.null(nu))
                stop("gss error: nu is missing for nbinomial family with 2 column response" )
        }
        else nu <- object$nu
        if (!ci) return(nu/exp(est))
        else {
            z <- est$fit+outer(est$se.fit*qnorm(1-(1-level)/2),c(0,1,-1))
            z <- nu/exp(z)
            return(list(fit=z[,1],lcl=z[,2],ucl=z[,3]))
        }
    }
    if (object$family=="weibull") {
        gg <- gamma(1+1/object$nu)
        if (!ci) return(gg*exp(est))
        else {
            z <- est$fit+outer(est$se.fit*qnorm(1-(1-level)/2),c(0,-1,1))
            z <- gg*exp(z)
            return(list(fit=z[,1],lcl=z[,2],ucl=z[,3]))
        }
    }
    if (object$family=="lognorm") {
        gg <- exp(1/2/object$nu^2)
        if (!ci) return(gg*exp(est))
        else {
            z <- est$fit+outer(est$se.fit*qnorm(1-(1-level)/2),c(0,-1,1))
            z <- gg*exp(z)
            return(list(fit=z[,1],lcl=z[,2],ucl=z[,3]))
        }
    }
    if (object$family=="loglogis") {
        gg <- pi/object$nu/sin(pi/object$nu)
        if (!ci) return(gg*exp(est))
        else {
            z <- est$fit+outer(est$se.fit*qnorm(1-(1-level)/2),c(0,-1,1))
            z <- gg*exp(z)
            return(list(fit=z[,1],lcl=z[,2],ucl=z[,3]))
        }
    }
    if (object$family=="polr") {
        P <- plogis(outer(est,c(0,cumsum(object$nu)),"+"))
        z <- P[,1]
        J <- dim(P)[2]
        for (i in 2:J) z <- cbind(z,P[,i]-P[,i-1])
        z <- cbind(z,1-P[,J])
        colnames(z) <- NULL
        return(z)
    }
}  
