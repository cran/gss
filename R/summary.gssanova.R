## Summarize gssanova objects
summary.gssanova <- function(obj,diagnostics=FALSE)
{
    y <- model.response(obj$mf,"numeric")
    wt <- model.weights(obj$mf)
    offset <- model.offset(obj$mf)
    if (!is.null(obj$alpha)) y <- cbind(y,obj$alpha)
    dev.resid <- switch(obj$family,
                        binomial=dev.resid.binomial(y,obj$eta,wt),
                        nbinomial=dev.resid.nbinomial(y,obj$eta,wt),
                        poisson=dev.resid.poisson(y,obj$eta,wt),
                        inverse.gaussian=dev.resid.inverse.gaussian(y,obj$eta,wt),
                        Gamma=dev.resid.Gamma(y,obj$eta,wt))
    dev.null <- switch(obj$family,
                       binomial=dev.null.binomial(y,wt,offset),
                       nbinomial=dev.null.nbinomial(y,wt,offset),
                       poisson=dev.null.poisson(y,wt,offset),
                       inverse.gaussian=dev.null.inverse.gaussian(y,wt,offset),
                       Gamma=dev.null.Gamma(y,wt,offset))
    w <- obj$w
    if (is.null(offset)) offset <- rep(0,length(obj$eta))
    ## Residuals
    res <- 10^obj$nlambda*obj$c 
    ## Fitted values
    fitted <- obj$eta
    fitted.off <- fitted-offset
    ## dispersion
    sigma2 <- obj$varht
    ## RSS, deviance
    rss <- sum(res^2)
    dev <- sum(dev.resid)
    ## Penalty associated with the fit
    penalty <- sum(obj$c*fitted.off*sqrt(w))
    penalty <- as.vector(10^obj$nlambda*penalty)
    ## Calculate the diagnostics
    if (diagnostics) {
        ## Obtain retrospective linear model
        comp <- NULL
        for (label in obj$terms$labels) {
            if (label=="1") next
            if (label=="offset") next
            comp <- cbind(comp,predict(obj,obj$mf,inc=label))
        }
        comp <- cbind(comp,yhat=fitted.off,y=fitted.off+res/sqrt(w),e=res/sqrt(w))
        term.label <- obj$terms$labels[obj$terms$labels!="1"]
        term.label <- term.label[term.label!="offset"]
        colnames(comp) <- c(term.label,"yhat","y","e")
        ## Sweep out constant
        comp <- sqrt(w)*comp - outer(sqrt(w),apply(w*comp,2,sum))/sum(w)
        ## Obtain pi
        comp1 <- comp[,c(term.label,"yhat")]
        decom <- t(comp1) %*% comp1[,"yhat"]
        names(decom) <- c(term.label,"yhat")
        decom <- decom[term.label]/decom["yhat"]
        ## Obtain kappa, norm, and cosines        
        corr <- t(comp)%*%comp
        corr <- t(corr/sqrt(diag(corr)))/sqrt(diag(corr))
        norm <- apply(comp,2,function(x){sqrt(sum(x^2))})
        cosines <- rbind(corr[c("y","e"),],norm)
        rownames(cosines) <- c("cos.y","cos.e","norm")
        corr <- corr[term.label,term.label,drop=FALSE]
        if (qr(corr)$rank<dim(corr)[2]) kappa <- rep(Inf,len=dim(corr)[2])
        else kappa <- as.numeric(sqrt(diag(solve(corr))))
        ## Obtain decomposition of penalty
        rough <- as.vector(10^obj$nlambda*t(comp[,term.label])%*%obj$c/penalty)
        names(kappa) <- names(rough) <- term.label
    }
    else decom <- kappa <- cosines <- rough <- NULL
    ## Return the summaries
    z <- list(call=obj$call,family=obj$family,method=obj$method,iter=obj$iter,
              fitted=fitted,dispersion=sigma2,residuals=res/sqrt(w),rss=rss,
              deviance=dev,dev.resid=sqrt(dev.resid)*sign(res),
              dev.null=dev.null,penalty=penalty,
              pi=decom,kappa=kappa,cosines=cosines,roughness=rough)
    class(z) <- "summary.gssanova"
    z
}
