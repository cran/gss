## Summarize ssanova objects
summary.ssanova <- function(obj,diagnostics=FALSE)
{
    y <- model.response(obj$mf,"numeric")
    w <- model.weights(obj$mf)
    offset <- model.offset(obj$mf)
    if (is.null(offset)) offset <- rep(0,length(obj$c))
    ## Residuals
    res <- 10^obj$nlambda*obj$c         
    if (!is.null(w)) res <- res/sqrt(w)
    ## Fitted values
    fitted <- as.numeric(y-res)
    fitted.off <- fitted-offset
    ## (estimated) sigma
    sigma <- sqrt(obj$varht)
    ## R^2
    if (!is.null(w)) {
        r.squared <- sum(w*(fitted-sum(w*fitted)/sum(w))^2)
        r.squared <- r.squared/sum(w*(y-sum(w*y)/sum(w))^2)
    }
    else r.squared <- var(fitted)/var(y)       
    ## Residual sum of squares
    if (is.null(w)) rss <- sum(res^2)
    else rss <- sum(w*res^2)
    ## Penalty associated with the fit
    if (is.null(w)) 
        penalty <- sum(obj$c*fitted.off)
    else penalty <- sum(obj$c*fitted.off*sqrt(w))
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
        comp <- cbind(comp,yhat=fitted.off,y=fitted.off+res,e=res)
        term.label <- obj$terms$labels[obj$terms$labels!="1"]
        term.label <- term.label[term.label!="offset"]
        if (any(outer(term.label,c("yhat","y","e"),"==")))
            warning("gss warning: avoid using yhat, y, or e as variable names")
        colnames(comp) <- c(term.label,"yhat","y","e")
        ## Sweep out constant
        if (!is.null(w))
            comp <- sqrt(w)*comp - outer(sqrt(w),apply(w*comp,2,sum))/sum(w)
        else comp <- sweep(comp,2,apply(comp,2,mean))
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
        if (qr(corr)$rank<dim(corr)[2])
            kappa <- rep(Inf,len=dim(corr)[2])
        else kappa <- as.numeric(sqrt(diag(solve(corr))))
        ## Obtain decomposition of penalty
        rough <- as.vector(10^obj$nlambda*t(comp[,term.label])%*%obj$c/penalty)
        names(kappa) <- names(rough) <- term.label
    }
    else decom <- kappa <- cosines <- rough <- NULL
    ## Return the summaries
    z <- list(call=obj$call,method=obj$method,fitted=fitted,residuals=res,
              sigma=sigma,r.squared=r.squared,rss=rss,penalty=penalty,
              pi=decom,kappa=kappa,cosines=cosines,roughness=rough)
    class(z) <- "summary.ssanova"
    z
}
