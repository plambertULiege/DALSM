predict.DALSM = function(object, data, newdata, probs, ...){
    if (missing(data)){
        message("The data frame <data> used to fit the DALSM model should be provided !")
        return(NA)
    }
    if (missing(newdata)){
        message("The data frame <newdata> for which 'predictions' are desired is missing !")
        return(NA)
    }
    npred = nrow(newdata)
    ##
    ## Values for <mu> for the <newdata>
    temp = DesignFormula(object$formula1, plyr::rbind.fill(data,newdata))
    if (!identical(object$regr1$knots,temp$knots)){
        message("Predictions requested for covariate values outside their range in the fitted model for location !")
        return(NA)
    }
    Xcal.1 = tail(temp$Xcal, npred)
    mu.pred = Xcal.1 %*% object$psi1
    ## Values for <sd> for the <newdata>
    temp = DesignFormula(object$formula2, plyr::rbind.fill(data,newdata))
    if (!identical(object$regr2$knots,temp$knots)){
        message("Predictions requested for covariate values outside their range in the fitted model for dispersion !")
        return(NA)
    }
    Xcal.2 = tail(temp$Xcal, npred)
    sd.pred = exp(Xcal.2 %*% object$psi2)
    ## Quantiles
    if (!missing(probs)){
        if (!(all(probs>=0) & all(probs<=1))){
            message("The requested quantiles should be in (0,1) !")
            return(NA)
        }
        nquant = length(probs)
        qval = numeric(nquant) ; names(qval) = probs
        quant.pred = matrix(nrow=npred,ncol=nquant)
        colnames(quant.pred) = probs
        low = min(object$derr$knots) ; up = max(object$derr$knots) ; mid = .5*(low+up)
        for (k in 1:nquant){
            qval[k] = optim(mid, fn = function(y) abs(object$derr$pdist(y)-probs[k]),
                            lower=low, upper=up,
                            method="L-BFGS-B")$par
            quant.pred[,k] = mu.pred + qval[k]*sd.pred
        }
    }
    ans = list(mu=c(mu.pred), sd=c(sd.pred))
    if (!missing(probs)) ans = c(ans, list(quant=quant.pred, qval=qval, probs=probs))
    return(ans)
}



