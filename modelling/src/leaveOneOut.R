leaveOneOut <- function(D,plotIt=FALSE)
  {
    ## Find the bandwidth giving the best balance between bias and variance
    ## of the model by leave one out.

    ## Plot the data
    if(plotIt)
      {
        dev.new()
        par(mfrow=c(1,2))
        plot(D$xk, D$x)
      }
    ## Make the vector of bandwidths which are fitted
    span <- c(seq(0.2, 1, by=0.1),2,4,10)
    ## Matrix for keeping the residuals
    R <- matrix(nrow=nrow(D), ncol=length(span))
    ## Do it for each bandwidth
    for(ii in 1:length(span))
      {
        print(paste("  Fitting for bandwidth",ii,"of",length(span)))
        ## Do the local 2-order polynomial regression one time for each point
        ## leaving the point out while fitting, and then predicting the point.
        for(i in 1:nrow(D))
          {
            R[i,ii] <- D[i,"x"] - predict(loess(x ~ xk, dat=D[-i,], span=span[ii]), D[i,])
          }
      }
    ## Find the best bandwidth
    RSSkAll <- apply(R, 2, function(x){sum(x^2,na.rm=TRUE)})
    ## Take the RRS for the best span value
    spanBest <- span[which.min(RSSkAll)]
    ## Calculate the RSSk
    RSSk <- sum((D$x - predict(loess(x ~ xk, dat=D, span=spanBest), D))^2,na.rm=TRUE)
    ## Plot the fitted function
    if(plotIt)
      {
        DT <- D[order(D$xk),]
        lines(DT$xk, predict(loess(x ~ xk, dat=D, span=spanBest), DT), col=2)
        ## Plot RSSk as a function of the bandwidth
        plot(span, RSSkAll)
      }
    ##
    return(RSSk)
  }
