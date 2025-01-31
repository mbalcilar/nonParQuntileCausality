# nonparametric quantile in causality test
# x = cause (independent) variable
# y = dependent variable
# q = vector of quatiles, default = 0.01,....,0.99
# type = "mean" for causality in mean test, 
#        "variance" for causality in variance test
# hm = numeric scalar, bandwidth, 
#      if NULL  optimal bandwith of Yu and Jones (1998) is used
#
# code only considers first lags of x and y
#
# Returns:
#    stat = vector t statistics for causality at each quantile
#    q = vector of quantiles
#
# Mehmet Balcilar, 2014-7-4

lrq.causality.test <- function(x, y, type=c("mean","variance"), q=NULL, hm=NULL) {

    if(is.null(q)) {
    qvec <-seq(0.01, 0.99, by = 0.01)
  }
  else {
    qvec <- q
  }
  
  type <- match.arg(type)
  
  nq <- length(qvec)
  
  tstatvec <- vector(length=nq, mode="numeric") # initilize the tstat vector
  
  tn <- length(y)-1
  
  yall <- embed(y,2)
  yuv <- yall[,-1]   # y(t-1)
  yur <- yall[,1] # y(t)
  
  xall <- embed(x,2)
  x <- xall[,-1]
  
  if (type == "variance") {
    y2 <- y^2
    x2 <- x^2
    cat("Causality in Variance Test\n")
  }
  else if(type == "mean") {
    y2 <- y
    x2 <- x
    cat("Causality in Mean Test\n")
  }
  else
    stop("Test type should be mean or variance")
  
  y2all <- embed(y2,2)
  y2r <- y2all[,1]  # y(t)
  y2r1 <- y2all[,-1]  # y(t-1)
  x2all <- embed(y2,2)
  x2r <- y2all[,1]  # x(t)
  x2r1 <- y2all[,-1]  # x(t-1)
  
  #print(head(cbind(yuv,yur,x)))
  
  if(is.null(hm)) {
    h <- dpill(yuv, yur, gridsize = tn) # calculate the optimal bandwith qrh which is 
    # based on the optimal bandwidth from mean 
    # regression, as in yu and jones 1998
    #h <- npcdistbw(ydat = y2r, xdat = yuv, tol = 0.01, ftol = 0.01)$xbw
  }
  else{
    h <- hm
  }
  
  cat("bandwidth = ", h, "\n")
  
  for  ( jj in 1:nq) {
    qj <- qvec[jj]
    qrh <- h*((qj*(1-qj)/(dnorm(qnorm(p=qj))^2))^(1/5))
    fit <- lprq2(x=y2r1, y=y2r, h=qrh, tau=qj, x0=yuv, type=type)
    #cat(mean(fit$fv),"\n")
    iftemp <- (y2r <= fit$fv) - qj
    ifvector <- data.matrix(iftemp)
    kk <- matrix(data = 0, nrow = tn, ncol = tn)
    ymatrix = kronecker(y[1: tn], t(vector(length= tn, mode="numeric")+1))-t(kronecker(y[1: tn], t(vector(length= tn, mode="numeric")+1)))
    wmatrix = kronecker(x[1: tn], t(vector(length= tn, mode="numeric")+1))-t(kronecker(x[1: tn], t(vector(length= tn, mode="numeric")+1)))
    kk=dnorm(ymatrix/qrh)*dnorm(wmatrix/(qrh/sd(y)*sd(x)))
    tstat <-  t(ifvector)%*%kk%*%ifvector*sqrt(tn/2/qj/(1-qj)/(tn-1)/sum(kk^2)) # Theorem 3.1, Song et al. (2012)
    tstatvec[jj] <- tstat
    cat("t-stat = ",tstat, "at qunatile", qj, "\n")
  }
  return(list(stat=tstatvec,q=qvec))
}


do.causality.figure <- function(obj,title="") {
  
  #browser()
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  lab_title <- title
  #title <- paste("Causality in", lab1, "for", clab)
  pdata <-  data.frame(cbind(q=obj$q,stats=obj$stat),CV=rep(1.96,length(obj$q)))
  
  nq <- dim(pdata)[1]
  
  p1 <- ggplot(data = pdata, aes(x = q)) +
    #geom_hline(aes(yintercept = 1.96,colour="CV (5%)"), color = "gray70", size = 1.5, linetype=1) +
    geom_line(aes(y=stats,colour="Statistic",
                  linetype="Statistic"), size=2) +
    #geom_line(aes(y=Volatility,colour="Volatility",
    #              linetype="Volatility"), size=2) +
    geom_line(aes(y=CV,colour="CV (5%)",
                  linetype="CV (5%)"), size=1) +
    ggtitle(lab_title) +
    scale_x_continuous(breaks=obj$q) + 
    xlab("Quantiles") +
    ylab("Test statistic") +
    
    
    theme_bw() +
    scale_size(range=c(1, 2, 2)) +
    theme(legend.position = "none") + 
    theme(plot.title = element_text(size = 14)) +
    theme(axis.text.x = element_text(angle = 0, size=12)) +
    theme(axis.text.y = element_text(angle = 0, size=12)) +
    theme(plot.margin = unit(c(5, 5, 0, 0), "mm")) + 
    theme(legend.spacing = unit(0, "cm")) +
    #theme(legend.text = element_text(size = 12)) +
    #geom_abline(intercept = 0, slope=0, colour = "#00CCCC", size = .5) #+ 
    scale_linetype_manual(values=c("solid", "solid")) +
    scale_color_manual(values=cbPalette[c(1,2)])
  
  return(p1)  
}


"lprq2" <- function(x, y, h, tau, x0, type) # modified from lprq, s.t. we can specify where to estimate quantiles
{       xx0 <- x0
        xx <- x
        fv <- xx
        dv <- xx
        for(i in 1:length(xx)) {
          z <- x - xx[i]
          z0 <- xx0-x0[i] 
          wx <- dnorm(z0/h)
          r <- rq(y~z, weights=wx, tau=tau, ci=FALSE)
          fv[i] <- r$coef[1.]
          dv[i] <- r$coef[2.]
        }
        list(xx = xx, fv = fv, dv = dv)
}

