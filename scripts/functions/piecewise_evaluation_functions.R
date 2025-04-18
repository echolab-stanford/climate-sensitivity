
# piecewise linear evaluation functions for ag regressions

# function to take in a matrix of estimates, where columns are coeff estiamtes in different decades
#   and rows are bootstrap replicates, and make a plot of coeff estimates and CI over time
# markeroff allows you to nudge plotted estimates to left or right
# laboff is offset for plotted labels on point + whiskers
plottimeeffects <- function(est,title="",xlim=c(0.5,n+0.5),ylab="change in log yield per degree day",
                            col="black",setylim="no",ylim="",laboff=1.15,newplot="yes",att=1:n,markeroff=0,axs=1,yrs="yes",cexlab=1.1) {
  estsum <- apply(est,2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
  rg=range(estsum)
  if (setylim=="no") {ylim=c(min(0,rg[1]*1.2),max(0,rg[2])*1.2) } #set dynamic yrange
  else {ylim=ylim} 
  n=dim(est)[2]
  if (newplot=="yes") {
    plot(1,type="n",xlim=xlim,ylim=ylim,axes=F,ylab=ylab,main="",xlab="")
    abline(h=0,lty=2)
    axis(2,las=1,cex.axis=axs)
    if (yrs=="yes") {text(att,estsum[1,]*laboff,paste0(seq(basedecade,2010,10),"s"),cex=cexlab,srt=90)}
    mtext(title,3,adj=0,font=2,cex=0.8)
  }
  segments(att+markeroff,estsum[1,],att+markeroff,estsum[3,],col=col)
  points(att+markeroff,estsum[2,],pch=19,col=col)
}


# function to plot piecewise 
plotresponsefnc_piece <- function(betas,climvar=c("gs_gdd","gs_edd"),ylim=c(-0.1,0.01),ylab="log yield",col=colz,title="",yrs="yes") {
  plot(1,type="n",xlim=c(8,40),ylim=ylim,axes=F,xlab="temperature (C)",ylab=ylab)
  axis(1); axis(2,las=1)
  abline(h=0,lty=2)
  for (j in 1:length(decades)) { #now plot response functions for other decades without CI
    y <- betas[,paste0("decade::",decades[j],":",climvar)]
    ypred <- c()
    for (i in 1:dim(y)[1]) {
      ypred <- cbind(ypred,gdd*y[i,1]+edd*y[i,2])
    }
    ci <- apply(ypred,1,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
    if (decades[j]==basedecade) { #plot CI for base decade
      polygon(c(xx,rev(xx)),c(ci[1,],rev(ci[3,])),col=alpha(col[1],0.3),border=NA)
    }
    lines(xx,ci[2,],lwd=2,col=col[j])
  }
  if (yrs=="yes") {text(10,seq(min(ylim)*0.7,min(ylim),length.out=length(decades)),paste0(decades,"s"),col=col,cex=0.8)}
  mtext(title,3,adj=0,font=2,cex=0.8)
  
}


# function for quick plotting of marginal effects of a variable over time, from a single regression
plottimemarginals <- function(model,varname,title="",xlim=c(0.5,n+0.5),decades=seq(1990,2010,10)) {
  vars <- grep(varname,names(coef(model)))
  coefs <- coef(model)[vars]
  se <- se(model)[vars]
  cilo= coefs - 1.96*se
  cihi= coefs + 1.96*se
  rg=range(cilo,cihi)
  ylim=c(min(0,rg[1]*1.2),max(0,rg[2])*1.2) #set dynamic yrange
  n=length(coefs)
  plot(1,type="n",xlim=xlim,ylim=ylim,axes=F,ylab="change in log yield per EDD",main="",xlab="")
  abline(h=0,lty=2)
  axis(2,las=1)
  segments(1:n,cilo,1:n,cihi)
  points(1:n,coefs,pch=19)
  text(1:n,cihi+.0005,paste0(decades,"s"),cex=0.8)
  mtext(title,3,adj=0,font=2,cex=0.8)
}

# function to compute median and pvalue over a matrix of bootstrapped coefficients
median_and_pval <- function(x) {
  med <- median(x)
  pvals <- 2*min(mean(x>=0),mean(x<=0)) #compute 2-tailed pval
  rbind(med,pvals)
}

