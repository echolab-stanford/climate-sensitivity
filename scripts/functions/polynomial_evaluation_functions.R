
# a bunch of polynomial evaluation functions that we use over and over for the models where we estimate
# 4th order polynomials across periods


# function to plot response functions by decade, assuming that we estimated a 4th order polynomial
#  functions are centered on a value, and by default the CI is plotted for the first period but not others
plotresponsefnc <- function(coef=coef,
                            xx,
                            xname="pw_tmean_poly",
                            ylim,ylab,
                            decades,
                            center=20,
                            addtext="no",
                            textattx="",
                            textatty="",
                            textlab="",
                            title="",
                            addhist="no",
                            histlb="",
                            histh="",
                            histbreaks=mids,
                            histcounts=pwdall,
                            colz=rev(met.brewer("Homer1",length(decades)))) {
  # coef = coef
  # xx=-15:40
  # center=20
  # xname="pw_tmean_poly"
  # ylim=c(-0.01,0.03)
  # ylab="log annual mortality rate"
  # decades=decades
  # addtext="yes"
  # textattx=-5
  # textatty=seq(0.002,0.0005,length.out=length(decades))
  # textlab=paste0(decades,"s")
  # addhist="yes"
  # histlb=-0.01
  # histh=0.005
  # histbreaks=us_mids
  # histcounts=us_pwdall
  # title="Colombia monthly mortality"
  # colz=rev(met.brewer("Homer1",length(decades)))
  
  xxr = t(data.frame(xx,xx^2,xx^3,xx^4))
  plot(1,type="n",xlim=range(xx),ylim=ylim,axes=F,xlab="temperature (C)",ylab=ylab,main=title,yaxs="i")
  axis(1); axis(2,las=1); abline(h=0,lty=2)
  
  df_list <- list() 
  for (i in 1:length(decades)) {
    # i <- 1
    coefs <- coef[,paste0("decade::",decades[i],":",xname,1:4)]
    yy = coefs%*%xxr #rows are bootstraps, columns are sensitivity at each temperature
    yy <- t(apply(yy, 1, function(row) {row - row[xx==center]}))
    ci <- apply(yy,2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
    if (i==1) {
      polygon(c(xx,rev(xx)),c(ci[1,],rev(ci[3,])),col=alpha(colz[i],0.3),border=NA)
    }
    lines(xx,ci[2,],col=colz[i],lwd=2)
    
    df <- as.data.frame(t(ci))
    colnames(df) <- c("pred_025", "pred_50", "pred_975")
    df$bins <- xx
    df$decade <- decades[i]
    
    df_list[[as.character(decades[i])]] <- df
  }
  
  if (addtext=="yes") {text(textattx,textatty,textlab,cex=0.8,col=colz)}
  if (addhist=="yes") {
    counts=histcounts/max(histcounts)
    rect(histbreaks,histlb,histbreaks+1,histlb+counts*histh,col="grey80")}
  
  # finally rbind decade response
  final_out <- do.call(rbind, df_list)
  return(final_out)
}

# version of above that can just be run on feols output - does not plot confidence interval
plotresponsefnc_simple <- function(coef=coef,xx,xname="pw_tmean_poly",ylim,ylab,decades,center=20,addtext="no",textattx="",textatty="",textlab="",title="",
                            addhist="no",histlb="",histh="",histbreaks=mids,histcounts=pwdall,colz=rev(met.brewer("Homer1",length(decades)))) {
  xxr = t(data.frame(xx,xx^2,xx^3,xx^4))
  plot(1,type="n",xlim=range(xx),ylim=ylim,axes=F,xlab="temperature (C)",ylab=ylab,main=title,yaxs="i")
  axis(1); axis(2,las=1); abline(h=0,lty=2)
  for (i in 1:length(decades)) {
    coefs <- coef[paste0("decade::",decades[i],":",xname,1:4)]
    yy=coefs%*%xxr 
    yy <- yy - yy[xx==center]
    lines(xx,yy,col=colz[i],lwd=2)
  }
  if (addtext=="yes") {text(textattx,textatty,textlab,cex=0.8,col=colz)}
  if (addhist=="yes") {
    counts=histcounts/max(histcounts)
    rect(histbreaks,histlb,histbreaks+1,histlb+counts*histh,col="grey80")}
}


# function to compute sensitivities at different parts of resposne curve, i.e. f(X) - f(base)
#   where base is some baseline temperature to compare to 
plotpointsensitivities <- function(coef=coef,base_mmt=20,xvals,decades,xname,ylim,colz=rev(met.brewer("Homer1",length(decades))),title="") {
  
  # coef=coef
  # base_mmt = us_decade_mmt
  # xvals=c(-10,0,30,35)
  # decades = seq(1970,2010,10)
  # xname="pw_tmean_poly"
  # ylim=c(0,0.007)
  # title="point sensitivities"
  # colz=rev(met.brewer("Homer1",length(decades)))
  
  # par(mar=c(3,5,2,2), mfrow=c(1,3))
  nn=length(xvals)
  plot(1,type="n",xlim=c(0.8,nn+0.5),ylim=ylim,axes=F,xlab="",ylab=paste0("effect of additional day, relative to period-specific C day"),main=title)
  axis(2,las=1); abline(h=0,lty=2)
  
  df_list <- list()
  
  for (k in 1:length(decades)) {
    # k <- 1
    base = base_mmt[[as.character(decades[k])]]
    xx=c(base,xvals) #where we want to evaluate sensitivities at period specific mmt
    xxr = t(data.frame(xx,xx^2,xx^3,xx^4))
    
    coefs <- coef[,paste0("decade::",decades[k],":",xname,1:4)]
    yy=coefs%*%xxr #rows are bootstraps, columns are sensitivity at each temperature
    ydiff =yy[,2:length(xx)] - yy[,1] #subtract off effect at base temperature
    toplot <- apply(ydiff,2,function(x) quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)) #compute quantiles across bootstraps for each temp
    segments(1:nn+k*0.1,toplot[1,1:nn],1:nn+k*0.1,toplot[3,1:nn],col=colz[k])
    points(1:nn+k*0.1,toplot[2,1:nn],pch=19,col=colz[k])
    
    df <- as.data.frame(t(toplot))
    colnames(df) <- c("p025", "p50", "p975")
    df$eval_temp <- paste0(xx[2:length(xx)],"C day")
    df$mmt <- base
    df$decade <- decades[k]
    
    df_list[[as.character(decades[k])]] <- df
  }
  
  text(1:nn+0.3,max(ylim),paste0(xx[2:length(xx)],"C day"),cex=0.8)
  abline(v=2:length(decades)-0.2,lwd=0.5,lty=2,col="grey")
  
  final_out <- do.call(rbind, df_list)
  return(final_out)
}

# function to compute and plot average derivative, given a set of polynomial coefficients (vector or matrix), 
#   points at which to eval derivative, and set of weights at those points
#   adj is factor by which to rescale effects - e.g. if wanting to go from days to years 
#   xloc is where to plot on x-axis, "add" gives additional offset, addyr plots the period/set of years as desired
plottotalsensitivity <- function(coef=coef,pwd=pwd,locations,xname,adj=1,add=0,col="black",decades=decades,srtyr=0,
                                 ylim,xlim=c(0.7,length(decades)+0.3),addyrs="no",yrs=paste0(decades,"s"),fixedresponse="no",fixedexposure="no",newplot="yes") {

  # coef=coef
  # pwd=pwd
  # locations = mids
  # xname = "pw_tmean_poly"
  # adj=1
  # add=0
  # col="black"
  # decades=decades
  # srtyr=0
  # ylim = c(-0.0005, -0.00018)
  # xlim=c(0.7,length(decades)+0.3)
  # addyrs="no"
  # yrs=paste0(decades,"s")
  # fixedresponse="no"
  # fixedexposure="no"
  # newplot="yes"
  
  xx=locations
  xxr = t(data.frame(1,2*xx,3*xx^2,4*xx^3)) #converted to polynomials, taking derivative
  if (newplot=="yes") {
    plot(1,type="n",xlim=xlim,ylim=ylim,axes=F,xlab="",main="total sensitivity to warming", ylab="E[dy/dT]")
    axis(2,las=1); abline(h=0,lty=2)
  }
  
  df_list <- list()
  
  for (k in 1:length(decades)) {
    # k <- 3
    coefs <- coef[,paste0("decade::",decades[k],":",xname,1:4)]
    wts <- pwd[pwd$decade==decades[k],2:dim(pwd)[2]]
    if (fixedresponse=="yes") {coefs <- coef[,paste0("decade::",decades[1],":",xname,1:4)]}
    if (fixedexposure=="yes") {wts <- pwd[pwd$decade==decades[1],2:dim(pwd)[2]]}
    yy=coefs%*%xxr #derivative for each decade evaluated at bin midpoints
    mdv <- apply(yy,1,function(x) weighted.mean(x,w=wts)) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
    toplot <- quantile(mdv,probs = c(0.025,0.5,0.975))*adj
    segments(k+add,toplot[1],k+add,toplot[3],col=col)
    points(k+add,toplot[2],pch=19,col=col)
    
    df <- as.data.frame(t(toplot))
    colnames(df) <- c("p025", "p50", "p975")
    df$decade <- decades[k]
    
    df_list[[as.character(decades[k])]] <- df
  }
  if (addyrs=="yes") {text(1:length(decades),max(ylim),yrs,cex=0.7,srt=srtyr)}
  
  final_out <- do.call(rbind, df_list)
  return(final_out)
}

# wrapper for heatmap p-value function when we have run polynomial regression - applying it to total sensitivity
# calculates total sensitivity for each decade in model across all bootstraps, then computes heatmap
#   of pvalues per mkheatmap function
# `good` is whether increases in the outcome are good (1) or bad (-1)
pvalpolyplot <- function(coef=coef,xx=mids,xname,pwd=pwd,title="",good=1) {
  #xx=mids
  xxr = t(data.frame(1,2*xx,3*xx^2,4*xx^3)) #converted to polynomials, taking derivative
  toeval <- c()
  for (k in 1:length(decades)) {
    coefs <- coef[,paste0("decade::",decades[k],":",xname,1:4)]
    wts <- pwd[pwd$decade==decades[k],2:dim(pwd)[2]]
    yy=coefs%*%xxr
    mdv <- apply(yy,1,function(x) weighted.mean(x,w=wts)) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
    toeval <- cbind(toeval,mdv)
  }
  mkheatmap(toeval,labels=T,cex=0.8,plotratio=T,labs=paste0(decades,"s"),main=title,good=good) #need to have loaded this function from scripts/functions/mkheatmap.R
}


pvalpointsensitivities <- function(coef=coef,
                                   base=20,
                                   xval,
                                   decades,
                                   xname,
                                   title="", 
                                   good = 1) {
  # coef=coef
  # base=20
  # xval=30
  # decades = decades
  # xname="pw_tmean_poly"
  # title=paste0(xval,"C day")
  
  toeval <- c()
  for (k in 1:length(decades)) {
    # k <- 3
    base_mmt <- base[[as.character(decades[k])]]
    xx=c(base_mmt,xval) #where we want to evaluate sensitivities
    xxr = t(data.frame(xx,xx^2,xx^3,xx^4))
    coefs <- coef[,paste0("decade::",decades[k],":",xname,1:4)]
    yy=coefs%*%xxr #rows are bootstraps, columns are sensitivity at each temperature
    ydiff=yy[,2:length(xx)] - yy[,1] #subtract off effect at base temperature
    toeval <- cbind(toeval,ydiff)
   }
  mkheatmap(toeval,labels=T,cex=0.8,plotratio=T,labs=paste0(decades,"s"),main=title, good=good) #need to have loaded this function from scripts/functions/mkheatmap.R
}

# function to compute median and pvalue over a matrix of bootstrapped coefficients
median_and_pval <- function(x) {
  med <- median(x)
  pvals <- 2*min(mean(x>=0),mean(x<=0)) #compute 2-tailed pval
  rbind(med,pvals)
}

# function to find mmt from estimated coefficients
find_mmt <- function(coefs, mmt_min, mmt_max, temp_range, grep_name, additional_grep = NULL) {
  
  # coefs = coefs
  # mmt_min = temp_quantiles[1]
  # mmt_max = temp_quantiles[4]
  # temp_range = -50:50
  # grep_name = "tmean_poly"
  # additional_grep = NULL
  
  # Create polynomial matrix for prediction
  pred_temps_poly <- poly(temp_range, degree = 4, raw = TRUE)
  
  # Extract temperature-related coefficients and compute predicted values
  if (is.null(additional_grep)) {
    temp_coefs <- coefs[grepl(grep_name, names(coefs))]
  } else {
    temp_coefs <- coefs[grepl(pattern = grep_name, x = names(coefs)) & grepl(pattern = additional_grep, x = names(coefs))]
  }
  
  orig_names <- names(temp_coefs)
  
  # add estimated coefficients to each polynomials
  poly1_coefs <- sum(temp_coefs[grepl(paste0(grep_name, 1), names(temp_coefs))])
  poly2_coefs <- sum(temp_coefs[grepl(paste0(grep_name, 2), names(temp_coefs))])
  poly3_coefs <- sum(temp_coefs[grepl(paste0(grep_name, 3), names(temp_coefs))])
  poly4_coefs <- sum(temp_coefs[grepl(paste0(grep_name, 4), names(temp_coefs))])
  
  # append new coefficients while preserving the original named vector format
  temp_coefs <- c(poly1_coefs, poly2_coefs, poly3_coefs, poly4_coefs)
  names(temp_coefs) <- orig_names
  
  pred_vector <- c(pred_temps_poly %*% temp_coefs)
  
  # Restrict MMT to be within the observed temperature distribution
  # mmt_min <- quantile(observed_temps, mmt_quantiles[1], na.rm = TRUE)
  # mmt_max <- quantile(observed_temps, mmt_quantiles[2], na.rm = TRUE)
  
  # Find potential MMT options and their predicted values
  mmt_options <- which(temp_range < mmt_max & temp_range > mmt_min)
  mmt_options_values <- pred_vector[mmt_options]
  
  # Find the temperature corresponding to the minimum predicted value
  mmt_index_in_options <- which.min(mmt_options_values)
  mmt <- temp_range[mmt_options][mmt_index_in_options]
  
  # # Scale predicted temperatures around MMT
  # pred_temps_mmt <- poly(mmt, degree = degree, raw = TRUE)
  # pred_temps_scaled <- scale(pred_temps_poly, center = pred_temps_mmt, scale = FALSE)
  
  return(mmt)
}
