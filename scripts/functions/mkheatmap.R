
# function to make a heatmap of pvalues from two vectors of estimates (e.g. from bootstrapping coefficients)
distdiff <- function(v1,v2) {
#  d = v1-v2
  d = v2-v1
  md = mean(d)
  #difference ratio only meaningful if coeff vectors have the same sign
  if (sign(mean(v1))==sign(mean(v2))) {diffratio = round(mean(v2)/mean(v1),2)} else {diffratio=""}
  pval = 2*min(mean(d>=0),mean(d<=0)) #two sided p-value
 # if (md>0) {pval=sum(d<0)/length(d)} else {pval=sum(d>0)/length(d)} #one sided pvalue
  return(list(md,pval,diffratio))
}

#function to make heatmap. mat is the matrix of values to compare, where comparisons made between columns
# `good` is whether increases in the outcome are good (1) or bad (-1)
mkheatmap <- function(mat,good,labels=FALSE,labs=1:n,cex=1,main="",plotratio=F) { 
  # mat <- toeval
  # labels=T
  # cex=0.8
  # plotratio=T
  # labs=paste0(decades,"s")
  # main=title=paste0(30,"C day")
  # good=1
  
  colpos=c(alpha("blue",c(0.8,0.5,0.3)),"grey80")  #colors when difference is positive
  colneg=c(alpha("red",c(0.8,0.5,0.3)),"grey80") #colors when difference is negative
  pvalbins=c(1,0.1,0.05,0.01,0) #bins of pvalues to use
  n=dim(mat)[2]
  plot(1,xlim=c(1,n),ylim=c(1,n),type="n",xlab="",ylab="",axes=F,main=main)
  for (i in 1:(n-1)) {
    # i <- 1
    for (j in (i+1):n) {
      # j = 3
      out <- distdiff(mat[,i],mat[,j]) #a positive number here means the effect is going up over time, given fnc above
      valbin <- cut(out[[2]],pvalbins,include.lowest = TRUE,labels = F,right=F) #categorize the pval
      # if (out[[1]]<0) {col=colpos[valbin]} else {col=colneg[valbin]}
      if (sign(out[[1]])==good) {col=colpos[valbin]} else {col=colneg[valbin]}
      rect(i,j-1,i+1,j,col=col,border=NA)
      if (plotratio==T) {text(i+0.5,j-0.5,out[[3]])}
    }
  }
  if (labels==TRUE) { #labels the heatmap with numbers corresponding to the period
    axis(1,at=1:(n-1)+0.5,labs[1:(n-1)],lwd=0,cex.axis=cex)
    axis(2,at=1:(n-1)+0.5,labs[2:n],lwd=0,cex.axis=cex,las=1)
  }
}