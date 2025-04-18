# Reproducing Figure 3: Changes in climate impacts on income and economic output over time

# part 1 (a-b): Temperature and growth in per capita GDP; output: figure_3_panels_a_b.pdf
# *part 2 (c-d): Temperature and personal income in the US; output: figure_3_panels_c_d.pdf*
# part 3 (e): Tropical cyclone winds and GDP growth; 
# part 4 (f): Precipitation and flood damages in the US

rm(list=ls())
# remove anything in environment and load library
unique_packages <- c(
  "tidyverse",
  "fixest",
  "MetBrewer",
  "DescTools",
  "plotrix",
  "arrow",
  "data.table",
  "sf",
  "magrittr",
  "haven",
  "cowplot"
)

# Identify packages that are not installed
not_installed <- setdiff(unique_packages, rownames(installed.packages()))

# Install missing packages
if (length(not_installed) > 0) {
  install.packages(not_installed)
}

# Load all packages
lapply(unique_packages, library, character.only = TRUE)
source('scripts/functions/mkheatmap.R')
source('scripts/functions/polynomial_evaluation_functions.R')
load('data/mortality/usa/USTemperatureQuantiles_popweighted_v2.RDA') # pop weighted temperature distributions for evaluation 

# part 2 (c-d): Temperature and personal income in the US; output:
pdf(file="fig/main/figure_3_panels_c_d.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,1,2))


# function to compute and plot average derivative, given a set of polynomial coefficients (vector or matrix), 
#   points at which to eval derivative, and set of weights at those points
#   adj is factor by which to rescale effects - e.g. if wanting to go from days to years 
#   xloc is where to plot on x-axis, "add" gives additional offset, addyr plots the period/set of years as desired
mkplot <- function(coef,locations,weights,adj=1,xloc=1,add=0,col="black",addyr="no",yr="") {
  xx=locations
  xxr = t(data.frame(1,2*xx,3*xx^2,4*xx^3)) #converted to polynomials, taking derivative
  yy=coef%*%xxr #derivative for each decade evaluated at bin midpoints
  mdv <- apply(yy,1,function(x) weighted.mean(x,w=weights)) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  toplot <- quantile(mdv,probs = c(0.025,0.5,0.975))*adj
  segments(xloc+add,toplot[1],xloc+add,toplot[3],col=col)
  points(xloc+add,toplot[2],pch=19,col=col)
  if (addyr=="yes") {text(xloc+0.2,toplot[1],yr,cex=0.7)}
}

# plot overall response function
coef <- read_rds(file='data/gdp/USIncome_decade_polynomial_bootstraps.rds')
# coefp <- read_rds(file='data/gdp/USIncome_pooled_polynomial_bootstraps.rds')
xx=-10:35
xxr = t(data.frame(xx,xx^2,xx^3,xx^4))
center=20
plot(1,type="n",xlim=c(-10,35),ylim=c(-0.001,0.0003),axes=F,xlab="",ylab="log income")
axis(1); axis(2,las=1)
abline(h=0,lty=2)
# plotting main response
# yy=coefp%*%xxr
# yy <- t(apply(yy,1,function(x) x-x[xx==center])) #recenter
# toplot <- apply(yy,2,function(x) quantile(x,probs = c(0.025,0.5,0.975))) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
# polygon(c(xx,rev(xx)),c(toplot[1,],rev(toplot[3,])),col="grey",border = NA)
# lines(xx,toplot[2,])
decades=seq(1970,2010,10)
colz = rev(met.brewer("Homer1",5))
for (d in 1:length(decades)) {
  nm = paste0("decade::",decades[d],":pw_tmean_poly",1:4)
  cf <- coef[,nm]
  yy=cf%*%xxr
  yy <- t(apply(yy,1,function(x) x-x[xx==center])) #recenter
  toplot <- apply(yy,2,function(x) quantile(x,probs = c(0.025,0.5,0.975))) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  if (d==1) {polygon(c(xx,rev(xx)),c(toplot[1,],rev(toplot[3,])),col=alpha(colz[d],0.3),border = NA)}
  lines(xx,toplot[2,],col=colz[d],lwd=2)
}
counts=us_pwdall/max(us_pwdall); breaks=us_mids
lb=-0.001; mh=0.00015
rect(breaks,lb,breaks+1,lb+counts*mh,col="grey80")
text(-10,seq(-0.0006,-0.0009,length.out=5),c("1968-79","1980-89","1990-99", "2000-09", "2010-19"),pos=4,cex=0.8,col=colz)


# plot pop-weighted sensitivity over time
plot(1,type="n",xlim=c(0.7,5.3),ylim=c(-0.012,0),axes=F,xlab="",main="average dy/dT",ylab="E[dy/dT]")
axis(2,las=1); abline(h=0,lty=2)
decades=seq(1970,2010,10)
for (d in 1:length(decades)) {
  nm = paste0("decade::",decades[d],":pw_tmean_poly",1:4)
  cf <- coef[,nm]
  wts <- us_pwd[us_pwd$decade==decades[d],2:dim(us_pwd)[2]]
  mkplot(cf,locations = us_mids,weights=wts,adj=365,xloc=d,addyr="yes",yr=paste0(decades[d],"s"))
}
dev.off()

# make -pvalue chart. need to have loaded coefs and estimated pop-weighted distributions
# function to compute weighted derivative for a set of coefficients
getvals <- function(coef,locations,weights) {
  xx=locations
  xxr = t(data.frame(1,2*xx,3*xx^2,4*xx^3)) #converted to polynomials, taking derivative
  yy=coef%*%xxr #derivative for each decade evaluated at bin midpoints
  mdv <- apply(yy,1,function(x) weighted.mean(x,w=weights)) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  return(mdv)
}
decades=seq(1970,2010,10)
out <- c()
for (d in 1:length(decades)) {
  nm = paste0("decade::",decades[d],":pw_tmean_poly",1:4)
  cf <- coef[,nm]
  wts <- us_pwd[us_pwd$decade==decades[d],2:dim(us_pwd)[2]]
  out <- cbind(out,getvals(cf,locations=us_mids,weights=wts))
}

pdf(file="fig/supplementary/figure_s12_c.pdf",width=3,height=3)
mkheatmap(out, good =1, labels=T,cex=0.8,plotratio=T,labs=c("'70s","'80s","'90s","'00s","'10s"),main="US income - temperature") #need to have loaded this function from scripts/functions/mkheatmap.R
dev.off()