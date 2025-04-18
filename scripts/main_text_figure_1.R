# Reproducing Figure 1: Sensitivity of agricultural yields to temperature over time

# part 1: Plot Response Functions; output: figure_1_part_a.pdf
# part 2: Total Sensitivities; output: figure_1_part_b.pdf

# remove anything in environment and load library
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

# Load necessary helper functions
source('scripts/functions/mkheatmap.R')
source('scripts/functions/piecewise_evaluation_functions.R')

# PLOT RESPONSE FUNCTIONS 
xx=8:40 #range over which to plot response functions for all agricultural plots
gdd=ifelse(xx<30,xx,30)
edd=ifelse(xx>=30,xx-30,0)

# first plot response functions
pdf(file=paste0("fig/main/figure_1_part_a.pdf"),width=12,height=12)
par(mfrow=c(3,3),mar=c(3,4,2,1))


# US corn, wheat, soybeans
vars=c("CORN","SOYBEANS")
basedecade=1950
decades=seq(basedecade,2010,10)
# colz <- allcolz[alldecades%in%decades]
colz = rev(met.brewer("Homer1",length(decades)))
for (k in vars) {
  betas <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_yield_bootstraps.rds'))
  plotresponsefnc_piece(betas,ylim=c(-0.1,0.01),yrs="no",title=paste0("US ",k))  # response function plot
}
k="WHEAT"
betas <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_yield_bootstraps.rds'))
plotresponsefnc_piece(betas,ylim=c(-0.2,0.01),climvar=c("gs_gdd2","gs_edd2"), title="US wheat",yrs="no")


# EU
basedecade=1990
decades=seq(basedecade,2010,10)
# colz <- allcolz[alldecades%in%decades]
colz = rev(met.brewer("Homer1",length(decades)))
for (k in c('wheat','maize')) {
  betas <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_yield_bootstraps.rds'))
  if (k=="wheat") {ylim=c(-0.1,0.01)} else {ylim=c(-0.07,0.01)}
  plotresponsefnc_piece(betas,ylim=ylim,climvar = c(paste0("cw_",k,"_dd_newx030"),paste0("cw_",k,"_dd_newx30")),title=paste0("EU ",k),yrs="no")
}


# BRAZIL
basedecade=1970
decades=seq(basedecade,2010,10)
# colz <- allcolz[alldecades%in%decades]
colz = rev(met.brewer("Homer1",length(decades)))
for (k in c('maize','soybean')) {
  betas <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_yield_bootstraps.rds'))
  if (k=="soybean") {ylim=c(-0.1,0.01)} else {ylim=c(-0.07,0.01)}
  plotresponsefnc_piece(betas,ylim=ylim,climvar = c(paste0("cw_",k,"_dd_newx030"),paste0("cw_",k,"_dd_newx30")),title=paste0("Brazil ",k),yrs="no")
}


# India
basedecade=1990
decades=seq(basedecade,2010,10)
# colz <- allcolz[alldecades%in%decades]
colz = rev(met.brewer("Homer1",length(decades)))
k="wheat"
betas <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_yield_bootstraps.rds'))
plotresponsefnc_piece(betas,ylim=c(-0.07,0.01),climvar = c(paste0("cw_",k,"_dd_newx030"),paste0("cw_",k,"_dd_newx30")),title=paste0("India ",k),yrs="no")


#Ag TFP
coef <- read_rds(file='data/agriculture/processed/TFP/AgTFPBootstraps_periodinteract.rds')
xx=5:30
center=20
decades=seq(1960,2010,10)
colz = rev(met.brewer("Homer1",length(decades)))
plot(1,type="n",xlim=c(5,30),ylim=c(-0.15,0.15),axes=F,xlab="temperature",ylab="ag TFP growth")
abline(h=0,lty=2); axis(1); axis(2,las=1)
for (i in 1:length(decades)) {
  yy=(coef[,i])%*%t(data.frame(xx))
  yy <- t(apply(yy,1,function(x) x-x[xx==center])) #recenter
  toplot <- apply(yy,2,function(x) quantile(x,probs = c(0.05,0.5,0.95))) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  if (i==1) {polygon(c(xx,rev(xx)),c(toplot[1,],rev(toplot[3,])),col=alpha(colz[i],0.3),border = NA)}
  lines(xx,toplot[2,],lwd=2,col=colz[i])
}

dev.off()


########################################################################
# total sensitivities inset plots - combining with above in illustrator
########################################################################
pdf(file="fig/main/figure_1_part_b.pdf",width=12,height=12)
ms=9
par(mfrow=c(3,3),mar=c(ms,ms,ms,ms))


# US
vars=c("CORN","SOYBEANS","WHEAT")
basedecade=1950
decades=seq(basedecade,2010,10)
colz = rev(met.brewer("Homer1",length(decades)))
for (k in vars) {
  betas <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_yield_bootstraps.rds'))
  bgdd <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_gdd_bootstraps.rds'))
  bedd <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_edd_bootstraps.rds'))
  if (k=="WHEAT") {add="2"} else {add=""} #toggle growing season for wheat
  totalsensitivity=betas[,paste0("decade::",decades,":gs_gdd",add)]*bgdd + betas[,paste0("decade::",decades,":gs_edd",add)]*bedd
  plottimeeffects(totalsensitivity,title="total sensitivity",ylab="d(logYield)/dC",col=colz)
}


# EU
basedecade=1990
decades=seq(basedecade,2010,10)
colz = rev(met.brewer("Homer1",length(decades)))
for (k in c('wheat','maize')) {
  betas <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_yield_bootstraps.rds'))
  bgdd <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_gdd_bootstraps.rds'))
  bedd <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_edd_bootstraps.rds'))
  totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
  plottimeeffects(totalsensitivity,title=paste0("EU ",k,", total sensitivity"),ylab="d(logYield)/dC",col=colz)
}


# BRAZIL
basedecade=1970
decades=seq(basedecade,2010,10)
colz = rev(met.brewer("Homer1",length(decades)))
for (k in c('maize','soybean')) {
  betas <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_yield_bootstraps.rds'))
  bgdd <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_gdd_bootstraps.rds'))
  bedd <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_edd_bootstraps.rds'))
  totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
  plottimeeffects(totalsensitivity,title=paste0("Brazil ",k,", total sensitivity"),ylab="d(logYield)/dC",col=colz)
}


# INDIA
basedecade=1990
decades=seq(basedecade,2010,10)
colz = rev(met.brewer("Homer1",length(decades)))
k="wheat"
betas <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_yield_bootstraps.rds'))
bgdd <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_gdd_bootstraps.rds'))
bedd <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_edd_bootstraps.rds'))
totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
plottimeeffects(totalsensitivity,title=paste0("India ",k,", total sensitivity"),ylab="d(logYield)/dC",col=colz,laboff = 1.1)


# Ag TFP
basedecade=1960
decades=seq(basedecade,2010,10)
colz = rev(met.brewer("Homer1",length(decades)))

coef <- read_rds(file='data/agriculture/processed/TFP/AgTFPBootstraps_periodinteract.rds')
plottimeeffects(coef,col=colz,title="",ylab="dGrowth/dC")  

dev.off()


# #########################################################
# # MAKE p-value plot for SI
# #########################################################
# 
# 
# pdf(file="plots/agriculture/Combined_YieldSensitivity_pvals.pdf",width=12,height=12)
# par(mfrow=c(3,3),mar=c(5,5,3,1))
# 
# 
# # US
# vars=c("CORN","SOYBEANS","WHEAT")
# basedecade=1950
# decades=seq(basedecade,2010,10)
# for (k in vars) {
#   betas <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_yield_bootstraps.rds'))
#   bgdd <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_gdd_bootstraps.rds'))
#   bedd <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_edd_bootstraps.rds'))
#   if (k=="WHEAT") {add="2"} else {add=""} #toggle growing season
#   totalsensitivity=betas[,paste0("decade::",decades,":gs_gdd",add)]*bgdd + betas[,paste0("decade::",decades,":gs_edd",add)]*bedd
#   mkheatmap(totalsensitivity,labels=T,cex=0.8,plotratio=T,labs=decades,main=paste0("US ",k)) 
# }
# 
# 
# # EU
# basedecade=1990
# decades=seq(basedecade,2010,10)
# for (k in c('wheat','maize')) {
#   betas <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_yield_bootstraps.rds'))
#   bgdd <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_gdd_bootstraps.rds'))
#   bedd <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_edd_bootstraps.rds'))
#   totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
#   mkheatmap(totalsensitivity,labels=T,cex=0.8,plotratio=T,labs=decades,main=paste0("EU ",k)) 
# }
# 
# 
# # brazil
# basedecade=1970
# decades=seq(basedecade,2010,10)
# for (k in c('maize','soybean')) {
#   betas <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_yield_bootstraps.rds'))
#   bgdd <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_gdd_bootstraps.rds'))
#   bedd <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_edd_bootstraps.rds'))
#   totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
#   mkheatmap(totalsensitivity,labels=T,cex=0.8,plotratio=T,labs=decades,main=paste0("Brazil ",k)) 
# }
# 
# 
# # India wheat
# basedecade=1990
# decades=seq(basedecade,2010,10)
# k="wheat"
# betas <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_yield_bootstraps.rds'))
# bgdd <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_gdd_bootstraps.rds'))
# bedd <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_edd_bootstraps.rds'))
# totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
# mkheatmap(totalsensitivity,labels=T,cex=0.8,plotratio=T,labs=decades,main="India wheat") 
# 
# 
# # Ag TFP
# basedecade=1960
# decades=seq(basedecade,2010,10)
# coef <- read_rds(file='data/agriculture/processed/TFP/AgTFPBootstraps_periodinteract.rds')
# mkheatmap(coef,labels=T,cex=0.8,plotratio=T,labs=decades,main="Global ag TFP") 
# 
# 
# dev.off()
# 
# 
# 
# 
# ##########################################################
# # DECOMPOSE CHANGING SENSITIVITIES
# ##########################################################
# 
# 
# pdf(file="plots/agriculture/Combined_YieldSensitivity_decomposition.pdf",width=8,height=8)
# par(mfrow=c(3,3),mar=c(4,4,2,1))
# 
# 
# # US
# basedecade=1950
# decades=seq(basedecade,2010,10)
# vars=c("CORN","SOYBEANS","WHEAT")
# for (k in vars) {
#   betas <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_yield_bootstraps.rds'))
#   bgdd <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_gdd_bootstraps.rds'))
#   bedd <- read_rds(file=paste0('data/agriculture/processed/usa/Bootstraps/',k,'_edd_bootstraps.rds'))
#   
#   if (k=="WHEAT") {add="2"} else {add=""} #toggle growing season
#   
#   
#   n=length(decades)
#   totalsensitivity=betas[,paste0("decade::",decades,":gs_gdd",add)]*bgdd + betas[,paste0("decade::",decades,":gs_edd",add)]*bedd
#   plottimeeffects(totalsensitivity,title="",ylab="dYield/dC",laboff=0.6,cexlab=0.7,xlim=c(0.5,n*3 + 3))
#   
#   # changing responses
#   chgresponse=betas[,paste0("decade::",decades,":gs_gdd",add)]*matrix(rep(bgdd[,1],n),ncol = n) + betas[,paste0("decade::",decades,":gs_edd",add)]*matrix(rep(bedd[,1],n),ncol = n)
#   plottimeeffects(chgresponse,newplot="no",col="grey30",att=(1:n)+n+1)
#   # changing exposure: base period betas, changing dEDD/dC and dGDD/dC
#   chgexposure=matrix(rep(betas[,paste0("decade::",basedecade,":gs_gdd",add)],n),ncol=n)*bgdd + matrix(rep(betas[,paste0("decade::",basedecade,":gs_edd",add)],n),ncol=n)*bedd
#   plottimeeffects(chgexposure,newplot="no",col="grey70",att=(1:n)+n*2+2)
#   abline(v=c(n+1,n*2+2),col="grey",lty=2,lwd=0.5)
# }
# 
# 
# # EU
# basedecade=1990
# decades=seq(basedecade,2010,10)
# n=length(decades)
# for (k in c('wheat','maize')) {
#   betas <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_yield_bootstraps.rds'))
#   bgdd <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_gdd_bootstraps.rds'))
#   bedd <- read_rds(file=paste0('data/agriculture/processed/eu/Bootstraps/',k,'_edd_bootstraps.rds'))
#   totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
#   plottimeeffects(totalsensitivity,ylab="dYield/dC",laboff=1.1,cexlab=0.7,xlim=c(0.5,n*3 + 3))
#   
#   # changing responses
#   chgresponse=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*matrix(rep(bgdd[,1],n),ncol = n) + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*matrix(rep(bedd[,1],n),ncol = n)
#   plottimeeffects(chgresponse,newplot="no",col="grey30",att=(1:n)+n+1)
#   # changing exposure: base period betas, changing dEDD/dC and dGDD/dC
#   chgexposure=matrix(rep(betas[,paste0("decade::",basedecade,":cw_",k,"_dd_newx030")],n),ncol=n)*bgdd + matrix(rep(betas[,paste0("decade::",basedecade,":cw_",k,"_dd_newx30")],n),ncol=n)*bedd
#   plottimeeffects(chgexposure,newplot="no",col="grey70",att=(1:n)+n*2+2)
#   abline(v=c(n+1,n*2+2),col="grey",lty=2,lwd=0.5)
# }
# 
# 
# # brazil
# basedecade=1970
# decades=seq(basedecade,2010,10)
# n=length(decades)
# for (k in c('maize','soybean')) {
#   betas <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_yield_bootstraps.rds'))
#   bgdd <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_gdd_bootstraps.rds'))
#   bedd <- read_rds(file=paste0('data/agriculture/processed/brazil/Bootstraps/',k,'_edd_bootstraps.rds'))
#   totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
#   plottimeeffects(totalsensitivity,ylab="dYield/dC",laboff=1.1,cexlab=0.7,xlim=c(0.5,n*3 + 3))
#   
#   # changing responses
#   chgresponse=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*matrix(rep(bgdd[,1],n),ncol = n) + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*matrix(rep(bedd[,1],n),ncol = n)
#   plottimeeffects(chgresponse,newplot="no",col="grey30",att=(1:n)+n+1)
#   # changing exposure: base period betas, changing dEDD/dC and dGDD/dC
#   chgexposure=matrix(rep(betas[,paste0("decade::",basedecade,":cw_",k,"_dd_newx030")],n),ncol=n)*bgdd + matrix(rep(betas[,paste0("decade::",basedecade,":cw_",k,"_dd_newx30")],n),ncol=n)*bedd
#   plottimeeffects(chgexposure,newplot="no",col="grey70",att=(1:n)+n*2+2)
#   abline(v=c(n+1,n*2+2),col="grey",lty=2,lwd=0.5)
# }
# 
# 
# # India wheat
# basedecade=1990
# decades=seq(basedecade,2010,10)
# n=length(decades)
# k="wheat"
# betas <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_yield_bootstraps.rds'))
# bgdd <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_gdd_bootstraps.rds'))
# bedd <- read_rds(file=paste0('data/agriculture/processed/india/Bootstraps/',k,'_edd_bootstraps.rds'))
# totalsensitivity=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*bgdd + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*bedd
# plottimeeffects(totalsensitivity,ylab="dYield/dC",laboff=1.1,cexlab=0.7,xlim=c(0.5,n*3 + 3))
# 
# 
# # changing responses
# chgresponse=betas[,paste0("decade::",decades,":cw_",k,"_dd_newx030")]*matrix(rep(bgdd[,1],n),ncol = n) + betas[,paste0("decade::",decades,":cw_",k,"_dd_newx30")]*matrix(rep(bedd[,1],n),ncol = n)
# plottimeeffects(chgresponse,newplot="no",col="grey30",att=(1:n)+n+1)
# # changing exposure: base period betas, changing dEDD/dC and dGDD/dC
# chgexposure=matrix(rep(betas[,paste0("decade::",basedecade,":cw_",k,"_dd_newx030")],n),ncol=n)*bgdd + matrix(rep(betas[,paste0("decade::",basedecade,":cw_",k,"_dd_newx30")],n),ncol=n)*bedd
# plottimeeffects(chgexposure,newplot="no",col="grey70",att=(1:n)+n*2+2)
# abline(v=c(n+1,n*2+2),col="grey",lty=2,lwd=0.5)
# 
# 
# 
# 
# dev.off()
