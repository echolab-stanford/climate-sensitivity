# Reproducing Figure 3: Changes in climate impacts on income and economic output over time

# part 1 (a-b): Temperature and growth in per capita GDP; output: figure_3_panels_a_b.pdf
# part 2 (c-d): Temperature and personal income in the US; output: figure_3_panels_c_d.pdf
# part 3 (e): Tropical cyclone winds and GDP growth; ouput: figure_3_panel_e.pdf
# *part 4 (f): Precipitation and flood damages in the US; output: figure_3_panel_f.pdf*


# replication of Davenport, Burke, Diffenbaugh 2021 PNAS


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

source('scripts/functions/predictpoly.R')
source('scripts/functions/mkheatmap.R')
source('scripts/functions/polynomial_evaluation_functions.R')


dt <- read_rds('data/floods/state_panel_data.Rds')


# main model
# feols(damage_value ~ monthly_precip_NORM | STATENAME^month + STATENAME^year ,dt)


# # run bootstraps
# Nboot = 1000
# options(warn = -1) 
# set.seed(94305)
# yperiod <- coefyr <- c()
# dts <- dt %>% drop_na(damage_value,monthly_precip_NORM)
# state = unique(dts$STATENAME)
# for (i in 1:Nboot) {
#   cl  <- data.frame(STATENAME=sample(state, size = length(state), replace=T))
#   dtboot <- inner_join(dts,cl,by="STATENAME")
#   mod <- feols(damage_value ~ i(decade,monthly_precip_NORM) | STATENAME^month + STATENAME^year ,dtboot) 
#   yperiod <- rbind(yperiod,coef(mod))
#   
#   # year model
#   mod <- feols(damage_value ~ monthly_precip_NORM*year | STATENAME^month + STATENAME^year ,dtboot) 
#   coefyr <- rbind(coefyr,coef(mod)["monthly_precip_NORM:year"])
#   print(i)
# }
# write_rds(yperiod,file='data/floods/processed/FloodPeriodModelBootstraps.rds')
# write_rds(coefyr,file='data/floods/processed/FloodBootstraps_yearmodel.rds')
# 
# 
# # calculate percentage change annually and write out
# coefyr <- read_rds('data/floods/processed/FloodBootstraps_yearmodel.rds')
# mod <- feols(damage_value ~ monthly_precip_NORM | STATENAME^month + STATENAME^year ,dt) #pooled model
# ychg = coefyr/abs(coef(mod)['monthly_precip_NORM'])*100 #percent change per year
# out <- data.frame(total=median_and_pval(ychg))
# write_rds(out,file='data/combined_plots/AnnualChanges_USFlooding.rds')


#make a fig
pdf(file="fig/main/figure_3_panel_f.pdf",width=4,height=4)
par(mar=c(4,4,2,1))
colz = rev(met.brewer("Homer1",3))
coef <- read_rds(file='data/floods/processed/FloodPeriodModelBootstraps.rds')
xx=seq(-2,4,0.1)
center=0
plot(1,type="n",xlim=c(-2,4),ylim=c(-4,6),axes=F,xlab="precip (sd)",ylab="ln(normalized flood damage)")
abline(h=0,lty=2); axis(1); axis(2,las=1)


for (i in 1:3) {
  yy=coef[,i]%*%t(data.frame(xx))
  yy <- t(apply(yy,1,function(x) x-x[xx==center])) #recenter
  toplot <- apply(yy,2,function(x) quantile(x,probs = c(0.025,0.5,0.975))) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  if (i==1) {polygon(c(xx,rev(xx)),c(toplot[1,],rev(toplot[3,])),col=alpha(colz[i],0.3),border = NA)}
  lines(xx,toplot[2,],lwd=2,col=colz[i])
}
ht <- hist(dt$monthly_precip_NORM,plot = F,breaks = 50)
bb <- length(ht$breaks)  
att <- -4
scl <- 1
rect(ht$breaks[1:(bb-1)],att,ht$breaks[2:bb],att+ht$counts/max(ht$counts)*scl,col="grey")
text(-2,seq(4,5.5,length.out=3),c("1988-1997","1998-2007","2008-2017"),col=colz,pos=4)
dev.off()




# heatmap of p values
pdf(file="fig/supplementary/figure_s12_e.pdf",width=3,height=3)
mkheatmap(coef, good =1, labels=T,cex=0.8,plotratio=T,labs=c("'88-97","'98-07","'08-19"),main="US flood damages") #need to have loaded this function from scripts/functions/mkheatmap.R
dev.off()












#########################################
# # OLDER 
# # bootstrap polynomial/spline
# Nboot = 1000
# options(warn = -1) 
# set.seed(94305)
# ypred <- yinter <- c()
# dts <- dt %>% drop_na(damage_value,monthly_precip_NORM)
# state = unique(dts$STATENAME)
# xx=seq(-2,6,0.1)
# for (i in 1:Nboot) {
#   cl  <- data.frame(STATENAME=sample(state, size = length(state), replace=T))
#   dtboot <- inner_join(dts,cl,by="STATENAME")
#   # mod <- feols(damage_value ~ poly(monthly_precip_NORM,4,raw=T) | STATENAME^month + STATENAME^year ,dtboot)
#   # ypred <- cbind(ypred,predictpoly(model=mod,xval=xx,fit="poly",order=4,recenter=0))
#   mod <- feols(damage_value ~ monthly_precip_NORM | STATENAME^month + STATENAME^year ,dtboot)
#   yy = coef(mod)['monthly_precip_NORM']*xx
#   ypred <- cbind(ypred,yy-yy[xx==0])
#   mod <- feols(damage_value ~ monthly_precip_NORM*year | STATENAME^month + STATENAME^year ,dtboot)
#   yinter <- rbind(yinter,c(coef(mod)['monthly_precip_NORM'],coef(mod)['monthly_precip_NORM:year']))
#   print(i)
# }
# write_rds(data.frame(xx,ypred),file='data/floods/processed/FloodPooledModelBootstraps.rds')
# write_rds(yinter,file="data/floods/processed/FloodBootstraps_linearinteract.rds")
# 
# 
# # make plot
# pdf(file="plots/floods/FloodImpact.pdf",width=9,height=3)
# par(mfrow=c(1,3),mar=c(4,4,1,2))
# 
# 
# ypred <- read_rds('data/floods/processed/FloodPooledModelBootstraps.rds')
# xx = ypred$xx
# xr=quantile(dt$monthly_precip_NORM,probs=c(0.01,0.995))
# plot(1,type="n",xlim=xr,ylim=c(-3,6),axes=F,xlab="precip (sd)",ylab="ln(normalized flood damage)")
# abline(h=0,lty=2)
# ci <- apply(ypred,1,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
# polygon(c(xx,rev(xx)),c(ci[1,],rev(ci[3,])),col="lightblue",border=NA)
# axis(1)
# axis(2,las=1)
# lines(xx,ci[2,],col="black",lwd=1.5)
# 
# 
# ht <- hist(dts$monthly_precip_NORM,plot = F,breaks = 50)
# bb <- length(ht$breaks)  
# att <- -3
# scl <- 1
# rect(ht$breaks[1:(bb-1)],att,ht$breaks[2:bb],att+ht$counts/max(ht$counts)*scl,col="grey")
# # plot 4th order
# mod <- feols(damage_value ~ poly(monthly_precip_NORM,4,raw=T) | STATENAME^month + STATENAME^year ,dt)
# yy=predictpoly(model=mod,xval=xx,fit="poly",order=4,recenter=0)
# lines(xx,yy,col="black",lwd=1.5,lty=2)
# 
# 
# 
# 
# # decade interactions    
# summary(mod <- feols(damage_value ~ i(decade,monthly_precip_NORM) | STATENAME^month + STATENAME^year ,dt))
# coefs <- coef(mod)
# se <- se(mod)
# cihi=coefs+1.96*se
# cilo=coefs-1.96*se
# plot(1:3,coefs,pch=19,axes=F,xlab="",ylab="change in log damage per sd",ylim = c(1,1.6),xlim=c(0.5,3.5))
# abline(h=0,lty=2)
# segments(1:3,cilo,1:3,cihi)
# axis(2,las=1)
# text(1:3,cihi+0.05,c("1988-1997","1998-2007","2008-2017"))
# 
# 
# #effect over time
# yinter <- read_rds('data/floods/processed/FloodBootstraps_linearinteract.rds')
# histy <- hist(yinter[,'monthly_precip_NORM:year']*100,breaks=30,plot=F)
# ymx = 180
# xloc <- median(histy$breaks)
# plot(histy,main="",las=1,xlab="change in marginal effect per year (x100)",ylab="bootstrap frequency",ylim=c(0,ymx))
# ci <- round(quantile(yinter[,'monthly_precip_NORM:year']*100,probs=c(0.025,0.975)),2)
# abline(v=0,lty=2,col="red")
# text(xloc,ymx,paste0("95% CI = [",ci[1],",",ci[2],"]"),pos=1)
# 
# 
# dev.off()
