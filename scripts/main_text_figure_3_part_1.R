# Reproducing Figure 3: Changes in climate impacts on income and economic output over time

# *part 1 (a-b): Temperature and growth in per capita GDP; output: figure_3_panels_a_b.pdf*
# part 2 (c-d): Temperature and personal income in the US; output: figure_3_panels_c_d.pdf
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

#using updated data from loss and damage
dt <- read_rds('data/gdp/temp_gdp_world_panel.rds')
dt <- dt %>% group_by(ISO3) %>% 
  mutate(growth = log(NY.GDP.PCAP.KD) - log(lag(NY.GDP.PCAP.KD,n=1)), year2=year*year) 
clim="era"
dt <- dt %>% mutate(temp = get(paste0(clim,"_mwtemp")),temp2=temp^2,prec=get(paste0(clim,"_mwprecip")),prec2=prec^2)
dt <- dt %>% filter(year>=1960)
dt <- dt %>% mutate(period=floor((year-1960)/20)+1)  #make 3 periods: 1960-1979, 1980-1999, 2000-2020
dt$period[dt$year==2020] <- 3
dt <- dt %>% mutate(gdp = NY.GDP.PCAP.KD*SP.POP.TOTL,pop=SP.POP.TOTL)


#main regression
#summary(mod<-feols(growth ~ i(period,temp) + i(period,temp2) + i(period,prec) + i(period,prec2) | ISO3 + year + ISO3[[year]] + ISO3[[year2]], dt))

# # bootstrap by period regressions
# Nboot=1000
# options(warn = -1) 
# set.seed(94305)
# dts <- dt %>% drop_na(growth,temp,prec)  
# iso = unique(dts$ISO3)
# coefperiod <- coefyr <- c()
# for (i in 1:Nboot) {
#   cl  <- data.frame(ISO3=sample(iso, size = length(iso), replace=T))
#   dtboot <- inner_join(dts,cl,by="ISO3")
#   mod <- feols(growth ~ i(period,temp) + i(period,temp2) + i(period,prec) + i(period,prec2) | ISO3 + year + ISO3[[year]] + ISO3[[year2]], data=dtboot )
#   cf <- coef(mod)[c(paste0("period::",1:3,":temp"),paste0("period::",1:3,":temp2"))]
#   coefperiod <- rbind(coefperiod,cf)
#   mod <- feols(growth ~ (temp + temp2 + prec + prec2)*year | ISO3 + year + ISO3[[year]] + ISO3[[year2]], data=dtboot )
#   cf <- coef(mod)[c("temp","temp2","temp:year","temp2:year")]  
#   coefyr <- rbind(coefyr,cf)
#   print(i)
# }
# write_rds(coefperiod,file="data/gdp/GDPtimeperiod_bootstraps_v2.rds")
# write_rds(coefyr,file="data/gdp/GDPbootstraps_yearmodel.rds")




# calculate % change over time for summary plot
forhist <- dt %>% dplyr::select(temp,gdp,pop) %>% ungroup() %>% drop_na()
temps <- Quantile(forhist$temp,weights = forhist$pop, probs = c(0.01,0.05,0.95,0.99),na.rm = T)

coefyr <- read_rds('data/gdp/GDPbootstraps_yearmodel.rds')
coefyr <- coefyr[,c("temp:year","temp2:year")]
# run base model to get point estimate of point sensitivity, we will measure % change relative to that
mod <- feols(growth ~ (temp + temp2 + prec + prec2)*year | ISO3 + year + ISO3[[year]] + ISO3[[year2]], data=dt )
coefpooled = coef(mod)[c("temp","temp2")]
base=20
xx=c(base,temps) #the 1st, 5th, 95th, and 99th percentile of daily exposure
xxr = t(data.frame(xx,xx^2))
yy=coefpooled%*%xxr 
ybase=yy[,2:length(xx)] - yy[,1] #pooled sensitivities at each xval
yy=coefyr%*%xxr
ydiff=yy[,2:length(xx)] - yy[,1]
ychg=sweep(ydiff,2,abs(ybase),"/")*100 # divide each bootstrap estimate of point sensitivities (rows) by pooled model estimates
# dividing by absolute value of pooled sensitivities to make sure changes have the correct sign (otherwise could divide by negative etc)
# out <- as.data.frame(apply(ychg,2,function(x) quantile(x,probs=c(0.025,0.05,0.5,0.95,0.975))))
out <- data.frame(apply(ychg,2,median_and_pval))
names(out) <- paste0("p",c(1,5,95,99)); rownames(out)=c("median","pval")


#total sensitivities
bins=-10:40 #temperature values to compute binned exposure at
hist_t <- weighted.hist(forhist$temp,forhist$pop,breaks=xx,plot = F)
pwdall <- hist_t$counts #popweighted exposure counts at each bin
#pop weighting maybe the most sensible default here since this is per cap growth
xx <- hist_t$mids
xxr = t(data.frame(1,2*xx)) #converted to polynomials, taking derivative
# first for pooled model
yy=coefpooled%*%xxr #derivative for each decade evaluated at bin midpoints
wts=pwdall #weighting by average exposure over the time period
mdv <- apply(yy,1,function(x) weighted.mean(x,w=wts)) #weighted average derivative
# now for bootstrapped year coefficients
yy=coefyr%*%xxr
mdvyr <- apply(yy,1,function(x) weighted.mean(x,w=wts)) #weighted avg change in derivative by year, per bootstrap
ychg = mdvyr/abs(mdv)*100 #percent change per year
out <- data.frame(out,total=median_and_pval(ychg))
# write_rds(out,file='data/combined_plots/AnnualChanges_GlobalGDP.rds')

#### VERSION 2 THAT LOOKS AT CHANGING SENSITIVITY
pdf(file="fig/main/figure_3_panels_a_b.pdf",width=8,height=4)
par(mfrow=c(1,2),mar=c(4,4,1,2))


#first plot effects over time with confidence interval
center=20
xx=1:30
coef <- read_rds("data/gdp/GDPtimeperiod_bootstraps_v2.rds")
est <- array(dim=c(dim(coef)[1],length(xx),3))
opt <- array(dim=c(dim(coef)[1],3))
for (j in 1:dim(coef)[1]) { #calculate predicted values
  for (k in 1:3) { #looping over periods
    b1=coef[j,paste0("period::",k,":temp")] ; b2=coef[j,paste0("period::",k,":temp2")]
    yy=b1*xx +b2*xx*xx
    est[j,,k] <- yy-yy[xx==center]
  }
}
colz = rev(met.brewer("Homer1",3))
plot(1,type="n",xlim=c(1,30),ylim=c(-0.3,0.05),axes=F,xlab="temperature (C)",ylab="growth")
axis(1)
axis(2,las=1)
# plotting CI for first period, then point est for all
for (k in 1:3) {
  estci <- apply(est[,,k],2,function(x) quantile(x,probs=c(0.025,0.5,0.975)))
  if (k==1) {polygon(c(xx,rev(xx)),c(estci[1,],rev(estci[3,])),col=alpha(colz[k],0.3),border="NA")}
  lines(xx,estci[2,],lwd=2,col=colz[k])
}
att=c(0.05,0.03,0.01)  #yaxis location for optima histograms
text(25,att,c("1961-79","1980-99","2000-20"),pos=4,cex=0.8,col=colz)


# plot GDP and temperature exposure histograms
forhist <- dt %>% dplyr::select(temp,gdp,pop) %>% ungroup() %>% drop_na()
ht <- weighted.hist(forhist$temp, forhist$gdp,plot=F,breaks = seq(-10,40,0.5))
lb=-0.25; mh=0.03
rect(ht$breaks,lb,ht$breaks+1,lb+(ht$counts/max(ht$counts)*mh),col="grey30",border = NA)
# pop weights
ht <- weighted.hist(forhist$temp, forhist$pop,plot=F,breaks = seq(-10,40,0.5))
lb=-0.3; mh=0.03
rect(ht$breaks,lb,ht$breaks+1,lb+(ht$counts/max(ht$counts)*mh),col="red",border = NA)


# now calculate changing sensitivity and its sources
# evaluate derivative of response function at different temperatures, weighting by the amt of
#   GDP at that temperature


mkplot <- function(cf,add=0,col="black",addyr="no") {
  xx=hist_t$mids #values to evaluate derivative at - bin midpoints
  xxr = t(data.frame(1,2*xx)) #converted to polynomials (derivative of quadratic)
  wts=hist_t$counts #counts of obs at each bin
  yy=cf%*%xxr #derivative for each decade evaluated at bin midpoints
  mdv <- apply(yy,1,function(x) weighted.mean(x,wts)) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  toplot <- quantile(mdv,probs = c(0.025,0.5,0.975))
  segments(k+add,toplot[1],k+add,toplot[3],col=col)
  points(k+add,toplot[2],pch=19,col=col)
  if (addyr=="yes") {text(k+0.2,toplot[1],periods[k],cex=0.7)}
  return(mdv)
}


bins=-10:40 #temperature values to compute binned exposure at
periods = c("1961-79","1980-99","2000-20")
plot(1,type="n",xlim=c(0.7,3.3),ylim=c(-0.01,0.007),axes=F,xlab="",main="average dy/dT",ylab="E[dy/dT]")
axis(2,las=1); abline(h=0,lty=2)


out1 <- c() #writing out vectors of derivatives to make pvalue plot
for (k in 1:length(periods)) { #first with GDP weights
  forhist <- dt %>% filter(period==k) %>% group_by(ISO3)  %>% summarize(across(c(temp,gdp,pop),mean,na.rm=T)) %>% drop_na()
  hist_t <- weighted.hist(forhist$temp,forhist$gdp,breaks=bins,plot = F) #GDP weighted exposure counts at each bin
  cf = coef[,c(paste0("period::",k,":temp"),paste0("period::",k,":temp2"))]
  out1 <- cbind(out1,mkplot(cf=cf,col="black",add=0,addyr="yes"))
}
out2 <- c()
for (k in 1:length(periods)) { #now with pop weights
  forhist <- dt %>% filter(period==k) %>% group_by(ISO3)  %>% summarize(across(c(temp,gdp,pop),mean,na.rm=T)) %>% drop_na()
  hist_t <- weighted.hist(forhist$temp,forhist$pop,breaks=bins,plot = F) #GDP weighted exposure counts at each bin
  cf = coef[,c(paste0("period::",k,":temp"),paste0("period::",k,":temp2"))]
  out2 <- cbind(out2,mkplot(cf=cf,col="red",add=0.2,addyr="no"))
}


dev.off()




# make p-value plot, separately for both sets of weights
pdf(file="fig/supplementary/figure_s12_a.pdf",width=3,height=3)
mkheatmap(out1, good =1, labels=T,cex=0.8,plotratio=T,labs=c("'60-79","'80-99","'00-19"),main="GDP-temperature:  GDP weights") #need to have loaded this function from scripts/functions/mkheatmap.R
dev.off()


pdf(file="fig/supplementary/figure_s12_b.pdf",width=3,height=3)
mkheatmap(out2,good=1,labels=T,cex=0.8,plotratio=T,labs=c("'60-79","'80-99","'00-19"),main="GDP-temperature:  pop weights") #need to have loaded this function from scripts/functions/mkheatmap.R
dev.off()


# ###. NOW TRY TO BREAK DOWN CHANGES IN TOTAL SENSITIVITY
# # first using GDP weights, then using pop weights
# 
# 
# forhist <- dt %>% group_by(ISO3,period)  %>% summarize(across(c(temp,gdp,pop),function(x) mean(x,na.rm=T))) 
# forhist <- forhist %>% pivot_wider(names_from = period,values_from=c(temp,gdp,pop))
# 
# 
# mkplot2 <- function(pc,pt,pw,add=0,col="black",addyr="no",at=1) { #pc = period for coefficients, pt=period for temperature, pw=period gdp or population
#   df = data.frame(forhist[paste0("temp_",pt)],forhist[paste0(wtvar,"_",pw)]) %>% drop_na()
#   hist_t <- weighted.hist(df[,1],df[,2],breaks=bins,plot = F) #GDP weighted exposure counts at each bin
#   cf = coef[,c(paste0("period::",pc,":temp"),paste0("period::",pc,":temp2"))]
#   xx=hist_t$mids #values to evaluate derivative at - bin midpoints
#   xxr = t(data.frame(1,2*xx)) #converted to polynomials (derivative of quadratic)
#   wts=hist_t$counts #counts of obs at each bin
#   yy=cf%*%xxr #derivative for each decade evaluated at bin midpoints
#   mdv <- apply(yy,1,function(x) weighted.mean(x,wts)) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
#   toplot <- quantile(mdv,probs = c(0.025,0.5,0.975))
#   segments(at+add,toplot[1],at+add,toplot[3],col=col)
#   points(at+add,toplot[2],pch=19,col=col)
#   if (addyr=="yes") {text(at+0.2,toplot[1],periods[at],cex=0.7)}
# }
# 
# 
# 
# 
# # first for GDP: plot total sensitivity, changing response, changing location, changing climate
# # total sensitivity
# pdf(file="plots/gdp/GDPImpact_decomposition.pdf",width=8,height=4)
# par(mfrow=c(1,2),mar=c(4,4,1,2))
# 
# 
# wtvar="gdp"
# col="black"
# plot(1,type="n",xlim=c(0.9,5),ylim=c(-0.01,0.007),axes=F,xlab="",main="average dy/dT - GDP weights",ylab="E[dy/dT]")
# axis(2,las=1); abline(h=0,lty=2)
# for (k in 1:length(periods)) { 
#   mkplot2(pc=k,pt=k,pw=k,col=col,addyr="no",at=1,add=0.2*k) #total sensitivity
#   mkplot2(pc=k,pt=1,pw=1,col=alpha(col,0.8),addyr="no",at=2,add=0.2*k) #changing response
#   mkplot2(pc=1,pt=1,pw=k,col=alpha(col,0.5),addyr="no",at=3,add=0.2*k) #changing location
#   mkplot2(pc=1,pt=k,pw=1,col=alpha(col,0.2),addyr="no",at=4,add=0.2*k) #changing climate
# }
# abline(v=1.9,lty=2)
# 
# 
# 
# 
# wtvar="pop"
# col="red"
# plot(1,type="n",xlim=c(0.9,5),ylim=c(-0.01,0.007),axes=F,xlab="",main="average dy/dT - population weights",ylab="E[dy/dT]")
# axis(2,las=1); abline(h=0,lty=2)
# for (k in 1:length(periods)) { 
#   mkplot2(pc=k,pt=k,pw=k,col=col,addyr="no",at=1,add=0.2*k) #total sensitivity
#   mkplot2(pc=k,pt=1,pw=1,col=alpha(col,0.8),addyr="no",at=2,add=0.2*k) #changing response
#   mkplot2(pc=1,pt=1,pw=k,col=alpha(col,0.5),addyr="no",at=3,add=0.2*k) #changing location
#   mkplot2(pc=1,pt=k,pw=1,col=alpha(col,0.2),addyr="no",at=4,add=0.2*k) #changing climate
# }
# abline(v=1.9,lty=2)
# dev.off()
