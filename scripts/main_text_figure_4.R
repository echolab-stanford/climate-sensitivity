# Reproducing Figure 4: Changes in climate impacts on conflict, violence, and injury

# combined conflict, crime, injury, suicide plot
# Uses bootstrap estimates from three scrips:
#   african conflict: conflict/MakeConflictFigure_v2.R
#   crime/USHomicideSuicideIngury_polynomial_v2
#   crime/USCrimeFigures_polynomial_v2


rm(list = ls())
gc()


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
load('data/mortality/usa/USTemperatureQuantiles_popweighted_v2.RDA') #pop weighted temperature distributions for evaluation

# temperature distribution for crime sample
dt <- read_rds('data/crime/US_crime_pwtemperature_poly_panel_1980_2019.rds')
dt <- dt %>% rename(fips=fips_state_county_code)
pwd_crime <- dt %>% dplyr::select(fips,pop,year,decade,pw_temp_below_neg20:pw_temp_above_45) %>% 
  drop_na(pop) %>% ungroup() %>% group_by(decade) %>% 
  summarise(across(starts_with("pw_temp"),~weighted.mean(.,w=pop,na.rm=T)))
pwdall_crime <- dt %>% dplyr::select(fips,pop,year,decade,pw_temp_below_neg20:pw_temp_above_45) %>% 
  drop_na(pop) %>% ungroup() %>% 
  summarise(across(starts_with("pw_temp"),~weighted.mean(.,w=pop,na.rm=T)))




pdf(file=paste0("fig/main/figure_4.pdf"),height=7,width=8)
par(mar=c(3,5,2,2),mfrow=c(2,2))


# african conflict
coef <- read_rds(file='data/conflict/processed/Conflict_Decade_Bootstraps.rds')
afr <- read_rds(file='data/conflict/Africa_yearly.rds')
decades=c(1990,2000,2020)
xx=10:35 
center=20
colz = rev(met.brewer("Homer1",3))
plot(1,type="n",xlim=c(10,35),ylim=c(-0.1,0.15),axes=F,xlab="temperature",ylab="conflict risk",main="African civil conflict",yaxs="i")
abline(h=0,lty=2); axis(1); axis(2,las=1)
for (i in 1:length(decades)) {
  yy=(coef[,i])%*%t(data.frame(xx))
  yy <- t(apply(yy,1,function(x) x-x[xx==center])) #recenter
  toplot <- apply(yy,2,function(x) quantile(x,probs = c(0.025,0.5,0.975))) #weighted average derivative for each bootstrap, where weights are gdp- or pop-weighted exposure counts at each bin
  if (i==1) {polygon(c(xx,rev(xx)),c(toplot[1,],rev(toplot[3,])),col=alpha(colz[i],0.3),border = NA)}
  lines(xx,toplot[2,],lwd=2,col=colz[i])
}
ht <- hist(afr$temp,plot=F,breaks = 50)
bb <- length(ht$breaks)  
att <- -0.1; scl <- 0.02
rect(ht$breaks[1:(bb-1)],att,ht$breaks[2:bb],att+ht$counts/max(ht$counts)*scl,col="grey")
text(10,seq(0.1,0.14,length.out=3),decades,col=colz,pos=4)


# violent crime
coef <- read_rds(file='data/crime/USViolentCrimeDecade_polynomial_bootstrap.rds')
decades = seq(1980,2010,10)
plotresponsefnc(coef=coef,xx=-14:40,ylim=c(-0.015,0.015),ylab="log monthly violent crime rate", decades=decades,
                addtext="yes",textattx=-10,textatty=seq(0.01,0.014,length.out=4),textlab=paste0(decades,"s"),
                addhist="yes",histlb=-0.015,histh=0.003,histbreaks=us_mids,histcounts=pwdall_crime,title="US violent crime")


decades = seq(1970,2010,10)
coef <- read_rds(file=paste0('data/crime/US_unintentional_injuries_Decade_polynomial_bootstrap.rds'))
plotresponsefnc(coef=coef,xx=-14:40,ylim=c(-0.005,0.02),ylab="log monthly injuries mortality", decades=decades,
                addtext="yes",textattx=-10,textatty=seq(0.01,0.015,length.out=5),textlab=paste0(decades,"s"),
                addhist="yes",histlb=-0.005,histh=0.002,histbreaks=us_mids,histcounts=us_pwdall, title="US injury mortality")
coef <- read_rds(file=paste0('data/crime/US_suicide_Decade_polynomial_bootstrap.rds'))
plotresponsefnc(coef=coef,xx=-14:40,ylim=c(-0.01,0.02),ylab="log monthly suicide mortality", decades=decades,
                addtext="yes",textattx=-10,textatty=seq(0.01,0.015,length.out=5),textlab=paste0(decades,"s"),
                addhist="yes",histlb=-0.01,histh=0.003,histbreaks=us_mids,histcounts=us_pwdall, title="US suicide")
dev.off()




# SI FIGURES SHOWING POINT SENSITIVITY, TOTAL SENSITIVITY, PVALUES
# assumes above stuff loaded
pdf(file=paste0("fig/supplementary/figure_s13.pdf"),height=10,width=9)
par(mar=c(3,5,2,3),mfrow=c(4,3))


# conflict
coef <- read_rds(file='data/conflict/processed/Conflict_Decade_Bootstraps.rds')
afr <- read_rds(file='data/conflict/Africa_yearly.rds')
periodmean = afr %>% group_by(decade) %>% summarise(baserate=mean(ged_binary,na.rm = T)) 
decades=c(1990,2000,2020)
plot(1,type="n",axes=F,xlab="",ylab="") #empty plot since we dont need a point sensitivity plot
# coefficients over time
coefs <- apply(coef,2,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
plot(1:3,coefs[2,]/periodmean$baserate*100,pch=19,axes=F,xlab="",ylab="% change in conflict per +1C",ylim = c(-0.1,0.55)*100,xlim=c(0.5,3.5),main="average dy/dT")
abline(h=0,lty=2)
segments(1:3,coefs[1,]/periodmean$baserate*100,1:3,coefs[3,]/periodmean$baserate*100)
axis(2,las=1)
text(1:3,coefs[3,]/periodmean$baserate*100+3,paste0(decades,"s"))


cf = t(apply(coef,1,function(x) x/periodmean$baserate*100)) #adjust coefficients by period specific base rates
mkheatmap(cf,labels=T,cex=0.8,plotratio=T,labs=decades,main="African civil conflict",good=-1) 


coef <- read_rds(file='data/crime/USViolentCrimeDecade_polynomial_bootstrap.rds')
decades = seq(1980,2010,10)
plotpointsensitivities(coef=coef,base_mmt=us_decade_mmt,xvals=c(-10,0,30,35),decades = decades,xname="pw_tmean_poly",ylim=c(-0.013,0.01),title="US violent crime")
plottotalsensitivity(coef=coef,pwd=pwd_crime,xname="pw_tmean_poly",locations=us_mids,decades=decades,ylim=c(0,0.0005),addyrs="yes")
pvalpolyplot(coef=coef,xx=us_mids,xname="pw_tmean_poly",pwd=pwd_crime,title="",good=-1)


decades = seq(1970,2010,10)
coef <- read_rds(file=paste0('data/crime/US_unintentional_injuries_Decade_polynomial_bootstrap.rds'))
plotpointsensitivities(coef=coef,base_mmt =us_decade_mmt,xvals=c(-10,0,30,35),decades = decades,xname="pw_tmean_poly",ylim=c(-0.005,0.023),title="US unintentional injury mortality")
plottotalsensitivity(coef=coef,pwd=us_pwd,xname="pw_tmean_poly",locations=us_mids,decades=decades,ylim=c(0,0.0004),addyrs="yes")
pvalpolyplot(coef=coef,xx=us_mids,xname="pw_tmean_poly",pwd=us_pwd,title="",good=-1)


coef <- read_rds(file=paste0('data/crime/US_suicide_Decade_polynomial_bootstrap.rds'))
plotpointsensitivities(coef=coef,base_mmt = us_decade_mmt,xvals=c(-10,0,30,35),decades = decades,xname="pw_tmean_poly",ylim=c(-0.005,0.01),title="US suicide")
plottotalsensitivity(coef=coef,pwd=us_pwd,xname="pw_tmean_poly",locations=us_mids,decades=decades,ylim=c(0,0.0004),addyrs="yes")
pvalpolyplot(coef=coef,xx=us_mids,xname="pw_tmean_poly",pwd=us_pwd,title="",good=-1)


dev.off()
