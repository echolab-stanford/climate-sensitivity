# Reproducing Figure 3: Changes in climate impacts on income and economic output over time

# part 1 (a-b): Temperature and growth in per capita GDP; output: figure_3_panels_a_b.pdf
# part 2 (c-d): Temperature and personal income in the US; output: figure_3_panels_c_d.pdf
# *part 3 (e): Tropical cyclone winds and GDP growth; ouput: figure_3_panel_e.pdf*
# part 4 (f): Precipitation and flood damages in the US

# Script to estimate cumulative impacts of tropical cyclones on growth

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


dt <- read_csv('data/cyclone/TC_growth_panel_1950-2022.csv')
dt <- dt %>% group_by(iso) %>% mutate(growth = log(gdppercap) - lag(log(gdppercap))) %>% filter(year<2020)
dt <- dt %>% mutate(period=case_when(year>=1965 & year<=1983 ~ 1,
                                     year>=1984 & year<=2001 ~ 2,
                                     year>=2002 & year<=2019 ~ 3)) 

# 
# # bootstrap by-period regression to get confidence intervals
# Nboot=1000
# options(warn = -1) 
# set.seed(94305)
# df <- dt %>% drop_na(growth,wind)
# #need to pre-generate lags, or bootstrap fails
# setDT(df); setkey(df, iso,year)
# df[, c(paste0("wind",0:15)) := shift(wind, 0:15, type = "lag"), by = iso] 
# # generate formulas for by-period model
# vp <- paste("i(period,",paste0("wind",0:15),")",sep="")
# fmla <- as.formula(paste0("growth ~ ",paste(vp,collapse=" + ",sep=' ')," | iso + year")) 
# # run once to get variable names - not strictly necessary but can get weird samples where vars get dropped
# mod <- feols(fmla,data=df)
# nms = names(coef(mod))
# #generate formulas for by-year model
# fmlayr <- as.formula(paste0("growth ~ (",paste0("wind",0:15,collapse=" + ",sep=' '),")*year | iso + year")) 
# nmsyr <- paste0("wind",0:15,":year")
# iso = unique(df$iso)
# coefperiod <- coefyr <- c()
# for (i in 1:Nboot) {
#   cl  <- data.frame(iso=sample(iso, size = length(iso), replace=T)) %>% mutate(ISOnew=row_number())
#   dtboot <- inner_join(df,cl,by="iso")
#   mod <- feols(fmla,data=dtboot)
#   coefperiod <- rbind(coefperiod,coef(mod)[nms])
#   mod <- feols(fmlayr,data=dtboot)
#   coefyr <- rbind(coefyr,coef(mod)[nmsyr])
#   print(i)
# }
# write_rds(coefperiod,file="data/cyclone/TCtimeperiod_bootstraps.rds")
# write_rds(coefyr,file='data/cyclone/TCbootstraps_yearmodel.rds')


# write out percent changes from annual model. can start here if above run
# coefyr <- read_rds('data/cyclone/TCbootstraps_yearmodel.rds')
# fmla <- as.formula(paste0("growth ~ ",paste0("wind",0:15,collapse=" + ",sep=' ')," | iso + year")) 
# mod <- feols(fmla,data=df)
# base = sum(coef(mod)[paste0("wind",0:15)])
# eff <- apply(coefyr,1,sum)
# ychg = eff/abs(base)*100 #percent change per year
# out <- data.frame(total=median_and_pval(ychg))
# write_rds(out,file='data/combined_plots/AnnualChanges_CyclonesGDP.rds')

# now make figure. can start here if above run
coef <- read_rds("data/cyclone/TCtimeperiod_bootstraps.rds")


colz = rev(met.brewer("Homer1",3))
pdf(file="fig/main/figure_3_panel_e.pdf",width=4,height=4)
par(mar=c(4,4,1,2))
plot(0,type="n",xlab="years since storm",ylab="cumulative growth",axes=F,ylim=c(-0.003,0.001),xlim=c(0,20))
abline(h=0,lty=2)
axis(1,at=seq(0,15,5),seq(0,15,5))
axis(2,las=1)
cumulatt=c(16,18,20) #where to plot cumulative effect on right side of plot
n=16 #number of lags + contemporaneous
out <- c()
for (i in 1:3) {
  eff <- coef[,paste0("period::",i,":wind",0:15)]
  eff <- t(apply(eff,1,cumsum)) 
  out <- cbind(out,eff[,paste0("period::",i,":wind15")])
  toplot <- apply(eff,2,function(x) quantile(x,probs = c(0.025,0.5,0.975)))
  if (i==1) {polygon(c(1:n,n:1)-1,c(toplot[1,],rev(toplot[3,])),col=alpha(colz[i],0.2),border = NA)}
  lines((1:n)-1,toplot[2,],lwd=2,col=colz[i])
  segments(cumulatt[i],toplot[1,16],cumulatt[i],toplot[3,16],col=colz[i])
  points(cumulatt[i],toplot[2,16],pch=19,col=colz[i])
}
dev.off()




# generate p-value heatmap
pdf(file="fig/supplementary/figure_s12_d.pdf",width=3,height=3)
mkheatmap(out, good = 1, labels=T,cex=0.8,plotratio=T,labs=c("'88-97","'98-07","'08-19"),main="TC - GDP growth") #need to have loaded this function from scripts/functions/mkheatmap.R
dev.off()


# #### old plot
# pdf(file="plots/cyclone/TCImpact.pdf",width=10,height=5)
# par(mfrow=c(1,2),mar=c(4,4,4,4))
# 
# 
# mod <- feols(growth~ l(wind,0:15) | iso + year,dt, panel.id = ~iso + year)
# cf <- coef(mod)
# n=length(cf)
# vcov <- vcov(mod)
# cumf <- cumsum(cf)
# sem <- c()
# for (i in 1:n) {
#   sem <- c(sem,sqrt(sum(vcov[1:i,1:i])))
# }
# cilo <- cumf-1.96*sem
# cihi <- cumf+1.96*sem
# plot(1:n,cumf,type="l",xlab="years since storm",ylab="cumulative growth",ylim=range(c(cilo,cihi)),axes=F)
# abline(h=0,lty=2)
# axis(1,at=seq(0,15,5)+1,seq(0,15,5))
# axis(2,las=1)
# polygon(c(1:n,n:1),c(cilo,rev(cihi)),col=alpha("lightblue",0.5),border = NA)
# lines(1:n,cumf)
# 
# 
# # by period
# mod <- feols(growth~ l(wind,0:15)*as.factor(period) | iso + year, dt,panel.id = ~iso + year)
# cf <- coef(mod)
# vcov <- vcov(mod)
# p2 <- grep("period)2",names(cf))
# p3 <- grep("period)3",names(cf))
# p1 <- which(1:length(cf)%in%c(p2,p3)==F)
# c1=sum(cf[p1]); se1=sqrt(sum(vcov[p1,p1]))
# c2=sum(cf[c(p1,p2)]); se2=sqrt(sum(vcov[c(p1,p2),c(p1,p2)]))
# c3=sum(cf[c(p1,p3)]); se3=sqrt(sum(vcov[c(p1,p3),c(p1,p3)]))
# cs=c(c1,c2,c3)
# cilo=cs-1.96*c(se1,se2,se3)
# cihi=cs+1.96*c(se1,se2,se3)
# 
# 
# plot(1:3,cs,ylim=c(-0.003,0),axes=F,xlab="",ylab="dGrowth/dWind",pch=19)
# abline(h=0,lty=2)
# segments(1:3,cilo,1:3,cihi)
# axis(2,las=1)
# text(1:3,cilo-0.0002,c("1965-83","1984-01","2002-19"))
# 
# 
# dev.off()
# 
# 
# 
# 
# ##. v2 plot with confidence intervals by period
# 
# 
# mod <- feols(growth~ l(wind,0:15)*as.factor(period) | iso + year, dt,panel.id = ~iso + year)
# cf <- coef(mod)
# vcov <- vcov(mod)
# p2 <- grep("period)2",names(cf))
# p3 <- grep("period)3",names(cf))
# p1 <- which(1:length(cf)%in%c(p2,p3)==F)
# c1=cumsum(cf[p1])
# c2=cumsum(cf[p1]+cf[p2])
# c3=cumsum(cf[p1]+cf[p3])
# 
# 
# n=length(p1)
# s1 <- s2 <- s3 <- c()
# for (i in 1:n) {
#   s1 <- c(s1,sqrt(sum(vcov[p1[1:i],p1[1:i]])))
#   s2 <- c(s2,sqrt(sum(vcov[c(p1[1:i],p2[1:i]),c(p1[1:i],p2[1:i])])))
#   s3 <- c(s3,sqrt(sum(vcov[c(p1[1:i],p3[1:i]),c(p1[1:i],p3[1:i])])))
# }
# cf <- list(c1,c2,c3); se <- list(s1,s2,s3)
# 
# 
# colz = rev(met.brewer("Homer1",3))
# pdf(file="plots/cyclone/TCImpact_v2.pdf",width=4,height=4)
# par(mar=c(4,4,1,2))
# plot(0,type="n",xlab="years since storm",ylab="cumulative growth",axes=F,ylim=c(-0.003,0.001),xlim=c(0,20))
# abline(h=0,lty=2)
# axis(1,at=seq(0,15,5),seq(0,15,5))
# axis(2,las=1)
# cumulatt=c(16,18,20) #where to plot cumulative effect on right side of plot
# for (i in 1:3) {
#   cilo=cf[[i]] -1.96*se[[i]]
#   cihi=cf[[i]] + 1.96*se[[i]]
#   if (i==1) {polygon(c(1:n,n:1)-1,c(cilo,rev(cihi)),col=alpha(colz[i],0.2),border = NA)}
#   lines(1:n-1,cf[[i]],lwd=2,col=colz[i])
#   segments(cumulatt[i],cilo[16],cumulatt[i],cihi[16],col=colz[i])
#   points(cumulatt[i],cf[[i]][16],pch=19,col=colz[i])
# }
# dev.off()
