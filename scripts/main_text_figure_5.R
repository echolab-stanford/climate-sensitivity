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


# v2 of summary plot - just showing total sensitivity in the main plot
#tedious list of outcomes: file name, name to print, whether more is good or bad, which type of sensitivity to estimate (for SI), period over which effects calculated
out <- c(
  "USAg_CORN","US maize",1,"ag","1950-2019","+1C growing season",
  "USAg_SOYBEANS","US soybeans",1,"ag","1950-2019","+1C growing season",
  "USAg_WHEAT","US wheat",1,"ag","1950-2019","+1C growing season",
  "EUAg_wheat","EU wheat",1,"ag","1990-2019","+1C growing season",
  "EUAg_maize", "EU maize",1,"ag","1990-2019","+1C growing season",
  "BrazilAg_soybean", "Brazil soy",1,"ag","1970-2019","+1C growing season",
  "BrazilAg_maize", "Brazil maize",1,"ag","1970-2019","+1C growing season",
  "IndiaAg_wheat", "India wheat",1,"ag","1990-2019","+1C growing season",
  "GlobalAgTFP","Global Ag TFP",1,"lin","1960-2019","+1C growing season",
  'USMortality',"US mortality - temperature",-1,"nonag","1968-2019","+1C monthly",
  'MEXMortality', "Mexico mortality - temperature", -1, "nonag", "1990-2019", "+1C monthly",
  'COMortality', 'Colombia mortality - temperature', -1, 'nonag', '1981-2019', '+1C monthly',
  "EUMortality", "EU mortality - temperature (annual)",-1,"nonag","1990-2019","+1C annual",
  "EUMortality_weeklydata", "EU mortality - temperature (weekly)",-1,"nonag","2000-2019","+1C weekly",
  "CyclonesMortality", "US mortality - cyclones",-1,"lin","1952-2015","+1 m/s wind speed",
  "GlobalGDP", "Global GDP - temperature",1,"nonag","1961-2019","+1C annual",
  "USIncome", "US income - temperature",1,"nonag","1968-2019","+1C annual",
  "CyclonesGDP", "Global GDP - cyclones",1,"lin","1965-2019","+1 m/s wind speed",
  "USFlooding", "US damages - floods",-1,"lin","1988-2017","+1sd monthly rainfall",
  "AfricanConflict", "African conflict",-1,"lin","1989-2019","+1C annual",
  "USViolentCrime", "US violent crime",-1,"nonag","1980-2019","+1C monthly",
  "USunintentional_injuries", "US injury mortality",-1,"nonag","1968-2019","+1C monthly",
  "USsuicide", "US suicide",-1,"nonag","1968-2019","+1C monthly"
)


toplot <- as.data.frame(matrix(out,nrow=23,byrow = T))
names(toplot) <- c("nms","nms1","posneg","outtype","years","exposure")
toplot$posneg <- as.numeric(toplot$posneg)
nn=dim(toplot)[1]




data <- list()
for (fl in 1:nn) {
  dir <- paste0('data/combined_plots/AnnualChanges_',toplot$nms[fl],".rds")
  data[[fl]] <- read_rds(dir)
}


# colors for rectangle shading
colpos=c(alpha("blue",c(0.6,0.4,0.2)),"grey90")  #colors when change is "good"
colneg=c(alpha("red",c(0.6,0.4,0.2)),"grey90") #colors when change is "bad"
pvalbins=c(1,0.1,0.05,0.01,0) #bins of pvalues to use




# now make plot, looping over outcomes
pdf(file='fig/main/figure_5.pdf',height=5,width=4)
par(mar=c(1,8,2,1))
plot(1,type="n",xlim=c(1,8),ylim=c(0,nn),xlab="",ylab="",axes=F)
for (i in 1:nn) {
  y = data[[i]][1,"total"]
  pval = data[[i]][2,"total"]
  valbin <- cut(pval,pvalbins,include.lowest = TRUE,labels = F,right=F) #categorize the pval
  if (sign(y)==toplot$posneg[i]) {
    col=colpos[valbin]
    plotval = -abs(y)
  } else {
    col=colneg[valbin]
    plotval = abs(y)
  }
  rect(7,nn-i,8,nn-i+1,col=col,border=NA)
  text(7.5,nn-i+0.5,sprintf("%.1f", round(plotval,1)),cex=0.7) #the `sprintf()` call makes in print same number of digits
}


text(1,nn:1-0.5,toplot$years,cex=0.7,pos=4) # add period
text(3,nn:1-0.5,toplot$exposure,cex=0.7,pos=4) # add exposure


# add y-axis labels
text(0,nn:1-0.5,toplot$nms1,adj=1,xpd=TRUE,cex=0.7)
dev.off()
