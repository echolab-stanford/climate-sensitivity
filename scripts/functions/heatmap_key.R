
# make a key for the heatmap plots
colpos=c(alpha("blue",c(0.8,0.5,0.3)),"grey80")  #colors when difference is positive
colneg=c(alpha("red",c(0.8,0.5,0.3)),"grey80") #colors when difference is negative

pdf(file="plots/combined/PvalKey.pdf",width=3,height=3)
plot(1,xlim=c(1,5),ylim=c(1,3),type="n",xlab="",ylab="",axes=F)
for (i in 1:4) {
  rect(i,1,i+1,2,col=colpos[i],border=NA)
  rect(i,2,i+1,3,col=colneg[i],border = NA)
}
dev.off()
