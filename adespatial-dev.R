##
### Changes to adespatial
##
## rm(list=ls())
##
### Move the data file:
if(FALSE) {
  load("Tiahura.rda")
  save(Tiahura,file="Tiahura.rda",version=2L)
  file.copy("Tiahura.rda","../adespatial/data",overwrite=TRUE)
  file.copy("Tiahura.R","../adespatial/R",overwrite=TRUE)
}
##
library(adespatial)
data("Tiahura", package = "adespatial")
##
tiah.jac <- dist.ldc(t(Tiahura$fish),method = "jaccard")
##
tiah.chclust <- constr.hclust(tiah.jac, coords=Tiahura$habitat[,"distance"],
                              chron=TRUE)
##
par(mfrow=c(3,1))
par(mar=c(3,6.5,2,2))
dst <- Tiahura$habitat[,"distance"]
plot(NA, xlim=range(dst), ylim=c(0.5,5.5), yaxt="n",
     ylab="Partitions\n\n", xlab="")
parts <- c(2,3,5,7,12)
cols <- c("turquoise", "orange", "chartreuse", "aquamarine", "blue",
          "violet", "pink", "cyan", "green", "red", "cornsilk", "purple")
for(i in 1L:length(parts)) {
   tiah.chclust$coords[,"y"] <- i
   plot(tiah.chclust, parts[i], link=TRUE, lwd=3, hybrids="none",
        lwd.pt=0.5, cex=3, pch=21, plot=FALSE,
        col=cols[round(seq(1,length(cols), length.out=parts[i]))])
   
}
axis(2, at=1:length(parts), labels=paste(parts,"groups"), las=1)
par(mar=c(4,6.5,1,2))
plot(x=dst, y=Tiahura$habitat[,"depth"],
     ylim=c(max(range(Tiahura$habitat[,"depth"])),-300),
     las=1, ylab="Depth\n(cm)\n", xlab="", type="l", lwd=2)
for(i in 1:nrow(Tiahura$reef)) {
   abline(v=Tiahura$reef[i,2], lty=3)
   abline(v=Tiahura$reef[i,3], lty=3)
   if((Tiahura$reef[i,3] - Tiahura$reef[i,2])<100) {
      text(x=(Tiahura$reef[i,2] + Tiahura$reef[i,3])/2, y=2350,
           labels=toupper(Tiahura$reef[i,1]),srt=90,adj=0)
   } else {
      text(x=(Tiahura$reef[i,2] + Tiahura$reef[i,3])/2, y=-150,
           labels=toupper(Tiahura$reef[i,1]))
   }
}
par(mar=c(5,6.5,0,2))
plot(NA,xlim=range(dst), ylim=c(0,1), las=1,
     ylab="Bottom composition\n(proportions)\n", xlab="Distance (m)")
bot <- cbind(0, Tiahura$habitat[,3:10])
for(i in 2:9) bot[,i] <- bot[,i] + bot[,i-1]
cols <- c("", "grey75", "brown", "grey25", "green", "purple",
          "lightgreen", "yellow", "white")
for(i in 2:9)
   polygon(x=c(dst, rev(dst)),y=c(bot[,i], rev(bot[,i-1]))/50,
           col=cols[i])
text(x=c(44, 365, 707, 538, 957, 111, 965),
     y=c(0.05, 0.47, 0.37, 0.58, 0.42, 0.80, 0.88),
     labels=colnames(bot)[2:8], xpd=TRUE)
##
