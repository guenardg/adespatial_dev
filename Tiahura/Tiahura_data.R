##
### Tiahura data
##
## rm(list=ls())
##
library(magrittr)
library(spdep)
##
Tiahura <-
    list(
        fish = read.table("Fish.txt"),
        species = scan("Species.txt",what=character(),sep="\n",quiet=TRUE),
        trait = read.table("Behav.txt"),
        habitat = read.table("Habitat.txt")
        
    )
##
Tiahura$fish %<>% t
rownames(Tiahura$fish) <-
    sprintf("Site_%d",1L:nrow(Tiahura$fish))
##
## head(Tiahura$trait)
colnames(Tiahura$trait) <-
    c("Feeding","Behavior","Adult","Egg","Activity")
Tiahura$trait[["Feeding"]] <-
    c("herbivorous","omnivorous","diurnal_grazer","carnivorous_type_1",
      "carnivorous_type_2","piscivorous",
      "zooplanktivorous")[Tiahura$trait[["Feeding"]]] %>% as.factor
Tiahura$trait[["Behavior"]] <-
    c("hiding","bottom","circling","above","swimmer","sub-surface",
      "pelagic")[Tiahura$trait[["Behavior"]]] %>% as.factor
Tiahura$trait[["Adult"]] <-
    c("0:0-15cm","1:16-30cm","2:31-60cm","3:61-120cm","4:121-240cm",
      "5:>240 cm")[Tiahura$trait[["Adult"]]] %>% as.ordered
Tiahura$trait[["Egg"]] <-
    c("pelagic","benthic","viviparous")[Tiahura$trait[["Egg"]]] %>% as.factor
Tiahura$trait[["Activity"]] <-
    c("diurnal","nocturnal","any")[Tiahura$trait[["Activity"]]] %>% as.factor
##
## head(Tiahura$habitat)
Tiahura$habitat %<>% t
rownames(Tiahura$habitat) <- rownames(Tiahura$fish)
colnames(Tiahura$habitat) <-
    c("distance","depth",
      "slab","sand","debris","turf","coral","algae","calcareous","other")
## rowSums(Tiahura$habitat[,3L:10L])
##
Tiahura$reef <-
    data.frame(
        zone=c("fringing reef","channel","barrier reef","ridge",
               "upper platform","outer\nslope"),
        from=c( 25, 250,  350,  815, 835,  900),
        to=c(  250, 350,  815,  835, 900, 1020)
    )
##
save(Tiahura, file="Tiahura.rda", version=2L)
system("cp -av ./Tiahura.rda ../../adespatial/data/Tiahura.rda")
system("cp -av ../Tiahura.R ../../adespatial/R/Tiahura.R")
##
### Package examples:
##
### Plotting the different sections of the transect:
plot(x=Tiahura$habitat[,"distance"], y=Tiahura$habitat[,"depth"],
     ylim=rev(range(Tiahura$habitat[,"depth"])), las=1,
     ylab="Depth\n(cm)\n", xlab="", type="l", lwd=2)
##
for(i in 1L:nrow(Tiahura$reef)) {
  abline(v=Tiahura$reef[i,2L],lty=3L)
  abline(v=Tiahura$reef[i,3L],lty=3L)
  if((Tiahura$reef[i,3L]-Tiahura$reef[i,2L])<100) {
    text(x=(Tiahura$reef[i,2L]+Tiahura$reef[i,3L])/2,y=2250,
         labels=toupper(Tiahura$reef[i,1L]),srt=90,adj=0,cex=0.75)
    } else {
      text(x=(Tiahura$reef[i,2L]+Tiahura$reef[i,3L])/2,y=2250,
           labels=toupper(Tiahura$reef[i,1L]),cex=0.75)
    }
}
abline(h=0)
##
## Plotting species presence in the sites:
plot_slice <- function(sl,split) {
  size <- ceiling(length(Tiahura$species)/split)
  sp_slice <- size*(sl - 1L) + (1L:size)
  image(z=t(as.matrix(Tiahura$fish[,sp_slice])),y=1:nrow(Tiahura$fish),
        x=1:length(sp_slice),zlim=c(0,1),col=c("white","black"),axes=FALSE,
        ylab="",xlab="")
  axis(1L,at=1:length(sp_slice),labels=Tiahura$species[sp_slice],las=2L)
  axis(2L,at=1:nrow(Tiahura$fish),label=rownames(Tiahura$fish),las=1L)
  invisible(NULL)
}
par(mar=c(15,5,2,2))
plot_slice(1L,5L)
## plot_slice(2L,5L)
## plot_slice(3L,5L)
## plot_slice(4L,5L)
## plot_slice(5L,5L)
##
