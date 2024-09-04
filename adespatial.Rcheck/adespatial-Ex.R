pkgname <- "adespatial"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "adespatial-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('adespatial')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Cperiodogram")
### * Cperiodogram

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Cperiodogram
### Title: Contingency periodogram
### Aliases: Cperiodogram

### ** Examples

# Data from the numerical example of Subsection 12.4.2 of Legendre and Legendre (2012).
test.vec <- c(1,1,2,3,3,2,1,2,3,2,1,1,2,3,3,1)
# Periodogram with tests using the chi-square distribution
res <- Cperiodogram(test.vec)
# Periodogram with permutation tests
res <- Cperiodogram(test.vec, nperm=2000, graph=FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Cperiodogram", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LCBD.comp")
### * LCBD.comp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LCBD.comp
### Title: Compute LCBD from any D matrix
### Aliases: LCBD.comp
### Keywords: spatial

### ** Examples


### Example 1
### Compute the Hellinger distance, then the LCBD indices.
if(require("vegan", quietly = TRUE)){
data(mite)
mite.hel = decostand(mite, "hellinger")
mite.D = dist(mite.hel)
out.mite.D = LCBD.comp(mite.D, sqrt.D=FALSE)
}

### Example 2
if(require("ade4", quietly = TRUE) & require("adegraphics", quietly = TRUE)){
data(doubs)
fish.sp = doubs$fish[-8,]   # Fish data; site 8 is removed because no fish were caught

out.comp = beta.div.comp(fish.sp, coef="S", quant=TRUE)

out.fish.D = LCBD.comp(out.comp$D, sqrt.D=TRUE)   # out.comp.D is not Euclidean
out.fish.D$beta
out.fish.Repl = LCBD.comp(out.comp$repl, sqrt.D=TRUE)   # out.comp$repl is not Euclidean
out.fish.Repl$beta
out.fish.AbDiff = LCBD.comp(out.comp$rich, sqrt.D=FALSE)   # out.comp$rich is Euclidean
out.fish.AbDiff$beta

### Plot maps of the LCBD indices
fish.xy = doubs$xy[-8,]   # Geographic coordinates; site 8 removed because no fish were caught

# Map of LCBD indices for %difference dissimilarity
s.value(fish.xy, out.fish.D$LCBD, method="size", symbol = "circle",
col = c("white", "brown"), main = "Doubs fish LCBD, %difference D")

# Map of LCBD indices for replacement component of %difference dissimilarity
s.value(fish.xy, out.fish.Repl$LCBD, method="size", symbol = "circle",
col = c("white", "brown"), main = "Doubs fish replacement LCBD")

# Map of LCBD indices for abundance difference component of %difference dissimilarity
s.value(fish.xy, out.fish.AbDiff$LCBD, method="size", symbol = "circle", 
col = c("white", "brown"), main = "Doubs fish abundance diff. LCBD")
}

## No test: 
if(require("ade4", quietly = TRUE) & require("betapart", quietly = TRUE)){
### Example 3
### This example requires packages \code{"betapart"} and \code{"ade4"} for data. 
### For the Baselga-family indices, the same partitioning results are obtained using
### (1) beta.div.comp or (2) beta.pair.abund() of \code{"betapart"} and LCBD.comp()

data(doubs)   # Data available in \code{"ade4"}
fish.sp = doubs$fish[-8,]   
# Fish data; site 8 is removed because no fish were caught
# We use abundance data in this example, not presence-absence data

# Partition into Baselga-family replacement and nestedness components 
# using \code{"beta.div.comp"} with the percentage difference index (aka Bray-Curtis)
out.comp = beta.div.comp(fish.sp, coef="BS", quant=TRUE)
out.comp$part

# Compute the D and component matrices using \code{"beta.pair.abund"}
out3 = beta.pair.abund(fish.sp, index.family = "bray")
summary(out3)

is.euclid(out3$beta.bray)    # D matrix out3$beta.bray is not Euclidean
out3.D = LCBD.comp(out3$beta.bray, sqrt.D=TRUE)
out3.D$beta
# Compare BDtotal here to BDtotal in out.comp$part (above)

out3.Repl = LCBD.comp(out3$beta.bray.bal, sqrt.D=TRUE)
out3.Repl$beta
# Compare BDtotal here to RichDiff in out.comp$part (above)

out3.AbDiff = LCBD.comp(out3$beta.bray.gra, sqrt.D=TRUE)
out3.AbDiff$beta
# Compare BDtotal here to RichDiff/Nes in out.comp$part (above)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LCBD.comp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ScotchWhiskey")
### * ScotchWhiskey

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ScotchWhiskey
### Title: Scotch Whiskey Data Set
### Aliases: ScotchWhiskey
### Keywords: Scotch Whiskey

### ** Examples

data(ScotchWhiskey)
lapply(ScotchWhiskey,ncol)
ScotchWhiskey$nbChar
ScotchWhiskey$listW  ## attr(ScotchWhiskey$listW,"class")
names(ScotchWhiskey)
names(ScotchWhiskey$dist)

plotWhiskey <- function(main) {
    plot(x=ScotchWhiskey$geo@coords[,1L]/1000,
         xlab="Eastings (km)",
         y=ScotchWhiskey$geo@coords[,2L]/1000,
         ylab="Northings (km)",
         main=main,
         type="n",asp=1)
    apply(
        ScotchWhiskey$neighbor@data,1L,
        function(X,coords) {
            segments(
                coords[X[1L],1L]/1000,
                coords[X[1L],2L]/1000,
                coords[X[2L],1L]/1000,
                coords[X[2L],2L]/1000
            )
        },
        coords=ScotchWhiskey$geo@coords
    )
    invisible(NULL)
}

plotWhiskey("Scotch whiskey: peat nose")
cols <- c("blue","orange")
points(ScotchWhiskey$geo@coords/1000,pch=21L,
       bg=cols[ScotchWhiskey$nose[,"peat"]+1L])
legend(x=50,y=1000,legend=c("Has a peat nose","Has no peat nose"),
       pch=21L,pt.bg=rev(cols))

plotWhiskey("Scotch whiskey: soft body")
cols <- c("red","green")
points(ScotchWhiskey$geo@coords/1000,pch=21L,
       bg=cols[ScotchWhiskey$body[,"soft"]+1L])
legend(x=50,y=1000,legend=c("Has a soft body","Has no soft body"),
       pch=21L,pt.bg=rev(cols))

plotWhiskey("Scotch whiskey: spicy palate")
cols <- c("red","green")
points(ScotchWhiskey$geo@coords/1000,pch=21L,
       bg=cols[ScotchWhiskey$palate[,"spice"]+1L])
legend(x=50,y=1000,legend=c("Has a spicy palate","Has no spicy palate"),
       pch=21L,pt.bg=rev(cols))

plotWhiskey("Scotch whiskey: sweet finish")
cols <- c("red","green")
points(ScotchWhiskey$geo@coords/1000,pch=21L,
       bg=cols[ScotchWhiskey$finish[,"sweet"]+1L])
legend(x=50,y=1000,legend=c("Has a sweet finish","Has no sweet finish"),
       pch=21L,pt.bg=rev(cols))

## To visualize (part of) the distance matrices:
as.matrix(ScotchWhiskey$dist$nose)[1:5,1:5]
as.matrix(ScotchWhiskey$dist$body)[1:5,1:5]
as.matrix(ScotchWhiskey$dist$palate)[1:5,1:5]
as.matrix(ScotchWhiskey$dist$finish)[1:5,1:5]

## The data tables:
ScotchWhiskey$colour
head(ScotchWhiskey$nose)
head(ScotchWhiskey$body)
head(ScotchWhiskey$palate)
head(ScotchWhiskey$finish)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ScotchWhiskey", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("TBI")
### * TBI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: TBI
### Title: TBI: Difference between multivariate observations at T1 and T2
### Aliases: TBI

### ** Examples

if(require("vegan", quietly = TRUE)) {

## Example 1 -

## Invertebrate communities subjected to insecticide treatment.

## As an example in their paper on Principal Response Curves (PRC method), van den
## Brink & ter Braak (1999) used observations on the abundances of 178 invertebrate
## species (macroinvertebrates and zooplankton) subjected to treatments in 12 mesocosms by
## the insecticide chlorpyrifos. The mesocosms were sampled at 11 occasions. The data,
## available in the {vegan} package, are log-transformed species abundances, ytranformed =
## log(10*y+1).

## The data of survey #4 will be compared to those of survey #11 in this example.
## Survey #4 was carried out one week after the insecticide treatment, whereas the fauna
## of the mesocosms was considered by the authors to have fully recovered from the
## insecticide treatment at survey #11.

data(pyrifos)

## The mesocosms had originally been attributed at random to the treatments. However,
## to facilitate presentation of the results, they will be listed here in order of
## increased insecticide doses: {0, 0, 0, 0, 0.1, 0.1, 0.9, 0.9, 6, 6, 44, 44} micro g/L.

## Select the 12 data rows of surveys 4 and 11 from the data file and reorder them

ord4 = c(38,39,41,47,37,44,40,46,43,48,42,45)

ord11 = c(122,123,125,131,121,128,124,130,127,132,126,129)

## Run the TBI function

res1 <- TBI(pyrifos[ord4,], pyrifos[ord11,], method = "%diff", nperm = 0, test.t.perm = FALSE)

res1$BCD.mat

## Example 2 -

## This example uses the mite data available in vegan. Let us pretend that sites 1-20
## represent T1 and sites 21-40 represent T2.


data(mite)

# Run the TBI function

res2 <- TBI(mite[1:20,], mite[21:40,], method = "%diff", nperm = 0, test.t.perm = FALSE)

summary(res2)

res2$BCD.mat

}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("TBI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Tiahura")
### * Tiahura

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Tiahura
### Title: Tiahura Transect Fish Data Set
### Aliases: Tiahura
### Keywords: Tiahura

### ** Examples

data(Tiahura)

## Compute dissimilary matrix from Jaccard's similarity coefficient:
tiah.jac <- dist.ldc(Tiahura$fish,method = "jaccard")

## Constrained clustering of the fish species:
tiah.chclust <- constr.hclust(tiah.jac, coords=Tiahura$habitat[,"distance"],
                              chron=TRUE)

## Plotting the results
par(mfrow=c(3,1))

## First graph: constrained clusters
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

## Second graph: transect profile
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

## Third graph: bottom composition
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

## Species presence graph set:
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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Tiahura", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("WRperiodogram")
### * WRperiodogram

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: WRperiodogram
### Title: Whittaker-Robinson periodogram
### Aliases: WRperiodogram plot.WRperio
### Keywords: multivariate

### ** Examples

 
### 1. Numerical example from Subsection 12.4.1 of Legendre and Legendre (2012)

test.vec <- c(2,2,4,7,10,5,2,5,8,4,1,2,5,9,6,3)

# Periodogram with permutation tests of significance
res <- WRperiodogram(test.vec)
plot(res) # Plot the periodogram

#####

### 2. Simulated data

# Generate a data series with periodic component using Legand's (1958) equation.
# Ref. Legendre and Legendre (2012, eq. 12.8, p. 753)
# x = time points, T = generated period, c = shift of curve, left (+) or right (-)

periodic.component <- function(x,T,c){cos((2*pi/T)*(x+c))}

n <- 500   # corresponds to 125 days with 4 observations per day
# Generate a lunar cycle, 29.5 days (T=118)
moon <- periodic.component(1:n, 118, 59)
# Generate a circadian cycle (T=4)
daily <- periodic.component(1:n, 4, 0)
# Generate an approximate tidal cycle (T=2)
# A real tidal signal would have a period of 12.42 h
tide <- periodic.component(1:n, 2, 0)

# Periodogram of the lunar component only 
res.moon.250 <- WRperiodogram(moon, nperm=0)  # T1=2, T2=n/2=250; no test
res.moon.130 <- WRperiodogram(moon, T2=130, nperm=499)
oldpar <- par(mfrow=c(1,2))
# Plot 2 moon cycles, n = 118*2 = 236 points
plot(moon[1:236], xlab="One time unit = 6 hours") 
plot(res.moon.130, prog=1) # Plot the periodogram

#####

# Add the daily and tidal components, plus a random normal error. With daily (T=4) and 
# tide (T=2), period 4 and its harmonics should have a higher W statistic than period 2
var1 <- daily + tide + rnorm(n, 0, 0.5)
# Plot a portion (40 points) of the data series
# Two periodic components identifiable. Tide (T=2) reinforces the daily signal (T=4)
par(mfrow=c(1,2))
plot(var1[1:40], pch=".", cex=1, xlab="One time unit = 6 hours")
lines(var1[1:40])
# Periodogram of 'var'
res.var1 <- WRperiodogram(var1, T2=40, nperm=499)
plot(res.var1, prog=3, line.col="blue") # Plot the periodogram
# The progressive correction for multiple tests (prog=3) was used in the periodogram.

#####

# Add the three components, plus a random normal error term
# to show that the WRperiodogram can test several periodic components at the same time.
# (5*moon) makes the lunar periods stronger than the daily and tidal periods
var2 <- 5*moon + daily + tide + rnorm(n, 0, 0.5)
# Plot a portion (150 points) of the data series
# The three periodic components are identifiable
par(mfrow=c(1,2))
plot(var2[1:150], pch=".", cex=1, xlab="One time unit = 6 hours")
lines(var2[1:150])

# Periodogram of 'var'
res.var2 <- WRperiodogram(var2, T2=130, nperm=499)
plot(res.var2, prog=1, line.col="blue") # Plot the periodogram
# Find the position of the maximum W statistic value in this periodogram
(which(res.var2[,2] == max(res.var2[,2])) -1)
# "-1" correction at the end of the previous line: the first computed period is T=2, 
# so period #118 is on line #117 of file res.var2

#####

# Illustration that the WR periodogram can handle missing values:
# Replace 10% of the 500 data by NA
select <- sort(sample(1:500)[1:50])
var.na <- var2
var.na[select] <- NA
res.var.na <- WRperiodogram(var.na, T2=130, nperm=499)
# Plot the periodogram with no correction for multiple tests
plot(res.var.na, prog=1)
# Plot periodogram again with progressive correction for multiple tests
plot(res.var.na, prog=3) 

#####

### 3. Data used in the examples of the documentation file of function afc() of {stats}
# Data file "ldeaths"; time series, 6 years x 12 months of deaths in UK hospitals
# First, examine the data file ldeaths. Then:
ld.res.perio <- WRperiodogram(ldeaths, nperm=499)
# Plot the periodogram with two types of corrections for multiple tests
par(mfrow=c(1,2))
plot(ld.res.perio, prog=1) # No correction for multiple testing
plot(ld.res.perio, prog=3) # Progressive correction for multiple tests
# The yearly cycle and harmonics are significant
# Compare the results of afc() to those of WRperiodogram above
acf(ldeaths)   # lag=1.0 is one year; see ?acf
par(oldpar)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("WRperiodogram", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("aem")
### * aem

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aem
### Title: Construct asymmetric eigenvector maps (AEM)
### Aliases: aem
### Keywords: spatial

### ** Examples

# Construction of object of class nb (spdep)
if(require("spdep", quietly = TRUE)){
nb <- cell2nb(5,5,"queen")

# Create fictitious geographical coordinates 
xy <- cbind(1:25,expand.grid(1:5,1:5))

# Build binary site-by-link matrix
bin.mat <- aem.build.binary(nb,xy)

# Construct AEM eigenfunctions from an object of class aem.build.binary
res <- aem(aem.build.binary=bin.mat,rm.link0=FALSE)
res$values

# Illustrate 4 AEM eigenfunctions using bubble plots
opal <- palette()
palette(c("black","white"))
oldpar <- par(mfrow=c(2,2))
symbols(x=xy[,2:3], circles=abs(res$vectors[,1]), inches=FALSE, asp=1,
 fg=ifelse(sign(-res$vectors[,1])+1>0,1,0), 
 bg=ifelse(sign(res$vectors[,1])+1>0,1,0), xlab="x", ylab="y")
title("AEM 1")
symbols(x=xy[,2:3], circles=abs(res$vectors[,2]), inches=FALSE, 
asp=1, fg=ifelse(sign(-res$vectors[,2])+1>0,1,0),
 bg=ifelse(sign(res$vectors[,2])+1>0,1,0), xlab="x", ylab="y")
title("AEM 2")
symbols(x=xy[,2:3], circles=abs(res$vectors[,3]), inches=FALSE, 
asp=1, fg=ifelse(sign(-res$vectors[,3])+1>0,1,0), 
bg=ifelse(sign(res$vectors[,3])+1>0,1,0), xlab="x", ylab="y")
title("AEM 3")
symbols(x=xy[,2:3], circles=abs(res$vectors[,4]), inches=FALSE, asp=1,
 fg=ifelse(sign(-res$vectors[,4])+1>0,1,0), 
 bg=ifelse(sign(res$vectors[,4])+1>0,1,0), xlab="x", ylab="y")
title("AEM 4")

# Construct AEM eigenfunctions using only a site-by-link matrix
res2 <- aem(binary.mat=bin.mat[[1]])
res2$values

# Illustrate 4 AEM eigenfunctions using bubble plots
par(mfrow=c(2,2))
symbols(x=xy[,2:3], circles=abs(res2$vectors[,1]), inches=FALSE, 
asp=1, fg=ifelse(sign(-res2$vectors[,1])+1>0,1,0), 
bg=ifelse(sign(res2$vectors[,1])+1>0,1,0), xlab="x", ylab="y")
title("AEM 1")
symbols(x=xy[,2:3], circles=abs(res2$vectors[,2]), inches=FALSE,
asp=1, fg=ifelse(sign(-res2$vectors[,2])+1>0,1,0), 
bg=ifelse(sign(res2$vectors[,2])+1>0,1,0), xlab="x", ylab="y")
title("AEM 2")
symbols(x=xy[,2:3], circles=abs(res2$vectors[,3]), inches=FALSE,
asp=1, fg=ifelse(sign(-res2$vectors[,3])+1>0,1,0), 
bg=ifelse(sign(res2$vectors[,3])+1>0,1,0), xlab="x", ylab="y")
title("AEM 3")
symbols(x=xy[,2:3], circles=abs(res2$vectors[,4]), inches=FALSE,asp=1,
 fg=ifelse(sign(-res2$vectors[,4])+1>0,1,0), 
 bg=ifelse(sign(res2$vectors[,4])+1>0,1,0), xlab="x", ylab="y")
title("AEM 4")

palette(opal)
par(oldpar)

# Construct AEM eigenfunctions with a function of the distance
# as weights to put on the links

# Construction of object of class nb (spdep)
nb<-cell2nb(5,5,"queen")

# Create fictitious geographical coordinates
xy <- cbind(1:25,expand.grid(1:5,1:5))

# Build binary site-by-link matrix
bin.mat <- aem.build.binary(nb,xy)

# Construct a matrix of distances
long.lien.mat<-as.matrix(dist(xy))

# Extract the edges, remove the ones directly linked to site 0
lien.b<-bin.mat$edges[-1:-5,]

# Construct a vector giving the length of each edge
long.lien<-vector(length=nrow(lien.b))

for(i in 1:nrow(lien.b)){
	long.lien[i]<-long.lien.mat[lien.b[i,1],lien.b[i,2]]
}

# Construct a vector of weights based on distance
weight.vec<-1-(long.lien/max(long.lien))^2

# Construct AEM eigenfunctions from an object of class aem.build.binary
res <- aem(aem.build.binary=bin.mat,weight=weight.vec,rm.link0=TRUE)
res

# Computing Moran's I for AEMs

# Building AEMs
xy <- cbind(1:25,expand.grid(1:5,1:5))
Wdist <- 1/as.matrix(dist(xy[,2:3]))

nb <- cell2nb(5,5,"queen")
bin.mat <- aem.build.binary(nb,xy)
linkBase <- bin.mat[[2]]
link <- linkBase[-which(linkBase[,1] == 0),]
weight <- numeric()

for(i in 1:nrow(link)){
   weight[i] <- Wdist[link[i,1],link[i,2]]
}

AEM <- aem(bin.mat, weight = weight, rm.link0 = TRUE)

# Constructing asymmetric matrix
matasym <- matrix(0,ncol=25, nrow=25)

for(i in 1:nrow(link)){
    matasym[link[i,1],link[i,2]]<- weight[i]
}

# Build a listw object from the asymmetric matrix
listwAsym <-  mat2listw(matasym, style = "B", zero.policy = TRUE)

# Calculate Moran's I for AEM
MoranIAEM <- moran.randtest(AEM$vectors, listwAsym)

}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aem", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("aem.build.binary")
### * aem.build.binary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aem.build.binary
### Title: Construct a site-by-edge binary matrix
### Aliases: aem.build.binary
### Keywords: spatial

### ** Examples


### Create an object of class nb (spdep)
if(require("spdep", quietly = TRUE)){
nb<-cell2nb(5,5,"queen")

### Create fictitious geographical coordinates 
xy <- cbind(1:25,expand.grid(1:5,1:5))

### Build a binary site-by-link matrix; remove the site which have identical Y coordinate
### (by default argument: rm.same.y = TRUE)
bin.mat <- aem.build.binary(nb,xy)
str(bin.mat)

### Build a binary site-by-link matrix using the argument link: remove the site which
### have identical Y coordinate (by default argument: rm.same.y = TRUE)
edges<-expand.grid(1,2:25)
bin.mat <- aem.build.binary(coords=xy,link=edges)
str(bin.mat)

### Build a binary site-by-link matrix, making the process affect the points at 
### an angle of 45 degrees
bin.mat.45 <- aem.build.binary(nb,xy, rot.angle=45)
str(bin.mat.45)

### Build a binary site-by-link matrix, making the process affect the points at
### an angle of pi/3 radians
bin.mat.pi3 <- aem.build.binary(nb,xy,unit.angle="radians", rot.angle=pi/3)
str(bin.mat.pi3)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aem.build.binary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aem.time")
### * aem.time

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aem.time
### Title: AEM for time series
### Aliases: aem.time
### Keywords: multivariate spatial

### ** Examples


# Time series containing 20 equispaced observations
out <- aem.time(20, moran = TRUE)

# Time series containing 20 observations with unequal spacing
# Generate (n-1) random interpoint distances
distances <- runif(19,1,5)

# Compute weights representing the ease of communication among points
w <- 1/(distances/max(distances))

# Compute the AEM eigenfunctions
out <- aem.time(20, w = w, moran = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aem.time", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aem.weight.edges")
### * aem.weight.edges

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aem.weight.edges
### Title: Weight edges when constructing AEM variables
### Aliases: aem.weight.edges aem.weight.time
### Keywords: spatial ts

### ** Examples


### Time serie example
### Example - 12 dates (days from January 1st of year 1) 
### in a 6-year study starting September 5, 2000
if(require("spdep", quietly = TRUE)){
dates <- as.Date(c(129,269,500,631,864,976,1228,1352,1606,1730,1957,2087),origin="2000/1/1")
autocor.limit <- 522  # Limit of autcorrelation in the correlogram

### Using aem.weight.time()
(wtime <- aem.weight.time(dates, alpha=2, max.d=autocor.limit))
### Using aem.weight.edges()
n <- length(dates)
nb <- cell2nb(1, n)
xy.dates <- cbind(1:n, rep(1, n), dates)
(wtime <- aem.weight.edges(nb, xy.dates, alpha=2, max.d=autocor.limit))

n <- length(dates)
nb <- cell2nb(1, n)
xy.dates <- cbind(1:n, dates, rep(1, n)) ## Note the inversion of 'dates' and 'rep(1,n)'
wtime <- aem.weight.edges(nb, xy.dates, alpha=2, 
max.d=autocor.limit,rot.angle=90) # Note that 'rot.angle=90' was used

### Spatial example using default d.max (notice the warning)
###########################################################################
nb<-cell2nb(5,5,"queen")
xy <- cbind(1:25,expand.grid(1:5,1:5))
(wspace <- aem.weight.edges(nb,xy))
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aem.weight.edges", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("beta.div")
### * beta.div

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: beta.div
### Title: Beta diversity computed as Var(Y)
### Aliases: beta.div

### ** Examples


if(require("vegan", quietly = TRUE) & require("adegraphics", quietly = TRUE)){
data(mite)
res = beta.div(mite, "hellinger", nperm=999)

# Plot a map of the LCBD indices using the Cartesian coordinates
data(mite.xy)
s.value(mite.xy, res$LCBD, symbol = "circle", col = c("white", "brown"), main="Map of mite LCBD")

### Example using the mite abundance data and the percentage difference dissimilarity
res = beta.div(mite, "percentdiff", nperm=999, clock=TRUE)

# Plot a map of the LCBD indices
signif = which(res$p.LCBD <= 0.05)	# Which are the significant LCBD indices?
nonsignif = which(res$p.LCBD > 0.05)	# Which are the non-significant LCBD indices?
g1 <- s.value(mite.xy[signif,], res$LCBD[signif], ppoint.alpha = 0.5, plegend.drawKey = FALSE,
 symbol = "circle", col = c("white", "red"), main="Map of mite LCBD (red = significant indices)")
g2 <- s.value(mite.xy[nonsignif,], res$LCBD[nonsignif], ppoint.alpha = 0.5,
symbol = "circle", col = c("white", "blue"))
g2+g1
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("beta.div", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("beta.div.comp")
### * beta.div.comp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: beta.div.comp
### Title: Decompose D in replacement and richness difference components
### Aliases: beta.div.comp

### ** Examples


if(require(ade4, quietly = TRUE)){
data(doubs)
fish.sp = doubs$fish[-8,]   # Fish data; site 8 is removed because no fish were caught

# Compute and partition a matrix of Jaccard indices (presence-absence data)
out1 = beta.div.comp(fish.sp, coef="J", quant=FALSE)
out1$part

# Compute and partition a matrix of percentage difference indices
# (quantitative form of Sorensen index)
out2 = beta.div.comp(fish.sp, coef="S", quant=TRUE)
out2$part
# In paragraph Value, see the description of the 5 elements of vector part. 
# Is the fish beta diversity dominated by replacement or richness/abundance difference?
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("beta.div.comp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("chooseCN")
### * chooseCN

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: chooseCN
### Title: Function to choose a connection network
### Aliases: chooseCN
### Keywords: spatial utilities

### ** Examples


if(require("ade4", quietly = TRUE)){
data(mafragh)

oldpar <- par(mfrow=c(2,2))
cn1 <- chooseCN(mafragh$xy,ask=FALSE,type=1)
cn2 <- chooseCN(mafragh$xy,ask=FALSE,type=2)
cn3 <- chooseCN(mafragh$xy,ask=FALSE,type=3)
cn4 <- chooseCN(mafragh$xy,ask=FALSE,type=4)
par(oldpar)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("chooseCN", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("constr.hclust")
### * constr.hclust

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: constr.hclust
### Title: Space- And Time-Constrained Clustering
### Aliases: constr.hclust

### ** Examples


## First example: Artificial map data from Legendre & Legendre
##                (2012, Fig. 13.26): n = 16

dat <- c(41,42,25,38,50,30,41,43,43,41,30,50,38,25,42,41)
coord.dat <- matrix(c(1,3,5,7,2,4,6,8,1,3,5,7,2,4,6,8,
                      4.4,4.4,4.4,4.4,3.3,3.3,3.3,3.3,
                      2.2,2.2,2.2,2.2,1.1,1.1,1.1,1.1),16,2)

## Obtaining a list of neighbours:
library(spdep)
listW <- nb2listw(tri2nb(coord.dat), style="B")
links.mat.dat <- listw2mat(listW)
neighbors <- listw2sn(listW)[,1:2]

## Calculating the (Euclidean) distance between points:
D.dat <- dist(dat)

## Display the points:
plot(coord.dat, type='n',asp=1)
title("Delaunay triangulation")
text(coord.dat, labels=as.character(as.matrix(dat)), pos=3)
for(i in 1:nrow(neighbors))
    lines(rbind(coord.dat[neighbors[i,1],],
          coord.dat[neighbors[i,2],]))

## Unconstrained clustring by hclust:
grpWD2_hclust <- hclust(D.dat, method="ward.D2")
plot(grpWD2_hclust, hang=-1)

## Clustering without a contiguity constraint;
## the result is represented as a dendrogram:
grpWD2_constr_hclust <- constr.hclust(D.dat, method="ward.D2")
plot(grpWD2_constr_hclust, hang=-1)

## Clustering with a contiguity constraint described by a list of
## links:
grpWD2cst_constr_hclust <-
    constr.hclust(
        D.dat, method="ward.D2",
        neighbors, coord.dat)

## To visualize using hclust's plotting method:
## stats:::plot.hclust(grpWD2cst_constr_hclust, hang=-1)

## Plot the results on a map with k=3 clusters:
plot(grpWD2cst_constr_hclust, k=3, links=TRUE, las=1, xlab="Eastings",
     ylab="Northings", cex=3, lwd=3)

## Generic functions from hclust can be used, for instance to obtain
## a list of members of each cluster:
cutree(grpWD2cst_constr_hclust, k=3)

## Now with k=5 clusters:
plot(grpWD2cst_constr_hclust, k=5, links=TRUE, las=1, xlab="Eastings",
     ylab="Northings", cex=3, lwd=3)
cutree(grpWD2cst_constr_hclust, k=5)

## End of the artificial map example


## Second example: Scotch Whiskey distilleries clustered using tasting
## scores (nose, body, palate, finish, and the four distances combined)
## constrained with respect to the distillery locations in Scotland.

## Documentation file about the Scotch Whiskey data: ?ScotchWhiskey

data(ScotchWhiskey)

## Cluster analyses for the nose, body, palate, and finish D
## matrices:

grpWD2cst_ScotchWhiskey <-
    lapply(
        ScotchWhiskey$dist,    ## A list of distance matrices
        constr.hclust,         ## The function called by function lapply
        links=ScotchWhiskey$neighbors@data,         ## The list of links
        coords=ScotchWhiskey$geo@coords/1000
    )

## The four D matrices (nose, body, palate, finish), represented as
## vectors in the ScotchWiskey data file, are combined as follows to
## produce a single distance matrix integrating all four types of
## tastes:

Dmat <- ScotchWhiskey$dist
ScotchWhiskey[["norm"]] <-
    sqrt(Dmat$nose^2 + Dmat$body^2 + Dmat$palate^2 + Dmat$finish^2)

## This example shows how to apply const.clust to a single D matrix when
## the data file contains several matrices.

grpWD2cst_ScotchWhiskey[["norm"]] <-
    constr.hclust(
        d=ScotchWhiskey[["norm"]],method="ward.D2",
        ScotchWhiskey$neighbors@data,
        coords=ScotchWhiskey$geo@coords/1000
    )

## A fonction to plot the Whiskey clustering results:

plotWhiskey <- function(wh, k) {
   par(fig=c(0,1,0,1))
   plot(grpWD2cst_ScotchWhiskey[[wh]], k=k, links=TRUE, las=1,
        xlab="Eastings (km)", ylab="Northings (km)", cex=0.1, lwd=3,
        main=sprintf("Feature: %s",wh))
   text(ScotchWhiskey$geo@coords/1000,labels=1:length(ScotchWhiskey$geo))
   legend(x=375, y=700, lty=1L, lwd=3, col=rainbow(1.2*k)[1L:k],
          legend=sprintf("Group %d",1:k), cex=1.25)
   SpeyZoom <- list(xlim=c(314.7,342.2), ylim=c(834.3,860.0))
   rect(xleft=SpeyZoom$xlim[1L], ybottom=SpeyZoom$ylim[1L],col="#E6E6E680",
        xright=SpeyZoom$xlim[2L], ytop=SpeyZoom$ylim[2L], lwd=2, lty=1L)
   par(fig=c(0.01,0.50,0.46,0.99), new=TRUE)
   plot(grpWD2cst_ScotchWhiskey[[wh]], xlim=SpeyZoom$xlim,
        ylim=SpeyZoom$ylim, k=k, links=TRUE, las=1, xlab="", ylab="",
        cex=0.1, lwd=3, axes=FALSE)
   text(ScotchWhiskey$geo@coords/1000,labels=1:length(ScotchWhiskey$geo))
   rect(xleft=SpeyZoom$xlim[1L], ybottom=SpeyZoom$ylim[1L],
        xright=SpeyZoom$xlim[2L], ytop=SpeyZoom$ylim[2L], lwd=2, lty=1L)
}

## Plot the clustering results on the map of Scotland for 5 groups.
## The inset map shows the Speyside distilleries in detail:
plotWhiskey("nose", 5L)
plotWhiskey("body", 5L)
plotWhiskey("palate", 5L)
plotWhiskey("finish", 5L)
plotWhiskey("norm", 5L)

## End of the Scotch Whiskey tasting data example

## Not run: 
##D 
##D ## Third example: Fish community composition along the Doubs River,
##D ## France. The sequence is analyzed as a case of chronological
##D ## clustering, substituting space for time.
##D 
##D if(require("ade4", quietly = TRUE)){
##D data(doubs, package="ade4")
##D Doubs.D <- dist.ldc(doubs$fish, method="hellinger")
##D grpWD2cst_fish <- constr.hclust(Doubs.D, method="ward.D2", chron=TRUE,
##D                                 coords=as.matrix(doubs$xy))
##D plot(grpWD2cst_fish, k=5, las=1, xlab="Eastings (km)",
##D      ylab="Northings (km)", cex=3, lwd=3)
##D 
##D ## Repeat the plot with other values of k (number of groups)
##D 
##D ## End of the Doubs River fish assemblages example
##D 
##D ## Example with 6 connected points, shown in Fig. 2 of Guénard & Legendre paper 
##D 
##D var = c(1.5, 0.2, 5.1, 3.0, 2.1, 1.4)
##D ex.Y = data.frame(var)
##D 
##D ## Site coordinates, matrix xy
##D x.coo = c(-1, -2, -0.5, 0.5, 2, 1)
##D y.coo = c(-2, -1, 0, 0, 1, 2)
##D ex.xy = data.frame(x.coo, y.coo)
##D 
##D ## Matrix of connecting edges E
##D from = c(1,1,2,3,4,3,4)
##D to = c(2,3,3,4,5,6,6)
##D ex.E = data.frame(from, to)
##D 
##D ## Carry out constrained clustering analysis
##D test.out <-
##D     constr.hclust(
##D         dist(ex.Y),       # Response dissimilarity matrix
##D         method="ward.D2", # Clustering method
##D         links=ex.E,       # File of link edges (constraint) E
##D         coords=ex.xy      # File of geographic coordinates
##D     )
##D 
##D par(mfrow=c(1,2))
##D ## Plot the map of the results for k = 3
##D plot(test.out, k=3)
##D ## Plot the dendrogram
##D stats:::plot.hclust(test.out, hang=-1)
##D }
##D 
##D ## Same example modified: disjoint clusters
##D ## Same ex.Y and ex.xy as in the previous example
##D var = c(1.5, 0.2, 5.1, 3.0, 2.1, 1.4)
##D ex.Y = data.frame(var)
##D 
##D ## Site coordinates, matrix xy
##D x.coo = c(-1, -2, -0.5, 0.5, 2, 1)
##D y.coo = c(-2, -1, 0, 0, 1, 2)
##D ex.xy = data.frame(x.coo, y.coo)
##D 
##D ## Matrix of connecting edges E2
##D from = c(1,1,2,4,4)
##D to = c(2,3,3,5,6)
##D ex.E2 = data.frame(from, to)
##D 
##D ## Carry out constrained clustering analysis
##D test.out2 <-
##D     constr.hclust(
##D         dist(ex.Y),       # Response dissimilarity matrix
##D         method="ward.D2", # Clustering method
##D         links=ex.E2,      # File of link edges (constraint) E
##D         coords=ex.xy      # File of geographic coordinates
##D     )
##D cutree(test.out2, k=2)
##D 
##D par(mfrow=c(1,2))
##D ## Plot the map of the results for k = 3
##D plot(test.out2, k=3)
##D ## Plot the dendrogram showing the disconnected groups
##D stats:::plot.hclust(test.out2, hang=-1)
##D axis(2,at=0:ceiling(max(test.out2$height,na.rm=TRUE)))
##D 
##D ## End of the disjoint clusters example
##D 
## End(Not run)
## End of examples




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("constr.hclust", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("create.dbMEM.model")
### * create.dbMEM.model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: create.dbMEM.model
### Title: Combine dbMEM matrices corresponding to groups of sites
### Aliases: create.dbMEM.model
### Keywords: spatial

### ** Examples

{
 # Generate random coordinates for 35 sites forming 6 distinct groups on the map
 Easting <- runif(35)+c(rep(0,6),rep(1.5,7),rep(3,6), rep(0,5),rep(1.5,5),rep(3,6))
 Northing<- runif(35)+c(rep(2.8,6),rep(2.3,7),rep(2.8,6), rep(0,5),rep(0.5,5),rep(0,6))
 cartesian <- cbind(Easting,Northing)
 rownames(cartesian) <- paste("S",1:nrow(cartesian),sep='')
 nsites.per.group <- c(6,7,6,5,5,6)

 result <- create.dbMEM.model(coord=cartesian, nsites=nsites.per.group)

 # Draw a map to check the coding of the sites into the groups
 site.codes <- unlist(apply(cbind(1:6),1,n=nsites.per.group,function(a,n) rep(a,n[a])))

 col.vec <- c("green3","gray99","orange2","gold1","brown3","gray70")
 plot(cartesian, pch=22, col="black", bg=col.vec[site.codes], cex=2, ylim=c(0,4),asp=1)
 text(cartesian,labels=rownames(cartesian), cex=0.5, pos=3)

 # Examine the staggered matrix of dbMEM eigenfunctions
 # Not run:
 result
}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("create.dbMEM.model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dbmem")
### * dbmem

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dbmem
### Title: dbMEM spatial eigenfunctions
### Aliases: dbmem

### ** Examples

if(require("ade4", quietly = TRUE) & require("adegraphics", quietly = TRUE)){

data(oribatid)
mite <- oribatid$fau      # 70 peat cores, 35 species
mite.xy <- oribatid$xy    # Geographic coordinates of the 70 cores

# Example 1: Compute the MEMs corresponding to all non-null eigenvalues
# thresh=1.012 is the value used in Borcard and Legendre (2002)
mite.dbmem1 <- dbmem(mite.xy, thresh=1.012, MEM.autocor = "non-null", silent = FALSE)
mite.dbmem1

# Print the (n-1) non-null eigenvalues
attributes(mite.dbmem1)$values
# or:  attr(mite.dbmem1, "values")

# Plot the associated spatial weighting matrix
s.label(mite.xy, nb = attr(mite.dbmem1, "listw"))

# Plot maps of the first 3 dbMEM eigenfunctions
s.value(mite.xy, mite.dbmem1[,1:3])

# Compute and test associated Moran's I values
# Eigenvalues are proportional to Moran's I

test <- moran.randtest(mite.dbmem1, nrepet = 99)
plot(test$obs, attr(mite.dbmem1, "values"), xlab = "Moran's I", ylab = "Eigenvalues")

# Decreasing values of Moran's I for the successive MEM.
# The red line is the expected value of Moran's I under H0.

plot(test$obs, xlab="MEM rank", ylab="Moran's I")
abline(h=-1/(nrow(mite.xy) - 1), col="red")

# Example 2: Compute only the MEMs with positive eigenvalues (and positive Moran's I)
mite.dbmem2 <- dbmem(mite.xy, thresh=1.012)
# or:  mite.dbmem2 <- dbmem(dist(mite.xy), thresh=1.012, silent=FALSE)
mite.dbmem2

# Examine the eigenvalues
attributes(mite.dbmem2)$values
# or:  attr(mite.dbmem2, "values")

# Examine (any portion of) the dbmem spatial eigenvectors
tmp <- as.matrix(mite.dbmem2)
tmp[1:10,1:6]
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dbmem", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("directional.response")
### * directional.response

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: directional.response
### Title: Directional indices of community change
### Aliases: directional.response

### ** Examples


# Artificial Example
art <- c(1,1,1,0,0,0,
         0,0,0,1,1,0,
         0,0,0,0,0,1)
art.data <- matrix(art, nrow=3, ncol=6, byrow=TRUE)

art.out <- directional.response(art.data, method="overlap",relativize=NULL)

# Real data example: the Doubs River fish data (Verneaux 1973), available in ade4.
# 30 sites, 27 species. No fish had been caught at site 8; remove that site
if(require("ade4", quietly = TRUE)) {

data(doubs)
dim(doubs$fish)   
fish <- doubs$fish[-8,] 
dim(fish)
doubs.out <- directional.response(fish, method="gain", relativize="S")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("directional.response", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dist.ldc")
### * dist.ldc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dist.ldc
### Title: Dissimilarity matrices for community composition data
### Aliases: dist.ldc

### ** Examples


if(require("vegan", quietly = TRUE)) {
data(mite)
mat1  = as.matrix(mite[1:10, 1:15])   # No column has a sum of 0
mat2 = as.matrix(mite[61:70, 1:15])   # 7 of the 15 columns have a sum of 0

#Example 1: compute Hellinger distance for mat1
D.out = dist.ldc(mat1,"hellinger")

#Example 2: compute chi-square distance for mat2
D.out = dist.ldc(mat2,"chisquare")

#Example 3: compute percentage difference dissimilarity for mat2
D.out = dist.ldc(mat2,"percentdiff")

}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dist.ldc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("envspace.test")
### * envspace.test

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: envspace.test
### Title: Perform a test of the shared space-environment fraction of a
###   variation partitioning using torus-translation (TT) or Moran Spectral
###   Randomisation (MSR)
### Aliases: envspace.test
### Keywords: spatial

### ** Examples

## No test: 
if(require(vegan)) { 
# Illustration of the test of the SSEF on the oribatid mite data
# (Borcard et al. 1992, 1994 for details on the dataset):
# Community data (response matrix):
data(mite)
# Hellinger-transformation of the community data (Legendre and Gallagher 2001):
Y <- decostand(mite, method = "hellinger")
# Environmental explanatory dataset:
data(mite.env)
# We only use two numerical explanatory variables:
env <- mite.env[, 1:2]
dim(Y)
dim(env)
# Coordinates of the 70 sites:
data(mite.xy)
coord <- mite.xy

### Building a list of candidate spatial weighting matrices (SWMs) for the 
### optimisation of the SWM selection, separately for 'Y' and 'env':
# We create five candidate SWMs: a connectivity matrix based on a Gabriel graphs, on
# a minimum spanning tree (i.e., two contrasted graph-based SWMs), either
# not weighted, or weighted by a linear function decreasing with the distance),
# and a distance-based SWM corresponding to the connectivity and weighting
# criteria of the original PCNM method:
candidates <- listw.candidates(coord, nb = c("gab", "mst", "pcnm"), weights = c("binary",
                                                                                "flin"))
### Optimisation of the selection of a SWM:
# SWM for 'Y' (based on the best forward-selected subset of MEM variables):
modsel.Y <- listw.select(Y, candidates, method = "FWD", MEM.autocor = "positive",
                         p.adjust = TRUE)
                         
names(candidates)[modsel.Y$best.id]                 # Best SWM selected
modsel.Y$candidates$Pvalue[modsel.Y$best.id]        # Adjusted p-value of the global model
modsel.Y$candidates$N.var[modsel.Y$best.id]         # Nb of forward-selected MEM variables
modsel.Y$candidates$R2Adj.select[modsel.Y$best.id]  # Adjusted R2 of the selected MEM var.

# SWM for 'env' (method = "global" for the optimisation, as all MEM variables are required
# to use MSR):
modsel.env <- listw.select(env, candidates, method = "global", MEM.autocor = "positive",
                           p.adjust = TRUE)

names(candidates)[modsel.env$best.id]                  # Best SWM selected
modsel.env$candidates$Pvalue[modsel.env$best.id]       # Adjusted p-value of the global model
modsel.env$candidates$N.var[modsel.env$best.id]        # Nb of forward-selected MEM variables
modsel.env$candidates$R2Adj.select[modsel.env$best.id] # Adjusted R2 of the selected MEM var.

### We perform the variation partitioning:
# Subset of selected MEM variables within the best SWM:
MEM.spe <- modsel.Y$best$MEM.select

VP <- varpart(Y, env, MEM.spe)
plot(VP)

# Test of the shared space-environment fraction (fraction [b]):
SSEF.test <- envspace.test(Y, env, coord, MEM.spe, 
                           listw.env = candidates[[modsel.env$best.id]], 
                           regular = FALSE, nperm = 999)
SSEF.test

# The SSEF is highly significant, indicating a potential induced spatial dependence.
}
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("envspace.test", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("forward.sel")
### * forward.sel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: forward.sel
### Title: Forward selection with multivariate Y using permutation under
###   reducel model
### Aliases: forward.sel
### Keywords: multivariate

### ** Examples


x <- matrix(rnorm(30),10,3)
y <- matrix(rnorm(50),10,5)
    
forward.sel(y,x,nperm=99, alpha = 0.5)
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("forward.sel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("forward.sel.par")
### * forward.sel.par

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: forward.sel.par
### Title: Parametric forward selection of explanatory variables in
###   regression and RDA
### Aliases: forward.sel.par
### Keywords: multivariate

### ** Examples


x <- matrix(rnorm(30),10,3)
y <- matrix(rnorm(50),10,5)
    
forward.sel.par(y,x, alpha = 0.5)
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("forward.sel.par", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("give.thresh")
### * give.thresh

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: give.thresh
### Title: Compute the maximum distance of the minimum spanning tree based
###   on a distance matrix
### Aliases: give.thresh

### ** Examples

xy <- matrix(rnorm(60),30,2)
dxy <- dist(xy)
th <- give.thresh(dxy)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("give.thresh", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("global.rtest")
### * global.rtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: global.rtest
### Title: Global and local tests
### Aliases: global.rtest local.rtest
### Keywords: multivariate spatial

### ** Examples



# wait for a generic dataset





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("global.rtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("listw.candidates")
### * listw.candidates

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: listw.candidates
### Title: Function to create a list of spatial weighting matrices
### Aliases: listw.candidates
### Keywords: spatial

### ** Examples

### Create 100 random sampling locations in a squared grid of 120 x 120:
xy <- matrix(nrow = 100, ncol = 2)
xy[, 1] <- sample(c(1:120), 100, replace = FALSE)
xy[, 2] <- sample(c(1:120), 100, replace = FALSE)
### The function listw.candidates is used to build the spatial weighting matrices that
### we want to test and compare (with the listw.select function). We test a Gabriel's graph, 
### a minimum spanning tree, and a distance-based connectivity defined by a threshold
### distance corresponding to the smallest distance keeping all sites connected (i.e., 
### the defaut value of d2). These connectivity matrices are then either not weighted 
### (binary weighting), or weighted by the linearly decreasing function:
candidates <- listw.candidates(coord = xy, nb = c("gab", "mst", "dnear"), 
                               weights = c("binary", "flin"))
names(candidates)                              
plot(candidates[[1]], xy)
plot(candidates[[3]], xy)
### Construction of a different list of spatial weighting matrices. This time, the
### connexions are defined by a distance-based criterion based on the same threshold
### value, but the connections are weighted by the concave-down function with a y parameter
### varying between 2 and 5, and a concave-up function with a y parametre of 0.2.
candidates2 <- listw.candidates(coord = xy, nb = "dnear", weights = c("fdown", "fup"),
                                y_fdown = 1:5, y_fup = 0.2)
### Number of spatial weighting matrices generated:
length(candidates2) 
### A single SWM can also easily be generated with listw.candidates:
lw <- listw.candidates(xy, nb = "gab", weights = "bin")
plot(lw[[1]], xy)

### Generating MEM variables from an object of listw.candidates with scores.listw:
MEM <- scores.listw(lw[[1]])
### See functions mem.select and listw.select for examples of how to use an object
### created by listw.candidates with these functions.





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("listw.candidates", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("listw.explore")
### * listw.explore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: listw.explore
### Title: Interactive tool to generate R code that creates a spatial
###   weighting matrix
### Aliases: listw.explore

### ** Examples

if(interactive()){
## a matrix or an object of class 'Spatial*' should be in the global environment
xy <- matrix(rnorm(50), 25)
listw.explore()
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("listw.explore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("listw.select")
### * listw.select

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: listw.select
### Title: Function to optimize the selection of a spatial weighting matrix
###   and select the best subset of eigenvectors (MEM, Moran's Eigenvector
###   Maps)
### Aliases: listw.select
### Keywords: spatial

### ** Examples

## No test: 
if(require(spdep)) {
### Create a grid of 15 x 15:
grid <- expand.grid(x = seq(1, 15, 1), y = seq(1, 15, 1))
### Generate a response variable Y structured at broad scale by linear combination of
### the first three MEM variables to which a normal noise is added:
nb <- cell2nb(nrow = 15, ncol = 15, "queen")
lw <- nb2listw(nb, style = "B")
MEM <- scores.listw(lw, MEM.autocor = "positive")
# Degree of spatial autocorrelation:
intensity <- 0.8
Y_space <- scale(MEM[, 1] + MEM[, 2] + MEM[, 3]) * intensity
Y_noise <- scale(rnorm(n = nrow(MEM), mean = 0, sd = 1)) * (1 - intensity)
Y <- Y_space + Y_noise
### Y is sampled in 100 randomly-chosen sites of the grid:
idx.sample <- sample(c(1:nrow(grid)), 100, replace = FALSE)
xy <- grid[idx.sample, ]
Y_sampled <- Y[idx.sample]
### The function listw.candidates is used to build the spatial weighting matrices that
### we want to test and compare (with the listw.select function). We test a Gabriel's graph,
### a minimum spanning tree, and a distance-based connectivity defined by a threshold
### distance corresponding to the smallest distance keeping all sites connected (i.e.,
### the defaut value of d2; see help of function listw.candidates).
### These connectivity matrices are then either not weighted (binary weighting), or
### weighted by the linearly decreasing function (see help of the function listw.candidates):
candidates <- listw.candidates(coord = xy, nb = c("gab", "mst"), weights = c("binary", "flin"))
### Number of candidate W matrices generated:
nbw <- length(candidates)
### Significance threshold value after p-value correction (Sidak correction):
1 - (1 - 0.05)^(1/nbw)
### Optimization of the selection of the SWM among the candidates generated above,
### using the corrected significance threshold calculated above for the global tests:
W_sel <- listw.select(Y_sampled, candidates, MEM.autocor = "positive", method = "FWD",
                    p.adjust = TRUE, nperm = 299)
### Some characteristics of the best spatial model:
# Best SWM:
W_sel$best.id
# Selected subset of spatial predictor within the best SWM:
W_sel$best$MEM.select
nrow(W_sel$best$summary)
# Corrected p-value of the global test of the best SWM:
W_sel$best$global.test$Pvalue
# Adjusted R2 of the subset of spatial predictors selected within the chosen SWM:
max(W_sel$best$summary$R2Adj)
# p-values of all the tested W matrices:
W_sel$candidates$Pvalue
# Adjusted R2 of the subset of spatial predictors selected for all the significant
# W matrices:
W_sel$candidates$R2Adj.select

# See Appendix S3 of Bauman et al. 2018 for more extensive examples and illustrations.
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("listw.select", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mem")
### * mem

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: scores.listw
### Title: Function to compute and manage Moran's Eigenvector Maps (MEM) of
###   a listw object
### Aliases: scores.listw mem orthobasis.listw [.orthobasisSp
### Keywords: spatial

### ** Examples


if(require("ade4", quietly = TRUE) & require("spdep", quietly = TRUE)){
data(oribatid)
nbtri <- tri2nb(as.matrix(oribatid$xy))
sc.tri <- scores.listw(nb2listw(nbtri, style = "B"))
summary(sc.tri)
}
if(require("adegraphics", quietly = TRUE)){
s.value(oribatid$xy,sc.tri[,1:9])
plot(sc.tri[,1:6], oribatid$xy, pSp.cex = 5, pSp.alpha = 0.5, pbackground.col = 'lightblue')
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mem", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mem.select")
### * mem.select

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mem.select
### Title: Selection of the best subset of spatial eigenvectors (MEM,
###   Moran's Eigenvector Maps)
### Aliases: mem.select
### Keywords: spatial

### ** Examples

if(require(vegan)){ 
# Illustration of the MIR selection on the oribatid mite data
# (Borcard et al. 1992, 1994 for details on the dataset):
# *******************************************************
# Community data (response matrix):
data(mite)
# We will compute the example on a single species:
spe <- mite[, 2]
# Environmental explanatory dataset:
data(mite.env)
# We only use two numerical explanatory variables:
env <- mite.env[, 1:2]
dim(env)
# Coordinates of the 70 sites:
data(mite.xy)
coord <- mite.xy
# We build the model we are interested in:
mod <- lm(spe ~ ., data = env)


# In order to avoid possible type I error rate inflation issues, we check 
# whether the model residuals are independent, and if they are spatially
# autocorrelated, we select a small subset of MEM variables to add to the
# model as covariables with the MIR selection:

# 1) We build a spatial weighting matrix based on Gabriel graph with a
# weighting function decreasing linearly with the distance:
w <- listw.candidates(coord, nb = "gab", weights = "flin")


# 2) We test the spatial autocorrelation of the model residuals and, if
# necessary, select a subset of spatial predictors:
y <- residuals(mod)
MEM <- mem.select(x = y, listw = w[[1]], method = "MIR", MEM.autocor = "positive",
         nperm = 999, alpha = 0.05)
dim(MEM$MEM.select)
# The residuals of the model presented spatial autocorrelation. The selection
# of MEM variables is thus performed to remove residual autocorrelation.

# 3) We can reconstruct our model adding the selected MEM variable as covariables:
env2 <- cbind(env, MEM$MEM.select)
mod_complete <- lm(spe ~ ., data = env2)
summary(mod_complete)$coefficient[, 1]   # Coefficient estimates
summary(mod_complete)$coefficient[, 2]   # Standard errors
}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mem.select", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mfpa")
### * mfpa

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mfpa
### Title: Multi-frequential periodogram analysis
### Aliases: mfpa plot.mfpa print.mfpa

### ** Examples


### Example 1

# Simulate data with frequencies 2.3 and 6.1 and a random component, n = 100. 
# No trend, no autocorrelated residuals.

y <- as.matrix(0.4*(sin(2.3*2*pi*(1:100)/100)) +
0.4*(sin(6.1*2*pi*(1:100)/100)) + 0.2*rnorm(100))

res <- mfpa(y, MaxNFreq = 2, MinFreq = 2, ntrend = 0, nlags = 0)

# Compute the periods associated with the two periodic components. Each
# frequency in element $frequencies is a number of cycles in the whole series.
# The periods are expressed in numbers of time intervals of the data series. In
# this example, if the data are measured every min, the periods are in min.

periods <- 100/res$frequencies$frequency 

# Draw the data series and the fitted (or predicted) values

plot(res)

### Example 2

# Generate hourly periodic data with tide signal (tide period T = 12.42 h)
# during 1 year, hence 24*365 = 8760 hourly data. See
# https://en.wikipedia.org/wiki/Tide.

# In this simulation, constant (c = 0) puts the maximum value of the cosine at
# midnight on the first day of the series.

periodic.component <- function(x, T, c) cos((2*pi/T)*(x+c))

tide.h <- periodic.component(1:8760, 12.42, 0)

# The number of tides in the series is: 8760/12.42 = 705.314 tidal cycles
# during one year.

# Sample the hourly data series once a day at 12:00 noon every day. The
# periodic signal to be detected has a period smaller then the interval between
# consecutive observations and its frequency is larger than (n-1). The sequence
# of sampling hours for the tide.h data is:

h.noon <- seq(12, 8760, 24) 
tide.data <- tide.h[h.noon]
length(tide.data)   

# The series contains 365 sampling units

# Compute Dutilleul's multi-frequential periodogram

res.noon <- mfpa(tide.data, MaxNFreq = 1, MinFreq = 2, ntrend = 1, nlags = 2)

# Examine the frequency detected by the periodogram, element
# res.noon$frequencies. This is a harmonic of the tide signal in the original
# data series tide.h.

# Compute the period of the signal in the data series sampled hourly:

period <- 365/res.noon$frequencies$frequency

# Draw the data series and the adjusted values

plot(res.noon)

# Repeat this analysis after addition of random noise to the tide data

tide.noise <- tide.data + rnorm(365, 0, 0.25)

res.noise <- mfpa(tide.noise, MaxNFreq = 1, MinFreq = 2, ntrend = 1, nlags = 2)

plot(res.noise)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mfpa", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("moran.bounds")
### * moran.bounds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: moran.bounds
### Title: Function to compute extreme values of Moran's I
### Aliases: moran.bounds
### Keywords: spatial

### ** Examples


if(require("ade4", quietly = TRUE)){
 if(require("spdep", quietly = TRUE)){
     data(oribatid)
     nbtri <- tri2nb(as.matrix(oribatid$xy))
     lwB <- nb2listw(nbtri, style = "B")
     lwW <- nb2listw(nbtri, style = "W")
     scB <- mem(lwB)
     scW <- mem(lwW)
     moran.bounds(lwB)
     moran.mc(scB[,1], lwB, 9)
     moran.mc(scB[,69], lwB, 9)
     moran.bounds(lwW)
     moran.mc(scW[,1], lwW, 9)
     moran.mc(scW[,69], lwW, 9)
 }
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("moran.bounds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("moran.randtest")
### * moran.randtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: moran.randtest
### Title: Function to compute Moran's index of spatial autocorrelation
### Aliases: moran.randtest
### Keywords: spatial

### ** Examples


if(require("ade4", quietly = TRUE)  & require("spdep", quiet = TRUE)){
data(mafragh)
tests <- moran.randtest(mafragh$env, nb2listw(mafragh$nb))
tests
plot(tests)

}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("moran.randtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("moranNP.randtest")
### * moranNP.randtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: moranNP.randtest
### Title: Function to compute positive and negative parts of Moran's index
###   of spatial autocorrelation
### Aliases: moranNP.randtest
### Keywords: spatial

### ** Examples

if(require("ade4", quietly = TRUE)  & require("spdep", quiet = TRUE)){
data(mafragh)
tests <- moranNP.randtest(mafragh$env[,1], nb2listw(mafragh$nb),
 alter = "two-sided", p.adjust.method = "holm")
tests
moran.randtest(mafragh$env[,1], nb2listw(mafragh$nb))$obs
sum(tests$obs)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("moranNP.randtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mspa")
### * mspa

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mspa
### Title: Multi-Scale Pattern Analysis
### Aliases: mspa print.mspa scatter.mspa
### Keywords: multivariate spatial

### ** Examples



####################################
### using oribatib mites dataset ###
####################################

if(require("ade4", quietly = TRUE)){
## load data
data(oribatid)

## get the list of spatial weights
cn <- chooseCN(oribatid$xy, res = "listw", ask = FALSE, type = 1)

## Hellinger transformation
hellTrans <- function(X){
  if (!( is.matrix(X) | is.data.frame(X) )) stop("Object is not a matrix.")  
  if (any(is.na(X))) stop("na entries in table.")
  
  sumRow <- apply(X,1,sum)
  Y <- X/sumRow
  Y <- sqrt(Y)
  
  return(Y)
}


## ENVIRONMENTAL VARIABLES ##
## Hill and Smith analysis for environmental variables
## (for a mixture of quantitative / qualitative variables)
hsEnv <- dudi.hillsmith(oribatid$envir,scannf=FALSE)

## detrending of the analysis (residuals of regression onto xy coordinates)
hsEnv.detr <- pcaivortho(hsEnv,oribatid$xy,scannf=FALSE)

## MSPA of the detrended analysis
mspaEnv <- mspa(hsEnv.detr,cn,scannf=FALSE,nf=2)
scatter(mspaEnv)



## SPECIES DATA ##
## PCA of species abundances, after Hellinger transformation
pcaFau <- dudi.pca(hellTrans(oribatid$fau),scale=FALSE,scannf=FALSE)

## detrending of this PCA
pcaFau.detr <- pcaivortho(pcaFau,oribatid$xy,scannf=FALSE)

# MSPA of the detrended analysis
mspaFau <- mspa(pcaFau.detr,cn,scannf=FALSE,nf=2)
scatter(mspaFau)



## CANONICAL MSPA ##
## RDA species ~ envir
## (species abundances predicted by environment)
## note: RDA = 'PCAIV' (PCA with Instrumental Variables)
rda1 <- pcaiv(dudi=pcaFau.detr, df=oribatid$envir,scannf=FALSE,nf=2)

## canonical MSPA (species predicted by environment)
mspaCan1 <- mspa(dudi=rda1, lw=cn, scannf=FALSE, nf=2)
scatter(mspaCan1)

## same analysis, using a non-parametric centring
mspaCan1NP <- mspa(dudi=rda1, lw=cn, scannf=FALSE, nf=2,cent="sim",nper=999)
scatter(mspaCan1NP) # basically no change



## PARTIAL CANONICAL MSPA ##
## partial RDA species ~ envir
## (species abundances not predicted by environment)
rda2 <- pcaivortho(dudi=pcaFau.detr,df=oribatid$envir,scannf=FALSE,nf=2)

## partial canonical MSPA
mspaCan2 <- mspa(dudi=rda2, lw=cn, scannf=FALSE, nf=2)
scatter(mspaCan2) # nothing left
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mspa", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("msr.4thcorner")
### * msr.4thcorner

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: msr.4thcorner
### Title: Moran spectral randomization for fourth-corner analysis
### Aliases: msr.4thcorner
### Keywords: spatial

### ** Examples


if(require("ade4", quietly = TRUE) & require("adephylo", quietly = TRUE)
& require("spdep", quietly = TRUE) & require("ape", quietly = TRUE)){
data(mafragh, package = "ade4")
fr1 <- fourthcorner(mafragh$env, mafragh$flo, mafragh$traits$tabQuantitative, nrepet = 49)
phy <- read.tree(text = mafragh$tre)
lw <- nb2listw(mafragh$nb)
fr1.msr <- msr(fr1, listwORorthobasis = lw, phyloORorthobasis = phy)

fr1
fr1.msr
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("msr.4thcorner", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("msr")
### * msr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: msr
### Title: Moran spectral randomization
### Aliases: msr msr.default
### Keywords: spatial

### ** Examples


library(spdep)
x1 <- matrix(rnorm(81*5), nrow = 81)
lw1 <- nb2listw(cell2nb(9, 9))

moran.mc(x1[,1], lw1, 2)$statistic

## singleton
x1.1 <- msr(x1[,1], lw1, nrepet = 9, method = "singleton")
apply(x1.1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## triplet
x1.2 <- msr(x1[,1], lw1, nrepet = 9, method = "triplet")
apply(x1.2, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## pair
x1.3 <- msr(x1[,1], lw1, nrepet = 9, method = "pair")
apply(x1.3, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## pair with cor.fixed
x1.4 <- msr(x1[,1], lw1, nrepet = 9, cor.fixed = 0.5)
apply(x1.4, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
cor(x1[,1], x1.4)

## pair preserving correlations for multivariate data
x1.5 <- msr(x1, lw1, nrepet = 9, cor.fixed = 0.5)
cor(x1)
lapply(x1.5, cor)

apply(x1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
apply(x1.5[[1]], 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)

## singleton preserving correlations for multivariate data
x1.6 <- msr(x1, lw1, nrepet = 9, method = "singleton")
cor(x1)
lapply(x1.6, cor)

apply(x1, 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)
apply(x1.6[[1]], 2, function(x) moran.mc(x, listw = lw1, nsim = 2)$statistic)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("msr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("msr.mantelrtest")
### * msr.mantelrtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: msr.mantelrtest
### Title: Moran spectral randomization for Mantel test
### Aliases: msr.mantelrtest
### Keywords: spatial

### ** Examples


if(require("ade4", quietly = TRUE)
& require("spdep", quietly = TRUE)){
data(mafragh, package = "ade4")

d1 <- dist(mafragh$env[,1:3])
d2 <- dist(mafragh$env[,7])
t1 <- mantel.randtest(d1,d2)
t1

lw <- nb2listw(mafragh$nb)
t2 <- msr(t1, listwORorthobasis = lw)
t2

}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("msr.mantelrtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("msr.varipart")
### * msr.varipart

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: msr.varipart
### Title: Moran spectral randomization for variation partitioning
### Aliases: msr.varipart

### ** Examples

library(ade4)
library(spdep)
data(mafragh)
## Performing standard variation partitioning
dudiY <- dudi.pca(mafragh$flo, scannf = FALSE, scale = FALSE)
mafragh.lw <- nb2listw(mafragh$nb)
me <- mem(mafragh.lw, MEM.autocor = "positive")
vprda <- varipart(dudiY, mafragh$env, me, type = "parametric")

## Adjust estimation and compute p-value by msr methods
vprda.msr <- msr(vprda, mafragh.lw, nrepet=99)
vprda.msr



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("msr.varipart", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("mst.nb")
### * mst.nb

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: mst.nb
### Title: Function to compute neighborhood based on the minimum spanning
###   tree
### Aliases: mst.nb
### Keywords: spatial

### ** Examples


xy <- matrix(rnorm(60),30,2)
dxy <- dist(xy)
th <- give.thresh(dxy)
nb1 <- mst.nb(dxy)
nb1
wh1 <- which(as.matrix(dxy)==th,arr.ind=TRUE)
plot(nb1,xy,pch=20,cex=2,lty=3)
lines(xy[wh1[1,],1],xy[wh1[1,],2],lwd=2)
title(main="Maximum distance of the minimum spanning tree in bold")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("mst.nb", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("multispati")
### * multispati

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: multispati
### Title: Multivariate spatial analysis
### Aliases: multispati plot.multispati summary.multispati print.multispati
### Keywords: multivariate spatial

### ** Examples



if (require(spdep, quiet = TRUE) & require(ade4, quiet = TRUE)) {
    data(mafragh)
    maf.xy <- mafragh$xy
    maf.flo <- mafragh$flo
    maf.listw <- nb2listw(mafragh$nb)
    if(adegraphicsLoaded()) {
      g1 <- s.label(maf.xy, nb = mafragh$nb, plab.cex = 0.75)
    } else {
      s.label(maf.xy, neig = mafragh$neig, clab = 0.75)
    }
    maf.coa <- dudi.coa(maf.flo,scannf = FALSE)
    maf.coa.ms <- multispati(maf.coa, maf.listw, scannf = FALSE, nfposi = 2, nfnega = 2)
    maf.coa.ms
    
    ### detail eigenvalues components
    fgraph <- function(obj){
      # use multispati summary
      sum.obj <- summary(obj)
      # compute Imin and Imax
      Ibounds <- moran.bounds(eval(as.list(obj$call)$listw))
      Imin <- Ibounds[1]
      Imax <- Ibounds[2]
      I0 <- -1/(nrow(obj$li)-1)
      # create labels
      labels <- lapply(1:length(obj$eig),function(i) bquote(lambda[.(i)]))
      # draw the plot
      xmax <- eval(as.list(obj$call)$dudi)$eig[1]*1.1
      oldpar <- par(las=1)
      var <- sum.obj[,2]
      moran <- sum.obj[,3]
      plot(x=var,y=moran,type='n',xlab='Inertia',ylab="Spatial autocorrelation (I)",
           xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n')
      text(x=var,y=moran,do.call(expression,labels))
      ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
      ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
      ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),
           ytlab[2:4],as.character(round(Imax,1)))
      axis(side=2,at=ytick,labels=ytlab)
      rect(0,Imin,xmax,Imax,lty=2)
      segments(0,I0,xmax,I0,lty=2)
      abline(v=0)
      title("Spatial and inertia components of the eigenvalues")
      par(oldpar)
    }
    fgraph(maf.coa.ms)
    ## end eigenvalues details


    if(adegraphicsLoaded()) {
      g2 <- s1d.barchart(maf.coa$eig, p1d.hori = FALSE, plot = FALSE)
      g3 <- s1d.barchart(maf.coa.ms$eig, p1d.hori = FALSE, plot = FALSE) 
      g4 <- s.corcircle(maf.coa.ms$as, plot = FALSE)
      G1 <- ADEgS(list(g2, g3, g4), layout = c(1, 3))
    } else {
      oldpar <- par(mfrow = c(1, 3))
      barplot(maf.coa$eig)
      barplot(maf.coa.ms$eig) 
      s.corcircle(maf.coa.ms$as)
      par(oldpar)
    }
 
 
    if(adegraphicsLoaded()) {
      g5 <- s.value(maf.xy, -maf.coa$li[, 1], plot = FALSE)
      g6 <- s.value(maf.xy, -maf.coa$li[, 2], plot = FALSE)
      g7 <- s.value(maf.xy, maf.coa.ms$li[, 1], plot = FALSE)
      g8 <- s.value(maf.xy, maf.coa.ms$li[, 2], plot = FALSE)
      G2 <- ADEgS(list(g5, g6, g7, g8), layout = c(2, 2))
    } else {
      oldpar <- par(mfrow = c(2, 2))
      s.value(maf.xy, -maf.coa$li[, 1])
      s.value(maf.xy, -maf.coa$li[, 2])
      s.value(maf.xy, maf.coa.ms$li[, 1])
      s.value(maf.xy, maf.coa.ms$li[, 2])
      par(oldpar)
    }


    w1 <- -maf.coa$li[, 1:2]
    w1m <- apply(w1, 2, lag.listw, x = maf.listw)
    w1.ms <- maf.coa.ms$li[, 1:2]
    w1.msm <- apply(w1.ms, 2, lag.listw, x = maf.listw)
    if(adegraphicsLoaded()) {
      g9 <- s.match(w1, w1m, plab.cex = 0.75, plot = FALSE)
      g10 <- s.match(w1.ms, w1.msm, plab.cex = 0.75, plot = FALSE)
      G3 <- cbindADEg(g9, g10, plot = TRUE)
    } else {
      oldpar <- par(mfrow = c(1,2))
      s.match(w1, w1m, clab = 0.75)
      s.match(w1.ms, w1.msm, clab = 0.75)
      par(oldpar)
    }

    maf.pca <- dudi.pca(mafragh$env, scannf = FALSE)
    multispati.randtest(maf.pca, maf.listw)
    maf.pca.ms <- multispati(maf.pca, maf.listw, scannf=FALSE)
    plot(maf.pca.ms)
}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("multispati", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ortho.AIC")
### * ortho.AIC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ortho.AIC
### Title: Compute AIC for models with orthonormal explanatory variables
### Aliases: ortho.AIC
### Keywords: models

### ** Examples


y <- matrix(rnorm(50),50,1)
x <- svd(scale(y %*% c(0.1,0.5,2,0,0.7)+matrix(rnorm(250),50,5)))$u
res <- ortho.AIC(y,x,ord.var=TRUE)
minAIC <- which.min(res$AICc)
nvar <- length(1:minAIC)+1 # number of orthogonal vectors + 1 for intercept
lm1 <- lm(y~x[,res$ord[1:minAIC]])
summary(lm1)$r.squared # R2
res$R2[minAIC] # the same
min(res$AICc) # corrected AIC
extractAIC(lm1) # classical AIC
min(res$AICc)-2*(nvar*(nvar+1))/(nrow(x)-nvar-1) # the same

lm2 <- lm(y~1)

res$AICc0 # corrected AIC for the null model
extractAIC(lm2) # classical AIC
res$AICc0-2*(1*(1+1))/(nrow(x)-1-1) # the same




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ortho.AIC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("orthobasis.poly")
### * orthobasis.poly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: orthobasis.poly
### Title: Function to compute polynomial of geographical coordinates
### Aliases: orthobasis.poly
### Keywords: spatial

### ** Examples


if(require("ade4", quietly = TRUE)){
data(mafragh, package = "ade4")
pol2 <- orthobasis.poly(mafragh$Spatial) 

if(require("adegraphics", quietly = TRUE)){
plot(pol2, mafragh$Spatial)
}
}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("orthobasis.poly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.TBI")
### * plot.TBI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.TBI
### Title: Plots of the outputs of a temporal beta diversity analysis
### Aliases: plot.TBI

### ** Examples


if(require("vegan", quietly = TRUE)) {

## Example 1 -

## Invertebrate communities subjected to insecticide treatment.

## As an example in their paper on Principal Response Curves (PRC method), van den 
## Brink & ter Braak (1999) used observations on the abundances of 178 invertebrate 
## species (macroinvertebrates and zooplankton) subjected to treatments in 12 mesocosms 
## by the insecticide chlorpyrifos. The mesocosms were sampled at 11 occasions. The 
## data, available in the {vegan} package, are log-transformed species abundances, 
## ytranformed = loge(10*y+1).

## The data of survey #4 will be compared to those of survey #11 in this example. 
## Survey #4 was carried out one week after the insecticide treatment, whereas the  
## fauna of the mesocosms was considered by the authors to have fully recovered from  
## the insecticide treatment at survey #11.

data(pyrifos)

## The mesocosms had originally been attributed at random to the treatments. However, 
## to facilitate presentation of the results, they will be listed here in order of 
## increased insecticide doses: {0, 0, 0, 0, 0.1, 0.1, 0.9, 0.9, 6, 6, 44, 44} 
## micro g/L.

## Select the 12 data rows of surveys 4 and 11 from the data file and reorder them

ord4 <- c(38,39,41,47,37,44,40,46,43,48,42,45)
ord11 <- c(122,123,125,131,121,128,124,130,127,132,126,129)

## Run the TBI function

res1 <- TBI(pyrifos[ord4,], pyrifos[ord11,], method = "%diff", nperm = 0, test.t.perm = FALSE)

res1$BCD.mat

## Draw BC plots

oldpar <- par(mfrow=c(1,2))

s.names <- paste("Surv",1:12,sep=".")

## In the 1st plot, the symbols have diameters proportional to the site TBI statistics 

plot(res1, s.names=s.names, col.bg="red", pch.loss=21, pch.gain=22, 
main="B-C plot, Pyrifos, surveys 4 & 11")

## In the 2nd plot, control the axes limit values by specifying xlim and ylim

plot(res1, s.names=1:12, col.bg="green", pch.loss=23, pch.gain=24, 
main="B-C plot, Pyrifos, surveys 4 & 11", xlim=c(0,0.5), ylim=c(0.1,0.6))

## In the 3rd plot, draw all symbols small and of the same size, using cex.symb=NULL

par(oldpar)

plot(res1, s.names=1:12, col.bg="gold", pch.loss=23, pch.gain=24, 
main="B-C plot, Pyrifos, surveys 4 & 11", cex.symb=NULL)

## Example 2 -

## This example uses the mite data available in vegan. Let us pretend that sites 1-20 
## represent a survey at time 1 (T1) and sites 21-40 a survey at time 2 (T2).

data(mite)

## Run the TBI function

res2 <- TBI(mite[1:20,],mite[21:40,],method="%diff",nperm=0,test.t.perm=FALSE)

res2$BCD.mat

## Draw BC plots

oldpar <- par(mfrow=c(1,2))

s.names=rownames(res2$BCD.mat)

## In the 1st plot, the symbols have diameters proportional to the site TBI statistics

plot(res2, s.names=s.names, col.bg="cadetblue2", pch.loss=21, pch.gain=22, 
main="B-C plot, Mite data")

# In the 2nd plot, control the axes limit values by specifying xlim and ylim

plot(res2, s.names=1:20, col.rim="coral2", pch.loss=19, pch.gain=15, 
main="B-C plot, Mite data", xlim=c(0,0.6), ylim=c(0,0.6))
par(oldpar)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.TBI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot.constr.hclust")
### * plot.constr.hclust

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.constr.hclust
### Title: Plotting Method For Space- And Time-Constrained Clustering
### Aliases: plot.constr.hclust

### ** Examples


## Artificial map data from Legendre & Legendre (2012, Fig. 13.26)
## n = 16

dat <- c(41,42,25,38,50,30,41,43,43,41,30,50,38,25,42,41)
coord.dat <- matrix(c(1,3,5,7,2,4,6,8,1,3,5,7,2,4,6,8,
                      4.4,4.4,4.4,4.4,3.3,3.3,3.3,3.3,
                      2.2,2.2,2.2,2.2,1.1,1.1,1.1,1.1),16,2)

## Obtaining a list of neighbours:
library(spdep)
listW <- nb2listw(tri2nb(coord.dat), style="B")
links.mat.dat <- listw2mat(listW)
neighbors <- listw2sn(listW)[,1:2]

## Calculating the (Euclidean) distance between points:
D.dat <- dist(dat)

## Display the points:
plot(coord.dat, type='n',asp=1)
title("Delaunay triangulation")
text(coord.dat, labels=as.character(as.matrix(dat)), pos=3)
for(i in 1:nrow(neighbors))
    lines(rbind(coord.dat[neighbors[i,1],],
          coord.dat[neighbors[i,2],]))

## Clustering with a contiguity constraint described by a list of
## links:
grpWD2cst_constr_hclust <-
    constr.hclust(
        D.dat, method="ward.D2",
        neighbors, coord.dat)

## Plot the results with k=5 clusters on a map:
plot(grpWD2cst_constr_hclust, k=5, links=TRUE, las=1,
     xlab="Eastings", ylab="Northings", cex=3, lwd=3)

## Repeat the plot with other values of k (number of groups)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.constr.hclust", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.orthobasisSp")
### * plot.orthobasisSp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.orthobasisSp
### Title: Function to display Moran's Eigenvector Maps (MEM) and other
###   spatial orthogonal bases
### Aliases: plot.orthobasisSp

### ** Examples

if(require("ade4", quietly = TRUE) & require("spdep", quietly = TRUE)){
data(mafragh)
me <- mem(nb2listw(mafragh$nb))

if(require("adegraphics", quietly = TRUE)){
plot(me[,1:6], mafragh$xy)
plot(me[,1:6], mafragh$Spatial) 
}
}
        



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.orthobasisSp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rotation")
### * rotation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rotation
### Title: Rotate a set of point by a certain angle
### Aliases: rotation
### Keywords: manip

### ** Examples


### Create a set of coordinates
coords<-cbind(runif(20),runif(20))

### Create a series of angles
rad<-seq(0,pi,l=20)

for(i in rad){
	coords.rot<-rotation(coords,i)
	plot(coords.rot)
}

### Rotate the coordinates by an angle of 90 degrees
coords.90<-rotation(coords,90*pi/180)
coords.90

plot(coords,xlim=range(rbind(coords.90,coords)[,1]),ylim=range(rbind(coords.90,coords)[,2]),asp=1)
points(coords.90,pch=19)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rotation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("scalogram")
### * scalogram

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: scalogram
### Title: Function to compute a scalogram
### Aliases: scalogram plot.scalogram
### Keywords: spatial

### ** Examples

if(require("ade4", quietly = TRUE) & require("spdep", quietly = TRUE)){
data(mafragh)
me <- mem(nb2listw(mafragh$nb))

if(require("adegraphics", quietly = TRUE)){
sc1 <- scalogram(mafragh$env$Conduc, me, nblocks = 10)
plot(sc1) 
}
}
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("scalogram", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stimodels")
### * stimodels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stimodels
### Title: Space-time interaction in ANOVA without replication
### Aliases: stimodels quicksti
### Keywords: multivariate spatial

### ** Examples

# The "trichoptera" data set contains the abundances of 56 Trichoptera species captured 
# during 100 consecutive days in 22 emergence traps along a stream. The 22 traps 
# (sites) form a regular transect, with geographic positions 1 to 22. The original 
# daily data collected at each site were pooled into 10 survey periods for the study 
# of Legendre et al. (2010) in order to reduce the very high proportion of zeros in the 
# original data matrix. Order of the observations in the data set: the 22 traps (sites) 
# are nested within the survey weeks, as required by the 'stimodels' and 'quicksti' 
# functions.

data(trichoptera)

# log-transform the species data, excluding the Site and Date colums

trich.log <- log1p(trichoptera[,-c(1,2)]) 

# A log-chord transformation (Legendre & Borcard 2018) would also be appropriate for 
# these data: trich.tr <- decostand(log1p(trichoptera[,-c(1,2)]), method="norm")

# Example 1. Compute the space-time interaction test using model 5. By specifying the  
# number of sites (traps), the sofware assumes that they form a regular transect with  
# equispaced sites. Important note to users – In real analyses, use more than 99 
# permutations.

out.1 <- stimodels(trich.log, S=22, Ti=10, model="5", nperm=99)

# The interaction is significant. Because of that, test results for the main effects, 
# space and time, obtained with model 5, cannot be interpreted. Tests of spatial 
# variation can be done for individual times using simple RDA against dbMEM.S  
# variables. Likewise, tests of temporal variation can be done for individual sites  
# using simple RDA against dbMEM.T variables. A global test of the hypothesis that none 
# of the times shows a significant spatial structure can be done with model 6a. For a  
# global test of temporal structure at the various sites, use model 6b.

## No test: 
     # Code not run during CRAN software tests

# Example 2. Run space-time analysis with global tests for main effects after testing
# the interaction, which is significant in this example 

out.2 <- quicksti(trich.log, S=22, Ti=10, nperm=999)

# Since the interaction is significant, function 'quicksti' will carry out the 
# tests of existence of a spatial (at least in one of the time periods) and temporal 
# (at least at one of the sites) structures using models 6a and 6b, respectively.

# 3. Run space-time analysis for two time blocks only, i.e. time periods 6 and 7,   
# then time periods 8 and 9.

# Example 3.1. Time periods 6 and 7. The interaction is not significant. In that case,  
# function 'quicksti' carries out the tests of the main effects using model 5.

out.3 <- quicksti(trich.log[111:154,], S=22, Ti=2, nperm=999)

# Example 3.2. Time periods 8 and 9. The interaction is significant. In that case,  
# 'quicksti' carries out the tests of the spatial effects using model 6a. Model 6b 
# cannot proceed with the test of the temporal effect because Ti=2. An explanation is 
# printed in the output list.

out.4 <- quicksti(trich.log[155:198,], S=22, Ti=2, nperm=999)

# 4. Illustrations of the use of 'COD.S' and 'COD.T' in STI analysis

# The following examples illustrate how to use other representations of the spatial or 
# temporal relationships among observations, through arguments 'COD.S' and 
# 'COD.T' of functions 'stimodels' and 'quicksti'. The trichoptera data 
# are used again.

# Example 4.1. Explicitly compute dbMEMs for the spatial structure along the regular 
# transect, using function 'dbmem' of adespatial, and provide it to 'stimodels' 
# or 'quicksti' as argument 'COD.S'. The dbMEMs must first be computed on the
# transect, then repeated (Ti-1) times to provide Ti repeats in total.

dbMEM.S1 <- as.matrix(dbmem(1:22))
dbMEM.S10 <- dbMEM.S1
for(j in 2:10) dbMEM.S10 <- rbind(dbMEM.S10, dbMEM.S1)
out.5 <- stimodels(trich.log, S=22, Ti=10, model="5", COD.S=dbMEM.S10, nperm=999)

# Results should be identical to those in output file out.1 of Example 1, except for 
# P-values which can vary slightly.

# Example 4.2. Assume now that the sampling sites have irregular positions, as 
# described by the following matrix of geographic coordinates 'xy.trich'. Provide 
# this matrix to argument S of 'stimodels'

xy.trich = matrix(c(1:5,11:15,21:25,31:35,41,42,rep(c(1,2),11)),22,2)
plot(xy.trich, asp=1)   # Plot a quick map of the site positions
out.6 <- stimodels(trich.log, S=xy.trich, Ti=10, model="5", nperm=999)

# Example 4.3. Compute a matrix of dbMEMs for the sites. The coding matrix provided to 
# argument 'COD.S' must contain repeated dbMEM.S codes because that matrix must have 
# the same number of rows as matrix Y. Construct coding matrix dbMEM.S10 containing the 
# dbMEM.S codes repeated 10 times. 

dbMEM.S1 <- as.matrix(dbmem(xy.trich))
dbMEM.S10 = dbMEM.S1
for(i in 1:9) dbMEM.S10 <- rbind(dbMEM.S10, dbMEM.S1)
out.7 <- stimodels(trich.log, S=22, Ti=10, model="5", COD.S=dbMEM.S10, nperm=999)

# Compare the results with those obtained in the output file out6, example 4.2.

# Careful: If an analysis requires a dbMEM coding matrix for 'COD.T', the dbMEM.T    
# codes must follow the required data arrangement: sites must be nested within times.
# The following function can be used to construct a dbMEM.T matrix.

MEM.T <- function(s, tt, coord.T=NULL)
  # Documentation of function MEM.T –
 # Generate a matrix of temporal eigenfunctions for input into stimodels, 
 # with sites nested within times.
 # Arguments –
 # s : number of space points (sites)
 # tt : number of time points
 # coord.T : optional matrix or vector giving the time point coordinates
  {
  n <- s*tt
  if(is.null(coord.T)) coord.T <- as.matrix(1:tt)
  MEM.TT <- as.matrix(dbmem(coord.T))
  dbMEM.T <- matrix(0,n,ncol(MEM.TT))    # Empty matrix to house dbMEM.T		
  beg.x <- seq(1, n, by=s)
  for(i in 1:tt) { # Fill tt blocks of rows with identical MEM.TT values
	  for(j in 1:s) dbMEM.T[(beg.x[i]+j-1),] <- MEM.TT[i,]
	  }
  dbMEM.T
  }

# Example of use of function MEM.T
 
dbMEM.T <- MEM.T(s=6, tt=5)
# Check the size of the dbMEM.T output matrix
dim(dbMEM.T)

## End(No test)
   # End of code not run during CRAN software tests




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stimodels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("test.W")
### * test.W

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: test.W
### Title: Function to compute and test eigenvectors of spatial weighting
###   matrices
### Aliases: test.W
### Keywords: spatial

### ** Examples


if(require(ade4) & require(spdep)){

data(oribatid)
# Hellinger transformation
fau <- sqrt(oribatid$fau / outer(apply(oribatid$fau, 1, sum), rep(1, ncol(oribatid$fau)), "*"))
# remove gradient effect
faudt <- resid(lm(as.matrix(fau) ~ as.matrix(oribatid$xy)))

# test a binary spatial weighting matrix
nbtri <- tri2nb(as.matrix(oribatid$xy))
tri.res <- test.W(faudt, nbtri)

maxi <- max(unlist(nbdists(nbtri, as.matrix(oribatid$xy))))

# test a simple spatial weighting function of the distance
f1 <- function(x) {1-(x)/(maxi)}
tri.f1 <- test.W(faudt, nbtri, f = f1, xy = as.matrix(oribatid$xy))

# test a spatial weighting function with various values of parameters
f2 <- function(x,dmax,y) {1-(x^y)/(dmax)^y}
tri.f2 <- test.W(faudt,nbtri, f = f2, y = 2:10, dmax = maxi, xy = as.matrix(oribatid$xy))
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("test.W", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tpaired.krandtest")
### * tpaired.krandtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tpaired.krandtest
### Title: Paired t-tests of differences between T1 and T2 for each species
### Aliases: tpaired.krandtest

### ** Examples


if(require("vegan", quietly = TRUE)) {

## Invertebrate communities subjected to insecticide treatment.

## As an example in their paper on Principal Response Curves (PRC), van den Brink & ter 
## Braak (1999) used observations on the abundances of 178 invertebrate species 
## (macroinvertebrates and zooplankton) subjected to treatments in 12 mesocosms by the 
## insecticide chlorpyrifos. The mesocosms were sampled at 11 occasions. The data, 
## available in the {vegan} package, are log-transformed species abundances, 
## y.tranformed = loge(10*y+1).

## The data of survey #4 will be compared to those of survey #11 in this example.  
## Survey #4 was carried out one week after the insecticide treatment, whereas the 
## fauna of the mesocosms was considered to have fully recovered from the treatment 
## at the time of survey #11.

data(pyrifos)

## The mesocosms had originally been attributed at random to the treatments. However,  
## to facilitate presentation of the results, they will be listed here in order of 
## increased insecticide doses: {0, 0, 0, 0, 0.1, 0.1, 0.9, 0.9, 6, 6, 44, 44} 
## micro g/L.

survey4.order = c(38,39,41,47,37,44,40,46,43,48,42,45)

survey11.order = c(122,123,125,131,121,128,124,130,127,132,126,129)

## Paired t-tests of differences between survey.4 and survey.11 for the p species

res <- tpaired.krandtest(pyrifos[survey4.order,],pyrifos[survey11.order,])

}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tpaired.krandtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tpaired.randtest")
### * tpaired.randtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tpaired.randtest
### Title: Permutational paired t-test
### Aliases: tpaired.randtest

### ** Examples


## Deer leg length, data from Zar (1999, p. 162).

deer <- matrix(c(142,140,144,144,142,146,149,150,142,148,138,136,147,139,143,141,143,
145,136,146),10,2)

rownames(deer) <- paste("Deer",1:10,sep=".")

colnames(deer) <- c('Hind.leg', 'Fore.leg')

res <- tpaired.randtest(deer[,1], deer[,2])   # Two-tailed test by default

## Compare the results to:  res2 = t.test(deer[,1], deer[,2], paired=TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tpaired.randtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("variogmultiv")
### * variogmultiv

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: variogmultiv
### Title: Function to compute multivariate empirical variogram
### Aliases: variogmultiv
### Keywords: spatial

### ** Examples


if(require(ade4)){
data(oribatid)
# Hellinger transformation
fau <- sqrt(oribatid$fau / outer(apply(oribatid$fau, 1, sum), rep(1, ncol(oribatid$fau)), "*"))
# Removing linear effect
faudt <- resid(lm(as.matrix(fau) ~ as.matrix(oribatid$xy))) 
mvspec <- variogmultiv(faudt, oribatid$xy, nclass = 20)
mvspec
plot(mvspec$d, mvspec$var,type = 'b', pch = 20, xlab = "Distance", ylab = "C(distance)")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("variogmultiv", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
