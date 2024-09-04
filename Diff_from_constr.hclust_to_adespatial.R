
## Guillaume Guénard, Université de Montréal, 2018-2024

## rm(list=ls())

constr.hclust_dir <- "~/git/constr.hclust"
adespatial_dir <- "~/git/adespatial"

files_to_transfer <- c(
  "src/constr.hclust.c",
  "src/constr.hclust.h",
  "R/constr.hclust.R",
  "R/constr.hclust-class.R",
  "R/plot.constr.hclust.R",
  "R/ScotchWhiskey.R"
)

## i <- 1L
for(i in 1L:length(files_to_transfer))
  system(
    sprintf(
      "meld %s %s",
      file.path(constr.hclust_dir,files_to_transfer[i]),
      file.path(adespatial_dir,files_to_transfer[i])
    )
  )

## Le fichier R/Tiahura.R ne se trouve pas dans constr.hclust mais directement
## dans adespatial_dev.
system(
  sprintf(
    "meld %s %s",
    "Tiahura.R",
    file.path(adespatial_dir,"R/Tiahura.R")
  )
)

## Note:
## 1. File adespatial/src/init.c may sometimes require editing.
##    For instance, function prototype line: extern void cclust(...);
##    and static const R_CMethodDef CEntries[] = {...} entry line:
##       ...{"cclust",          (DL_FUNC) &cclust,          11},
##    will require editing if the function arguments are being altered, and
##    any new C function will require to be registered and maintained.
