

library(bigstatsr)
library(data.table)
library(devtools)
library(doParallel)
library(fst)
library(ggplot2)
library(igraph) # library graph theory
library(magrittr)
library(maps)
library(parallel)
library(plyr)
library(RColorBrewer)
library(Rcpp)
library(RcppAlgos)
library(rgdal)
library(rPref)
library(rprojroot)
library(sf)
library(sp)
library(stringr)
library(tictoc)
library(zoom)

#Get expand.grid version that won't run out of memory when > 25 dams
devtools::source_gist("https://gist.github.com/r2evans/e5531cbab8cf421d14ed", filename = "lazyExpandGrid.R") 