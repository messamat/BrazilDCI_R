#############         Script 4 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ###############################
#############              Manuscript in preparation - PNAS                       ##############################

####################################################################################################################
###########        Sampling-based PRIORITIZATION ANALYSIS - BRAZIL         ########################################
###########                         (PARETO-FRONT Analysis)                #######################################
###########                        Sample future portifolios               ######################################
###########            (National-level multi-objective optimization)       #####################################
###############################################################################################################

## Choose what type of DCI function to run (DCIp or DCIi)
DCIfunc <- DCIp_opti3
# DCIfunc <- DCId_opti

## Choose the number of scenarios to sample
numbScen <- 100

## Packages
require(tictoc)
require(plyr)
require(bigstatsr)
require(parallel)
require(doParallel)
require(Rcpp)
library(fst)
devtools::source_gist("https://gist.github.com/r2evans/e5531cbab8cf421d14ed", filename = "lazyExpandGrid.R") #Get expand.grid version that won't run out of memory when > 25 dams


# # Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')
outdir_permut = file.path(resdir, 'outpermut_basins')

# # Import network and dams dataset (alternative)
# rootdir <- find_root(has_dir("PythonOutputs"))
# datadir <- file.path(rootdir, "PythonOutputs")
# dcigdb <- file.path(datadir, 'dci.gdb')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

## Create a vector with unique basin IDs
basinList <- DamAttributes[ESTAGIO_1 == "Planned", .N, by=DAMBAS_ID08ext] %>% 
  setorder(N) %>% 
  .[,DAMBAS_ID08ext]

## Make a vector of basins that are currently free of hydropower
resuDCI <- as.data.frame(fread("DCI_Brazil_L8_DCIp.csv", header = T))
HydroFreeBasins <- as.character(resuDCI$HYBAS_ID[resuDCI$All_curr == 100])


#Compute DCI for Allcurrent scenario
DCI_L8_current <- NetworkBRAZIL[,
                                list(DCI = DCIfunc(
                                  d2= DamAttributes[DAMBAS_ID08ext == HYBAS_ID08ext,
                                                    list(
                                                      id1 = DownSeg,
                                                      id2 = UpSeg,
                                                      pass = Allcurrent)],
                                  d3 = .SD[, list(id=as.character(SEGID),
                                                  l=Shape_Length)],
                                  print = F)),
                                by=.(HYBAS_ID08ext)]


## Get the total number of future dams
MaxFutDams <- sum(DamAttributes$ESTAGIO_1 == "Planned")

## Get the total number of possible portfolios for each basin
basinpermutlen <- DamAttributes[ESTAGIO_1 == "Planned", 2^(.N), by=DAMBAS_ID08ext]
sum(log10(basinpermutlen$V1)) #Total number of possible portfolios in brazil 10^682

#Import scenarios with 0 dams
outdir_permut = file.path(resdir, 'outpermut_basins')
allscens <- lapply(file.path(outdir_permut, list.files(outdir_permut)), 
                   function(x) read.fst(x, columns = c("scenbasin", "DAMBAS_ID08ext", "ndams", "POT_KWbasin", "SHPnum", "LHPnum",
                                                       "prevfree", "DCI"))) %>%
  do.call(rbind, .) %>%
  setDT

#Import scenarios with 0 dams
outdir_permut = file.path(resdir, 'outpermut_basins2')
allscens2 <- lapply(file.path(outdir_permut, list.files(outdir_permut)), 
                   function(x) read.fst(x, columns = c("scenbasin", "DAMBAS_ID08ext", "ndams", "POT_KWbasin", "SHPnum", "LHPnum",
                                                       "prevfree","DCI"))) %>%
  do.call(rbind, .) %>%
  setDT

allscens <- allscens[, .(scenbasin, DAMBAS_ID08ext, DCI)][
  allscens2[, .(scenbasin, DAMBAS_ID08ext, ndams, POT_KWbasin, SHPnum, LHPnum, prevfree), ],
  on=.(scenbasin, DAMBAS_ID08ext)]

for (j in unique(allscens$DAMBAS_ID08ext)) {
  print(j)
  write.fst(allscens[DAMBAS_ID08ext == j,], file.path(outdir_permut, paste0('permut_', j, '.fst')))
}


#Add scenario with no new dam
nbasins <- allscens[, length(unique(DAMBAS_ID08ext))]

#Generate 100,000 basin-specific scenarios and then compute statistics over them
nsamp=10000
tic()
check <- replicate(n=10, simplify=FALSE,
                   exp = 
                     DCIportfolio <- allscens[, .SD[sample(.N, nsamp, replace=T)], by=DAMBAS_ID08ext] %>%
                     .[, brazilscen := rep(1:nsamp, nbasins)] %>% 
                     .[, list(NatAverageDCI = mean(DCI),
                              AddCapacity = sum(POT_KWbasin),
                              NFutDams = sum(ndams), 
                              NFutSHP = sum(SHPnum),
                              NFutLHP = sum(LHPnum),
                              NFreeDammed = sum(prevfree, na.rm=T)
                     ), by=brazilscen]) %>%
  do.call(rbind, .) %>%
  setDT
toc()
  














#Generate dummy dt
dummyscens <- DamAttributes[ESTAGIO_1 == "Planned", .(scenbasin=seq(1,2^(.N)), 
                                                      DCI=runif(2^(.N)),
                                                      ndams=.N/2, #nnew dams
                                                      POT_KWbasin = .N/2*1000,
                                                      SHPnum = .N/3,
                                                      LHPnum = .N*2/3,
                                                      prevfree = 1),
                            by=DAMBAS_ID08ext]
dummyscens




#####################################################################################
#Get list of all scenarios
Scens <- lapply(sample(x = 1:MaxFutDams, size = numbScen, replace=T), function(N) {
  sample(unlist(ListIDs), N, FALSE, NULL)
})

#Initialize parallel analysis
cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
on.exit(stopCluster(cl))
doParallel::registerDoParallel(cl)

#Launch DCI analysis for all scenarios in parallel
DCIscens <- foreach(i=1:numbScen, 
                    .packages = c("data.table", "Rcpp", "rccpcomb"), 
                    .noexport = c('combi2inds')) %dopar% {
  scen_DAMIDlist <- Scens[[i]]
  
  ## Get the number of future dams in that scenario
  ScenNdams <- length(scen_DAMIDlist)

  ## Get the number of basins that will be no longer free-flowing in this given scenario
  UniqueBasin <- unique(DamAttributes[DAMID %in% scen_DAMIDlist, DAMBAS_ID08ext])  
  NFreeDammed <- sum(UniqueBasin %in% HydroFreeBasins)
   
  #Subset dam attributes to only process DCI in basins with changed dam composition
  DamAttributesScen <- DamAttributes[DAMBAS_ID08ext %in% UniqueBasin,]

  ## Change the permeability value for 0.1 for the sampled dams from that scenario (plus keep all dams from All Current scenario)
  DamAttributesScen[DAMID %in% scen_DAMIDlist, Allcurrent := 0.1]

  #### rbind DCI Analysis for all basins inside this scenario with all other basins
  DCI_L8_ScenPlanned <- rbind(DCI_L8_current[!(HYBAS_ID08ext %in% UniqueBasin),],
                              NetworkBRAZIL[HYBAS_ID08ext %in% UniqueBasin,
                                            list(DCI = DCIfunc(
                                              DamAttributesScen[DAMBAS_ID08ext == HYBAS_ID08ext,
                                                                list(
                                                                  id1 = DownSeg,
                                                                  id2 = UpSeg,
                                                                  pass = Allcurrent
                                                                )],
                                              .SD[, list(id=as.character(SEGID),
                                                         l=Shape_Length)],
                                              print = F)),
                                            by=.(HYBAS_ID08ext)]
  )

  return(list(
    mean(DCI_L8_ScenPlanned$DCI, na.rm = T),
    sum(DamAttributes$POT_KW[DamAttributes$DAMID %in% scen_DAMIDlist]/1000),
    ScenNdams,
    sum(DamAttributes$Tipo_1[DamAttributes$DAMID %in% scen_DAMIDlist] == "SHP"),
    sum(DamAttributes$Tipo_1[DamAttributes$DAMID %in% scen_DAMIDlist] == "LHP"),
    NFreeDammed,
    toString(scen_DAMIDlist))
  )
}
toc()

# Summarized national-level data in data.frame format (average DCIs, capacity gain)
NationalScen <- as.data.frame(do.call("rbind", lapply(DCIscens, unlist)))
colnames(NationalScen) <- c("NatAverageDCI", "AddCapacity", "NFutDams", "NFutSHP", "NFutLHP",  "NFreeDammed", "DamIDs")

#Save as csv
fwrite(NationalScen, file = "NationalScen_DCIp.csv")
#fwrite.csv(NationalScen, file = "NationalScen_DCIi.csv")

