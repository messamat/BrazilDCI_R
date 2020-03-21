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
require(rccpcomb)

# # Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')

# # Import network and dams dataset (alternative)
# rootdir <- find_root(has_dir("PythonOutputs"))
# datadir <- file.path(rootdir, "PythonOutputs")
# dcigdb <- file.path(datadir, 'dci.gdb')
NetworkBRAZILCrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='networkattributes'))
DamAttributesCrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='damattributes'))

## Remove the dams that have NAs from the dataset
DamAttributesCrude <- DamAttributesCrude[!is.na(DamAttributesCrude$HYBAS_ID08),]

## Create a vector with unique basin IDs
basinList <- as.character(unique(DamAttributes$HYBAS_ID08ext))

## Make a vector of basins that are currently free of hydropower
resuDCI <- as.data.frame(fread("DCI_Brazil_L8_DCIp.csv", header = T))
HydroFreeBasins <- as.character(resuDCI$HYBAS_ID[resuDCI$All_curr == 100])


#Compute DCI for Allcurrent scenario
DCI_L8_current <- NetworkBRAZIL[,
                                list(DCI = DCIp_opti(
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

## Create a list IDs of planned dams to be sampled
ListIDs <- DamAttributes[ESTAGIO_1 == "Planned", 'DAMID']
MaxNListIDs <- nrow(ListIDs)

tic()
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

