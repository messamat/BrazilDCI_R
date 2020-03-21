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
DCIfunc <- DCIp_opti
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

#------ DEAL WITH COASTAL DRAINAGES ------
#Keep only network within basins that contain dams, excluding coastal segments without dams
#(only keep segments that are either in UpSeg or DownSeg of a dam and within a basin with a dam)
Basin08WithDams <- DamAttributesCrude[!is.na(HYBAS_ID08), unique(HYBAS_ID08)] #vector of all basins that contain a dam
NetworkBRAZIL <- NetworkBRAZILCrude[which((paste(NetworkBRAZILCrude$SEGID, NetworkBRAZILCrude$HYBAS_ID08) %in% 
                                             paste(DamAttributesCrude$UpSeg, DamAttributesCrude$HYBAS_ID08)) |
                                            (paste(NetworkBRAZILCrude$SEGID, NetworkBRAZILCrude$HYBAS_ID08) %in% 
                                               paste(DamAttributesCrude$DownSeg, DamAttributesCrude$HYBAS_ID08)) &
                                            NetworkBRAZILCrude$HYBAS_ID08 %in% Basin08WithDams),]
NetworkBRAZIL$SEGIDBAS <- with(NetworkBRAZIL, paste(SEGID, HYBAS_ID08, sep='_'))


#Create a new ID for dams that combines segment IDs and basin IDs
DamAttributesCrude$DownSegBAS <- with(DamAttributesCrude, paste(DownSeg, HYBAS_ID08, sep='_'))
DamAttributesCrude$UpSegBAS <- with(DamAttributesCrude, paste(UpSeg, HYBAS_ID08, sep='_'))
DamAttributesCrude$SEGIDBAS <- with(DamAttributesCrude, paste(SEGID, HYBAS_ID08, sep='_'))

#Identify basins with multiple networks that have dams 
#(Find those basins that have multiple segments with no downstream dam i.e. multiple outlets)
multibas <- as.data.frame(table(NetworkBRAZIL[(NetworkBRAZIL$SEGID %in% DamAttributesCrude$DownSeg) & 
                                                !(NetworkBRAZIL$SEGID %in% DamAttributesCrude$UpSeg), 'HYBAS_ID08']))

#For those, assign an extra digit to basins to account for those that have multiple networks with dams
coastdamSEG <- NetworkBRAZIL[(NetworkBRAZIL$SEGID %in% DamAttributesCrude$DownSeg) & 
                               !(NetworkBRAZIL$SEGID %in% DamAttributesCrude$UpSeg) &
                               NetworkBRAZIL$HYBAS_ID08 %in% multibas[multibas$Freq > 1, 'Var1'],]
coastdamSEG <- ddply(coastdamSEG, .(HYBAS_ID08), mutate, 
                     HYBAS_ID08ext = paste0(HYBAS_ID08, seq_along(SEGID)))

#Find connected segments within that basin to identify them as part of that subbasin
seg_list <- c()
for (segnum in seq_along(coastdamSEG$HYBAS_ID08)) {
  seg_mouth <- coastdamSEG[segnum,]
  seg_iter <- seg_mouth
  while (nrow(seg_iter)>0) {
    seg_list <- c(seg_list, paste(seg_iter$SEGID, seg_mouth$HYBAS_ID08ext, sep='_'))
    seg_iter <- DamAttributesCrude[DamAttributesCrude$DownSegBAS %in% seg_iter$SEGIDBAS,] 
  }
}

#Assign new subbasin ID too all river segments (simply HYBAS_ID08 + '1' if not multiple networks)
NetworkBRAZIL <- merge(NetworkBRAZIL, 
                       data.frame(HYBAS_ID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                                  SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1),
                                  stringsAsFactors = F), 
                       by='SEGIDBAS', all.x=T)
NetworkBRAZIL[is.na(NetworkBRAZIL$HYBAS_ID08ext), "HYBAS_ID08ext"] <-  NetworkBRAZIL[is.na(NetworkBRAZIL$HYBAS_ID08ext), 
                                                                                     paste0(HYBAS_ID08, 1)]

#Assign new subbasin ID too all dams (simply HYBAS_ID08 + '1' if not multiple networks)
DamAttributesCrude <- merge(DamAttributesCrude, 
                            data.frame(HYBAS_ID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                                       SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1),
                                       stringsAsFactors = FALSE), 
                            by='SEGIDBAS', all.x=T)
DamAttributesCrude[is.na(HYBAS_ID08ext), "HYBAS_ID08ext"] <-  DamAttributesCrude[is.na(HYBAS_ID08ext),
                                                                                 paste0(HYBAS_ID08, 1)]

## Remove a basin with problems (Probably the dam is too close to the upstream edge)
DamAttributesCrude <- DamAttributesCrude[-which(DamAttributesCrude$HYBAS_ID08 == 6080595090),]
NetworkBRAZIL <- NetworkBRAZIL[-which(NetworkBRAZIL$HYBAS_ID08 == 6080595090),]

#------ FORMAT DATA ------
## Final Dataframes to run DCI analyses
DamAttributes <- DamAttributesCrude
NetworkBRAZIL <- NetworkBRAZIL                     


# Organize the matrix based on type and stage of each dam
DamAttributes[, ESTAGIO_1 := ifelse(ESTAGIO_1 != "OperaÃ§Ã£o", 'Planned', 'Operation')]
DamAttributes[, Tipo_1 := ifelse(Tipo_1 != "UHE", "SHP", "LHP")]


## Fix mistakes in the dataset
DamAttributes$Tipo_1[which(DamAttributes$Tipo_1 == "SHP" & DamAttributes$POT_KW > 30000)] <- "LHP"    #UHE Buritizal(Das Mortes), UHE Resplendor(Doce), UHE Pouso Alto (Sucuriú)
DamAttributes$Tipo_1[which(DamAttributes$Tipo_1 == "LHP" & DamAttributes$POT_KW < 30000 &
                             DamAttributes$AREA_NA_MA < 13.0 & DamAttributes$ESTAGIO_1 == "Planned")] <- "SHP"   #Keep old dams as UHEs and new ones as SHPs


## Create a vector with unique basin IDs
basinList <- as.character(unique(DamAttributes$HYBAS_ID08ext))

## Make a vector of basins that are currently free of hydropower
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIp.csv", header = T))
HydroFreeBasins <- as.character(resuDCI$HYBAS_ID[resuDCI$All_curr == 100])

## Create a vector of NAs to be filled out with DCIs of each basin
DCIbasins <- rep(NA, times = length(basinList))

#Convert data frmes to data tables to speed up analysis
DamAttributes <- as.data.table(DamAttributes)
NetworkBRAZIL <- as.data.table(NetworkBRAZIL)

DamAttributes[, `:=`(Allcurrent = ifelse(DamAttributes$ESTAGIO_1 == 'Operation', 0.1, 1),
                     Allfuture = 1,
                     SHPcurrent = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                            DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHPcurrent = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                            DamAttributes$Tipo_1 == 'LHP', 0.1, 1),
                     SHPfuture = ifelse(DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHPfuture = ifelse(DamAttributes$Tipo_1 == 'LHP', 0.1, 1)
)]

setnames(DamAttributes, 'HYBAS_ID08ext', 'DAMBAS_ID08ext')

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

# ## Save as csv
# write.csv(NationalScen, file = "NationalScen_DCIp.csv")
# # write.csv(NationalScen, file = "NationalScen_DCIi.csv")

###########################################################################################################
################ Plot PARETO and select OPTIMAL and NON-OPTIMAL scenarios #################################
###########################################################################################################

## Loand packages
require(rPref)

## Import data on scenarios
NationalScenarios <- read.csv("NationalScen_DCIp.csv", header = T)
#NationalScenarios <- read.csv("NationalScen_DCIi.csv", header = T)

Best <- psel(NationalScenarios, pref = high(NationalScenarios$AddCapacity) * high(NationalScenarios$NatAverageDCI))
Worst <- psel(NationalScenarios, pref = low(NationalScenarios$AddCapacity) * low(NationalScenarios$NatAverageDCI))


## Plot 6 (DCI loss Vs Capacity)
tiff(filename = "Figure6.tiff", height = 2396, width = 3500, res = 300, compression = c("lzw"))
par(oma = c(8, 7.5, 7, 2), mar = c(0.5, 0, 0, 0), bty = "n")

## Create a matrix for the layout function
mat <- matrix(c(1, 1, 2,
                1, 1, 1,
                1, 1, 1), ncol = 3, nrow = 3, byrow = T)

layout(mat = mat, widths = c(3.5, 3.5, 6.0), heights = c(rep(3.5, 0.3, 0.3)))



## Plot Individual DCIs Vs Capacity
plot(NationalScenarios$AddCapacity, NationalScenarios$NatAverageDCI, type = "n", 
     ylim = c(60, 85), xlim = c(0, 65), ylab = "", xlab = "", xaxt = "n", yaxt = "n")
     
axis(side = 2, at = c(60, 65, 70, 75, 80, 85), cex.axis = 3.2, line = 0)
axis(side = 1, at = c(0, 10, 20, 30, 40, 50, 60), cex.axis = 3.2, line = 0, mgp = c(3, 1.9, 0))
     
mtext("Nation-wide river connectivity (DCI)", side = 2, cex = 2.7, line = 4.9)
mtext("Nation-wide capacity gain (gigawatts)", side = 1, cex = 2.7, line = 6.1)

## Plot Max and min future energy demands
segments(x0 = 15, x1 = 15, y0 = 0, y1 = 85, col = "gray60", lty = 3, lwd = 4)
segments(x0 = 36, x1 = 36, y0 = 0, y1 = 73, col = "gray60", lty = 3, lwd = 4)


## Plot average estimates
points(NationalScenarios$AddCapacity/1000, NationalScenarios$NatAverageDCI, pch = 21, col = "#42424220", 
       bg = "#9400D340", cex = 1.6)

## Plot best and worse scenarios
points(Worst$AddCapacity/1000, Worst$NatAverageDCI, pch = 21, col = "#424242", bg = "white", cex = 2.1)
points(Best$AddCapacity/1000, Best$NatAverageDCI, pch = 21, col = "black", bg = "#303030", cex = 2.1)



## Select the best and worst scenarios that fall inside thresholds of future demands

## Define the window of future demands
CurrentGen_GW <- sum(DamAttributes$POT_KW[which(DamAttributes$ESTAGIO_1 == "Operation")])/1000/1000

DemandLow <- 120 - CurrentGen_GW
DemandHigh <- 141 - CurrentGen_GW


BestDemand <- Best[Best$AddCapacity/1000 >= DemandLow & Best$AddCapacity/1000 <= DemandHigh, 1:6]
WorstDemand <- Worst[Worst$AddCapacity/1000 >= DemandLow & Worst$AddCapacity/1000 <= DemandHigh, 1:6]

## Organize data for boxplot
boxplot(BestDemand$NFutSHP, BestDemand$NFutLHP, WorstDemand$NFutSHP, WorstDemand$NFutLHP,
        col = c("#8B000099", "#4876FF99", "#8B000099", "#4876FF99"),
        ylim = c(0, 1500), at = c(1, 2, 3.3, 4.3), xaxt = "n", yaxt = "n")


axis(side = 2, at = c(0, 250, 500, 750, 1000, 1250, 1250), cex.axis = 2.4, line = 0)
axis(side = 1, at = c(1, 2), labels = c("SHP", "LHP"), cex.axis = 2.4, line = 0, mgp = c(3, 1.6, 0))
axis(side = 1, at = c(3.3, 4.3), labels = c("SHP", "LHP"), cex.axis = 2.4, line = 0, mgp = c(3, 1.6, 0))


mtext("Number of projects", side = 2, cex = 1.9, line = 4)
mtext("Optimal", side = 3, cex = 2.0, line = 1.5, at = 1.5)
mtext("Non-optimal", side = 3, cex = 2.0, line = 1.5, at = 3.9)


dev.off()




##############################################################################################################
##############################################################################################################

## Calculate some statistics
round(mean(BestDemand$NFutSHP), digits = 0)
round(mean(BestDemand$NFutLHP), digits = 0)

round(sd(BestDemand$NFutSHP), digits = 0)
round(sd(BestDemand$NFutLHP), digits = 0)

round(mean(WorstDemand$NFutSHP), digits = 0)
round(mean(WorstDemand$NFutLHP), digits = 0)

round(sd(WorstDemand$NFutSHP), digits = 0)
round(sd(WorstDemand$NFutLHP), digits = 0)

## Find an example to illustrate added capacity
round(BestDemand[order(BestDemand$AddCapacity),], digits = 0) 
## 25949 GW      616 All     553 SHPs     63 LHPs         DCI 75

round(WorstDemand[order(WorstDemand$AddCapacity),], digits = 0)
## 25726 GW     1271 All    1149 SHPs    122 LHPs     DCI 69


###############################################################################################################
## Plot DCI Vs N dams
plot(NationalScenarios$NFutDams, NationalScenarios$NatAverageDCI, ylab = "Average DCI", xlab = "Number of new dams")


## Plot Capacity Vs N dams
plot(NationalScenarios$NFutDams, NationalScenarios$AddCapacity, ylab = "Capacity gain in megawatts", xlab = "Number of new dams")
sky <- psel(NationalScenarios, pref = low(NationalScenarios$NFutDams) * high(NationalScenarios$AddCapacity))
#plot_front(NationalScenarios, pref = ChoicePref1, col = "blue")
points(sky$NFutDams, sky$AddCapacity, lwd = 3, col = "blue")
