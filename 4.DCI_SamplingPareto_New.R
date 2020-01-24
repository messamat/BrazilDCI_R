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
numbScen <- 5

## Packages
require(tictoc)
require(plyr)

# ## Import network and dams dataset
# DamAttributesCrude <- read.csv("damattributes.txt", header = T)
# NetworkBRAZILCrude <- read.csv("networkattributes.txt", header = T)

## Import network and dams dataset (alternative)
rootdir <- find_root(has_dir("PythonOutputs"))
datadir <- file.path(rootdir, "PythonOutputs")
dcigdb <- file.path(datadir, 'dci.gdb')
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
                                  SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1)), 
                       by='SEGIDBAS', all.x=T)
NetworkBRAZIL[is.na(NetworkBRAZIL$HYBAS_ID08ext), "HYBAS_ID08ext"] <-  NetworkBRAZIL[is.na(NetworkBRAZIL$HYBAS_ID08ext), 
                                                                                     paste0(HYBAS_ID08, 1)]

#Assign new subbasin ID too all dams (simply HYBAS_ID08 + '1' if not multiple networks)
DamAttributesCrude <- merge(DamAttributesCrude, 
                            data.frame(HYBAS_ID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                                       SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1)), 
                            by='SEGIDBAS', all.x=T)
DamAttributesCrude[is.na(HYBAS_ID08ext), "HYBAS_ID08ext"] <-  DamAttributesCrude[is.na(HYBAS_ID08ext),
                                                                                 paste0(HYBAS_ID08, 1)]

## Remove a basin with problems (Probably the dam is too close to the upstream edge)
DamAttributesCrude <- DamAttributesCrude[-which(DamAttributesCrude$HYBAS_ID08 == 6080595090),]
NetworkBRAZIL <- NetworkBRAZIL[-which(NetworkBRAZIL$HYBAS_ID08 == 6080595090),]

## Final Dataframes to run DCI analyses
DamAttributes <- DamAttributesCrude
NetworkBRAZIL <- NetworkBRAZIL                     


# Organize the matrix based on type and stage of each dam
levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) != "OperaÃ§Ã£o"] <- "Planned"
levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) == "OperaÃ§Ã£o"] <- "Operation"
levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) != "UHE"] <- "SHP"
levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) == "UHE"] <- "LHP"

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

## Create a matrix to be filled with the summarized national-level data (average DCIs, capacity gain)
NationalMat <- matrix(NA, ncol = 7,  nrow = numbScen)
colnames(NationalMat) <- c("NatAverageDCI", "AddCapacity", "NFutDams", "NFutSHP", "NFutLHP",  "NFreeDammed", "DamIDs")
NationalScen <- as.data.frame(NationalMat)

## Create a vector to be filled out with permeability values (operation = 0.1 and Planned = 1)
Permeability <- rep(NA, times = dim(DamAttributes)[1])
Permeability[DamAttributes$ESTAGIO_1 == "Operation"] <- 0.1
Permeability[DamAttributes$ESTAGIO_1 == "Planned"] <- 1

## Get the total number of future dams
MaxFutDams <- sum(DamAttributes$ESTAGIO_1 == "Planned")

## Create a list IDs of planned dams to be sampled
DamID <- 1:dim(DamAttributes)[1]
ListIDs <- DamID[DamAttributes$ESTAGIO_1 == "Planned"]
MaxNListIDs <- length(ListIDs)


tic("total")
## Loop over the scenarios
for (n in 1: numbScen){
  
  ## Get the number of future dams in that scenario
  ScenNdams <- sample(x = 1:MaxFutDams, size = 1)
  
  ## Sample the planned dams IDs that will get built in that scenario
  SampledDamsIDs <- sample(x = 1:MaxNListIDs, size = ScenNdams)
  SampledIDs <- ListIDs[SampledDamsIDs]
  
  ## Change the permeability value for 0.1 for the sampled dams from that scenario
  ScenPermeability <- Permeability
  ScenPermeability[which(DamID %in% SampledIDs)] <- 0.1
  
  ## Get the number of basins that will be no longer free-flowing in this given scenario
  UniqueBasin <- as.character(unique(DamAttributes$HYBAS_ID08ext[which(DamID %in% SampledIDs)]))
  UniqueBasinIDs <- substr(UniqueBasin, 1, nchar(basinList)-1)
  NFreeDammed <- sum(UniqueBasinIDs %in% HydroFreeBasins)
  
  #### DCI Analysis for all basins inside this scenario
  ## Loop over basins
  for (j in 1: length(basinList)){
    
    ## filter attributes of the basin j         
    BasinX <- NetworkBRAZIL[NetworkBRAZIL$HYBAS_ID08ext == basinList[j], ]
    DamX <- DamAttributes[DamAttributes$HYBAS_ID08ext == basinList[j], ]
    PermeabilityX <- ScenPermeability[DamAttributes$HYBAS_ID08ext == basinList[j]]
    
    # Create a sequence ranging from 1 to the maximum number of segments
    # Lsegments <- min(unique(BasinX$SEGID))
    Numbsegments <- length(unique(BasinX$SEGID))
    
    ## Determine that basins with no dams have DCI = 100 (for current scenarios or for just one type of hydropower)
    if (Numbsegments == 1){
      DCI <- 100
    }
    
    ## For basins with dams, DCI is calculated below
    if (Numbsegments > 1){
      
      ## Create vectors of segids and their total lenght
      listSeg <- paste("seg", BasinX$SEGID, sep ="")
      listLengths <- BasinX$Shape_Length
      
      
      ## Compile the data on the edges for the DCI analysis for the basin j
      DowSeg <- paste("seg", DamX$DownSeg, sep ="")
      UpSeg <- paste("seg", DamX$UpSeg, sep ="")
      DamPermeability <- PermeabilityX
      
      # EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin, Capacity, DamPermeability)
      EdgesData <- data.frame(DowSeg, UpSeg, DamPermeability)
      
      
      ### Run the DCI analysis for the basin j
      # attributes of edges; links between nodes
      d2 = data.frame (id1 = EdgesData$DowSeg, id2 = EdgesData$UpSeg, pass = EdgesData$DamPermeability)
      
      # attributes of nodes; node sizes
      d3 = data.frame(id = listSeg, l = listLengths)
      
      # Run DCI analysis for the basin
      DCI <- DCIfunc (d2, d3, print = F)
      
      # Fill out with DCI values
      DCIbasins[j] <- DCI
      
    }
    
    ## Print indexes to follow the calculations
    print(c(n, j))
    
  }
  
  NationalScen[n, 1] <- mean(DCIbasins, na.rm = T)
  NationalScen[n, 2] <- sum(DamAttributes$POT_KW[DamID %in% SampledIDs]/1000)
  NationalScen[n, 3] <- ScenNdams
  NationalScen[n, 4] <- sum(DamAttributes$Tipo_1[DamID %in% SampledIDs] == "SHP")
  NationalScen[n, 5] <- sum(DamAttributes$Tipo_1[DamID %in% SampledIDs] == "LHP")
  NationalScen[n, 6] <- NFreeDammed
  NationalScen[n, 7] <- toString(SampledIDs)
  
}

toc()


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
