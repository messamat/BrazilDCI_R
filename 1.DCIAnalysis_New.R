#############         Script 1 to run DCI analysis - September 2019                ###############################
#############       Based on distances calculated in Arcgis - PythonCode          ###############################
#############             Manuscript in preparation - PNAS                        ##############################


###############################################################################################################
#########     ESTIMATE DCI FOR ALL THE SCENARIOS - BRAZIL    #################################################
##########         (CURRENT VS FUTURE) - All, SHPs, LHPs    #################################################
###########               Level 8 All Scenarios            ################################################
#########################################################################################################

## Choose what type of DCI function to run (DCIp or DCIi)
DCIfunc <- DCIp_opti
# DCIfunc <- DCId_opti

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
                                                  
## Organize the matrix based on type and stage of each dam
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

## Create a matrix of basins and DCIs               
DCI_BRAZIL_L8 <- matrix(NA, nrow = length(basinList), ncol = 7)
colnames(DCI_BRAZIL_L8) <- c("HYBAS_ID", "All_curr", "All_fut", "SHP_curr", "SHP_fut", "LHP_curr", "LHP_fut")
DCI_BRAZIL_L8 [ , 1] <- as.numeric(substr(basinList, 1, nchar(basinList)-1))

## Create two matrixes to fill with dam atributes of each basin (N dams, MW production)
basinNDams_L8 <- matrix(NA, nrow = length(basinList), ncol = 7)
colnames(basinNDams_L8) <- c("HYBAS_IDext", "N_All_curr", "N_All_fut", "N_SHP_curr", "N_SHP_fut", "N_LHP_curr", "N_LHP_fut")
basinNDams_L8 [ , 1] <- as.numeric(substr(basinList, 1, nchar(basinList)-1))

basinMWDams_L8 <- matrix(NA, nrow = length(basinList), ncol = 7)
colnames(basinMWDams_L8) <- c("HYBAS_IDext", "MW_All_curr", "MW_All_fut", "MW_SHP_curr", "MW_SHP_fut", "MW_LHP_curr", "MW_LHP_fut")
basinMWDams_L8 [ , 1] <- as.numeric(substr(basinList, 1, nchar(basinList)-1))

## Possible Scenarios of time-period
Current <- "Operation"
Future <- c("Operation", "Planned")

## Possible Scenarios of hydro type
Small <- "SHP"
Large <- "LHP"
All <- c("SHP", "LHP")

## Create vectors with possible combination of scenarios
ScenarioSituation <- c("Current", "Future", "Current", "Future", "Current", "Future")
ScenarioType <- c("All", "All", "Small", "Small", "Large", "Large")

tic("total")
## Loop over scenarios
for (x in 1: 6){
  
  ## Define types and situation according to the scenario x
  SituationScen <- get(ScenarioSituation[x])
  TypeScen <- get(ScenarioType[x])
  
  
  ## Loop over basins
  for (j in 1: length(basinList)){
    
    ## filter attributes of the basin j         
    BasinX <- NetworkBRAZIL[NetworkBRAZIL$HYBAS_ID08ext == basinList[j], ]
    DamX <- DamAttributes[DamAttributes$HYBAS_ID08ext == basinList[j], ]
    
    # Create a sequence ranging from 1 to the maximum number of segments
    #Lsegments <- min(unique(BasinX$batNetID))
    Numbsegments <- length(BasinX$SEGID)
    
    ## If there is no dam in the network in that scenario, DCI = 100
    if (Numbsegments == 1){
      DCI <- 100
    }
    
    ## For basins with dams in that scenario, DCI is calculated below
    if (Numbsegments > 1){
      
      # Paste the name seg before each segment number and make a vector
      # Create a vector of fragment lenghts
      listSeg <- rep (NA, times = Numbsegments)
      listLengths <- rep (NA, times = Numbsegments)
      
      ## Loop over segments
      # Fill these two vectors
      for (i in 1: Numbsegments){
        listSeg[i] <- paste("seg", BasinX$SEGID[i], sep ="")
        listLengths[i] <- BasinX$Shape_Length [i]
      }  
      
      ## Compile the data on the edges (dam attributes)
      DowSeg <- paste("seg", DamX$DownSeg, sep ="")
      UpSeg <- paste("seg", DamX$UpSeg, sep ="")
      Type <- DamX$Tipo_1
      Situation <- DamX$ESTAGIO_1
      ID_number <- DamX$TARGET_FID
      ID_name <- DamX$NOME
      Basin <- DamX$HYBAS_ID08ext
      
      ## Put these data in a data frame with all necessary information about the edges
      EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin)
      
      ## Create a passibility vector based on the number of segments of the basin
      ## Assign 100% permeability for all dams with a vector
      passVec <- rep(1, times = dim(EdgesData)[1])
      
      ## Find the position of dams that should be present according to the scenarios (Type and Sitiation = T)
      PresentDams <- which(EdgesData$Type %in% TypeScen * EdgesData$Situation %in% SituationScen == 1)
      
      ## Assign a permeability value of 0.1 for the dams of the scenario based on their position
      passVec[PresentDams] <- 0.1
      
      # List of edges and passability
      d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[, 2], pass = passVec)
      
      # attributes of nodes; node sizes
      d3 = data.frame(id = listSeg, l = listLengths)
      
      # Run DCI analyses for the basin
      DCI <- DCIfunc (d2, d3, print = F)
      #DCIi <- DCIi(d2, d3, print = F) #Downseg as the reference
      
      
    }
    ## Compute the DCI of the Basin in the matrix
    DCI_BRAZIL_L8[j, x + 1] <- DCI
    
    ### Attributes of the basin (N of dams, cumulative MW)
    basinNDams_L8[j, x + 1] <- dim(DamX)[1]
    basinMWDams_L8[j, x + 1] <- sum(DamX$POT_KW, na.rm = T)/1000 #kw to MW
    
    ## Print Scenarios and basin numbers
    print(c(x,j))
    
  }  
 
}
toc()

## Create a csv file with the output of DCI analysis for all the scenarios
write.csv(DCI_BRAZIL_L8, file = "DCI_Brazil_L8_DCIp.csv")
write.csv(basinNDams_L8, file = "basinNDams_L8_DCIp.csv")
write.csv(basinMWDams_L8, file = "basinMWDams_L8_DCIp.csv")

# write.csv(DCI_BRAZIL_L8, file = "DCI_Brazil_L8_DCIi.csv")
# write.csv(basinNDams_L8, file = "basinNDams_L8_DCIi.csv")
# write.csv(basinMWDams_L8, file = "basinMWDams_L8_DCIi.csv")


#############################################################################################################
#############################################################################################################
#############################################################################################################
#####################################     PLOTS FIGURE 2     ################################################

## Import data
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIp.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8_DCIp.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8_DCIp.csv", header = T))

# ## Import data
# resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIi.csv", header = T))
# resuNDams <- as.data.frame(read.csv("basinNDams_L8_DCIi.csv", header = T))
# resuMWDams <- as.data.frame(read.csv("basinMWDams_L8_DCIi.csv", header = T))


## Filter results of basins with both SHPs and LHPs in the future
resuDCI_Both <- resuDCI[resuDCI$SHP_fut != 100 & resuDCI$LHP_fut != 100, ]

## Filter results just for basins that are currently free of hydropower
resuDCI_Free <- resuDCI[resuDCI$All_curr == 100, ]

## Filter results just for basins that are currently regulated
resuDCI_Regul <- resuDCI[resuDCI$All_curr != 100, ]


## Calculate differences
DCILoss_All <- resuDCI$All_fut - resuDCI$All_curr
DCILoss_SHP <- resuDCI$SHP_fut - resuDCI$SHP_curr
DCILoss_LHP <- resuDCI$LHP_fut - resuDCI$LHP_curr

DCILoss_Both_All <- resuDCI_Both$All_fut - resuDCI_Both$All_curr
DCILoss_Both_SHP <- resuDCI_Both$SHP_fut - resuDCI_Both$SHP_curr
DCILoss_Both_LHP <- resuDCI_Both$LHP_fut - resuDCI_Both$LHP_curr

DCILoss_Free_All <- resuDCI_Free$All_fut - resuDCI_Free$All_curr
DCILoss_Free_SHP <- resuDCI_Free$SHP_fut - resuDCI_Free$SHP_curr
DCILoss_Free_LHP <- resuDCI_Free$LHP_fut - resuDCI_Free$LHP_curr

DCILoss_Regul_All <- resuDCI_Regul$All_fut - resuDCI_Regul$All_curr
DCILoss_Regul_SHP <- resuDCI_Regul$SHP_fut - resuDCI_Regul$SHP_curr
DCILoss_Regul_LHP <- resuDCI_Regul$LHP_fut - resuDCI_Regul$LHP_curr


## Calculate % change (new - old / old * 100)
PercLoss_All <- DCILoss_All/resuDCI$All_curr * 100
PercLoss_SHP <- DCILoss_SHP/resuDCI$SHP_curr * 100
PercLoss_LHP <- DCILoss_LHP/resuDCI$LHP_curr * 100

PercLoss_Both_All <- DCILoss_Both_All/resuDCI_Both$All_curr * 100
PercLoss_Both_SHP <- DCILoss_Both_SHP/resuDCI_Both$SHP_curr * 100
PercLoss_Both_LHP <- DCILoss_Both_LHP/resuDCI_Both$LHP_curr * 100

PercLoss_Free_All <- DCILoss_Free_All/resuDCI_Free$All_curr * 100
PercLoss_Free_SHP <- DCILoss_Free_SHP/resuDCI_Free$SHP_curr * 100
PercLoss_Free_LHP <- DCILoss_Free_LHP/resuDCI_Free$LHP_curr * 100

PercLoss_Regul_All <- DCILoss_Regul_All/resuDCI_Regul$All_curr * 100
PercLoss_Regul_SHP <- DCILoss_Regul_SHP/resuDCI_Regul$SHP_curr * 100
PercLoss_Regul_LHP <- DCILoss_Regul_LHP/resuDCI_Regul$LHP_curr * 100

## Export csv with the percentage change
# PercChangeData <- data.frame(resuDCI$HYBAS_ID, PercLoss_All, PercLoss_SHP, PercLoss_LHP)
# colnames(PercChangeData) <- c("HYBAS_ID", "PercLoss_All", "PercLoss_SHP", "PercLoss_LHP")
# write.csv(PercChangeData, file = "DCI_Brazil_L8_PercentageDCIp.csv")

# # Export csv with the percentage change
# PercChangeData <- data.frame(resuDCI$HYBAS_ID, PercLoss_All, PercLoss_SHP, PercLoss_LHP)
# colnames(PercChangeData) <- c("HYBAS_ID", "PercLoss_All", "PercLoss_SHP", "PercLoss_LHP")
# write.csv(PercChangeData, file = "DCI_Brazil_L8_PercentageDCIi.csv")


###############################################################################################################
##### Plot % change
tiff(filename = "Figure2.tiff", height = 1656, width = 3200, res = 300, compression = c("lzw"))

mat <- matrix(c(rep(1, times = 2), rep(2, times = 2), rep(3, times = 2)), ncol = 3, nrow = 2, byrow = F)
layout(mat = mat, widths = c(4, 4, 4), heights = c(rep(0.5, times = 3)))

# Par
par (oma = c(5.5, 6, 6.5, 4.0), mar = c(0.5, 0, 0, 0), bty = "n")

par(bty = "n")

### Plot data for all basins

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI)[1])

plot(x = treatment, y = PercLoss_All, ylim = c(-100, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = PercLoss_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4.0, pch = 21)
points(x = treatment * 2, y = PercLoss_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4.0, pch = 21)
points(x = treatment * 3, y = PercLoss_All, bg = "#9400D305", col = c("#42424220"), cex = 4.0, pch = 21)

# Plot 95% CI
UppCI_SHP <- quantile(PercLoss_SHP, 0.975, na.rm = T)
LowCI_SHP <- quantile(PercLoss_SHP, 0.025, na.rm = T)
UppCI_LHP <- quantile(PercLoss_LHP, 0.975, na.rm = T)
LowCI_LHP <- quantile(PercLoss_LHP, 0.025, na.rm = T)
UppCI_All <- quantile(PercLoss_All, 0.975, na.rm = T)
LowCI_All <- quantile(PercLoss_All, 0.025, na.rm = T)

segments(y0 = UppCI_SHP, x0 = treatment * 1, 
         y1 = LowCI_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_LHP, x0 = treatment * 2, 
         y1 = LowCI_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_All, x0 = treatment * 3, 
         y1 = LowCI_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(PercLoss_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 2, at = c(0, -20, -40, -60, -80, -100), cex.axis = 2.5, line = - 3)
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.5, line = 1, mgp = c(3, 1.7, 0))
mtext("Change in river connectivity (%)", side = 2, cex = 2.2, line = 1.7)
#mtext("Predicted % change in river connectivity (DCI)", side = 2, cex = 2.0, line = 1.7)

mtext("All basins", side = 3, cex = 2.0, line = 1)
mtext("A", side = 3, cex = 2.6, line = 2.9, at = 0.2)


### Plot data just for hydropower-free basin
# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Free)[1])

plot(x = treatment, y = PercLoss_Free_All, ylim = c(-100, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = PercLoss_Free_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 2, y = PercLoss_Free_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 3, y = PercLoss_Free_All, bg = "#9400D305", col = c("#42424220"), cex = 4, pch = 21)

# Plot 95% CI
UppCI_Free_SHP <- quantile(PercLoss_Free_SHP, 0.975, na.rm = T)
LowCI_Free_SHP <- quantile(PercLoss_Free_SHP, 0.025, na.rm = T)
UppCI_Free_LHP <- quantile(PercLoss_Free_LHP, 0.975, na.rm = T)
LowCI_Free_LHP <- quantile(PercLoss_Free_LHP, 0.025, na.rm = T)
UppCI_Free_All <- quantile(PercLoss_Free_All, 0.975, na.rm = T)
LowCI_Free_All <- quantile(PercLoss_Free_All, 0.025, na.rm = T)

segments(y0 = UppCI_Free_SHP, x0 = treatment * 1, 
         y1 = LowCI_Free_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_Free_LHP, x0 = treatment * 2, 
         y1 = LowCI_Free_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Free_All, x0 = treatment * 3, 
         y1 = LowCI_Free_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(PercLoss_Free_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_Free_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_Free_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.5, line = 1, mgp = c(3, 1.7, 0))

mtext("Hydro-free basins", side = 3, cex = 1.6, line = 1, at = treatment * 3 - 0.42)
mtext("B", side = 3, cex = 2.6, line = 2.9, at = 0.2)


### Plot data just for basins with both SHPs and LHPs
# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Regul)[1])

plot(x = treatment, y = PercLoss_Regul_All, ylim = c(-100, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = PercLoss_Regul_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 2, y = PercLoss_Regul_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 3, y = PercLoss_Regul_All, bg = "#9400D305", col = c("#42424220"), cex = 4, pch = 21)

# Plot 95% CI
UppCI_Regul_SHP <- quantile(PercLoss_Regul_SHP, 0.975, na.rm = T)
LowCI_Regul_SHP <- quantile(PercLoss_Regul_SHP, 0.025, na.rm = T)
UppCI_Regul_LHP <- quantile(PercLoss_Regul_LHP, 0.975, na.rm = T)
LowCI_Regul_LHP <- quantile(PercLoss_Regul_LHP, 0.025, na.rm = T)
UppCI_Regul_All <- quantile(PercLoss_Regul_All, 0.975, na.rm = T)
LowCI_Regul_All <- quantile(PercLoss_Regul_All, 0.025, na.rm = T)

segments(y0 = UppCI_Regul_SHP, x0 = treatment * 1, 
         y1 = LowCI_Regul_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_Regul_LHP, x0 = treatment * 2, 
         y1 = LowCI_Regul_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Regul_All, x0 = treatment * 3, 
         y1 = LowCI_Regul_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(PercLoss_Regul_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_Regul_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_Regul_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.5, line = 1, mgp = c(3, 1.7, 0))

mtext("Regulated basins", side = 3, cex = 1.6, line = 1, at = treatment * 2 - 0.07)
mtext("C", side = 3, cex = 2.6, line = 2.9, at = 0.2)


dev.off()









######################################### OLD PLOT - Basins occupied by both ###################################
##### Plot % change
tiff(filename = "Figure2.tiff", height = 1656, width = 3200, res = 300, compression = c("lzw"))

mat <- matrix(c(rep(1, times = 2), rep(2, times = 2), rep(3, times = 2)), ncol = 3, nrow = 2, byrow = F)
layout(mat = mat, widths = c(4, 4, 4), heights = c(rep(0.5, times = 3)))

# Par
par (oma = c(5.5, 6, 6.5, 4.0), mar = c(0.5, 0, 0, 0), bty = "n")

par(bty = "n")

### Plot data for all basins

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI)[1])

plot(x = treatment, y = PercLoss_All, ylim = c(-100, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = PercLoss_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4.0, pch = 21)
points(x = treatment * 2, y = PercLoss_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4.0, pch = 21)
points(x = treatment * 3, y = PercLoss_All, bg = "#9400D305", col = c("#42424220"), cex = 4.0, pch = 21)

# Plot 95% CI
UppCI_SHP <- quantile(PercLoss_SHP, 0.975, na.rm = T)
LowCI_SHP <- quantile(PercLoss_SHP, 0.025, na.rm = T)
UppCI_LHP <- quantile(PercLoss_LHP, 0.975, na.rm = T)
LowCI_LHP <- quantile(PercLoss_LHP, 0.025, na.rm = T)
UppCI_All <- quantile(PercLoss_All, 0.975, na.rm = T)
LowCI_All <- quantile(PercLoss_All, 0.025, na.rm = T)

segments(y0 = UppCI_SHP, x0 = treatment * 1, 
         y1 = LowCI_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_LHP, x0 = treatment * 2, 
         y1 = LowCI_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_All, x0 = treatment * 3, 
         y1 = LowCI_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(PercLoss_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 2, at = c(0, -20, -40, -60, -80, -100), cex.axis = 2.5, line = - 3)
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.5, line = 1, mgp = c(3, 1.7, 0))
mtext("% change in river connectivity", side = 2, cex = 2.2, line = 1.7)
#mtext("Predicted % change in river connectivity (DCI)", side = 2, cex = 2.0, line = 1.7)

mtext("All basins", side = 3, cex = 2.0, line = 1)
mtext("A", side = 3, cex = 2.6, line = 2.9, at = 0.2)


### Plot data just for hydropower-free basin
# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Free)[1])

plot(x = treatment, y = PercLoss_Free_All, ylim = c(-100, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = PercLoss_Free_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 2, y = PercLoss_Free_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 3, y = PercLoss_Free_All, bg = "#9400D305", col = c("#42424220"), cex = 4, pch = 21)

# Plot 95% CI
UppCI_Free_SHP <- quantile(PercLoss_Free_SHP, 0.975, na.rm = T)
LowCI_Free_SHP <- quantile(PercLoss_Free_SHP, 0.025, na.rm = T)
UppCI_Free_LHP <- quantile(PercLoss_Free_LHP, 0.975, na.rm = T)
LowCI_Free_LHP <- quantile(PercLoss_Free_LHP, 0.025, na.rm = T)
UppCI_Free_All <- quantile(PercLoss_Free_All, 0.975, na.rm = T)
LowCI_Free_All <- quantile(PercLoss_Free_All, 0.025, na.rm = T)

segments(y0 = UppCI_Free_SHP, x0 = treatment * 1, 
         y1 = LowCI_Free_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_Free_LHP, x0 = treatment * 2, 
         y1 = LowCI_Free_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Free_All, x0 = treatment * 3, 
         y1 = LowCI_Free_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(PercLoss_Free_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_Free_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_Free_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.5, line = 1, mgp = c(3, 1.7, 0))

mtext("Hydropower-free basins", side = 3, cex = 1.6, line = 1, at = treatment * 3 - 0.42)
mtext("B", side = 3, cex = 2.6, line = 2.9, at = 0.2)


### Plot data just for basins with both SHPs and LHPs
# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Both)[1])

plot(x = treatment, y = PercLoss_Both_All, ylim = c(-100, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = PercLoss_Both_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 2, y = PercLoss_Both_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 3, y = PercLoss_Both_All, bg = "#9400D305", col = c("#42424220"), cex = 4, pch = 21)

# Plot 95% CI
UppCI_Both_SHP <- quantile(PercLoss_Both_SHP, 0.975, na.rm = T)
LowCI_Both_SHP <- quantile(PercLoss_Both_SHP, 0.025, na.rm = T)
UppCI_Both_LHP <- quantile(PercLoss_Both_LHP, 0.975, na.rm = T)
LowCI_Both_LHP <- quantile(PercLoss_Both_LHP, 0.025, na.rm = T)
UppCI_Both_All <- quantile(PercLoss_Both_All, 0.975, na.rm = T)
LowCI_Both_All <- quantile(PercLoss_Both_All, 0.025, na.rm = T)

segments(y0 = UppCI_Both_SHP, x0 = treatment * 1, 
         y1 = LowCI_Both_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_Both_LHP, x0 = treatment * 2, 
         y1 = LowCI_Both_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Both_All, x0 = treatment * 3, 
         y1 = LowCI_Both_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(PercLoss_Both_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(PercLoss_Both_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(PercLoss_Both_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.5, line = 1, mgp = c(3, 1.7, 0))

mtext("Basins with both", side = 3, cex = 1.6, line = 1, at = treatment * 2 - 0.07)
mtext("C", side = 3, cex = 2.6, line = 2.9, at = 0.2)


dev.off()


######################################################################################################
########################### Claculate some statistics for % change ###################################
######################################################################################################

## Absolute change in overall DCI
mean(100 - resuDCI$All_curr)
sd(100 - resuDCI$All_curr)

## Absolute change in SHP DCI
mean(100 - resuDCI$SHP_curr)
sd(100 - resuDCI$SHP_curr)

## Absolute change in LHP DCI
mean(100 - resuDCI$LHP_curr)
sd(100 - resuDCI$LHP_curr)

## Proportion of the current scenario between SHPs and LHPs
round(mean(100 - resuDCI$SHP_curr)/mean(100 - resuDCI$LHP_curr), digits = 0)



## Overall future % change values
mean(PercLoss_All)
sd(PercLoss_All)
min(PercLoss_All)
max(PercLoss_All)

## SHPs future % change values
mean(PercLoss_SHP)
sd(PercLoss_SHP)
min(PercLoss_SHP)
max(PercLoss_SHP)

## SHPs future % change values
mean(PercLoss_LHP)
sd(PercLoss_LHP)
min(PercLoss_LHP)
max(PercLoss_LHP)



## Free-flowing All future % change values
mean(PercLoss_Free_All)
sd(PercLoss_Free_All)
min(PercLoss_Free_All)
max(PercLoss_Free_All)

## Free-flowing SHPs future % change values
mean(PercLoss_Free_SHP)
sd(PercLoss_Free_SHP)
min(PercLoss_Free_SHP)
max(PercLoss_Free_SHP)

## Free-flowing LHPs future % change values
mean(PercLoss_Free_LHP)
sd(PercLoss_Free_LHP)
min(PercLoss_Free_LHP)
max(PercLoss_Free_LHP)



## Basins with both SHPs future % change values
mean(PercLoss_Both_All)
sd(PercLoss_Both_All)
min(PercLoss_Both_All)
max(PercLoss_Both_All)

## Free-flowing SHPs future % change values
mean(PercLoss_Both_SHP)
sd(PercLoss_Both_SHP)
min(PercLoss_Both_SHP)
max(PercLoss_Both_SHP)

## Free-flowing SHPs future % change values
mean(PercLoss_Both_LHP)
sd(PercLoss_Both_LHP)
min(PercLoss_Both_LHP)
max(PercLoss_Both_LHP)






####################################################################################################################################
####################################################################################################################################
####################################### DETERMINE AVERAGE LINEAR DISTANCES OF RIVER BASINS #########################################
######################################  AND NUMBER OF POSSIBLE PORTIFOLIOS OF DAMS  ################################################
######################################        (INFORM METHOD SECTION)                   ############################################

## Create an empty matrix
LinearDistKm<- rep(NA, times = length(basinList))
NPortifolios<- rep(NA, times = length(basinList))

## Loop over basins
for (j in 1: length(basinList)){
  
  ## filter attributes of the basin j         
  BasinX <- NetworkBRAZIL[NetworkBRAZIL$HYBAS_ID08ext == basinList[j], ]
  DamX <- DamAttributes[DamAttributes$HYBAS_ID08ext == basinList[j], ]

  ## Linear distance
  LinearDist <- sum(BasinX$Shape_Length)
  LinearDistKm[j] <- LinearDist/1000
  
  ## Possible dam portifolios per basin
  NFutDams <- sum(DamX$ESTAGIO_1 == "Planned")
  NPortifolios[j] <- 2^NFutDams
  
}

mean(LinearDistKm)
sum(NPortifolios)

##################################################################################
############################ Plot difference #####################################

tiff(filename = "Figure2b.tiff", height = 1396, width = 3200, res = 300, compression = c("lzw"))

# mat <- matrix(c(rep(1, times = 2), rep(2, times = 2)), ncol = 2, nrow = 2, byrow = F)
# layout(mat = mat, widths = c(4, 4), heights = c(rep(0.5, times = 2)))

mat <- matrix(c(rep(1, times = 2), rep(2, times = 2), rep(3, times = 2)), ncol = 3, nrow = 2, byrow = F)
layout(mat = mat, widths = c(4, 4, 4), heights = c(rep(0.5, times = 3)))


# Par
par (oma = c(5.5, 6, 4, 4.0), mar = c(0.5, 0, 0, 0), bty = "n")

par(bty = "n")

### Plot data for all basins

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI)[1])

plot(x = treatment, y = DCILoss_All, ylim = c(-90, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = DCILoss_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4.0, pch = 21)
points(x = treatment * 2, y = DCILoss_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4.0, pch = 21)
points(x = treatment * 3, y = DCILoss_All, bg = "#9400D305", col = c("#42424220"), cex = 4.0, pch = 21)

# Plot 95% CI
UppCI_SHP <- quantile(DCILoss_SHP, 0.975, na.rm = T)
LowCI_SHP <- quantile(DCILoss_SHP, 0.025, na.rm = T)
UppCI_LHP <- quantile(DCILoss_LHP, 0.975, na.rm = T)
LowCI_LHP <- quantile(DCILoss_LHP, 0.025, na.rm = T)
UppCI_All <- quantile(DCILoss_All, 0.975, na.rm = T)
LowCI_All <- quantile(DCILoss_All, 0.025, na.rm = T)

segments(y0 = UppCI_SHP, x0 = treatment * 1, 
         y1 = LowCI_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_LHP, x0 = treatment * 2, 
         y1 = LowCI_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_All, x0 = treatment * 3, 
         y1 = LowCI_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(DCILoss_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(DCILoss_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(DCILoss_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 2, at = c(0, -20, -40, -60, -80), cex.axis = 2.3, line = - 3)
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.3, line = 1, mgp = c(3, 1.7, 0))
mtext("DCI loss", side = 2, cex = 2.5, line = 1.5)


mtext("All basins", side = 3, cex = 2.0, line = 1)



### Plot data just for basins with both SHPs and LHPs

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Both)[1])

plot(x = treatment, y = DCILoss_Both_All, ylim = c(-90, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = DCILoss_Both_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 2, y = DCILoss_Both_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 3, y = DCILoss_Both_All, bg = "#9400D305", col = c("#42424220"), cex = 4, pch = 21)

# Plot 95% CI
UppCI_Both_SHP <- quantile(DCILoss_Both_SHP, 0.975, na.rm = T)
LowCI_Both_SHP <- quantile(DCILoss_Both_SHP, 0.025, na.rm = T)
UppCI_Both_LHP <- quantile(DCILoss_Both_LHP, 0.975, na.rm = T)
LowCI_Both_LHP <- quantile(DCILoss_Both_LHP, 0.025, na.rm = T)
UppCI_Both_All <- quantile(DCILoss_Both_All, 0.975, na.rm = T)
LowCI_Both_All <- quantile(DCILoss_Both_All, 0.025, na.rm = T)

segments(y0 = UppCI_Both_SHP, x0 = treatment * 1, 
         y1 = LowCI_Both_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_Both_LHP, x0 = treatment * 2, 
         y1 = LowCI_Both_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Both_All, x0 = treatment * 3, 
         y1 = LowCI_Both_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(DCILoss_Both_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(DCILoss_Both_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(DCILoss_Both_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.3, line = 1, mgp = c(3, 1.7, 0))

mtext("With SHPs and LHPs", side = 3, cex = 1.6, line = 1, at = treatment * 2 - 0.07)



### Plot data just for hydropower-free basin

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Free)[1])

plot(x = treatment, y = DCILoss_Free_All, ylim = c(-90, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment * 1, y = DCILoss_Free_SHP, bg = "#8B000005", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 2, y = DCILoss_Free_LHP, bg = "#4876FF05", col = c("#42424220"), cex = 4, pch = 21)
points(x = treatment * 3, y = DCILoss_Free_All, bg = "#9400D305", col = c("#42424220"), cex = 4, pch = 21)

# Plot 95% CI
UppCI_Free_SHP <- quantile(DCILoss_Free_SHP, 0.975, na.rm = T)
LowCI_Free_SHP <- quantile(DCILoss_Free_SHP, 0.025, na.rm = T)
UppCI_Free_LHP <- quantile(DCILoss_Free_LHP, 0.975, na.rm = T)
LowCI_Free_LHP <- quantile(DCILoss_Free_LHP, 0.025, na.rm = T)
UppCI_Free_All <- quantile(DCILoss_Free_All, 0.975, na.rm = T)
LowCI_Free_All <- quantile(DCILoss_Free_All, 0.025, na.rm = T)

segments(y0 = UppCI_Free_SHP, x0 = treatment * 1, 
         y1 = LowCI_Free_SHP, x1 = treatment * 1, lwd = 4.3, col = "black")
segments(y0 = UppCI_Free_LHP, x0 = treatment * 2, 
         y1 = LowCI_Free_LHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Free_All, x0 = treatment * 3, 
         y1 = LowCI_Free_All, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot means
points(x = treatment[1] * 1, y = mean(DCILoss_Free_SHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 2, y = mean(DCILoss_Free_LHP), pch = "-", col = "black", cex = 13.0)
points(x = treatment[1] * 3, y = mean(DCILoss_Free_All), pch = "-", col = "black", cex = 13.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("SHP", "LHP", "All"),
     cex.axis = 2.3, line = 1, mgp = c(3, 1.7, 0))

mtext("Currently free of hydropower", side = 3, cex = 1.6, line = 1, at = treatment * 3 - 0.42)


dev.off()


