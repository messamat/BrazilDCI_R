#############         Script 2 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ##############################
#############              Manuscript in preparation - PNAS                       ##############################

###################################################################################################################
###########      FRAGMENTATION OVER TIME ANALYSIS - BRAZIL   #####################################################
###########           DCI - current dams over years          ###################################################
##############################################################################################################

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

# Organize dates and place back as years in the dataframe 
DamYr <- as.Date(DamAttributes$INIC_OPER, format = "%m/%d/%Y")
DamYear <- as.numeric(substring(DamYr, first = 1, last = 4))
DamAttributes$INIC_OPER <- DamYear

## Vector with the years of the analysis
YearVec <- seq(from = 1899, to = 2018)

## Create a matrix of basins and DCIs per year             
DCI_TIME_L8 <- matrix(NA, nrow = length(basinList), ncol = length(YearVec))
colnames(DCI_TIME_L8) <- YearVec
rownames(DCI_TIME_L8) <- as.numeric(substr(basinList, 1, nchar(basinList)-1))


## Possible Scenarios of time-period
Current <- "Operation"

## Possible Scenarios of hydro type
Small <- "SHP"
Large <- "LHP"
All <- c("SHP", "LHP")

## Create vectors with possible combination of scenarios
ScenarioSituation <- c("Current", "Current", "Current")
ScenarioType <- c("All", "Small", "Large")

tic("total")
## Loop over scenarios
for (x in 1: 3){
  
  ## Define types and situation according to the scenario x
  SituationScen <- get(ScenarioSituation[x])
  TypeScen <- get(ScenarioType[x])
  
  
  ## Loop over years
  for(y in 1:length(YearVec)){
    
    ## Create a matrix to be filled with the DCIs of a given year
    DCI_yearY_L8 <- rep(NA, times = length(basinList))
    
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
        Basin <- DamX$HYBAS_ID08
        Year <- DamX$INIC_OPER
        
        ## Put these data in a data frame with all necessary information about the edges
        EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin, Year)
        
        ## Create a passibility vector based on the number of segments of the basin
        ## Assign 100% permeability for all dams with a vector
        passVec <- rep(1, times = dim(EdgesData)[1])
        
        ## Find the position of dams that should be present according to the scenarios (Type and Sitiation = T)
        PresentDams <- which(EdgesData$Type %in% TypeScen * EdgesData$Situation %in% SituationScen == 1)
        
        ## Assign a permeability value of 0.1 for the dams of the scenario based on their position
        passVec[PresentDams] <- 0.1
        
        ## Replace the permeability of dams not constructed yet according to year y for 100%
        passVec[which(EdgesData$Year >= YearVec[y])] <- 1
        
        # List of edges and passability
        d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[, 2], pass = passVec)
        
        # attributes of nodes; node sizes
        d3 = data.frame(id = listSeg, l = listLengths)
        
        # Run DCI analyses for the basin
        DCI <- DCIfunc (d2, d3, print = F)
        #DCIi <- DCIi(d2, d3, print = F) #Downseg as the reference
        
      }
      
      # Fill out the yearly vector of DCIs for all basins
      DCI_yearY_L8[j] <- DCI
      
      ## Compute the DCI of the Basin and year in the matrix
      DCI_TIME_L8[, y] <- DCI_yearY_L8
      
      ## Print indexes to follow the calculations
      print(c(x, y, j))
      
    }
   
    ## Create a csv file with the output of DCI analysis for all the scenarios
    Timefilename <- paste("DCI_TIME_L8_DCIp_", ScenarioType[x], ".csv", sep ="")
    write.csv(DCI_TIME_L8, file = Timefilename)
     
    # ## Create a csv file with the output of DCI analysis for all the scenarios
    # Timefilename <- paste("DCI_TIME_L8_DCIi_", ScenarioType[x], ".csv", sep ="")
    # write.csv(DCI_TIME_L8, file = Timefilename)
    
  }
  
}
  
toc()


############################################################################################################
############################         PLOT FIGURE 1          ################################################
############################################################################################################

## read table with the outputs
TimeTable_All <- as.data.frame(read.csv("DCI_TIME_L8_DCIp_All.csv"))
TimeTable_SHP <- as.data.frame(read.csv("DCI_TIME_L8_DCIp_Small.csv"))
TimeTable_LHP <- as.data.frame(read.csv("DCI_TIME_L8_DCIp_Large.csv"))

# Read table for the future projections per basin
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIp.csv", header = T))

# ## read table with the outputs
# TimeTable_All <- as.data.frame(read.csv("DCI_TIME_L8_DCIi_All.csv"))
# TimeTable_SHP <- as.data.frame(read.csv("DCI_TIME_L8_DCIi_Small.csv"))
# TimeTable_LHP <- as.data.frame(read.csv("DCI_TIME_L8_DCIi_Large.csv"))
# 
# # Read table for the future projections per basin
# resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIi.csv", header = T))

## Create a vector of years
YearVec <- seq(from = 1899, to = 2018)

## Transform imported data.frames in a matrices
## SHPs
maxBasinSHP <- dim(TimeTable_SHP)[2]
MatrixSHP <- as.matrix(TimeTable_SHP[, 2: maxBasinSHP])
colnames(MatrixSHP) <- YearVec
rownames(MatrixSHP) <- TimeTable_SHP$X

## LHPs
maxBasinLHP <- dim(TimeTable_LHP)[2]
MatrixLHP <- as.matrix(TimeTable_LHP[, 2: maxBasinLHP])
colnames(MatrixLHP) <- YearVec
rownames(MatrixLHP) <- TimeTable_LHP$X

## All
maxBasinAll <- dim(TimeTable_All)[2]
MatrixAll <- as.matrix(TimeTable_All[, 2: maxBasinAll])
colnames(MatrixAll) <- YearVec
rownames(MatrixAll) <- TimeTable_All$X


## Plot figure
tiff(filename = "Figure1.tiff", height = 2396, width = 4400, res = 300, compression = c("lzw"))

par (oma = c(5, 5, 3, 4), mar = c(0.5, 0, 0, 0), bty = "n")
plot(YearVec, MatrixSHP[1,], ylim = c(0, 100), xlim = c(1899, 2043), type = "n", ylab = "", xlab = "", 
     xaxt = "n", yaxt = "n")

# Plot axis
axis(side = 2, at = c(0, 25, 50, 75, 100), cex.axis = 2.6, line = - 2)
axis(side = 1, at = c(1900, 1920, 1940, 1960, 1980, 2000, 2018), cex.axis = 2.6, line = -0.8, mgp = c(3, 1.45, 0))
axis(side = 1, at = c(2025, 2042), labels = c("",""), cex.axis = 2.3, line = -0.8, mgp = c(3, 1.25, 0))
#axis(side = 1, at = c(2031), labels = c("Future"), cex.axis = 2.3, line = -1)
mtext("Future", side = 1, line = 0.6, cex = 2.6, at = 2034)
mtext("River connectivity (DCI)", side = 2, cex = 3.2, line = 2.1)
mtext("Year", side = 1, cex = 3.2, line = 3.5)

## Plot a line for each basin (rows)
# SHPs
for (i in 1: dim(MatrixSHP)[1]) {
  lines(YearVec, MatrixSHP[i,], col = "#8B000015", lwd = 3.8, lty = 1)
}
  
# LHPs
for (i in 1: dim(MatrixLHP)[1]) {
  lines(YearVec, MatrixLHP[i,], col = "#4876FF15", lwd = 3.8, lty = 1)
}
 

# Plot the average lines SHP
averageDCISHP <- colSums(MatrixSHP)/dim(MatrixSHP)[1]
lines(YearVec, averageDCISHP, col = "black", lwd = 4.5, lty = 3)

# Plot the average lines LHP
averageDCILHP <- colSums(MatrixLHP)/dim(MatrixLHP)[1]
lines(YearVec, averageDCILHP, col = "black", lwd = 4.5, lty = 5)

# Plot the average lines All
averageDCIAll <- colSums(MatrixAll)/dim(MatrixAll)[1]
lines(YearVec, averageDCIAll, col = "black", lwd = 4.5, lty = 1)

## Plot future lines
# SHPs
for (i in 1: dim(resuDCI)[1]){
  lines(c(2025, 2042), c(resuDCI$SHP_curr[i], resuDCI$SHP_fut[i]), col = "#8B000015", lwd = 3.8, lty = 1)
}

# LHPs
for (i in 1: dim(resuDCI)[1]){
  lines(c(2025, 2042), c(resuDCI$LHP_curr[i], resuDCI$LHP_fut[i]), col = "#4876FF15", lwd = 3.8, lty = 1)
}

# Calculate average lines for current and future dams
averageFutAll <- c(mean(resuDCI$All_curr), mean(resuDCI$All_fut))
averageFutSHP <- c(mean(resuDCI$SHP_curr), mean(resuDCI$SHP_fut))
averageFutLHP <- c(mean(resuDCI$LHP_curr), mean(resuDCI$LHP_fut))

# Plot average lines
lines(c(2025, 2042), averageFutSHP, col = "black", lwd = 4.5, lty = 3)
lines(c(2025, 2042), averageFutLHP, col = "black", lwd = 4.5, lty = 5)
lines(c(2025, 2042), averageFutAll, col = "black", lwd = 4.5, lty = 1)

# Plot average points
points(2043, averageFutSHP[2], pch = 21, bg = "#8B0000", cex = 2.8)
points(2043, averageFutLHP[2], pch = 21, bg = "#4876FF", cex = 2.8)
points(2043, averageFutAll[2], pch = 21, bg = "#9400D3", cex = 2.8)

points(2019, averageFutSHP[1], pch = 21, bg = "#8B0000", cex = 2.8)
points(2019, averageFutLHP[1], pch = 21, bg = "#4876FF", cex = 2.8)
points(2019, averageFutAll[1], pch = 21, bg = "#9400D3", cex = 2.8)

mtext("SHPs", side = 1, cex = 2.0, line = -20.5, at = 2050)
mtext("LHPs", side = 1, cex = 2.0, line = -28, at = 2050)
mtext("All", side = 1, cex = 2.0, line = -18.5, at = 2048)

dev.off()


###############################################################################################
################## Compute some statistics of fragmentation over time #########################
###############################################################################################

## read table with the outputs
TimeTable_All <- as.data.frame(read.csv("DCI_TIME_L8_DCIp_All.csv"))
TimeTable_SHP <- as.data.frame(read.csv("DCI_TIME_L8_DCIp_Small.csv"))
TimeTable_LHP <- as.data.frame(read.csv("DCI_TIME_L8_DCIp_Large.csv"))

# Read table for the future projections per basin
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIp.csv", header = T))


## Calculate the average decrease in absolute values of DCI
# All current
initialStage <- rep(100, times = length(resuDCI$All_curr))
mean(initialStage - resuDCI$All_curr)

# All future
initialStage <- rep(100, times = length(resuDCI$All_fut))
mean(initialStage - resuDCI$All_fut)

# SHP current
initialStage <- rep(100, times = length(resuDCI$SHP_curr))
mean(initialStage - resuDCI$SHP_curr)

# SHP future
initialStage <- rep(100, times = length(resuDCI$SHP_fut))
mean(initialStage - resuDCI$SHP_fut)

# LHP current
initialStage <- rep(100, times = length(resuDCI$LHP_curr))
mean(initialStage - resuDCI$LHP_curr)

# LHP future
initialStage <- rep(100, times = length(resuDCI$LHP_fut))
mean(initialStage - resuDCI$LHP_fut)


## Calculate the standard deviation for decrease in DCI
# All current
initialStage <- rep(100, times = length(resuDCI$All_curr))
sd(initialStage - resuDCI$All_curr)

# All future
initialStage <- rep(100, times = length(resuDCI$All_fut))
sd(initialStage - resuDCI$All_fut)

# SHP current
initialStage <- rep(100, times = length(resuDCI$SHP_curr))
sd(initialStage - resuDCI$SHP_curr)

# SHP future
initialStage <- rep(100, times = length(resuDCI$SHP_fut))
sd(initialStage - resuDCI$SHP_fut)

# LHP current
initialStage <- rep(100, times = length(resuDCI$LHP_curr))
sd(initialStage - resuDCI$LHP_curr)

# LHP future
initialStage <- rep(100, times = length(resuDCI$LHP_fut))
sd(initialStage - resuDCI$LHP_fut)


## Calculate the range of the decrease
# All current
100 - max(resuDCI$All_curr)
100 - min(resuDCI$All_curr)

# All future
100 - max(resuDCI$All_fut)
100 - min(resuDCI$All_fut)

# SHP
100 - max(resuDCI$SHP_curr)
100 - min(resuDCI$SHP_curr)

# LHP
100 - max(resuDCI$LHP_curr)
100 - min(resuDCI$LHP_curr)


## Calculate rates of decrease per year (I did the first round in excel)

## All
YearMeans_All <- round(colMeans(TimeTable_All), digits = 1)
YearRates_All <- rep(NA, times = length(YearMeans_All) - 2)

for (i in 3: length(YearMeans_All)){
  YearRates_All[i-2] <- YearMeans_All[i] - YearMeans_All[i - 1]
}

YearRatesMat <- cbind(1900:2018, YearRates_All)

mean(YearRatesMat[61:70,2]) ##60s
mean(YearRatesMat[71:80,2]) ##70s
mean(YearRatesMat[81:90,2]) ##80s
mean(YearRatesMat[91:100,2]) ##90s
mean(YearRatesMat[101:110,2]) ##2000s
mean(YearRatesMat[111:119,2]) ##2010s

round(mean(YearRatesMat[61:100,2]), digits = 2) ## 1960 -1999
round(mean(YearRatesMat[101:119,2]), digits = 2) ##2000 - 2018



## SHPs
YearMeans_SHP <- round(colMeans(TimeTable_SHP), digits = 1)
YearRates_SHP <- rep(NA, times = length(YearMeans_SHP) - 2)

for (i in 3: length(YearMeans_SHP)){
  YearRates_SHP[i-2] <- YearMeans_SHP[i] - YearMeans_SHP[i - 1]
}

YearRatesMatSHP <- cbind(1900:2018, YearRates_SHP)

mean(YearRatesMatSHP[61:70,2]) ##60s
mean(YearRatesMatSHP[71:80,2]) ##70s
mean(YearRatesMatSHP[81:90,2]) ##80s
mean(YearRatesMatSHP[91:100,2]) ##90s
mean(YearRatesMatSHP[101:110,2]) ##2000s
mean(YearRatesMatSHP[111:119,2]) ##2010s

round(mean(YearRatesMatSHP[61:100,2]), digits = 2) ## 1960 -1999
round(mean(YearRatesMatSHP[101:119,2]), digits = 2) ##2000 - 2018
