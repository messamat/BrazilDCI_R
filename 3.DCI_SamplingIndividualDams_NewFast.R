#############         Script 3 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ##############################
#############              Manuscript in preparation - PNAS                       ##############################

################################################################################################################
###########       FRAGMENTATION BY PLANNED PROJECTS - BRAZIL      #############################################
###########                         (INDIVIDUALLY)                ############################################
###########       RUN DCI for all possible future portifolios     ###########################################
###########            (Assign it back to each planned dam)       ##########################################
###########################################################################################################

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


## Remove the basins that will have more than 20 future dams
## To get fast outputs as a test
PlanDamsData <- DamAttributes[DamAttributes$ESTAGIO_1 == "Planned",]
BasinDamCounts <- table(PlanDamsData$HYBAS_ID08)
RemoBasinNames <- names(BasinDamCounts[BasinDamCounts > 20])
basinList <- basinList[-which(substr(basinList, 1, nchar(basinList)-1) %in% RemoBasinNames)]


## Create a list to add the outputs of each dam
DamList = list()

tic("total")
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
    
    ## Create vectors of segids and their total lenght
    listSeg <- paste("seg", BasinX$SEGID, sep ="")
    listLengths <- BasinX$Shape_Length
    
    ## Compile the data on the edges (dam attributes)
    DowSeg <- paste("seg", DamX$DownSeg, sep ="")
    UpSeg <- paste("seg", DamX$UpSeg, sep ="")
    Type <- DamX$Tipo_1
    Situation <- DamX$ESTAGIO_1
    ID_number <- DamX$TARGET_FID
    ID_name <- DamX$NOME
    Basin <- DamX$HYBAS_ID08
    Capacity <- DamX$POT_KW/1000 #MW
    
    ## Put these data in a data frame with all necessary information about the edges
    EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin, Capacity)
    
    # Determine which dams are operating and each ones are planned
    OpperatingDams <- EdgesData[EdgesData$Situation == "Operation", ]
    PlannedDams <- EdgesData[EdgesData$Situation == "Planned", ]
    
    
    ## Construct a matrix of possible values of permeability just for the new dams
    # n is the number of all new possible dams
    n <- dim(PlannedDams)[1]
    
    # Generate all the permeability options just for the future dams
    l <- rep(list(c(0.1, 1)), n)
    OutcomePossib <- expand.grid(l)
    
    # Exclude the scenario with zero new dams (i.e. last row of OutcomePossib)
    #futurePossib <- OutcomePossib
    futurePossib <- OutcomePossib[-dim(OutcomePossib)[1],]
    #futurePossibCurrDams <- OutcomePossib[dim(OutcomePossib)[1],]
    
    # ## Optional - just to gain time
    # ## Subselect possible scenarios for basins with a huge number of dams (computational reasons)
    # if (n > 10){
    #   
    #   futurePossib <- futurePossib[sample(nrow(futurePossib), size = 1000, replace = F),]
    #   
    # }
    
    ## Create a vector of possible DCIs per basin
    if(n > 1){
      DCIpossib <- rep(NA, times = dim(futurePossib)[1])
    }
    
    if(n == 1){
      DCIpossib <- NA
    }
    
    
    ## Create a matrix of dams to be filled out with average, upper and lower DCIs
    priorizationDCI <- matrix(NA, nrow = dim(OutcomePossib)[2], ncol = 11)
    colnames(priorizationDCI) <- c("DCIMeanDiff", "DCIUppLim", "DCIDownLim", "DCIUpCI", "DCIDownCI", "Type", "Situation",
                                   "Capacity", "ID", "Name", "Basin")
    rownames(priorizationDCI) <- names(futurePossib)
    
    priorizationDCI <- data.frame(priorizationDCI)
    
    
    ## Run the possible scenarios for each future dam when the n of future dams is bigger than 1
    if (n > 1){
      
      ## Loop for over all the future dams
      for (s in 1: dim(priorizationDCI)[1]){
        
        ## Select just the future scenarios where the future dam s is present
        DamPresentScen <- futurePossib[futurePossib[,s] == 0.1, ]
        
        ## Create a scenario in which just the dam s is absent to calculate the differential DCI
        DamAbsentScen <- DamPresentScen
        DamAbsentScen[,s] <- 1
        
        ## Create a vector to be filled out with all DCI estimates for each dam
        DCIpossib <- rep(NA, times = dim(DamPresentScen)[1])
        
        
        ## Loop over possible scenarios for that given dam to calculate DCIs
        for (z in 1: dim(DamPresentScen)[1]){
          
          passVecCur <- rep(0.1, times = dim(OpperatingDams)[1])
          passVecFut <- as.numeric(DamPresentScen[z,])
          passVecFutWithout <- as.numeric(DamAbsentScen[z,])
          
          PermeabilityWith <- c(passVecCur, passVecFut)
          PermeabilityWithout <- c(passVecCur, passVecFutWithout)
          
          # Combine all Edge data to run DCI
          mergedEdges <- rbind(OpperatingDams, PlannedDams)
          RunDCI <- cbind(mergedEdges, PermeabilityWith, PermeabilityWithout)
          
          # attributes of edges; links between nodes
          d2 = data.frame (id1 = RunDCI$DowSeg, id2 = RunDCI$UpSeg, pass = RunDCI$PermeabilityWith)
          d2_without = data.frame (id1 = RunDCI$DowSeg, id2 = RunDCI$UpSeg, pass = RunDCI$PermeabilityWithout)
          
          # attributes of nodes; node sizes
          d3 = data.frame(id = listSeg, l = listLengths)
          
          # Run DCI analysis for the basin
          DCI <- DCIfunc (d2, d3, print = F)
          DCI_without <- DCIfunc (d2_without, d3, print = F)
          
          # Fill out with the contribution of the dam s to scenario z (transformed in negative)
          DCIpossib[z] <- DCI - DCI_without
          
        }
        
        ## Calculate average, lower and upper DCI for each dam
        priorizationDCI[s, 1] <- mean(DCIpossib)
        priorizationDCI[s, 2] <- max(DCIpossib)
        priorizationDCI[s, 3] <- min(DCIpossib)
        priorizationDCI[s, 4] <- quantile(DCIpossib, 0.975)
        priorizationDCI[s, 5] <- quantile(DCIpossib, 0.025)
        priorizationDCI[s, 6] <- as.character(PlannedDams$Type[s])
        priorizationDCI[s, 7] <- as.character(PlannedDams$Situation[s])
        priorizationDCI[s, 8] <- PlannedDams$Capacity[s]
        priorizationDCI[s, 9] <- PlannedDams$ID_number[s]
        priorizationDCI[s, 10] <- as.character(PlannedDams$ID_name[s])
        priorizationDCI[s, 11] <- PlannedDams$Basin[s]
        
      }
      
    }
    
    ## Run the possible scenarios for each future dam when the n of future dams is equal to 1
    if (n == 1){
      
      ## Create a vector to be filled out with all DCI estimates for each dam
      DCIpossib <- NA
      
      passVecCur <- rep(0.1, times = dim(OpperatingDams)[1])
      passVecFut <- 0.1
      passVecFutWithout <- 1
      
      PermeabilityWith <- c(passVecCur, passVecFut)
      PermeabilityWithout <- c(passVecCur, passVecFutWithout)
      
      # Combine all Edge data to run DCI
      mergedEdges <- rbind(OpperatingDams, PlannedDams)
      RunDCI <- cbind(mergedEdges, PermeabilityWith, PermeabilityWithout)
      
      # attributes of edges; links between nodes
      d2 = data.frame (id1 = RunDCI$DowSeg, id2 = RunDCI$UpSeg, pass = RunDCI$PermeabilityWith)
      d2_without = data.frame (id1 = RunDCI$DowSeg, id2 = RunDCI$UpSeg, pass = RunDCI$PermeabilityWithout)
      
      # attributes of nodes; node sizes
      d3 = data.frame(id = listSeg, l = listLengths)
      
      # Run DCI analysis for the basin
      DCI <- DCIfunc (d2, d3, print = F)
      DCI_without <- DCIfunc (d2_without, d3, print = F)
      
      # Fill out with the contribution of the dam s to scenario z
      DCIpossib <- DCI - DCI_without
      
      # When there is just one future dam, s = 1
      s <- 1
      
      ## Calculate average, lower and upper DCI for each dam
      priorizationDCI[s, 1] <- mean(DCIpossib)
      priorizationDCI[s, 2] <- max(DCIpossib)
      priorizationDCI[s, 3] <- min(DCIpossib)
      priorizationDCI[s, 4] <- quantile(DCIpossib, 0.975)
      priorizationDCI[s, 5] <- quantile(DCIpossib, 0.025)
      priorizationDCI[s, 6] <- as.character(PlannedDams$Type[s])
      priorizationDCI[s, 7] <- as.character(PlannedDams$Situation[s])
      priorizationDCI[s, 8] <- PlannedDams$Capacity[s]
      priorizationDCI[s, 9] <- PlannedDams$ID_number[s]
      priorizationDCI[s, 10] <- as.character(PlannedDams$ID_name[s])
      priorizationDCI[s, 11] <- PlannedDams$Basin[s]
      
    }
      
    ## Add results to a list
    DamList[[j]] <- priorizationDCI
    
    ## Print indexes to follow the calculations
    print(c(j, s))
    
  }
  
}

toc()  


### Organize the data to plot
# Dam rank
big_data = do.call(rbind, DamList)
RankedDams <- big_data[order(big_data[,1]),]
DamRank <- seq(from = 1, to = dim(RankedDams)[1])
PriorityAnalysis <- cbind(RankedDams, DamRank)

## Export a csv of the prioritization analysis
# write.csv(PriorityAnalysis, file = "IndividualDam_DCIp.csv")
# write.csv(PriorityAnalysis, file = "IndividualDam_DCIi.csv")

###############################################################################################################
####################       PLOT DCI Loss individual future dams VS Capacity       #############################
###############################################################################################################

# ## Import the csv with the prioritization results
# PriorityAnalysis <- read.csv("IndividualDam_DCIp.csv", header = T)

## Clumsy version - Multiple runs (2,215 dams)
## Import the csv with the prioritization results
PriorityAnalysis1 <- read.csv("IndividualDam_DCIp_Forwards1_1038.csv", header = T)
PriorityAnalysis2 <- read.csv("IndividualDam_DCIp_Forwards1040_1048.csv", header = T)
PriorityAnalysis3 <- read.csv("IndividualDam_DCIp_Backwards1063_1049.csv", header = T)
PriorityAnalysis4 <- read.csv("IndividualDam_DCIp_Backwards1215_1065.csv", header = T)

## Merge all of them in one
PriorityAnalysisFull <- rbind(PriorityAnalysis1, PriorityAnalysis2, PriorityAnalysis3, PriorityAnalysis4)
PriorityAnalysisFull$DamRank <- 1:dim(PriorityAnalysisFull)[1]
PriorityAnalysis <- PriorityAnalysisFull

## Plot 5 (DCI loss Vs Capacity)
tiff(filename = "Figure5.tiff", height = 2396, width = 3700, res = 300, compression = c("lzw"))
par(oma = c(6, 9, 2, 0.5), mar = c(0.5, 0, 0, 0), bty = "n")


## Create a color vector to differentiate SHP and LHP
TypeFac <- factor(x = PriorityAnalysis$Type)
TypeColor <- factor(x = TypeFac, levels = levels(TypeFac), 
                    labels = c("#4876FF50", "#8B000050"))

## Plot figure #xaxt = "n", log(PriorityAnalysis$Capacity)
plot(PriorityAnalysis$Capacity, PriorityAnalysis$DCIMeanDiff, type = "n", ylim = c(-60, 0), xlim = c(0.1, 8000),
     ylab = "", xlab = "",  yaxt = "n", xaxt = "n", log = "x") 

options(scipen=5)
axis(side = 2, at = c(-60, -40, -20, 0), cex.axis = 2.6, line = 0)
axis(side = 1, at = c(1, 10, 30, 100, 500, 4000), cex.axis = 2.6, line = 0, mgp = c(3, 1.8, 0))
axis(side = 1, at = c(0.1, 1), labels = c(0.1, " "), cex.axis = 2.6, line = 0, mgp = c(3, 1.8, 0))

mtext("Change in basin-level", side = 2, cex = 3.2, line = 6.8)
mtext("river connectivity (DCI)", side = 2, cex = 3.2, line = 4.0)
mtext("Generation capacity (megawatts)", side = 1, cex = 3.2, line = 4.9)


## Plot range log(PriorityAnalysis$Capacity)
segments(y0 = PriorityAnalysis$DCIUppLim, x0 = PriorityAnalysis$Capacity,
         y1 = PriorityAnalysis$DCIDownLim, x1 = PriorityAnalysis$Capacity, lwd = 1.6, col = "#42424220")

# ## Plot 95% confidence intervals
# segments(y0 = PriorityAnalysis$DCIUppLim, x0 = log(PriorityAnalysis$Capacity), 
#          y1 = PriorityAnalysis$DCIDownCI, x1 = log(PriorityAnalysis$Capacity), lwd = 1.3, col = "#42424220")

## Plot average log(PriorityAnalysis$Capacity)
points(PriorityAnalysis$Capacity, PriorityAnalysis$DCIMeanDiff, pch = 21, col = "white", 
       bg = "white", cex = 2.5)
points(PriorityAnalysis$Capacity, PriorityAnalysis$DCIMeanDiff, pch = 21, col = "#42424220", 
       bg = as.vector(TypeColor), cex = 2.5)

## Plot legend
legend(x = 200, y = -47, pch = c(21), legend = c("Planned SHP", "Planned LHP"), cex = 2.5,
       pt.cex = 3.3, pt.bg = c("#8B000060","#4876FF60"), xpd = T, bty = "n")

dev.off()



################################################################################################################
#####################################  Some Statistics  ########################################################

##
PriorityAnalysis$Capacity
PriorityAnalysis$DCIMeanDiff

## Linear models effect on DCI Vs Capacity
plot(PriorityAnalysis$Capacity, PriorityAnalysis$DCIMeanDiff)

## All
overallLm <- lm(log(-PriorityAnalysis$DCIMeanDiff) ~ log(PriorityAnalysis$Capacity))
summary(overallLm)
hist(overallLm$resid, main="Histogram of Residuals",
     ylab="Residuals")
qqnorm(overallLm$resid)
qqline(overallLm$resid)
plot(log(PriorityAnalysis$Capacity), log(-PriorityAnalysis$DCIMeanDiff))
abline(overallLm)

cor.test(x=log(PriorityAnalysis$Capacity), y=log(-PriorityAnalysis$DCIMeanDiff), method = 'pearson')


## Just SHPs
SHPLm <- lm(log(-PriorityAnalysis$DCIMeanDiff[PriorityAnalysis$Type == "SHP"]) ~ log(PriorityAnalysis$Capacity[PriorityAnalysis$Type == "SHP"]))
summary(SHPLm)
hist(SHPLm$resid, main="Histogram of Residuals",
     ylab="Residuals")
qqnorm(SHPLm$resid)
qqline(SHPLm$resid)
plot(log(PriorityAnalysis$Capacity[PriorityAnalysis$Type == "SHP"]), log(-PriorityAnalysis$DCIMeanDiff[PriorityAnalysis$Type == "SHP"]))
abline(SHPLm)

cor.test(x=log(PriorityAnalysis$Capacity[PriorityAnalysis$Type == "SHP"]), y=log(-PriorityAnalysis$DCIMeanDiff[PriorityAnalysis$Type == "SHP"]), method = 'pearson')

## Just LHPs
LHPLm <- lm(log(-PriorityAnalysis$DCIMeanDiff[PriorityAnalysis$Type == "LHP"]) ~ log(PriorityAnalysis$Capacity[PriorityAnalysis$Type == "LHP"]))
summary(LHPLm)
hist(LHPLm$resid, main="Histogram of Residuals",
     ylab="Residuals")
qqnorm(LHPLm$resid)
qqline(LHPLm$resid)
plot(log(PriorityAnalysis$Capacity[PriorityAnalysis$Type == "LHP"]), log(-PriorityAnalysis$DCIMeanDiff[PriorityAnalysis$Type == "LHP"]))
abline(LHPLm)

cor.test(x=log(PriorityAnalysis$Capacity[PriorityAnalysis$Type == "LHP"]), y=log(-PriorityAnalysis$DCIMeanDiff[PriorityAnalysis$Type == "LHP"]), method = 'pearson')



## Create a csv of future dams rank basend on mean DCI

## Add columns with lat-long
## Created in ArcGIS two columns with Lat and long of the dams (Decimal degree, WGS 1984)
FullDamAttrubutes <- read.csv("DamAttributesCoordinates.txt", header = T)

Lat <- vector()
Long <- vector()

for (i in 1: dim(OrderDams)[1]){
  
  IDName_Position <- which(FullDamAttrubutes$TARGET_FID == PriorityAnalysis$ID[i] & FullDamAttrubutes$NOME %in% PriorityAnalysis$Name[i])
  DamLatLong <- FullDamAttrubutes[IDName_Position, 20:21]
  Lat <- c(Lat, FullDamAttrubutes$Lat[IDName_Position])
  Long <- c(Long, FullDamAttrubutes$Long[IDName_Position])
  
 }

## Organize order and headers
OrderDams <- data.frame(PriorityAnalysis$DamRank, PriorityAnalysis$Type, PriorityAnalysis$ID, PriorityAnalysis$Name, round(PriorityAnalysis$Capacity, digits = 1),
      round(PriorityAnalysis$DCIMeanDiff, digits = 1), round(PriorityAnalysis$DCIDownLim, digits = 1), round(PriorityAnalysis$DCIUppLim, digits = 1),
      Lat, Long)
colnames(OrderDams) <- c("Rank", "Type",  "DamID", "Name", "Capacity(MW)", "Mean effect on basin's DCI", "Lower limit", "Upper limit", "Latitude", "Longitude")



write.csv(OrderDams, file = "Supplement_FutureDamRank.csv")
