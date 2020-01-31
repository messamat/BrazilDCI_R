#############         Script 3 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ##############################
#############              Manuscript in preparation - PNAS                       ##############################

################################################################################################################
###########       FRAGMENTATION BY PLANNED PROJECTS - BRAZIL      #############################################
###########                         (INDIVIDUALLY)                ############################################
###########       RUN DCI for all possible future portifolios     ###########################################
###########            (Assign it back to each planned dam)       ##########################################
###########################################################################################################
#If run into memoery errors: see https://stackoverflow.com/questions/51248293/error-vector-memory-exhausted-limit-reached-r-3-5-0-macos


## Choose what type of DCI function to run (DCIp or DCIi)
DCIfunc <- DCIp_opti
# DCIfunc <- DCId_opti

## Packages
require(tictoc)
require(plyr)
require(bigstatsr)
require(parallel)
require(doParallel)
require(Rcpp)
require(rccpcomb)
require(devtools)
devtools::source_gist("https://gist.github.com/r2evans/e5531cbab8cf421d14ed", filename = "lazyExpandGrid.R") #Get expand.grid version that won't run out of memory when > 25 dams
require(ggplot2)

# # Import network and dams dataset (Mathis folder structure)
# rootdir <- find_root(has_dir("src"))
# resdir <- file.path(rootdir, "results")
# dcigdb <- file.path(resdir, 'dci.gdb')

## Import network and dams dataset (alternative)
rootdir <- find_root(has_dir("PythonOutputs"))
datadir <- file.path(rootdir, "PythonOutputs")
dcigdb <- file.path(datadir, 'dci.gdb')
NetworkBRAZILCrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='networkattributes'))
DamAttributesCrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='damattributes'))

#Source code just in case
sourceCpp(file.path(rootdir, 'src/BrazilDCI_R/combi2inds_rcpp.cpp'))

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

## Final Dataframes to run DCI analyses
DamAttributes <- DamAttributesCrude
NetworkBRAZIL <- NetworkBRAZIL                     

# Organize the matrix based on type and stage of each dam
DamAttributes[, ESTAGIO_1 := ifelse(ESTAGIO_1 != "Operação", 'Planned', 'Operation')]
DamAttributes[, Tipo_1 := ifelse(Tipo_1 != "UHE", "SHP", "LHP")]

## Fix mistakes in the dataset
DamAttributes$Tipo_1[which(DamAttributes$Tipo_1 == "SHP" & DamAttributes$POT_KW > 30000)] <- "LHP"    #UHE Buritizal(Das Mortes), UHE Resplendor(Doce), UHE Pouso Alto (Sucuri?)
DamAttributes$Tipo_1[which(DamAttributes$Tipo_1 == "LHP" & DamAttributes$POT_KW < 30000 &
                             DamAttributes$AREA_NA_MA < 13.0 & DamAttributes$ESTAGIO_1 == "Planned")] <- "SHP"   #Keep old dams as UHEs and new ones as SHPs

setnames(DamAttributes, 'HYBAS_ID08ext', 'DAMBAS_ID08ext')

#Estimate number of iterations required to process all dams
permuttable <- setorder(DamAttributes[ESTAGIO_1 == "Planned", 
                                      list(ndams = .N, 
                                           npermut =  do.call(lazyExpandGrid, rep(list(c(0.1, 1)), .N))$n), 
                                      by=DAMBAS_ID08ext], 'ndams') %>%
  .[,list(.N, sum=sum(npermut)), by=ndams] %>%
  .[, list(ndams = ndams, nbasins=N, cumpermut = cumsum(sum))]


ggplot(permuttable, aes(x=ndams, y=cumpermut)) + 
  geom_bar(stat='identity') +
  scale_y_log10()

## Create a vector with unique basin IDs
basinList <- as.character(unique(DamAttributes$DAMBAS_ID08ext))

# #Remove the basins that will have more than X future dams if needed
# X <- 10
# basinList <- DamAttributes[ESTAGIO_1 == "Planned", .N, by=DAMBAS_ID08ext][N<X, DAMBAS_ID08ext]

#### Parameters to run permutations-based DCI ####
#minimum number of dams in basin at which the function starts sampling permutation scenarios (for speed/memory sake)
#22 could be handled with 16GB of RAM
minlazy <- 22 
#number of seed scenarios to sample from lazy permutation. The actual number of scenario pairs that will actually 
#be evluated will often be at least >5-10 times more than that number
nsamples <- 10000 

############## LAUNCH ANALYSIS ####################
tic("total")
#Initialize parallel analysis
cl <- parallel::makeCluster(bigstatsr::nb_cores()) #make cluster based on recommended number of cores
on.exit(stopCluster(cl))
doParallel::registerDoParallel(cl)

DCIscens <- foreach(j=basinList, ## Loop over basins
                    .packages = c("data.table", "Rcpp", "rccpcomb", "magrittr", 'plyr','tcltk'), 
                    .noexport = c('combi2inds')) %dopar% {
                      
                      ## filter attributes of the basin j         
                      NetX <- NetworkBRAZIL[HYBAS_ID08ext == j, ]
                      DamX <- DamAttributes[DAMBAS_ID08ext == j, ]       
                      DAMIDvec <- as.character(DamX[ESTAGIO_1 == "Planned", DAMID])
                      
                      #Create progress bar for long processes
                      if (length(DAMIDvec) > 15) {
                      mypb <- tkProgressBar(title = "R progress bar", label = "",
                                            min = 0, max = 1, initial = 0, width = 300)
                      setTkProgressBar(mypb, 1, title = j, label = NULL)
                      }
                      
                      
                      #Get all permutation scenarios
                      if(length(DAMIDvec) < minlazy) {
                        #Create all combinations of future dams permeabilities
                        permut <- DamX[ESTAGIO_1 == "Planned", do.call(CJ, rep(list(c(0.1, 1)), .N))] #The equivalent of expand.grid but faster
                        colnames(permut) <- DAMIDvec
                        
                        #Add all existing dams
                        permut[, as.character(DamX[ESTAGIO_1 == "Operation", as.character(DAMID)]) := 0.1]
                        
                      } else { #If > 22 dams, the number of possible permutations leads to memory errors
                        permutlazy <- DamX[ESTAGIO_1 == "Planned", do.call(lazyExpandGrid, rep(list(c(0.1, 1)), .N))] #Get lazy permutation
                        permutsample <- as.data.table(permutlazy$nextItems(sample(check$n, size=nsamples))) #Sample 10000 permutations
                        colnames(permutsample) <- DAMIDvec
                        #For each dam, select all scenarios where dam is absent and create paired scenario where all other dams are unchanged, but dam is present
                        permut <- ldply(DAMIDvec, function(dam) {
                          rbind(permutsample[get(dam) == 1,],permutsample[get(dam) == 1,]) %>%
                            .[seq((.N/2+1),(.N)), (dam) := 0.1]
                        }) %>% 
                          #Remove duplicate scenarios
                          unique
                        
                        #Add all existing dams
                        permut[, as.character(DamX[ESTAGIO_1 == "Operation", as.character(DAMID)]) := 0.1]
                      }

                      #Melt and assign all combinations to a given basin-scenario ID to run with data.table
                      scenarios_reg <- melt(cbind(permut, 
                                                  DAMBAS_ID08ext_scenario = paste(j, seq_len(nrow(permut)), sep='_')), 
                                            id.var="DAMBAS_ID08ext_scenario") %>%
                        setnames('variable', 'DAMID') %>%
                        merge(DamX[, .(DAMID, DownSeg, UpSeg)], by='DAMID', allow.cartesian=F)
                    
                      #Run DCI for each scenario
                      DCIX<- scenarios_reg[, 
                                           list(DCI = DCIp_opti(
                                             d3=NetX[, list(id=as.character(SEGID),
                                                        l=Shape_Length)],
                                             d2=.SD[, list(id1 = DownSeg,
                                                        id2 = UpSeg,
                                                        pass = value)],
                                             print = F)),
                                           by=.(DAMBAS_ID08ext_scenario)]

                      #For each dam, get DCI diff for all scenarios with and without it (could replace ldply with a data.table solution -TBD)
                      #This just orders rows for each dam's permeability value in the scenario grid and then divides the output array in two (by rows)
                      #The first half of the rows have passability = 0.1 for the ordered dam, the second half has pasability =1
                      #By computing the row-wise difference between the two halves, it gives us the difference for all scenarios at once
                      DCIXpermut <- cbind(DCIX, permut)
                      
                      DCIdiffstats <- ldply(1:length(DAMIDvec), function(i) {
                        colorder <- c(i, (1:length(DAMIDvec))[-i])
                        setorderv(DCIXpermut, cols = DAMIDvec[colorder], na.last=FALSE)
                        DCIdiff <- DCIXpermut[(1+nrow(DCIXpermut)/2):(nrow(DCIXpermut)), DCI] - DCIXpermut[1:(nrow(DCIXpermut)/2), DCI] 
                        return(data.table(
                          DAMIDvec[i],
                          mean(DCIdiff), 
                          max(DCIdiff),
                          min(DCIdiff), 
                          quantile(DCIdiff, 0.975), 
                          quantile(DCIdiff, 0.025),
                          length(DCIdiff)*2))
                      }) %>% #Add dam attributes
                        cbind(DamX[ESTAGIO_1 == "Planned", .(Tipo_1, 
                                                             ESTAGIO_1, 
                                                             POT_KW/1000, 
                                                             NOME, 
                                                             DAMBAS_ID08ext)])
                  
                      colnames(DCIdiffstats) <- c("DAMID", "DCIMeanDiff", "DCIUppLim", "DCIDownLim", "DCIUpCI", "DCIDownCI", "Nscenarios",
                                                  "Type", "Situation", "Capacity","Name", "Basin")
                      
                      #Close progress bar
                      if (length(DAMIDvec) > 15) {close(mypb)}
                      
                      return(DCIdiffstats)
}
toc()
big_data <- as.data.table(do.call("rbind", DCIscens))

### Organize the data to plot
# Dam rank
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
