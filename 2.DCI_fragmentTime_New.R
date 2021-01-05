#############         Script 2 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ##############################
#############              Manuscript in preparation - PNAS                       ##############################

###################################################################################################################
###########      FRAGMENTATION OVER TIME ANALYSIS - BRAZIL   #####################################################
###########           DCI - current dams over years          ###################################################
##############################################################################################################
#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT


damfragment_time <- function(DamAttributes, NetworkBRAZIL, DCIfunc, DCIname) {

  ## Create a vector with unique basin IDs
  basinList <- as.character(unique(DamAttributes$DAMBAS_ID08ext))
  
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
    
    Timefilename <- file.path(resdir, 
                              paste0("DCI_TIME_L8_", DCIname, "_", 
                                     ScenarioType[x], ".csv"))
    
    if (!file.exists(Timefilename)) {
      ## Loop over years
      for(y in 1:length(YearVec)){
        
        ## Create a matrix to be filled with the DCIs of a given year
        DCI_yearY_L8 <- rep(NA, times = length(basinList))
        
        ## Loop over basins
        for (j in 1: length(basinList)){
          
          ## filter attributes of the basin j         
          BasinX <- NetworkBRAZIL[NetworkBRAZIL$HYBAS_ID08ext == basinList[j], ]
          DamX <- DamAttributes[DamAttributes$DAMBAS_ID08ext == basinList[j], ]
          
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
            ID_number <- DamX$DAMID
            ID_name <- DamX$NOME
            Basin <- DamX$DAMBAS_ID08
            Year <- DamX$INIC_OPER
            
            ## Put these data in a data frame with all necessary information about the edges
            EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, 
                                    ID_number, ID_name, Basin, Year)
            
            ## Create a passibility vector based on the number of segments of the basin
            ## Assign 100% permeability for all dams with a vector
            passVec <- rep(1, times = dim(EdgesData)[1])
            
            ## Find the position of dams that should be present according to the scenarios (Type and Sitiation = T)
            PresentDams <- which(EdgesData$Type %in% 
                                 TypeScen * EdgesData$Situation %in% 
                                 SituationScen == 1)
            
            ## Assign a permeability value of 0.1 for the dams of the scenario based on their position
            passVec[PresentDams] <- 0.1
            
            ## Replace the permeability of dams not constructed yet according to year y for 100%
            passVec[which(EdgesData$Year >= YearVec[y])] <- 1
            
            # List of edges and passability
            d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[, 2],
                             pass = passVec) %>%
              setDT
            
            # attributes of nodes; node sizes
            d3 = data.frame(id = listSeg, l = listLengths) %>%
              setDT
            
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
        write.csv(DCI_TIME_L8, file = Timefilename)
      }
    } else {
      print(paste0(Timefilename, ' already exists...'))
    }
  }
  
  toc()
  
  ##############################################################################
  ############################         PLOT FIGURE 1          ##################
  ##############################################################################
  ## read table with the outputs
  TimeTable_All <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_All.csv"))))
  TimeTable_SHP <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_Small.csv"))))
  TimeTable_LHP <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_Large.csv"))))
  
  # Read table for the future projections per basin
  resuDCI <- as.data.frame(read.csv(
    file.path(resdir,
              paste0("DCI_Brazil_L8_", DCIname, ".csv")), header = T))
  
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
  fig1func <- function() {
    par (oma = c(5, 5, 3, 4), mar = c(0.5, 0, 0, 0), bty = "n")
    plot(YearVec, MatrixSHP[1,], ylim = c(0, 100), xlim = c(1899, 2043),
         type = "n", ylab = "", xlab = "", 
         xaxt = "n", yaxt = "n")
    
    # Plot axis
    axis(side = 2, at = c(0, 25, 50, 75, 100), cex.axis = 2.6, line = - 2)
    axis(side = 1, at = c(1900, 1920, 1940, 1960, 1980, 2000, 2018), 
         cex.axis = 2.6, line = -0.8, mgp = c(3, 1.45, 0))
    axis(side = 1, at = c(2025, 2042), labels = c("",""), 
         cex.axis = 2.3, line = -0.8, mgp = c(3, 1.25, 0))
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
      lines(c(2025, 2042), c(resuDCI$SHP_curr[i], resuDCI$SHP_fut[i]),
            col = "#8B000015", lwd = 3.8, lty = 1)
    }
    
    # LHPs
    for (i in 1: dim(resuDCI)[1]){
      lines(c(2025, 2042), c(resuDCI$LHP_curr[i], resuDCI$LHP_fut[i]), 
            col = "#4876FF15", lwd = 3.8, lty = 1)
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
    
    
    text("SHPs", cex = 2.0, y=averageFutSHP[2]+1, x=2050, xpd=NA)
    text("LHPs", cex = 2.0, y=averageFutLHP[2]+1, x=2050, xpd=NA)
    text("All", cex = 2.0, y=averageFutAll[2]+1, x=2047, xpd=NA)
  }
  
  outfigname = file.path(figdir, paste0("Figure1_", DCIname, "_",
                                        format(Sys.Date(), '%Y%m%d')))
  jpeg(filename = paste0(outfigname, '.jpg'),
       height = 2396, width = 4400, res = 300)
  fig1func()
  dev.off()
  
  pdf(file = paste0(outfigname, '.pdf'),
      height = 2396/300, width = 4400/300)
  fig1func()
  dev.off()
}


damfragment_time(DamAttributes = DamAttributes,
                 NetworkBRAZIL = NetworkBRAZIL,
                 DCIfunc = DCIp_opti5,
                 DCIname = 'DCIp')


damfragment_time(DamAttributes = DamAttributes,
                 NetworkBRAZIL = NetworkBRAZIL,
                 DCIfunc = DCIi_opti,
                 DCIname = 'DCIi')


  ###############################################################################################
  ################## Compute some statistics of fragmentation over time #########################
  ###############################################################################################
  DCIname <- 'DCIp'
  
  ## read table with the outputs
  TimeTable_All <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_All.csv"))))
  TimeTable_SHP <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_Small.csv"))))
  TimeTable_LHP <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_Large.csv"))))
  
  # Read table for the future projections per basin
  resuDCI <- as.data.frame(read.csv(
    file.path(resdir,
              paste0("DCI_Brazil_L8_", DCIname, ".csv")), header = T)) %>%
    setDT
  
  ## Create a vector of years
  YearVec <- seq(from = 1899, to = 2018)
  
  #Number of future dams
  DamAttributes[DamAttributes$ESTAGIO_1 == 'Planned', .N]
  
  ## Calculate the average decrease in absolute values of DCI
  # All current 
  resuDCI[, mean(100 - All_curr)] #DCI
  # All future
  resuDCI[, mean(All_curr - All_fut)] #DCI
  resuDCI[, 100*mean((All_curr - All_fut)/All_curr)] #% DCI
  
  # SHP current
  resuDCI[, mean(100 - SHP_curr)] #DCI
  
  # SHP future
  resuDCI[, 100*mean((SHP_curr - SHP_fut)/SHP_curr)] #DCI
  
  # LHP current
  resuDCI[, mean(100 - LHP_curr)]
  
  # LHP future
  resuDCI[, 100*mean((LHP_curr - LHP_fut)/LHP_curr)]
  
  
  ## Calculate the standard deviation for decrease in DCI
  # All current
  resuDCI[, sd(100 - All_curr)]
  
  # All future
  resuDCI[, sd(All_curr - All_fut)]
  
  # SHP current
  resuDCI[, sd(100 - SHP_curr)] #DCI
  
  # SHP future
  resuDCI[, sd(100 - SHP_fut)] #DCI
  
  # LHP current
  resuDCI[, sd(100 - LHP_curr)] #DCI
  
  # LHP future
  resuDCI[, sd(100 - LHP_fut)] #DCI
  
  
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
  
  #Max decrease in % DCI
  resuDCI[, max((All_curr - All_fut)/All_curr)]
  
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
  
################################################################################
############################### DCIi ANALYSIS ##################################
################################################################################
  
  DCIname <- 'DCIi'
  
  ## read table with the outputs
  TimeTable_All <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_All.csv"))))
  TimeTable_SHP <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_Small.csv"))))
  TimeTable_LHP <- as.data.frame(read.csv(
    file.path(resdir, 
              paste0("DCI_TIME_L8_", DCIname, "_Large.csv"))))
  
  # Read table for the future projections per basin
  resuDCI <- as.data.frame(read.csv(
    file.path(resdir,
              paste0("DCI_Brazil_L8_", DCIname, ".csv")), header = T))
  
  ## Create a vector of years
  YearVec <- seq(from = 1899, to = 2018)
  
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

