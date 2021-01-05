#############         Script 1 to run DCI analysis - September 2019                ###############################
#############       Based on distances calculated in Arcgis - PythonCode          ###############################
#############             Manuscript in preparation - PNAS                        ##############################


###############################################################################################################
#########     ESTIMATE DCI FOR ALL THE SCENARIOS - BRAZIL    #################################################
##########         (CURRENT VS FUTURE) - All, SHPs, LHPs    #################################################
###########               Level 8 All Scenarios            ################################################
#########################################################################################################
#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

compute_DCIscenarios <- function(DamAttributes, NetworkBRAZIL, DCIfunc, DCIname) {
  #Base + for-loop way to compute for scenarios
  tic()
  ## Create a vector with unique basin IDs
  basinList <- as.character(unique(DamAttributes$DAMBAS_ID08ext))
  
  ## Create a matrix of basins and DCIs               
  DCI_BRAZIL_L8 <- matrix(NA, nrow = length(basinList), ncol = 7)
  colnames(DCI_BRAZIL_L8) <- c("DAMBAS_ID", "All_curr", "All_fut", "SHP_curr", "SHP_fut", "LHP_curr", "LHP_fut")
  DCI_BRAZIL_L8 [ , 1] <- as.numeric(substr(basinList, 1, nchar(basinList)-1))
  
  ## Create two matrixes to fill with dam atributes of each basin (N dams, MW production)
  basinNDams_L8 <- matrix(NA, nrow = length(basinList), ncol = 7)
  colnames(basinNDams_L8) <- c("DAMBAS_IDext", "N_All_curr", "N_All_fut", "N_SHP_curr", "N_SHP_fut", "N_LHP_curr", "N_LHP_fut")
  basinNDams_L8 [ , 1] <- as.numeric(substr(basinList, 1, nchar(basinList)-1))
  
  basinMWDams_L8 <- matrix(NA, nrow = length(basinList), ncol = 7)
  colnames(basinMWDams_L8) <- c("DAMBAS_IDext", "MW_All_curr", "MW_All_fut", "MW_SHP_curr", "MW_SHP_fut", "MW_LHP_curr", "MW_LHP_fut")
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
  
  #tic("total")
  ## Loop over scenarios
  for (x in 1: 6){
    
    ## Define types and situation according to the scenario x
    SituationScen <- get(ScenarioSituation[x])
    TypeScen <- get(ScenarioType[x])
    
    
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
        Basin <- DamX$DAMBAS_ID08ext
        
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
        d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[,2], pass = passVec) %>%
          setDT
        
        # attributes of nodes; node sizes
        d3 = data.frame(id = listSeg, l = listLengths) %>%
          setDT
        
        # Run DCI analyses for the basin
        DCI <- DCIfunc(d2, d3, print = F)
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
  write.csv(DCI_BRAZIL_L8, 
            file = file.path(resdir, paste0("DCI_Brazil_L8_", DCIname, ".csv")))
  write.csv(basinNDams_L8, 
            file = file.path(resdir, paste0("basinNDams_L8_", DCIname, ".csv")))
  write.csv(basinMWDams_L8, 
            file = file.path(resdir, paste0("basinMWDams_L8_", DCIname, ".csv")))
}


compute_DCIscenarios(DamAttributes = DamAttributes,
                     NetworkBRAZIL = NetworkBRAZIL,
                     DCIfunc = DCIp_opti5,
                     DCIname = "DCIp")

compute_DCIscenarios(DamAttributes = DamAttributes,
                     NetworkBRAZIL = NetworkBRAZIL,
                     DCIfunc = DCIi_opti,
                     DCIname = "DCIi")


######################## EXTRA STUFF ###########################################
#New way to compute for a scenario
# tic()
# DCI_L8_current <- NetworkBRAZIL[,
#                           list(DCI = DCIfunc(
#                             d2 = DamAttributes[DAMBAS_ID08ext == HYBAS_ID08ext,
#                                                list(
#                                                  id1 = DownSeg,
#                                                  id2 = UpSeg,
#                                                  pass = Allcurrent
#                                                )],
#                             d3 = .SD[, list(id=as.character(SEGID),
#                                             l=Shape_Length)],
#                             print = F)),
#                           by=.(HYBAS_ID08ext)]
# toc()