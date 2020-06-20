#############         Script 5 to run DCI analysis - September 2019                ###############################
#############        Based on distances calculated in Arcgis - PythonCode         ##############################
#############              Manuscript in preparation - PNAS                       ##############################

#############                    MIGRATORY SPECIES ANALYSIS                    ################################
########    Based on the output:  "DCIAnalysis.R" -> "DCI_Brazil_L8.csv"        ###############################

#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

###############################################################################################################
##################   FIRST STEP - Organize fish species point dataset   #######################################
##################            (merge dataset and filter species)        #######################################
###############################################################################################################

if (!file.exists("MigratorySpp_BasinID_L8.txt")) {
  ## Import raw data - Merge tables from the diffent sources
  
  MigFish_CEPTA <- as.data.frame(read.csv("MigratoryFish_RawData_CEPTA.txt", header = T))
  MigFish_CEPAM <- as.data.frame(read.csv("MigratoryFish_RawData_CEPAM.csv", header = T))
  MigFish_RedList <- as.data.frame(read.csv("MigratoryFish_RawData_RedList.txt", header = T))
  
  ##### (standardize names, remove uncesessary columns, and standardize colnames)
  
  ## Remove columns and standardize names of tables
  SimpMigFish_CEPTA <- cbind(MigFish_CEPTA[,5:10], MigFish_CEPTA$Datum)
  colnames(SimpMigFish_CEPTA) <- c("Order", "Family", "Genus", "Species", "Latitude", "Longitude", "Datum")
  
  SimpMigFish_CEPAM <- MigFish_CEPAM
  colnames(SimpMigFish_CEPAM) <- c("Order", "Family", "Genus", "Species", "Latitude", "Longitude", "Datum")
  
  GeneraRedList <- word(MigFish_RedList$taxon, 1, sep =" ")
  SimpMigFish_RedList <- cbind(MigFish_RedList[,6:7], GeneraRedList, MigFish_RedList$taxon, 
                               MigFish_RedList$decimalLat, MigFish_RedList$decimalLon,
                               MigFish_RedList$geodeticDa)
  colnames(SimpMigFish_RedList) <- c("Order", "Family", "Genus", "Species", "Latitude", "Longitude", "Datum")
  
  
  ## Keep just the two first names of the species (remove authors and years)
  FishCEPTA_names <- word(SimpMigFish_CEPTA$Species, 1,2, sep =" ")
  SimpMigFish_CEPTA$Species <- FishCEPTA_names
  
  FishCEPAM_names <- word(SimpMigFish_CEPAM$Species, 1,2, sep =" ")
  SimpMigFish_CEPAM$Species <- FishCEPAM_names
  
  ## Merge all the three tables (no need to remove) and a vector of FID
  ThreeMigFish <- rbind(SimpMigFish_RedList, SimpMigFish_CEPTA, SimpMigFish_CEPAM) 
  
  
  ## Update families names
  
  # Add levels with the new familily names and correct taxonomic errors
  levels(ThreeMigFish$Family) <- c(levels(ThreeMigFish$Family), "Triportheidae")
  levels(ThreeMigFish$Family) <- c(levels(ThreeMigFish$Family), "Bryconidae")
  
  ThreeMigFish$Family[ThreeMigFish$Genus == "Triportheus"] <- "Triportheidae"
  ThreeMigFish$Family[ThreeMigFish$Genus == "Brycon" | ThreeMigFish$Genus == "Salminus"] <- "Bryconidae"  
  
  ThreeMigFish$Family[ThreeMigFish$Species == "Bagropsis reinhardti"] <- "Pimelodidae"
  ThreeMigFish$Order[ThreeMigFish$Species == "Bagropsis reinhardti"] <- "Siluriformes"
  ThreeMigFish$Family <- droplevels(ThreeMigFish$Family, exclude = "Rivulidae")
  ThreeMigFish$Order <- droplevels(ThreeMigFish$Order, exclude = "Cyprinodontiformes")
  
  
  ## Filter species according to ichthyologists opinion
  
  # Keep the families that are entrely migratory
  famig <- ThreeMigFish[ThreeMigFish$Family == "Anostomidae" | ThreeMigFish$Family == "Bryconidae" | 
                          ThreeMigFish$Family == "Curimatidae" | ThreeMigFish$Family == "Hemiodontidae"|
                          ThreeMigFish$Family == "Prochilodontidae" | ThreeMigFish$Family == "Triportheidae" |
                          ThreeMigFish$Family == "Pimelodidae", ]
  
  # Keep jus a few Genus for some families
  cyno <- ThreeMigFish[ThreeMigFish$Genus == "Rhaphiodon", ]
  serr <- ThreeMigFish[ThreeMigFish$Genus == "Mylossoma" | ThreeMigFish$Genus == "Colossoma" |
                         ThreeMigFish$Genus == "Piaractus" | ThreeMigFish$Genus == "Myloplus", ]
  dora <- ThreeMigFish[ThreeMigFish$Genus == "Pterodoras", ]
  hept <- ThreeMigFish[ThreeMigFish$Genus == "Rhamdia", ]
  lori <- ThreeMigFish[ThreeMigFish$Genus == "Rhinelepis", ]
  
  # Combine all taxa in a new data frame
  MigFishGen <- rbind(famig, cyno, serr, dora, hept, lori)
  MigFishGen <- MigFishGen[order(MigFishGen$Species), ]
  MigFishGen <- MigFishGen[order(MigFishGen$Family), ]
  MigFishGen <- MigFishGen[order(MigFishGen$Order), ]
  
  
  ## Create an FID column
  FID <- 1: dim(MigFishGen)[1]
  
  
  ## Creat columns for red-listed and top fisheries species
  # Create vectors filled with zero
  redSpp <- levels(SimpMigFish_RedList$Species)
  valuedSpp <- c("Brachyplatystoma filamentosum", "Brachyplatystoma rousseauxii", "Brachyplatystoma vaillantii", 
                 "Colossoma macropomum",
                 "Hypophthalmus edentatus", "Hypophthalmus fimbriatus", "Hypophthalmus marginatus",
                 "Leporinus obtusidens",
                 "Myloplus asterias", "Myloplus rubripinnis", "Myloplus torquatus",
                 "Mylossoma aureum", "Mylossoma duriventre",
                 "Prochilodus argenteus", "Prochilodus brevis", "Prochilodus costatus", "Prochilodus hartii",
                 "Prochilodus lacustris", "Prochilodus lineatus", "Prochilodus nigricans", "Prochilodus rubrotaeniatus",
                 "Salminus brasiliensis",
                 "Semaprochilodus brama", "Semaprochilodus insignis", "Semaprochilodus taeniurus",
                 "Phractocephalus hemioliopterus", "Piaractus mesopotamicus", "Pinirampus pirinampu",
                 "Pseudoplatystoma corruscans", "Pseudoplatystoma punctifer", "Pseudoplatystoma reticulatum")
  
  
  RedListed <- rep(0, times = dim(MigFishGen)[1])
  Fisheries <- rep(0, times = dim(MigFishGen)[1])
  
  # Insert 1s for the the target species in the vector
  RedListed[MigFishGen$Species %in% redSpp] <- 1
  Fisheries[MigFishGen$Species %in% valuedSpp] <- 1
  
  
  ### Create the final data frame and export it as a csv
  MigFishGeneralData <- cbind(FID, MigFishGen, RedListed, Fisheries)
  write.csv(MigFishGeneralData, file = "MigSpeciesForGIS.csv")
  
  
  ###############################################################################################################
  ##################   SECOND STEP - Import "MigSpeciesForGIS.csv" into GIS   ###################################
  ##################    (See the document "AnalysisPLAN_MigratorySpp.docx")   ###################################
  ###############################################################################################################
}

###############################################################################################################
################      THIRD STEP - Prepare plot based on GIS outputs              #################################
###############    (Import species coordinates - "MigratorySpp_BasinID_L8.txt")   ###################################
###############################################################################################################

assess_migratoryfragment <- function(DCIname) {
  SppDatatab = "MigratorySpp_BasinID_L8.txt"
  resuDCItab = file.path(resdir, paste0("DCI_Brazil_L8_", DCIname, ".csv"))
  
  ## Import data for the plot (After processing in ArcGIS)
  SppData_L8 <- as.data.frame(read.csv(SppDatatab, header = T))
  resuDCI_L8 <- as.data.frame(read.csv(resuDCItab, header = T))
  
  ## List the species in a vector
  spList <- levels(SppData_L8$Species)
  
  ## Create a matrix to fill in with DCI averages and CI per species and scenario
  SpMeanDCI_curr <- matrix(NA, nrow = length(spList), ncol = 8)
  rownames(SpMeanDCI_curr) <- spList
  colnames(SpMeanDCI_curr) <- c("Curr_SHP_mean", "Curr_SHP_Upp", "Curr_SHP_Low", "Curr_LHP_mean", "Curr_LHP_Upp", "Curr_LHP_Low", "RedListed", "TopFishery")
  
  SpMeanDCI_fut <- matrix(NA, nrow = length(spList), ncol = 8)
  rownames(SpMeanDCI_fut) <- spList
  colnames(SpMeanDCI_fut) <- c("Fut_SHP_mean", "Fut_SHP_Upp", "Fut_SHP_Low", "Fut_LHP_mean", "Fut_LHP_Upp", "Fut_LHP_Low", "RedListed", "TopFishery")
  
  
  ## Find the average DCI and CI for basins where the species i occur
  for (i in 1: length(spList)){
    
    # Select a species and get basin data
    spID <- levels(SppData_L8$Species)[i]
    spIDvec <- SppData_L8[SppData_L8$Species == spID & SppData_L8$HYBAS_ID > 0, ]
    
    # Filter unique basins (some spp. have multiple occurrence records at the same basin)
    SpBasins <- unique(spIDvec$HYBAS_ID)
    
    # Get the DCI's for such basins
    SpDCIDF <- resuDCI_L8[resuDCI_L8$DAMBAS_ID %in% SpBasins, ]
    
    
    ## Fill out a table with average DCI and CI for current and future scenarios
    # mean current
    SpMeanDCI_curr[i, 1] <- mean(SpDCIDF$SHP_curr)
    SpMeanDCI_curr[i, 4] <- mean(SpDCIDF$LHP_curr)
    
    # Upper CI values current
    SpMeanDCI_curr[i, 2] <- quantile(SpDCIDF$SHP_curr, 0.975, na.rm = T)
    SpMeanDCI_curr[i, 5] <- quantile(SpDCIDF$LHP_curr, 0.975, na.rm = T)
    
    # Lower CI values current
    SpMeanDCI_curr[i, 3] <- quantile(SpDCIDF$SHP_curr, 0.025, na.rm = T)
    SpMeanDCI_curr[i, 6] <- quantile(SpDCIDF$LHP_curr, 0.025, na.rm = T)
    
    # mean future
    SpMeanDCI_fut[i, 1] <- mean(SpDCIDF$SHP_fut, na.rm = T)
    SpMeanDCI_fut[i, 4] <- mean(SpDCIDF$LHP_fut, na.rm = T)
    
    # Upper CI values future
    SpMeanDCI_fut[i, 2] <- quantile(SpDCIDF$SHP_fut, 0.975, na.rm = T)
    SpMeanDCI_fut[i, 5] <- quantile(SpDCIDF$LHP_fut, 0.975, na.rm = T)
    
    # Lower CI values future
    SpMeanDCI_fut[i, 3] <- quantile(SpDCIDF$SHP_fut, 0.025, na.rm = T)
    SpMeanDCI_fut[i, 6] <- quantile(SpDCIDF$LHP_fut, 0.025, na.rm = T)
    
    ## Indentify Red Listed or Top Fisheries Species 
    SpMeanDCI_curr[i, 7] <- max(SppData_L8$RedListed[SppData_L8$Species == spID])
    SpMeanDCI_curr[i, 8] <- max(SppData_L8$Fisheries[SppData_L8$Species == spID])
    
    SpMeanDCI_fut[i, 7] <- max(SppData_L8$RedListed[SppData_L8$Species == spID])
    SpMeanDCI_fut[i, 8] <- max(SppData_L8$Fisheries[SppData_L8$Species == spID])
    
  }
  
  ## Assign Brycon opalinus as red-listed
  SpMeanDCI_curr[rownames(SpMeanDCI_curr) == "Brycon opalinus", 7] <- 1
  SpMeanDCI_fut[rownames(SpMeanDCI_fut) == "Brycon opalinus", 7] <- 1
  
  
  ########################################################################################################
  ################################### PLOT FIGURE ########################################################
  ########################################################################################################
  
  ## Plot figure
  outfigname = file.path(figdir, paste0("Figure4_", DCIname, "_",
                                        format(Sys.Date(), '%Y%m%d'), ".tiff"))
  tiff(filename = outfigname, 
       height = 2596, width = 4300, res = 300, compression = c("lzw"))
  
  ## Plot figure
  mat <- matrix(c(1, 3, 2), ncol = 3, nrow = 1, byrow = T)
  layout(mat = mat, widths = c(2, 0.15, 2), heights = c(2, 2, 2))
  
  par(oma = c(9, 11.5, 9, 4), mar = c(0.5, 0, 0, 0), bty = "n")
  
  # Panel A - current
  plot(SpMeanDCI_curr[,4], SpMeanDCI_curr[,1], xlim = c(0, 100), ylim = c(0, 100), type = "n", 
       ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  
  axis(side = 2, at = c(0, 20, 40, 60, 80, 100), cex.axis = 3.5, line = 0)
  axis(side = 1, at = c(0, 20, 40, 60, 80, 100), cex.axis = 3.5, line = 0, mgp = c(3, 2.1, 0))
  
  mtext("SHP river connectivity (DCI)", side = 2, cex = 3.5, line = 5.8)
  mtext("LHP river connectivity (DCI)", side = 1, cex = 3.5, line = 7.5, at = 100)
  
  
  ## Plot CI lines
  segments(y0 = SpMeanDCI_curr[,2], x0 = SpMeanDCI_curr[,4],
           y1 = SpMeanDCI_curr[,3], x1 = SpMeanDCI_curr[,4], lwd = 3.0, col = "#CCCCCC60")
  
  segments(y0 = SpMeanDCI_curr[,1], x0 = SpMeanDCI_curr[,5],
           y1 = SpMeanDCI_curr[,1], x1 = SpMeanDCI_curr[,6], lwd = 3.0, col = "#CCCCCC60")
  
  ## Clean space over the lines for points
  points(SpMeanDCI_curr[,4], SpMeanDCI_curr[, 1], pch = 21, cex = 4.5, bg = "white", col = "white")
  
  ## Plot 1:1 ratio line
  abline(a = 0, b = 1, lwd = 6, col = "gray30", lty = 3)
  
  ## Plot fish illustrations
  require(jpeg)
  img <- readJPEG("VerticalComposition.jpg")
  rasterImage(img,0,0,27,84)
  
  ## Plot points
  # points(SpMeanDCI_curr[,4], SpMeanDCI_curr[, 1], pch = 21, cex = 4.5, bg = "#91BFDB95", col = "#66666699")
  points(SpMeanDCI_curr[,4], SpMeanDCI_curr[, 1], pch = 21, cex = 4.5, bg = "#73737395", col = "#4F4F4F")
  points(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 4], SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 1], 
         pch = 21, cex = 4.5, bg = "#FC8D5999", col = "#4F4F4F")
  points(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 4], SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 1],
         pch = 21, cex = 4.5, bg = "#FFF68F99", col = "#4F4F4F")
  
  # ## Plot CI red-listed and valued
  # CI_currSHP_Low_Red <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 4], 0.025, na.rm = T)
  # CI_currSHP_Up_Red <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 4], 0.975, na.rm = T)
  # CI_currLHP_Low_Red <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 1], 0.025, na.rm = T)
  # CI_currLHP_Up_Red <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 1], 0.975, na.rm = T)
  # 
  # segments(y0 = CI_currLHP_Low_Red, x0 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 4], na.rm = T),
  #          y1 = CI_currLHP_Up_Red, x1 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 4], na.rm = T), lwd = 3.0, col = "#FC8D59")
  # 
  # segments(y0 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 1], na.rm = T), x0 = CI_currSHP_Low_Red,
  #          y1 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 1], na.rm = T), x1 = CI_currSHP_Up_Red, lwd = 3.0, col = "#FC8D59")
  # 
  # 
  # CI_currSHP_Low_Val <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 4], 0.025, na.rm = T)
  # CI_currSHP_Up_Val <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 4], 0.975, na.rm = T)
  # CI_currLHP_Low_Val <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 1], 0.025, na.rm = T)
  # CI_currLHP_Up_Val <- quantile(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 1], 0.975, na.rm = T)
  # 
  # segments(y0 = CI_currLHP_Low_Val, x0 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 4], na.rm = T),
  #          y1 = CI_currLHP_Up_Val, x1 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 4], na.rm = T), lwd = 3.0, col = "#FFF68F")
  # 
  # segments(y0 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 1], na.rm = T), x0 = CI_currSHP_Low_Val,
  #          y1 = mean(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 1], na.rm = T), x1 = CI_currSHP_Up_Val, lwd = 3.0, col = "#FFF68F")
  # 
  # ## Average red-listed and valued
  # points(mean(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 4], na.rm = T), mean(SpMeanDCI_curr[SpMeanDCI_curr[,7] == 1, 1], na.rm = T),
  #        pch = 23, cex = 4, bg = "white", col = "#FC8D59")
  # points(mean(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 4], na.rm = T), mean(SpMeanDCI_curr[SpMeanDCI_curr[,8] == 1, 1], na.rm = T),
  #        pch = 23, cex = 4, bg = "white", col = "#FFF68F")
  
  ## Plot overall CI
  CI_currSHP_Low <- quantile(SpMeanDCI_curr[,4], 0.025, na.rm = T)
  CI_currSHP_Up <- quantile(SpMeanDCI_curr[,4], 0.975, na.rm = T)
  CI_currLHP_Low <- quantile(SpMeanDCI_curr[,1], 0.025, na.rm = T)
  CI_currLHP_Up <- quantile(SpMeanDCI_curr[,1], 0.975, na.rm = T)
  
  segments(y0 = CI_currLHP_Low, x0 = mean(SpMeanDCI_curr[,4], na.rm = T),
           y1 = CI_currLHP_Up, x1 = mean(SpMeanDCI_curr[,4], na.rm = T), lwd = 3.0, col = "black")
  
  segments(y0 = mean(SpMeanDCI_curr[,1], na.rm = T), x0 = CI_currSHP_Low,
           y1 = mean(SpMeanDCI_curr[,1], na.rm = T), x1 = CI_currSHP_Up, lwd = 3.0, col = "black")
  
  ## Plot overall average
  points(mean(SpMeanDCI_curr[,4], na.rm = T), mean(SpMeanDCI_curr[, 1], na.rm = T),
         pch = 25, cex = 5, bg = "black", col = "black")
  
  
  ## Title
  mtext("Present", side = 3, cex = 3.7, line = 2)
  mtext("A", side = 3, cex = 4.5, line = 3, at = 2)
  
  # Panel B - future
  plot(SpMeanDCI_fut[,4], SpMeanDCI_fut[,1], xlim = c(0, 100), ylim = c(0, 100), type = "n", 
       ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  
  axis(side = 2, at = c(0, 20, 40, 60, 80, 100), labels = c("", "", "", "", "", ""), cex.axis = 3.5, line = 0)
  axis(side = 1, at = c(0, 20, 40, 60, 80, 100), cex.axis = 3.5, line = 0, mgp = c(3, 2.1, 0))
  
  ## Plot lines
  segments(y0 = SpMeanDCI_fut[,2], x0 = SpMeanDCI_fut[,4],
           y1 = SpMeanDCI_fut[,3], x1 = SpMeanDCI_fut[,4], lwd = 3.0, col = "#CCCCCC60")
  
  segments(y0 = SpMeanDCI_fut[,1], x0 = SpMeanDCI_fut[,5],
           y1 = SpMeanDCI_fut[,1], x1 = SpMeanDCI_fut[,6], lwd = 3.0, col = "#CCCCCC60")
  
  ## Clean space over the lines
  points(SpMeanDCI_fut[,4], SpMeanDCI_fut[, 1], pch = 21, cex = 4.5, bg = "white", col = "white")
  
  ## Plot 1:1 ratio line
  abline(a = 0, b = 1, lwd = 6, col = "gray30", lty = 3)
  
  # Plot points
  #points(SpMeanDCI_fut[,4], SpMeanDCI_fut[, 1], pch = 21, cex = 4.5, bg = "#91BFDB95", col = "#66666699")
  points(SpMeanDCI_fut[,4], SpMeanDCI_fut[, 1], pch = 21, cex = 4.5, bg = "#73737395", col = "#4F4F4F")
  points(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 4], SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 1], 
         pch = 21, cex = 4.5, bg = "#FC8D5999", col = "#4F4F4F")
  points(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 4], SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 1], 
         pch = 21, cex = 4.5, bg = "#FFF68F99", col = "#4F4F4F")
  
  # ## Plot CI red-listed and valued
  # CI_futSHP_Low_Red <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 4], 0.025, na.rm = T)
  # CI_futSHP_Up_Red <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 4], 0.975, na.rm = T)
  # CI_futLHP_Low_Red <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 1], 0.025, na.rm = T)
  # CI_futLHP_Up_Red <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 1], 0.975, na.rm = T)
  # 
  # segments(y0 = CI_futLHP_Low_Red, x0 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 4], na.rm = T),
  #          y1 = CI_futLHP_Up_Red, x1 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 4], na.rm = T), lwd = 3.0, col = "#FC8D59")
  # 
  # segments(y0 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 1], na.rm = T), x0 = CI_futSHP_Low_Red,
  #          y1 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 1], na.rm = T), x1 = CI_futSHP_Up_Red, lwd = 3.0, col = "#FC8D59")
  # 
  # 
  # CI_futSHP_Low_Val <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 4], 0.025, na.rm = T)
  # CI_futSHP_Up_Val <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 4], 0.975, na.rm = T)
  # CI_futLHP_Low_Val <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 1], 0.025, na.rm = T)
  # CI_futLHP_Up_Val <- quantile(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 1], 0.975, na.rm = T)
  # 
  # segments(y0 = CI_futLHP_Low_Val, x0 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 4], na.rm = T),
  #          y1 = CI_futLHP_Up_Val, x1 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 4], na.rm = T), lwd = 3.0, col = "#FFF68F")
  # 
  # segments(y0 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 1], na.rm = T), x0 = CI_futSHP_Low_Val,
  #          y1 = mean(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 1], na.rm = T), x1 = CI_futSHP_Up_Val, lwd = 3.0, col = "#FFF68F")
  # 
  # ## Average red-listed and valued
  # points(mean(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 4], na.rm = T), mean(SpMeanDCI_fut[SpMeanDCI_fut[,7] == 1, 1], na.rm = T),
  #        pch = 23, cex = 4, bg = "white", col = "#FC8D59")
  # points(mean(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 4], na.rm = T), mean(SpMeanDCI_fut[SpMeanDCI_fut[,8] == 1, 1], na.rm = T),
  #        pch = 24, cex = 4, bg = "white", col = "#FFF68F")
  
  # Plot overall CI
  CI_futSHP_Low <- quantile(SpMeanDCI_fut[, 4], 0.025, na.rm = T)
  CI_futSHP_Up <- quantile(SpMeanDCI_fut[, 4], 0.975, na.rm = T)
  CI_futLHP_Low <- quantile(SpMeanDCI_fut[, 1], 0.025, na.rm = T)
  CI_futLHP_Up <- quantile(SpMeanDCI_fut[, 1], 0.975, na.rm = T)
  
  segments(y0 = CI_futLHP_Low, x0 = mean(SpMeanDCI_fut[,4], na.rm = T),
           y1 = CI_futLHP_Up, x1 = mean(SpMeanDCI_fut[,4], na.rm = T), lwd = 3.0, col = "black")
  
  segments(y0 = mean(SpMeanDCI_fut[, 1], na.rm = T), x0 = CI_futSHP_Low,
           y1 = mean(SpMeanDCI_fut[, 1], na.rm = T), x1 = CI_futSHP_Up, lwd = 3.0, col = "black")
  
  ## Overal average
  points(mean(SpMeanDCI_fut[,4], na.rm = T), mean(SpMeanDCI_fut[, 1], na.rm = T),
         pch = 25, cex = 5, bg = "black", col = "black")
  
  
  ## Plot legend
  polygon(x = c(1, 1, 37, 37), y = c(0, 26, 26, 0), col = "white", border = NA)
  legend(x = -3, y = 30, pch = c(21, 21, 21, 25), cex = 3.0, pt.cex = c(5, 5, 5, 4), bg = "white",
         legend = c("Red-listed", "Highly valuable", "Other migratory", "Overall average"),
         pt.bg = c("#FC8D5999","#FFF68F99", "#73737395", "black"), xpd = T, bty = "n")   #Blue "#91BFDB95"
  
  ## Title
  mtext("Future", side = 3, cex = 3.7, line = 2)
  mtext("B", side = 3, cex = 4.5, line = 3, at = 2)
  
  
  dev.off()
  
  
  
  ###########################################################
  ########### Calculate some statistics #####################
  ###########################################################
  
  ## Average current
  mean(SpMeanDCI_curr[,1], na.rm = T)
  mean(SpMeanDCI_curr[,4], na.rm = T)
  sd(SpMeanDCI_curr[,1], na.rm = T)
  sd(SpMeanDCI_curr[,4], na.rm = T)
  
  ## Average future
  mean(SpMeanDCI_fut[,1], na.rm = T)
  mean(SpMeanDCI_fut[,4], na.rm = T)
  sd(SpMeanDCI_fut[,1], na.rm = T)
  sd(SpMeanDCI_fut[,4], na.rm = T)
  
  
  # Get the number of species belonging to each category
  AnalysSppFut <- SpMeanDCI_fut[!is.na(SpMeanDCI_fut[,2]),]
  AnalysSppCurr <- SpMeanDCI_curr[!is.na(SpMeanDCI_curr[,2]),]
  
  FinalSpNumbFut <- dim(AnalysSppFut)[1]
  FinalSpNumbCurr <- dim(AnalysSppCurr)[1]
  
  RedListNumbSp <- sum(AnalysSppFut[,7])
  RedListNumbSp <- sum(AnalysSppCurr[,7])
  
  FisheryNumbSp <- sum(AnalysSppFut[,8])
  FisheryNumbSp <- sum(AnalysSppCurr[,8])
  
  
  ## Number of species above and below the 1:1 line
  # Current raw
  sum(SpMeanDCI_curr[,1]/SpMeanDCI_curr[,4] > 1, na.rm = T)
  sum(SpMeanDCI_curr[,1]/SpMeanDCI_curr[,4] < 1, na.rm = T)
  sum(SpMeanDCI_curr[,1]/SpMeanDCI_curr[,4] == 1, na.rm = T)
  
  # Current percentage
  sum(SpMeanDCI_curr[,1]/SpMeanDCI_curr[,4] > 1, na.rm = T)/FinalSpNumbFut * 100
  sum(SpMeanDCI_curr[,1]/SpMeanDCI_curr[,4] < 1, na.rm = T)/FinalSpNumbFut * 100
  sum(SpMeanDCI_curr[,1]/SpMeanDCI_curr[,4] == 1, na.rm = T)/FinalSpNumbFut * 100
  
  # Future raw
  sum(SpMeanDCI_fut[,1]/SpMeanDCI_fut[,4] > 1, na.rm = T)
  sum(SpMeanDCI_fut[,1]/SpMeanDCI_fut[,4] < 1, na.rm = T)
  sum(SpMeanDCI_fut[,1]/SpMeanDCI_fut[,4] == 1, na.rm = T)
  
  # Future percentage
  sum(SpMeanDCI_fut[,1]/SpMeanDCI_fut[,4] > 1, na.rm = T)/FinalSpNumbFut * 100
  sum(SpMeanDCI_fut[,1]/SpMeanDCI_fut[,4] < 1, na.rm = T)/FinalSpNumbFut * 100
  sum(SpMeanDCI_fut[,1]/SpMeanDCI_fut[,4] == 1, na.rm = T)/FinalSpNumbFut * 100
  
  
  ############################################################
  ## Number of red-listed species above and below the 1:1 line
  # Get the name of the red-listed species
  RedListedMean_curr <- SpMeanDCI_curr[which(SpMeanDCI_curr[,7] == 1),]
  RedListedMean_fut <- SpMeanDCI_fut[which(SpMeanDCI_fut[,7] == 1),]
  
  NRedList <- RedListedMean_curr[!is.na(RedListedMean_curr[,2]) == T,]
  
  
  ## Number of species above and below the 1:1 line -  SHP ratio < 1
  # Current raw
  sum(RedListedMean_curr[,1]/RedListedMean_curr[,4] > 1, na.rm = T) 
  sum(RedListedMean_curr[,1]/RedListedMean_curr[,4] < 1, na.rm = T) 
  sum(RedListedMean_curr[,1]/RedListedMean_curr[,4] == 1, na.rm = T)
  
  # Current percentage
  sum(RedListedMean_curr[,1]/RedListedMean_curr[,4] > 1, na.rm = T)/dim(NRedList)[1] * 100
  sum(RedListedMean_curr[,1]/RedListedMean_curr[,4] < 1, na.rm = T)/dim(NRedList)[1] * 100
  sum(RedListedMean_curr[,1]/RedListedMean_curr[,4] == 1, na.rm = T)/dim(NRedList)[1] * 100
  
  
  # Future raw
  sum(RedListedMean_fut[,1]/RedListedMean_fut[,4] > 1, na.rm = T)
  sum(RedListedMean_fut[,1]/RedListedMean_fut[,4] < 1, na.rm = T)
  sum(RedListedMean_fut[,1]/RedListedMean_fut[,4] == 1, na.rm = T)
  
  # Future percentage
  sum(RedListedMean_fut[,1]/RedListedMean_fut[,4] > 1, na.rm = T)/dim(NRedList)[1] * 100
  sum(RedListedMean_fut[,1]/RedListedMean_fut[,4] < 1, na.rm = T)/dim(NRedList)[1] * 100
  sum(RedListedMean_fut[,1]/RedListedMean_fut[,4] == 1, na.rm = T)/dim(NRedList)[1] * 100
  
  
  ############################################################
  ## Number of Valuable species above and below the 1:1 line
  # Get the name of the valuable species
  FisheryMean_curr <- SpMeanDCI_curr[which(SpMeanDCI_curr[,8] == 1),]
  FisheryMean_fut <- SpMeanDCI_fut[which(SpMeanDCI_fut[,8] == 1),]
  
  NFishery <- FisheryMean_curr[!is.na(FisheryMean_curr[,2]) == T,]
  
  
  ## Number of species above and below the 1:1 line
  # Current raw
  sum(FisheryMean_curr[,1]/FisheryMean_curr[,4] > 1, na.rm = T)
  sum(FisheryMean_curr[,1]/FisheryMean_curr[,4] < 1, na.rm = T)
  sum(FisheryMean_curr[,1]/FisheryMean_curr[,4] == 1, na.rm = T)
  
  # Current percentage
  sum(FisheryMean_curr[,1]/FisheryMean_curr[,4] > 1, na.rm = T)/dim(NFishery)[1] * 100
  sum(FisheryMean_curr[,1]/FisheryMean_curr[,4] < 1, na.rm = T)/dim(NFishery)[1] * 100
  sum(FisheryMean_curr[,1]/FisheryMean_curr[,4] == 1, na.rm = T)/dim(NFishery)[1] * 100
  
  
  # Future raw
  sum(FisheryMean_fut[,1]/FisheryMean_fut[,4] > 1, na.rm = T)
  sum(FisheryMean_fut[,1]/FisheryMean_fut[,4] < 1, na.rm = T)
  sum(FisheryMean_fut[,1]/FisheryMean_fut[,4] == 1, na.rm = T)
  
  # Future percentage
  sum(FisheryMean_fut[,1]/FisheryMean_fut[,4] > 1, na.rm = T)/dim(NFishery)[1] * 100
  sum(FisheryMean_fut[,1]/FisheryMean_fut[,4] < 1, na.rm = T)/dim(NFishery)[1] * 100
  sum(FisheryMean_fut[,1]/FisheryMean_fut[,4] == 1, na.rm = T)/dim(NFishery)[1] * 100
  
  
  
  #################################################################################################
  ######################### Rank percentage change ################################################
  
  ## Ordered table by Raw DCI
  OrderedSppTable <- as.data.frame(AnalysSppFut[order(AnalysSppFut[,1]),])
  OrderedSppTable[OrderedSppTable$RedListed == 1, ]
  OrderedSppTable[OrderedSppTable$TopFishery == 1, ]
  
  ## Claculate percentage change
  SppPercSHP <- AnalysSppFut[,1] - AnalysSppCurr[,1] / AnalysSppCurr[,1] * 100
  SppPercLHP <- AnalysSppFut[,4] - AnalysSppCurr[,4] / AnalysSppCurr[,4] * 100
  
  PercChange <- as.data.frame(cbind(SppPercSHP, AnalysSppFut))
  
  PercDataFrame <- PercChange[order(PercChange$SppPercSHP), ]
  PercDataFrame[PercDataFrame$RedListed == 1, ]
  PercDataFrame[PercDataFrame$TopFishery == 1, ]
  
  ## Write a csv for the supplement
  CombinedTableSpp <- cbind(SppPercSHP, SppPercLHP, AnalysSppCurr[,1:6], AnalysSppFut)
  # write.csv(CombinedTableSpp, file = "SupplementS1_SpeciesDCI.csv")
  
  
  ## Get the number of species that will experience a decrease in DCI of more than 40%
  sum(SppPercSHP < -40)
}

assess_migratoryfragment(DCIname = 'DCIp')
assess_migratoryfragment(DCIname = 'DCIi')
