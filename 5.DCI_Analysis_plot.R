#############################################################################################################
#############################################################################################################
#############################################################################################################
#####################################     PLOTS FIGURE 2     ################################################

#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT


DCIplot2 <- function(DCIname) {
  ## Import data
  resuDCI <- as.data.frame(read.csv(
    file.path(resdir, paste0("DCI_Brazil_L8_", DCIname, ".csv")), header = T))
  resuNDams <- as.data.frame(read.csv(
    file.path(resdir, paste0("basinNDams_L8_", DCIname, ".csv")), header = T))
  resuMWDams <- as.data.frame(read.csv(
    file.path(resdir, paste0("basinMWDams_L8_", DCIname, ".csv")), header = T))

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
  
  # Export csv with the percentage change
  PercChangeData <- data.frame(resuDCI$DAMBAS_ID, PercLoss_All, PercLoss_SHP, PercLoss_LHP)
  colnames(PercChangeData) <- c("HYBAS_ID", "PercLoss_All", "PercLoss_SHP", "PercLoss_LHP")
  write.csv(PercChangeData, file = 
              file.path(resdir, paste0("DCI_Brazil_L8_Percentage", DCIname,".csv")))

  ###############################################################################################################
  ##### Plot % change
  tiff(filename = file.path(figdir, paste0("Figure2", DCIname, ".tiff")), 
       height = 1656, width = 3200, res = 300, compression = c("lzw"))
  
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
  
  ######################################################################################################
  ########################### Claculate some statistics for % change ###################################
  ######################################################################################################
  
  message("Absolute change in overall DCI (mean, sd)")
  message(mean(100 - resuDCI$All_curr))
  message(sd(100 - resuDCI$All_curr))
  
  message("Absolute change in SHP DCI")
  message(mean(100 - resuDCI$SHP_curr))
  message(sd(100 - resuDCI$SHP_curr))
  
  message("Absolute change in LHP DCI")
  message(mean(100 - resuDCI$LHP_curr))
  message(sd(100 - resuDCI$LHP_curr))
  
  message("Proportion of the current scenario between SHPs and LHPs")
  message(round(mean(100 - resuDCI$SHP_curr)/mean(100 - resuDCI$LHP_curr), digits = 0))
  
  
  message("Overall future % change values (mean, sd, min, max)")
  message(mean(PercLoss_All))
  message(sd(PercLoss_All))
  message(min(PercLoss_All))
  message(max(PercLoss_All))
  
  message("SHPs future % change values")
  message( mean(PercLoss_SHP))
  message(sd(PercLoss_SHP))
  message(min(PercLoss_SHP))
  message(max(PercLoss_SHP))
  
  message("SHPs future % change values")
  message(mean(PercLoss_LHP))
  message(sd(PercLoss_LHP))
  message(min(PercLoss_LHP))
  message(max(PercLoss_LHP))
  
  message("Free-flowing All future % change values")
  message(mean(PercLoss_Free_All))
  message(sd(PercLoss_Free_All))
  message(min(PercLoss_Free_All))
  message(max(PercLoss_Free_All))
  
  message("Free-flowing SHPs future % change values")
  message(mean(PercLoss_Free_SHP))
  message(sd(PercLoss_Free_SHP))
  message(min(PercLoss_Free_SHP))
  message(max(PercLoss_Free_SHP))
  
  message("Free-flowing LHPs future % change values")
  message(mean(PercLoss_Free_LHP))
  message(sd(PercLoss_Free_LHP))
  message(min(PercLoss_Free_LHP))
  message(max(PercLoss_Free_LHP))
  
  message("Basins with both SHPs future % change values")
  message(mean(PercLoss_Both_All))
  message(sd(PercLoss_Both_All))
  message(min(PercLoss_Both_All))
  message(max(PercLoss_Both_All))
  
  message("Free-flowing SHPs future % change values")
  message(mean(PercLoss_Both_SHP))
  message(sd(PercLoss_Both_SHP))
  message(min(PercLoss_Both_SHP))
  message(max(PercLoss_Both_SHP))
  
  message("Free-flowing SHPs future % change values")
  message(mean(PercLoss_Both_LHP))
  message(sd(PercLoss_Both_LHP))
  message(min(PercLoss_Both_LHP))
  message(max(PercLoss_Both_LHP))
  
  
  ####################################################################################################################################
  ####################################################################################################################################
  ####################################### DETERMINE AVERAGE LINEAR DISTANCES OF RIVER BASINS #########################################
  ######################################  AND NUMBER OF POSSIBLE PORTIFOLIOS OF DAMS  ################################################
  ######################################        (INFORM METHOD SECTION)                   ############################################
  
  basinList <- unique(DamAttributes$DAMBAS_ID08ext)
  ## Create an empty matrix
  LinearDistKm<- rep(NA, times = length(basinList))
  NPortifolios<- rep(NA, times = length(basinList))
  
  ## Loop over basins
  for (j in 1: length(basinList)){
    
    ## filter attributes of the basin j         
    BasinX <- NetworkBRAZIL[NetworkBRAZIL$HYBAS_ID08ext == basinList[j], ]
    DamX <- DamAttributes[DamAttributes$DAMBAS_ID08ext == basinList[j], ]
    
    ## Linear distance
    LinearDist <- sum(BasinX$Shape_Length)
    LinearDistKm[j] <- LinearDist/1000
    
    ## Possible dam portifolios per basin
    NFutDams <- sum(DamX$ESTAGIO_1 == "Planned")
    NPortifolios[j] <- 2^NFutDams
    
  }
  
  message("mean(LinearDistKm)")
  message(mean(LinearDistKm))
  message("sum(NPortifolios)")
  message(sum(NPortifolios))
  
  ##################################################################################
  ############################ Plot difference #####################################
  
  tiff(filename = file.path(figdir, paste0("Figure2b", DCIname, ".tiff")), 
       height = 1396, width = 3200, res = 300, compression = c("lzw"))
  
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
}

DCIplot2('DCIp')
DCIplot2('DCIi')