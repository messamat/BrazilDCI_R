#############                   Script 7 Supplements - September 2019              ###############################
#############                             (DCI i)                                  ###############################
#############          Based on DCIs from script 1 and folder "Shapes"             ###############################
#############              Manuscript in preparation - PNAS                        ###############################
##################################################################################################################


#############################################################################################################
#####################################     PLOTS FIGURE 2     ################################################

## Import data
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIi.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8_DCIi.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8_DCIi.csv", header = T))


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


## Plot % change
tiff(filename = "SI_Appendix_Figure2.tiff", height = 1656, width = 3200, res = 300, compression = c("lzw"))

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


##################################################################################################################
################################# Supplement Figure 2 ############################################################
##################################################################################################################

## Import data
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIi.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8_DCIi.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8_DCIi.csv", header = T))


## Filter results of basins with both SHPs and LHPs in the future
resuDCI_Both <- resuDCI[resuDCI$SHP_fut != 100 & resuDCI$LHP_fut != 100, ]

## Filter results just for basins that are currently free of hydropower
resuDCI_Free <- resuDCI[resuDCI$All_curr == 100, ]

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


##### Plot % change
tiff(filename = "SI_Appendix_Figure2.tiff", height = 1656, width = 3200, res = 300, compression = c("lzw"))

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



##################################################################################################################
################################# Supplement Figure 3 ############################################################
##################################################################################################################

## Upload packages
library(rprojroot)
require(rgdal)
require(maps)
require(RColorBrewer)
require(sp)
require(zoom)

## Import network and dams dataset (alternative)
rootdir <- find_root(has_dir("Shapes"))
datadir <- file.path(rootdir, "Shapes")
SouthAmerDir <- file.path(datadir, "South America")
WatershedsDir <- file.path(datadir, "Watersheds ANA")
BasinsDir <- file.path(datadir, "Python R Outputs")

# Read table for the future projections per basin (absolut values and percentage change)
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8_DCIi.csv", header = T))
PercChangeData <- as.data.frame(read.csv("DCI_Brazil_L8_PercentageDCIi.csv"))

## Assign breaks to values of DCI (100 - DCI)
CurSHPCut <- cut(100 - resuDCI$SHP_curr, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
CurLHPCut <- cut(100 - resuDCI$LHP_curr, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))

FutSHPCut <- cut(100 - resuDCI$SHP_fut, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
FutLHPCut <- cut(100 - resuDCI$LHP_fut, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))

PercSHPCut <- cut(-1 * PercChangeData$PercLoss_SHP, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
PercLHPCut <- cut(-1 * PercChangeData$PercLoss_LHP, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))


## Create a factor with color names
CurrColors <- factor(paste(CurSHPCut, CurLHPCut, sep = ""), 
                     labels = c("#756BB1", "#A96B8C", "#DD6B66", "#DE2D26",
                                "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                "#3182BD", "#9ECAE1", "#DEEBF7", "white"))

FutuColors <- factor(paste(FutSHPCut, FutLHPCut, sep = ""), 
                     labels = c("#756BB1", "#A96B8C", "#DD6B66", "#DE2D26",
                                "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                "#3182BD", "#9ECAE1", "#DEEBF7", "white"))

PercColors <- factor(paste(PercSHPCut, PercLHPCut, sep = ""), 
                     labels = c("#756BB1", "#A96B8C", "#DD6B66", "#DE2D26",
                                "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                "#3182BD", "#9ECAE1", "#DEEBF7", "white"))



## Create a dataframe with the colors for all the scenarios for each basin
DCIColors <- data.frame(resuDCI$HYBAS_ID, CurrColors, FutuColors, PercColors)
colnames(DCIColors) <- c("HYBAS_ID", "CurrColors", "FutuColors", "PercColors")


## Import shapefiles
## Use the function "View" to view the columns of the shapefiles
SouthAmerShape <- readOGR(dsn = SouthAmerDir, layer = "South_America")
WatershedsShape <- readOGR(dsn = WatershedsDir, layer = "Regiões_Hidrográficas")
BasinsShape <- readOGR(dsn = BasinsDir, layer = "BasinAtlas_level08_Brazil_Clipped")

## Filter basins that have or will have at least on dam
BasinDams <- BasinsShape[which(BasinsShape$HYBAS_ID %in% resuDCI$HYBAS_ID),]

## Merge DCI colors in the Basin Shapefile
DCImapBasins <- merge(BasinDams, DCIColors, by = "HYBAS_ID", duplicateGeoms = TRUE)



## Plot figure - use zm() to zoom in the map
tiff(filename = "SI_Appendix_Figure3.tiff", height = 1926, width = 4400, res = 300, compression = c("lzw"))

par (mfrow = c(1,3), oma = c(0, 2, 0, 2), mar = c(0, 0, 0, 0), bty = "n")

#plot 1
BrazilLayer <- SouthAmerShape$COUNTRY == "Brazil"
plot(SouthAmerShape[BrazilLayer, ], col = "white", border = F, add = F)
mtext("Present", side = 3, line = -5, cex = 2.8)
mtext("A", side = 3, line = -4.5, at = -73, cex = 3.4)
plot(DCImapBasins, col = as.character(DCImapBasins$CurrColors), add = T, border= F)
plot(WatershedsShape, col = "#FFFFFF00", border= T, lwd = 1.5, add = T)

## Plot South America map
require(jpeg)
img <- readJPEG("SA.jpeg")
rasterImage(img,xleft = -73, xright = -60, ybottom = -34, ytop = -19)
mtext("South America", side = 1, cex = 1.4, line = -7.0, at = -66.6)

#plot 2
BrazilLayer <- SouthAmerShape$COUNTRY == "Brazil"
plot(SouthAmerShape[BrazilLayer, ], col = "white", border = F, add = F)
mtext("Future", side = 3, line = -5, cex = 2.8)
mtext("B", side = 3, line = -4.5, at = -73, cex = 3.4)
plot(DCImapBasins, col = as.character(DCImapBasins$FutuColors), add = T, border= F)
plot(WatershedsShape, col = "#FFFFFF00", border= T, lwd = 1.5, add = T)

## Plot scale bar
map.scale(x = -73, y = -32, ratio = FALSE, metric = T, relwidth = 0.20, cex = 1.8)  

#plot 3
BrazilLayer <- SouthAmerShape$COUNTRY == "Brazil"
plot(SouthAmerShape[BrazilLayer, ], col = "white", border = F, add = F)
mtext("% change", side = 3, line = -5, cex = 2.8)
mtext("C", side = 3, line = -4.5, at = -73, cex = 3.4)
plot(DCImapBasins, col = as.character(DCImapBasins$PercColors), add = T, border= F)
plot(WatershedsShape, col = "#FFFFFF00", border= T, lwd = 1.5, add = T)

## Plot North Arrow
require(jpeg)
img <- readJPEG("NorthArrow.jpg")
rasterImage(img,xleft = -42, xright = -39, ybottom = -32, ytop = -25)

## Plot color scale
require(jpeg)
img <- readJPEG("ColorScale.JPEG")
rasterImage(img,xleft = -76, xright = -61, ybottom = -34, ytop = -19)

## Insert the numbers
mtext("0", side = 1, cex = 1.4, line = -6.5, at = -76.2)
mtext("1", side = 1, cex = 1.4, line = -6.5, at = -72)
mtext("25", side = 1, cex = 1.4, line = -6.5, at = -68.5)
mtext("50", side = 1, cex = 1.4, line = -6.5, at = -64.8)
mtext("100", side = 1, cex = 1.4, line = -6.5, at = -60.5)

mtext("1", side = 2, cex = 1.4, line = 0.2, at = -30.4)
mtext("25", side = 2, cex = 1.4, line = 0.2, at = -26.85)
mtext("50", side = 2, cex = 1.4, line = 0.2, at = -22.9)
mtext("100", side = 2, cex = 1.4, line = 0.2, at = -19.0)


dev.off()




##################################################################################################################
################################# Supplement Figure 4 ############################################################
##################################################################################################################

## Import data for the plot (After processing in ArcGIS)
SppData_L8 <- as.data.frame(read.csv("MigratorySpp_BasinID_L8.txt", header = T))
resuDCI_L8 <- as.data.frame(read.csv("DCI_Brazil_L8_DCIi.csv", header = T))

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
  SpDCIDF <- resuDCI_L8[resuDCI_L8$HYBAS_ID %in% SpBasins, ]
  
  
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


## Plot figure
tiff(filename = "SI_Appendix_Figure4.tiff", height = 2596, width = 4300, res = 300, compression = c("lzw"))

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

# ## Plot fish illustrations
# require(jpeg)
# img <- readJPEG("VerticalComposition.jpg")
# rasterImage(img,0,0,27,84)

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



##################################################################################################################
################################# Supplement Figure 5 ############################################################
##################################################################################################################

## Import the csv with the prioritization results
PriorityAnalysis <- read.csv("IndividualDam_DCIi.csv", header = T)

## Plot 5 (DCI loss Vs Capacity)
tiff(filename = "SI_Appendix_Figure5.tiff", height = 2396, width = 3700, res = 300, compression = c("lzw"))
par(oma = c(6, 9, 2, 0.5), mar = c(0.5, 0, 0, 0), bty = "n")


## Create a color vector to differentiate SHP and LHP
TypeFac <- factor(x = PriorityAnalysis$Type)
TypeColor <- factor(x = TypeFac, levels = levels(TypeFac), 
                    labels = c("#4876FF50", "#8B000050"))

## Plot figure #xaxt = "n", log(PriorityAnalysis$Capacity)
plot(PriorityAnalysis$Capacity, PriorityAnalysis$DCIMeanDiff, type = "n", ylim = c(-120, 0), xlim = c(0.1, 8000),
     ylab = "", xlab = "",  yaxt = "n", xaxt = "n", log = "x") 

options(scipen=5)
axis(side = 2, at = c(-100, -80, -60, -40, -20, 0), cex.axis = 2.6, line = 0)
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
legend(x = 200, y = -95, pch = c(21), legend = c("Planned SHP", "Planned LHP"), cex = 2.5,
       pt.cex = 3.3, pt.bg = c("#8B000060","#4876FF60"), xpd = T, bty = "n")

dev.off()



## Some Statistics  ##

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





##################################################################################################################
################################# Supplement Figure 6 ############################################################
##################################################################################################################

## Loand packages
require(rPref)

## Import data on scenarios
NationalScenarios1 <- read.csv("NationalScen_DCIi_MultiRun1.csv", header = T)
NationalScenarios2 <- read.csv("NationalScen_DCIi_MultiRun2.csv", header = T)
NationalScenarios3 <- read.csv("NationalScen_DCIi_MultiRun3.csv", header = T)
NationalScenarios4 <- read.csv("NationalScen_DCIi_MultiRun4.csv", header = T)
NationalScenarios5 <- read.csv("NationalScen_DCIi_MultiRun5.csv", header = T)

## Merge all different runs in one single dataset
NationalScenarios <- rbind(NationalScenarios1, NationalScenarios2, NationalScenarios3, NationalScenarios4,
                           NationalScenarios5)

Best <- psel(NationalScenarios, pref = high(NationalScenarios$AddCapacity) * high(NationalScenarios$NatAverageDCI))
Worst <- psel(NationalScenarios, pref = low(NationalScenarios$AddCapacity) * low(NationalScenarios$NatAverageDCI))


## Plot 6 (DCI loss Vs Capacity)
tiff(filename = "SI_Appendix_Figure6.tiff", height = 2396, width = 3800, res = 300, compression = c("lzw"))
par(oma = c(8, 8.5, 7, 2), mar = c(0.5, 0, 0, 0), bty = "n")

## Create a matrix for the layout function
mat <- matrix(c(1, 1, 2, 2,
                1, 1, 0, 0,
                1, 1, 3, 3,
                1, 1, 1, 0), ncol = 4, nrow = 4, byrow = T)

layout(mat = mat, widths = c(3.5, 3.5, 2.5, 2), heights = c(0.8, 0.35, 0.7, 0.25))


## Plot Individual DCIs Vs Capacity
plot(NationalScenarios$AddCapacity, NationalScenarios$NatAverageDCI, type = "n", 
     ylim = c(40, 80), xlim = c(0, 65), ylab = "", xlab = "", xaxt = "n", yaxt = "n")

axis(side = 2, at = c(40, 45, 50, 55, 60, 65, 70, 75, 80), cex.axis = 3.2, line = 0)
axis(side = 1, at = c(0, 10, 20, 30, 40, 50, 60), cex.axis = 3.2, line = 0, mgp = c(3, 1.9, 0))

mtext("Nation-wide river connectivity (DCI)", side = 2, cex = 2.7, line = 4.9)
mtext("Nation-wide capacity gain (gigawatts)", side = 1, cex = 2.7, line = 6.1)

## Plot Panel letters
mtext("A", side = 3, cex = 3.8, line = 1, at = 2)
mtext("B", side = 3, cex = 3.8, line = 1, at = 40)
mtext("C", side = 3, cex = 3.8, line = -26, at = 55)


## Plot Max and min future energy demands
segments(x0 = 15, x1 = 15, y0 = 0, y1 = 80, col = "gray60", lty = 3, lwd = 4)
segments(x0 = 36, x1 = 36, y0 = 0, y1 = 80, col = "gray60", lty = 3, lwd = 4)

## Plot years of demands
mtext("2030", side = 2, cex = 1.7, line = -17.5, at = 78)
mtext("2040", side = 2, cex = 1.7, line = -38.5, at = 78)


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

BestDemand <- Best[Best$AddCapacity/1000 >= DemandLow & Best$AddCapacity/1000 <= DemandHigh, 1:7]
WorstDemand <- Worst[Worst$AddCapacity/1000 >= DemandLow & Worst$AddCapacity/1000 <= DemandHigh, 1:7]

## Organize data for boxplot 1 (numbers)
boxplot(BestDemand$NFutSHP, BestDemand$NFutLHP, WorstDemand$NFutSHP, WorstDemand$NFutLHP,
        col = c("#8B000099", "#4876FF99", "#8B000099", "#4876FF99"),
        ylim = c(0, 1500), at = c(1, 2, 3.3, 4.3), xaxt = "n", yaxt = "n")


axis(side = 2, at = c(0, 500, 1000, 1500), cex.axis = 2.4, line = 0)
axis(side = 1, at = c(1, 2), labels = c("SHP", "LHP"), cex.axis = 2.4, line = 0, mgp = c(3, 1.6, 0))
axis(side = 1, at = c(3.3, 4.3), labels = c("SHP", "LHP"), cex.axis = 2.4, line = 0, mgp = c(3, 1.6, 0))


mtext("# new projects", side = 2, cex = 1.9, line = 4)
mtext("Optimal", side = 3, cex = 2.4, line = 1.5, at = 1.4)
mtext("Least-optimal", side = 3, cex = 2.4, line = 1.5, at = 3.8)


## Organize data for boxplot 2 (free-flowing basins loss)
boxplot(NA, BestDemand$NFreeDammed, WorstDemand$NFreeDammed, NA,
        col = c("gray60", "gray60"),
        ylim = c(0, 500), at = c(1, 1.9, 3.2, 4), xaxt = "n", yaxt = "n")

axis(side = 4, at = c(0, 250, 500), cex.axis = 2.4, line = -7, mgp = c(3, 1.7, 0))
mtext("# hydro-free basins lost", side = 4, cex = 1.9, line = -2, at = 150)


dev.off()
