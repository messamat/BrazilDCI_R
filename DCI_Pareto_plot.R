###########################################################################################################
################ Plot PARETO and select OPTIMAL and NON-OPTIMAL scenarios #################################
###########################################################################################################

## Load packages
library(magrittr)
library(rPref)
library(data.table)
library(fst)
library(rprojroot)

# # Import network and dams dataset (Mathis folder structure)
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')
outdir_permut = file.path(resdir, 'outpermut_basins')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT

#Import scenarios
nationalfiles <- list.files(path=resdir, pattern="NationalScen_DCIp_.*[.]fst")
NationalScenarios <- read.fst(file.path(resdir, nationalfiles[length(nationalfiles)])) %>%
  setDT %>%
  unique #Only keep unique scenarios (mostly removes the multiple scenarios with all dams or no dams)

Best <- psel(NationalScenarios, pref = high(NationalScenarios$AddCapacity) * high(NationalScenarios$NatAverageDCI))
Worst <- psel(NationalScenarios, pref = low(NationalScenarios$AddCapacity) * low(NationalScenarios$NatAverageDCI))

## Plot 6 (DCI loss Vs Capacity)
tiff(filename = file.path(resdir, "Figure6.tiff"), height = 2396, width = 3500, res = 300, compression = c("lzw"))
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
