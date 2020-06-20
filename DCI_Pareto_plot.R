###########################################################################################################
################ Plot PARETO and select OPTIMAL and NON-OPTIMAL scenarios #################################
###########################################################################################################

#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT

################################################################################
#Get Pareto calculation output data
getParetodat <- function(DCIname) {
  #Import scenarios
  nationalfiles <- list.files(path=resdir,
                              pattern=paste0("NationalScen_", DCIname, "_.*[.]fst"))
  NationalScenarios <- read.fst(file.path(resdir, nationalfiles[length(nationalfiles)])) %>%
    setDT %>%
    unique 
  return(NationalScenarios)
}

#Make plot of outputs + write optimal scenarios to table
plottabPareto <- function(DCIname, in_scenarios) {
  NationalScenarios <- in_scenarios
  
  Best <- psel(NationalScenarios, 
               pref = high(NationalScenarios$AddCapacity) * high(NationalScenarios$NatAverageDCI))
  print('Number and stats on optimal scenarios:')
  nrow(Best)
  summary(Best)
  
  Worst <- psel(NationalScenarios, 
                pref = low(NationalScenarios$AddCapacity) * low(NationalScenarios$NatAverageDCI))
  print('Number and stats on least-optimal scenarios:')
  nrow(Worst)
  summary(Worst)
  
  ## Define the window of future demands
  CurrentGen_GW <- sum(DamAttributes$POT_KW[which(DamAttributes$ESTAGIO_1 == "Operation")])/1000/1000
  
  DemandLow <- 120 - CurrentGen_GW
  DemandHigh <- 141 - CurrentGen_GW
  
  BestDemand <- Best[Best$AddCapacity/1000 >= DemandLow & Best$AddCapacity/1000 <= DemandHigh, 1:7]
  WorstDemand <- Worst[Worst$AddCapacity/1000 >= DemandLow & Worst$AddCapacity/1000 <= DemandHigh, 1:7]
  
  # Plot 6 (DCI loss Vs Capacity)
  tiff(filename = file.path(figdir, paste0("Figure6_", DCIname, Sys.Date(), ".tiff")),
       height = 2396, width = 3800, res = 300, compression = c("lzw"))
  par(oma = c(8, 8.5, 7, 2), mar = c(0.5, 0, 0, 0), bty = "n")
  
  ## Create a matrix for the layout function
  mat <- matrix(c(1, 1, 2, 2,
                  1, 1, 0, 0,
                  1, 1, 3, 3,
                  1, 1, 1, 0), ncol = 4, nrow = 4, byrow = T)
  
  layout(mat = mat, 
         widths = c(3.5, 3.5, 2.5, 2), 
         heights = c(0.8, 0.35, 0.7, 0.25))
  
  
  ## Plot Individual DCIs Vs Capacity
  mincapfloor <- 5*floor(min(NationalScenarios$NatAverageDCI)/5)
  maxcapceil <- 5*ceiling(max(NationalScenarios$NatAverageDCI)/5)
  
  plot(NationalScenarios$AddCapacity, 
       NationalScenarios$NatAverageDCI, 
       type = "n", 
       ylim = c(mincapfloor, maxcapceil), xlim = c(0, 65), 
       ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  
  if (DCIname == 'DCI')
    c()
  
  axis(side = 2, at = seq(mincapfloor, maxcapceil, 5),
       cex.axis = 3.2, line = 0)
  axis(side = 1, at = c(0, 10, 20, 30, 40, 50, 60),
       cex.axis = 3.2, line = 0, mgp = c(3, 1.9, 0))
  
  mtext("Nation-wide river connectivity (DCI)", side = 2, cex = 2.7, line = 4.9)
  mtext("Nation-wide capacity gain (gigawatts)", side = 1, cex = 2.7, line = 6.1)
  
  ## Plot Panel letters
  mtext("A", side = 3, cex = 3.8, line = 1, at = 2)
  mtext("B", side = 3, cex = 3.8, line = 1, at = 40)
  mtext("C", side = 3, cex = 3.8, line = -26, at = 55)
  
  
  ## Plot Max and min future energy demands
  segments(x0 = 15, x1 = 15, y0 = 0, y1 = 85, col = "gray60", lty = 3, lwd = 4)
  segments(x0 = 36, x1 = 36, y0 = 0, y1 = 85, col = "gray60", lty = 3, lwd = 4)
  
  
  
  ## Plot average estimates
  points(NationalScenarios$AddCapacity/1000, NationalScenarios$NatAverageDCI,
         pch = 21, col = "#42424210", bg = "#9400D310", cex = 1.6)
  
  ## Plot best and worse scenarios
  points(Worst$AddCapacity/1000, Worst$NatAverageDCI, 
         pch = 21, col = "#424242", bg = "white", cex = 2.1)
  points(Best$AddCapacity/1000, Best$NatAverageDCI, 
         pch = 21, col = "black", bg = "#303030", cex = 2.1)
  
  
  ## Select the best and worst scenarios that fall inside thresholds of future demands
  
  ## Organize data for boxplot 1 (numbers)
  boxplot(BestDemand$NFutSHP, BestDemand$NFutLHP, WorstDemand$NFutSHP, WorstDemand$NFutLHP,
          col = c("#8B000099", "#4876FF99", "#8B000099", "#4876FF99"),
          ylim = c(0, 100*ceiling(max(WorstDemand$NFutSHP)/100)), 
          at = c(1, 2, 3.3, 4.3), xaxt = "n", yaxt = "n")
  
  
  axis(side = 2, at = c(0, 500, 1000, 1500), cex.axis = 2.4, line = 0)
  axis(side = 1, at = c(1, 2), labels = c("SHP", "LHP"), cex.axis = 2.4, line = 0, mgp = c(3, 1.6, 0))
  axis(side = 1, at = c(3.3, 4.3), labels = c("SHP", "LHP"), cex.axis = 2.4, line = 0, mgp = c(3, 1.6, 0))
  
  
  mtext("# new projects", side = 2, cex = 1.9, line = 4)
  mtext("Optimal", side = 3, cex = 2.4, line = 1.5, at = 1.4)
  mtext("Least-optimal", side = 3, cex = 2.4, line = 1.5, at = 3.8)
  
  
  ## Organize data for boxplot 2 (free-flowing basins loss)
  boxplot(NA, BestDemand$NFreeDammed, WorstDemand$NFreeDammed, NA,
          col = c("gray60", "gray60"),
          ylim = c(0, 100*ceiling(max(WorstDemand$NFreeDammed)/100)),
          at = c(1, 1.9, 3.2, 4), xaxt = "n", yaxt = "n")
  
  axis(side = 4, at = c(0, 250, 500), cex.axis = 2.4, line = -7, mgp = c(3, 1.7, 0))
  mtext("# hydro-free basins lost", side = 4, cex = 1.9, line = -2, at = 150)
  
  
  dev.off()
  
  #---------------- Write CSV with the dataset -----------------------------------
  #Import scenarios with damIDs and sort dam names
  outdir_permut = file.path(resdir, paste0('outpermut_basins', DCIname))
  allscens <- lapply(file.path(outdir_permut, list.files(outdir_permut)), read.fst) %>% #Read all scenarios for each basin
    do.call(rbind, .) %>% #Join them together
    setDT %>%
    .[!is.na(DamIDs), DamIDs := toString(sort(str_split(string=DamIDs, 
                                                        pattern='[,]\\s*',
                                                        simplify=T))), #Order dam names
      by=.(scenbasin, DAMBAS_ID08ext)] %>%
    setkey(DamIDs) #Set key
  
  #Retrieve dam IDs for each basin-specific scenario
  BestScens <-  Best[, tstrsplit(scenIDs, ',')] %>%
    setnames(DamAttributes[, sort(as.character(unique(DAMBAS_ID08ext)))]) %>%
    .[, optimizedID := .I]
  
  BestScens_damIDs <-  melt(BestScens, id.vars='optimizedID', 
                            variable.name = "DAMBAS_ID08ext", 
                            value.name = "scenbasin") %>%
    .[, scenbasin := as.integer(scenbasin)] %>%
    merge(., allscens[, .(DAMBAS_ID08ext, scenbasin, DamIDs)],
          by=c('DAMBAS_ID08ext', 'scenbasin'), all.y=F) %>%
    .[!is.na(DamIDs), list(DamIDs = toString(DamIDs)),
      by = optimizedID] %>%
    setorder(optimizedID) 
  
  Best_damIDs <- cbind(Best[BestScens_damIDs$optimizedID, -'scenIDs', with=FALSE], 
                       BestScens_damIDs)
  
  ## Write a csv with the dataset
  BestDemandPrint <- Best_damIDs[Best_damIDs$AddCapacity/1000 >= DemandLow & 
                                   Best_damIDs$AddCapacity/1000 <= DemandHigh,]
  OptimalDamPrint <- BestDemandPrint[,-'optimizedID', with=F] %>%
    setorder('AddCapacity')
  colnames(OptimalDamPrint) <- c("National average DCI", "Capacity gain (GW)", "N future dams",
                                 "N future SHPs", "N future LHPs", "N free basins lost", "DamIDs")
  write.csv(OptimalDamPrint, 
            file = file.path(resdir, 
                             paste0("Supplement_OptimalPortifolios", DCIname, ".csv")))
  
}

#Get statistics on optimal and least-optimal scenarios
getParetostats <- function(DCIname, in_scenarios) {
  NationalScenarios <- in_scenarios
  
  Best <- psel(NationalScenarios, 
               pref = high(NationalScenarios$AddCapacity) * high(NationalScenarios$NatAverageDCI))
  message('Number and stats on optimal scenarios:')
  message(nrow(Best))
  message(summary(Best))
  
  Worst <- psel(NationalScenarios, 
                pref = low(NationalScenarios$AddCapacity) * low(NationalScenarios$NatAverageDCI))
  message(('Number and stats on least-optimal scenarios:'))
  message(nrow(Worst))
  message(summary(Worst))
  
  ## Define the window of future demands
  CurrentGen_GW <- sum(DamAttributes$POT_KW[which(DamAttributes$ESTAGIO_1 == "Operation")])/1000/1000
  
  DemandLow <- 120 - CurrentGen_GW
  DemandHigh <- 141 - CurrentGen_GW
  
  BestDemand <- Best[Best$AddCapacity/1000 >= DemandLow & Best$AddCapacity/1000 <= DemandHigh, 1:7]
  WorstDemand <- Worst[Worst$AddCapacity/1000 >= DemandLow & Worst$AddCapacity/1000 <= DemandHigh, 1:7]
  
  ## Calculate some statistics
  message(('Mean number of SHP - optimal:'))
  message(round(mean(BestDemand$NFutSHP), digits = 0))
  message(('Mean number of LHP - optimal:'))
  message(round(mean(BestDemand$NFutLHP), digits = 0))
  
  message(('SD number of SHP - optimal:'))
  message(round(sd(BestDemand$NFutSHP), digits = 0))
  message(('SD number of LHP - optimal:'))
  message(round(sd(BestDemand$NFutLHP), digits = 0))
  
  message(('Mean number of SHP - least-optimal:'))
  message(round(mean(WorstDemand$NFutSHP), digits = 0))
  message(('Mean number of LHP - least-optimal:'))
  message(round(mean(WorstDemand$NFutLHP), digits = 0))
  
  message(('SD number of SHP - least-optimal:'))
  message(round(sd(WorstDemand$NFutSHP), digits = 0))
  message(('SD number of LHP - least-optimal:'))
  message(round(sd(WorstDemand$NFutLHP), digits = 0))
  
  ## Find an example to illustrate added capacity
  print(round(BestDemand[order(BestDemand$AddCapacity),-'scenIDs'], digits = 0)) 
  ## 75 DCI      28071 GW      547 SHPs      74 LHPs         233 NFreeDammed
  
  print(round(WorstDemand[order(WorstDemand$AddCapacity),-'scenIDs'], digits = 0))
  # 67 DCI       27966 GW    1346 SHPs    146 LHPs        441 NFreeDammed
  
  # (1346-547)/1346
  # (146-74)/146
  # (441-233)/441
  
  ## Free-flowing statistics
  message('mean(BestDemand$NFreeDammed)')
  message(round(mean(BestDemand$NFreeDammed), digits = 0))
  message('mean(WorstDemand$NFreeDammed)')
  message(round(mean(WorstDemand$NFreeDammed), digits = 0))
  
  message('sd for best and worst NFreeDammed, respectively')
  message(round(sd(BestDemand$NFreeDammed), digits = 0))
  message(round(sd(WorstDemand$NFreeDammed), digits = 0))
  
}

#Check whether Pareto frontier approximation is good
checkPareto <- function(DCIname, in_scenarios) {
  NationalScenarios <- in_scenarios
  
  #Assess whether reaching Pareto frontier
  nscen_seq <- c(seq(1, nrow(NationalScenarios), 100), nrow(NationalScenarios))
  growfrontier_dt <- lapply(nscen_seq,
                            function(nscen) {
                              print(nscen)
                              subscens <- NationalScenarios[sample(.N, nscen), -'scenIDs', with=F]
                              sub_best <- psel(subscens,
                                               pref = high(subscens$AddCapacity) *
                                                 high(subscens$NatAverageDCI)) %>%
                                .[, nscen := nscen]
                              return(sub_best)
                            }
  ) %>%
    rbindlist
  growfrontier_dt[, AddCapacity_int := round(AddCapacity)]
  
  growfrontierdiff <- growfrontier_dt[
    , DCIdiff := (max(NatAverageDCI) - NatAverageDCI)/max(NatAverageDCI),
    by=AddCapacity_int] %>%
    .[, list(meandcidiff = mean(DCIdiff)), by=nscen]
  
  ggplot(growfrontierdiff, aes(x=nscen, y=meandcidiff)) +
    geom_point()
  
}

for (i in c('DCIp', 'DCIi')) {
  scenarios <- getParetodat(i)
  #plottabPareto(i, scenarios)
  getParetostats(i, scenarios)
  #checkPareto(i, scenarios)
}