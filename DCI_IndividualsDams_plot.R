#Purpose: PLOT DCI Loss individual future dams VS Capacity 


#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% 
  setDT 


plotindivdams_util <- function(RankedDams, DCIname, prodmeasure) {
  tiff(filename = file.path(figdir, 
                            paste0("Figure5", DCIname, "_", 
                                   prodmeasure, "_", 
                                   format(Sys.Date(), "%Y%m%d"),  ".tiff")),
       height = 2396, width = 3700, res = 300, compression = c("lzw"))
  par(oma = c(6, 9, 2, 0.5), mar = c(0.5, 0, 0, 0), bty = "n")
  RankedDams <- as.data.frame(RankedDams)
  
  ## Create a color vector to differentiate SHP and LHP
  TypeFac <- factor(x = RankedDams$Type)
  TypeColor <- factor(x = TypeFac, levels = levels(TypeFac), 
                      labels = c("#4876FF50", "#8B000050"))
  
  maxdiffceil <- -20*ceiling(max(RankedDams$DCIMeanDiff)/20)
  maxcapaceil <- 1000*ceiling(max(RankedDams[, prodmeasure])/1000)
  ## Plot figure #xaxt = "n", log(RankedDams[, prodmeasure])
  plot(RankedDams[, prodmeasure], -RankedDams$DCIMeanDiff, type = "n", 
       ylim = c(maxdiffceil, 0), 
       xlim = c(0.1, maxcapaceil),
       ylab = "", xlab = "",  yaxt = "n", xaxt = "n", log = "x") 
  
  options(scipen=5)
  axis(side = 2, at = seq(maxdiffceil, 0, 20), 
       cex.axis = 2.6, line = 0)
  axis(side = 1, at = c(1, 10, 30, 100, 500, maxcapaceil), 
       cex.axis = 2.6, line = 0, mgp = c(3, 1.8, 0))
  axis(side = 1, at = c(0.1, 1), labels = c(0.1, " "), 
       cex.axis = 2.6, line = 0, mgp = c(3, 1.8, 0))
  
  mtext("Change in basin-level", side = 2, cex = 3.2, line = 6.8)
  mtext("river connectivity (DCI)", side = 2, cex = 3.2, line = 4.0)
  if (prodmeasure == 'Capacity') {
    mtext("Generation capacity (megawatts)", side = 1, cex = 3.2, line = 4.9)
  } else if (prodmeasure == 'Guarantee') {
    mtext("Physical guarantee (megawatts)", side = 1, cex = 3.2, line = 4.9)
  }
  
  
  ## Plot range log(RankedDams[, prodmeasure])
  segments(y0 = -RankedDams$DCIUppLim, x0 = RankedDams[, prodmeasure],
           y1 = -RankedDams$DCIDownLim, x1 = RankedDams[, prodmeasure], lwd = 1.6, col = "#42424220")
  
  # ## Plot 95% confidence intervals
  # segments(y0 = RankedDams$DCIUppLim, x0 = log(RankedDams[, prodmeasure]), 
  #          y1 = RankedDams$DCIDownCI, x1 = log(RankedDams[, prodmeasure]), lwd = 1.3, col = "#42424220")
  
  ## Plot average log(RankedDams[, prodmeasure])
  points(RankedDams[, prodmeasure], -RankedDams$DCIMeanDiff, pch = 21, col = "white", 
         bg = "white", cex = 2.5)
  points(RankedDams[, prodmeasure], -RankedDams$DCIMeanDiff, pch = 21, col = "#42424220", 
         bg = as.vector(TypeColor), cex = 2.5)
  
  ## Plot legend
  legend(x = 0.1, y = 0.8*maxdiffceil, pch = c(21), legend = c("Planned SHP", "Planned LHP"), cex = 2.5,
         pt.cex = 3.3, pt.bg = c("#8B000060","#4876FF60"), xpd = T, bty = "n")
  
  dev.off()
}

plotindivdams <- function(DamAttributes, DCIname, prodmeasure) {
  #Import data
  indivfiles <- list.files(path=resdir,
                           pattern=paste0('SamplingIndividualDams_results_', 
                                          DCIname, '.*[.]fst'))
  RankedDams <- read.fst(file.path(resdir, indivfiles[length(indivfiles)])) %>%
    setDT %>%
    setorder('DCIMeanDiff') %>%
    .[, DamRank := .I] #Organize the data to plot

  
  #Correct a small glitch for dams in basins for which they were sampled
  RankedDams[DCIDownLim<0, DCIDownLim := 0]
  RankedDams[DCIDownCI<0, DCIDownCI := 0]
  
  
  ## Plot 5 (DCI loss Vs Capacity)
  plotindivdams_util(RankedDams, DCIname, prodmeasure)
  
  ################################################################################################################
  #####################################  Some Statistics  ########################################################
  message("RankedDams[, prodmeasure]")
  message(RankedDams[, get(prodmeasure)])
  message("RankedDams$DCIMeanDiff")
  message(RankedDams$DCIMeanDiff)
  
  ## Linear models effect on DCI Vs Capacity
  plot(RankedDams[, get(prodmeasure)], RankedDams$DCIMeanDiff)
  
  #Add 0.001 to dam with 0 impact on DCI and log it
  RankedDams[DCIMeanDiff == 0, DCIMeanDiff := DCIMeanDiff+0.001]
  
  ## All
  overallLm <- lm(log(RankedDams$DCIMeanDiff) ~ log(RankedDams[, get(prodmeasure)]))
  message(summary(overallLm))
  hist(overallLm$resid, main="Histogram of Residuals",
       ylab="Residuals")
  qqnorm(overallLm$resid)
  qqline(overallLm$resid)
  plot(RankedDams[, get(prodmeasure)], RankedDams$DCIMeanDiff, log='xy')
  abline(overallLm)
  
  message("cor.test(x=log(RankedDams[, get(prodmeasure)]), y=log(RankedDams$DCIMeanDiff), method = 'pearson')")
  message(cor.test(x=log(RankedDams[, get(prodmeasure)]), y=log(RankedDams$DCIMeanDiff), method = 'pearson'))
  
  ## Just SHPs
  SHPLm <- lm(log(RankedDams$DCIMeanDiff[RankedDams$Type == "SHP"]) ~ log(RankedDams[, get(prodmeasure)][RankedDams$Type == "SHP"]))
  message(summary(SHPLm))
  hist(SHPLm$resid, main="Histogram of Residuals",
       ylab="Residuals")
  qqnorm(SHPLm$resid)
  qqline(SHPLm$resid)
  plot(log(RankedDams[, get(prodmeasure)][RankedDams$Type == "SHP"]), log(RankedDams$DCIMeanDiff[RankedDams$Type == "SHP"]))
  abline(SHPLm)
  
  message("cor.test(x=log(RankedDams[, get(prodmeasure)][RankedDams$Type == 'SHP']), y=log(RankedDams$DCIMeanDiff[RankedDams$Type == 'SHP']), method = 'pearson')")
  message(cor.test(x=log(RankedDams[, get(prodmeasure)][RankedDams$Type == 'SHP']), y=log(RankedDams$DCIMeanDiff[RankedDams$Type == 'SHP']), method = 'pearson'))
  
  ## Just LHPs
  LHPLm <- lm(log(RankedDams$DCIMeanDiff[RankedDams$Type == "LHP"]) ~ log(RankedDams[, get(prodmeasure)][RankedDams$Type == "LHP"]))
  message(summary(LHPLm))
  hist(LHPLm$resid, main="Histogram of Residuals",
       ylab="Residuals")
  qqnorm(LHPLm$resid)
  qqline(LHPLm$resid)
  plot(log(RankedDams[, get(prodmeasure)][RankedDams$Type == "LHP"]), log(RankedDams$DCIMeanDiff[RankedDams$Type == "LHP"]))
  abline(LHPLm)
  
  message("cor.test(x=log(RankedDams[, get(prodmeasure)][RankedDams$Type == 'LHP']), y=log(RankedDams$DCIMeanDiff[RankedDams$Type == 'LHP']), method = 'pearson')")
  message(cor.test(x=log(RankedDams[, get(prodmeasure)][RankedDams$Type == 'LHP']), y=log(RankedDams$DCIMeanDiff[RankedDams$Type == 'LHP']), method = 'pearson'))
  
  ## Create a csv of future dams rank basend on mean DCI
  RankedDams_format <- RankedDams[DamAttributes[, .(DAMID, POINT_X, POINT_Y)],
                                  on='DAMID']
  
  ## Organize order and headers
  numcols <- names(RankedDams_format)[sapply(RankedDams_format, is.numeric)]
  OrderDams <- RankedDams_format[
    , (numcols) := lapply(.SD, function(x) round(x, digits = 2)), 
    .SDcols=numcols] %>%
    .[, .(DamRank, Type, Name, DAMID, Capacity, Guarantee, 
          DCIMeanDiff, DCIUppLim, DCIDownLim, 
          POINT_X, POINT_Y)] %>%
    setorder(DCIMeanDiff) %>%
    setnames(c("Rank", "Type", "Name",  "DamID", 
               "Capacity", "Guarantee",
               "Mean effect on basin's DCI", 
               "Upper limit", "Lower limit",
               "Lat_WGS84", "Lon_WGS84"))
  write.csv(OrderDams[!is.na(Rank),], 
            file = file.path(resdir,
                             paste0("Supplement_FutureDamRank", DCIname, ".csv")))
}


plotindivdams(DamAttributes = DamAttributes,
              DCIname = 'DCIp',
              prodmeasure = 'Capacity')

plotindivdams(DamAttributes = DamAttributes,
              DCIname = 'DCIp',
              prodmeasure = 'Guarantee')

plotindivdams(DamAttributes = DamAttributes,
              DCIname = 'DCIi',
              prodmeasure = 'Capacity')

plotindivdams(DamAttributes = DamAttributes,
              DCIname = 'DCIi',
              prodmeasure = 'Guarantee')