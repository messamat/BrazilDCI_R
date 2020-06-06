#Purpose: PLOT DCI Loss individual future dams VS Capacity 

#Set folder structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
dcigdb <- file.path(resdir, 'dci.gdb')

#Import data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT

RankedDams <- read.fst(file.path(resdir, 
                                 list.files(path=resdir, pattern='SamplingIndividualDams_results.*[.]fst'))) %>%
  setDT %>%
  setorder('DCIMeanDiff') %>%
  .[, DamRank := .I] #Organize the data to plot

#Correct a small glitch
RankedDams[DCIDownLim<0, DCIDownLim := 0]
RankedDams[DCIDownCI<0, DCIDownCI := 0]

## Plot 5 (DCI loss Vs Capacity)
tiff(filename = file.path(resdir, "Figure5.tiff"), height = 2396, width = 3700, res = 300, compression = c("lzw"))
par(oma = c(6, 9, 2, 0.5), mar = c(0.5, 0, 0, 0), bty = "n")


## Create a color vector to differentiate SHP and LHP
TypeFac <- factor(x = RankedDams$Type)
TypeColor <- factor(x = TypeFac, levels = levels(TypeFac), 
                    labels = c("#4876FF50", "#8B000050"))

## Plot figure #xaxt = "n", log(RankedDams$Capacity)
plot(RankedDams$Capacity, -RankedDams$DCIMeanDiff, type = "n", ylim = c(-60, 0), xlim = c(0.1, 8000),
     ylab = "", xlab = "",  yaxt = "n", xaxt = "n", log = "x") 

options(scipen=5)
axis(side = 2, at = c(-60, -40, -20, 0), cex.axis = 2.6, line = 0)
axis(side = 1, at = c(1, 10, 30, 100, 500, 4000), cex.axis = 2.6, line = 0, mgp = c(3, 1.8, 0))
axis(side = 1, at = c(0.1, 1), labels = c(0.1, " "), cex.axis = 2.6, line = 0, mgp = c(3, 1.8, 0))

mtext("Change in basin-level", side = 2, cex = 3.2, line = 6.8)
mtext("river connectivity (DCI)", side = 2, cex = 3.2, line = 4.0)
mtext("Generation capacity (megawatts)", side = 1, cex = 3.2, line = 4.9)


## Plot range log(RankedDams$Capacity)
segments(y0 = -RankedDams$DCIUppLim, x0 = RankedDams$Capacity,
         y1 = -RankedDams$DCIDownLim, x1 = RankedDams$Capacity, lwd = 1.6, col = "#42424220")

# ## Plot 95% confidence intervals
# segments(y0 = RankedDams$DCIUppLim, x0 = log(RankedDams$Capacity), 
#          y1 = RankedDams$DCIDownCI, x1 = log(RankedDams$Capacity), lwd = 1.3, col = "#42424220")

## Plot average log(RankedDams$Capacity)
points(RankedDams$Capacity, -RankedDams$DCIMeanDiff, pch = 21, col = "white", 
       bg = "white", cex = 2.5)
points(RankedDams$Capacity, -RankedDams$DCIMeanDiff, pch = 21, col = "#42424220", 
       bg = as.vector(TypeColor), cex = 2.5)

## Plot legend
legend(x = 200, y = -47, pch = c(21), legend = c("Planned SHP", "Planned LHP"), cex = 2.5,
       pt.cex = 3.3, pt.bg = c("#8B000060","#4876FF60"), xpd = T, bty = "n")

dev.off()


################################################################################################################
#####################################  Some Statistics  ########################################################
RankedDams$Capacity
RankedDams$DCIMeanDiff


## Linear models effect on DCI Vs Capacity
plot(RankedDams$Capacity, RankedDams$DCIMeanDiff)

#Add 0.001 to dam with 0 impact on DCI and log it
RankedDams[DCIMeanDiff == 0, DCIMeanDiff := DCIMeanDiff+0.001]

## All
overallLm <- lm(log(RankedDams$DCIMeanDiff) ~ log(RankedDams$Capacity))
summary(overallLm)
hist(overallLm$resid, main="Histogram of Residuals",
     ylab="Residuals")
qqnorm(overallLm$resid)
qqline(overallLm$resid)
plot(RankedDams$Capacity, RankedDams$DCIMeanDiff, log='xy')
abline(overallLm)

cor.test(x=log(RankedDams$Capacity), y=log(RankedDams$DCIMeanDiff), method = 'pearson')

## Just SHPs
SHPLm <- lm(log(RankedDams$DCIMeanDiff[RankedDams$Type == "SHP"]) ~ log(RankedDams$Capacity[RankedDams$Type == "SHP"]))
summary(SHPLm)
hist(SHPLm$resid, main="Histogram of Residuals",
     ylab="Residuals")
qqnorm(SHPLm$resid)
qqline(SHPLm$resid)
plot(log(RankedDams$Capacity[RankedDams$Type == "SHP"]), log(RankedDams$DCIMeanDiff[RankedDams$Type == "SHP"]))
abline(SHPLm)

cor.test(x=log(RankedDams$Capacity[RankedDams$Type == "SHP"]), y=log(RankedDams$DCIMeanDiff[RankedDams$Type == "SHP"]), method = 'pearson')

## Just LHPs
LHPLm <- lm(log(RankedDams$DCIMeanDiff[RankedDams$Type == "LHP"]) ~ log(RankedDams$Capacity[RankedDams$Type == "LHP"]))
summary(LHPLm)
hist(LHPLm$resid, main="Histogram of Residuals",
     ylab="Residuals")
qqnorm(LHPLm$resid)
qqline(LHPLm$resid)
plot(log(RankedDams$Capacity[RankedDams$Type == "LHP"]), log(RankedDams$DCIMeanDiff[RankedDams$Type == "LHP"]))
abline(LHPLm)

cor.test(x=log(RankedDams$Capacity[RankedDams$Type == "LHP"]), y=log(RankedDams$DCIMeanDiff[RankedDams$Type == "LHP"]), method = 'pearson')

# ## Create a csv of future dams rank basend on mean DCI
# ## Add columns with lat-long
# ## Created in ArcGIS two columns with Lat and long of the dams (Decimal degree, WGS 1984)
# FullDamAttrubutes <- read.csv("DamAttributesCoordinates.txt", header = T)
# 
# Lat <- vector()
# Long <- vector()
# 
# for (i in 1: dim(OrderDams)[1]){
#   IDName_Position <- which(FullDamAttrubutes$TARGET_FID == RankedDams$ID[i] & FullDamAttrubutes$NOME %in% RankedDams$Name[i])
#   DamLatLong <- FullDamAttrubutes[IDName_Position, 20:21]
#   Lat <- c(Lat, FullDamAttrubutes$Lat[IDName_Position])
#   Long <- c(Long, FullDamAttrubutes$Long[IDName_Position])
# 
# }
# 
## Organize order and headers
numcols <- names(RankedDams)[sapply(RankedDams, is.numeric)]
OrderDams <- RankedDams[, (numcols) := lapply(.SD, function(x) round(x, digits = 2)), 
                        .SDcols=numcols] %>%
  .[, .(DamRank, Type, Name, DAMID, Capacity, DCIMeanDiff, DCIUppLim, DCIDownLim)] %>%
  merge(DamAttributes[,.(DAMID, NEAR_X, NEAR_Y)], by='DAMID', all.y=F) %>%
  setorder(DCIMeanDiff) %>%
  setnames(c("Rank", "Type", "Name",  "DamID", "Capacity(MW)",
             "Mean effect on basin's DCI", 
             "Upper limit", "Lower limit",
             "X_Sirgas", "Y_Sirgas"))
write.csv(OrderDams, file = file.path(resdir, "Supplement_FutureDamRank.csv"))
