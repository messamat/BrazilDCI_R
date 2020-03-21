
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