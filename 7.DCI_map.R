#############         Script 6 to map DCIs in Brazil - September 2019              ###############################
#############       Based on distances calculated in Arcgis - PythonCode           ###############################
#############          Based on DCIs from script 1 and folder "Shapes"             ###############################
#############              Manuscript in preparation - PNAS                        ###############################


###############################################################################################################
#########                 MAP DCIs - BRAZIL                  #################################################
##########      (CURRENT VS FUTURE) -  SHPs, LHPs, blend    #################################################
###########               Level 8 All Scenarios             ################################################
#########################################################################################################

#Import packages
source('00.DCI_packages.R')
source('00.DCI_functions.R')

## Import network and dams dataset (alternative)
datadir <- file.path(rootdir, "data")
SouthAmerDir <- file.path(datadir, "South America")
WatershedsDir <- file.path(datadir, "Watersheds ANA")
BasinsDir <- file.path(datadir, "Python R Outputs")

## Merge DCI colors in the Basin Shapefile
BasinsShape <- readOGR(dsn = dcigdb, layer = "BasinATLAS_lev08_v10_br")


mapDCI <- function(DCIname, in_bas) {
  # Read table for the future projections per basin (absolut values and percentage change)
  resuDCI <- as.data.frame(
    read.csv(file.path(resdir,
                       paste0("DCI_Brazil_L8_", DCIname, ".csv")), header = T))
  PercChangeData <- as.data.frame(
    read.csv(file.path(resdir, 
                       paste0("DCI_Brazil_L8_Percentage", DCIname, ".csv"))))
  
  ## Assign breaks to values of DCI (100 - DCI)
  CurSHPCut <- cut(100 - resuDCI$SHP_curr, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
  CurLHPCut <- cut(100 - resuDCI$LHP_curr, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
  
  FutSHPCut <- cut(100 - resuDCI$SHP_fut, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
  FutLHPCut <- cut(100 - resuDCI$LHP_fut, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
  
  PercSHPCut <- cut(-1 * PercChangeData$PercLoss_SHP, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
  PercLHPCut <- cut(-1 * PercChangeData$PercLoss_LHP, breaks = c(-1, 0, 25, 50, 100), labels = c("d", "c", "b", "a"))
  
  
  ## Create a factor with color names
  if (DCIname == 'DCIp') {
    CurrColors <- factor(paste(CurSHPCut, CurLHPCut, sep = ""), 
                         labels = c("#756BB1", "#A96B8C", "#DD6B66",
                                    "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                    "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                    "#3182BD", "#9ECAE1", "#DEEBF7", "white"))
  } else if (DCIname == 'DCIi') {
    CurrColors <- factor(paste(CurSHPCut, CurLHPCut, sep = ""), 
                         labels = c("#756BB1", "#A96B8C", "#DD6B66", "#DE2D26",
                                    "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                    "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                    "#3182BD", "#9ECAE1", "#DEEBF7", "white"))
  }

  
  FutuColors <- factor(paste(FutSHPCut, FutLHPCut, sep = ""), 
                       labels = c("#756BB1", "#A96B8C", "#DD6B66", "#DE2D26",
                                  "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                  "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                  "#3182BD", "#9ECAE1", "#DEEBF7"))
  
  PercColors <- factor(paste(PercSHPCut, PercLHPCut, sep = ""), 
                       labels = c("#756BB1", "#A96B8C", "#DD6B66", "#DE2D26",
                                  "#9198C8", "#BCBDDC", "#DBA8A6", "#FC9272",
                                  "#6FA3CE", "#ADC4DE", "#EFEDF5", "#FEE0D2",
                                  "#3182BD", "#9ECAE1", "#DEEBF7", "white"))
  
  Colvec <- c("#DE2D26", "#DD6B66", "#A96B8C", "#756BB1",
              "#FC9272", "#DBA8A6", "#BCBDDC", "#9198C8",
              "#FEE0D2", "#EFEDF5", "#ADC4DE", "#6FA3CE",
              "white", "#DEEBF7", "#9ECAE1", "#3182BD")
  
  ## Create a dataframe with the colors for all the scenarios for each basin
  DCIColors <- data.frame(resuDCI$DAMBAS_ID, CurrColors, FutuColors, PercColors)
  colnames(DCIColors) <- c("HYBAS_ID", "CurrColors", "FutuColors", "PercColors")
  
  
  ## Import shapefiles
  ## Use the function "View" to view the columns of the shapefiles
  SouthAmerShape <- readOGR(dsn = SouthAmerDir, layer = "South_America")
  WatershedsShape <- readOGR(dsn = WatershedsDir, layer = "Regiões_Hidrográficas")
  #BasinsShape <- readOGR(dsn = dcigdb, layer = "BasinATLAS_lev08_v10_br")
  
  ## Filter basins that have or will have at least on dam
  BasinDams <- in_bas[which(in_bas$HYBAS_ID %in% resuDCI$DAMBAS_ID),]
  DCImapBasins <- merge(BasinDams, DCIColors, by="HYBAS_ID", duplicateGeoms=TRUE)
  
  ## Plot figure - use zm() to zoom in the map
  tiff(filename = file.path(figdir, paste0("Figure3", DCIname, ".tiff")),
       height = 1926, width = 4400, res = 300, compression = c("lzw"))
  
  par (mfrow = c(1,3), oma = c(0, 2, 0, 2), mar = c(0, 0, 0, 0), bty = "n")
  
  #plot 1
  BrazilLayer <- SouthAmerShape$COUNTRY == "Brazil"
  plot(SouthAmerShape[BrazilLayer, ], col = "white", border = F, add = F)
  mtext("Present", side = 3, line = -5, cex = 2.8)
  mtext("A", side = 3, line = -4.5, at = -73, cex = 3.4)
  plot(DCImapBasins, col = as.character(DCImapBasins$CurrColors), add = T, border= F)
  plot(WatershedsShape, col = "#FFFFFF00", border= T, lwd = 1.5, add = T)
  
  ## Plot South America map
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
  mtext("From LHPs", side = 1, cex = 1.4, line = -4.5, at = -68.5)
  
  mtext("1", side = 2, cex = 1.4, line = 0.2, at = -30.4)
  mtext("25", side = 2, cex = 1.4, line = 0.2, at = -26.85)
  mtext("50", side = 2, cex = 1.4, line = 0.2, at = -22.9)
  mtext("100", side = 2, cex = 1.4, line = 0.2, at = -19.0)
  mtext("From SHPs", side = 2, cex = 1.4, line = 2.2, at = -26.85)
  
  
  dev.off()
  
  
  #######################################  Create a jpeg of the color scale   #########################################
  
  ## Create a jpeg of the color scale (1000 x 1000)
  par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), bty = "n")
  ## Create a matrix of colors in their positions
  matColVec <- matrix(Colvec, nrow = 4, ncol = 4, byrow = T)
  ## Plot an empty space
  plot(0:4, -4:0, xaxt = "n", yaxt = "n", type = "n", ylab = "", xlab = "", bty = "n", xlim = c(0,4), ylim = c(-4,0))
  ## Fill the plot with 4 x 4 rectangles
  # Loop over rows
  for( i in 1:4) {
    # Loop over columns
    for( j in 1:4) {
      rect(xleft = i - 1, xright = i, ybottom = - (j - 1), ytop = -j, col = matColVec[j,i]) 
    }
  }
  
  
  ###################################### Create a jpeg for major watershed codes ###########################################
  #plot
  ## Create a jpeg of the color scale (1200 x 1100)
  jpeg(filename = "SA.jpeg", height = 1500, width = 1300, quality = 100)
  
  par(mfrow = c(1,1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), bty = "n", bg = "white")
  BrazilLayer <- SouthAmerShape$COUNTRY == "Brazil"
  plot(SouthAmerShape, col = "gray20", border = T, xlim=c(-95.349, -26.585), ylim=c(-54.614, 8.705), add = F)
  plot(SouthAmerShape[BrazilLayer, ], col = "white", lwd = 4, border = T, add = T)
  plot(WatershedsShape, col = "#FFFFFF00", border= T, lwd = 2.5, add = T)
  
  dev.off()
  
  ########################################### FIND Color Pallets ##########################################################
  ## ColBrewer pallet's choice:
  ## Palettes names: Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd
  
  ## Create a color pallet for red and blue
  par(mfrow = c(1,1))
  brewer.pal(n = 3, name = "Reds")
  display.brewer.pal(n = 3, name = "Blues")
  display.brewer.pal(n = 3, name = "Reds")
  display.brewer.all()
  
  colpallete1 <- brewer.pal(n = 3, name = "Reds")
  colpallete2 <- brewer.pal(n = 3, name = "Blues")
  colpallete3 <- brewer.pal(n = 3, name = "Purples")
  
  ## Plot the legend
  par(mfrow = c(4,4), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), bty = "o")
  
  for (i in 1: length(Colvec)){
    
    plot(0, 0, type = "n", xaxt = "n", yaxt = "n")
    box(which = "plot", lwd = 0.3)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
         col = Colvec[i])
  }
  
  
  ##############################################################################################################
  ##########################   TAPAJOS and JURUENA  ############################################################
  ##############################################################################################################
  
  ## Upload the data
  PercChangeData <- as.data.frame(read.csv(
    file.path(resdir, paste0("DCI_Brazil_L8_Percentage", DCIname, ".csv"))))
  
  #Get HYBAS data level 08 and merge with DCI results
  basinatlas_l8merge <- merge(
    st_read(dsn = dcigdb, layer = 'BasinATLAS_lev08_v10_brproj')[, c('HYBAS_ID', 'PFAF_ID')],
    PercChangeData,
    by='HYBAS_ID',
    all.x=T
  ) %>%
    setDT
  
  #Subset to keep only data for Tapajos and Juruena basins
  #Tapajos HydroBASINS level 4: HYBAS_ID 6040262220 PFAF_ID 6224
  #Juruena HydroBASINS level 5: HYBAS_ID 6050405310 PFAF_ID 62244 
  tapajosbasins <- basinatlas_l8merge[substr(as.character(PFAF_ID), 1, 4)=='6224',]
  juruenabasins <- basinatlas_l8merge[substr(as.character(PFAF_ID), 1, 5)=='62244',] 
  buritibasin <- basinatlas_l8merge[substr(as.character(PFAF_ID), 1, 4)=='6224',]
  cuparileste <- basinatlas_l8merge[PFAF_ID=='62241404',] 

  
  ## Get the number of basins without future dams
  JuruenaFree <- juruenabasins[is.na(PercLoss_All), .N]
  TapajosFree <- tapajosbasins[is.na(PercLoss_All), .N]
  
  ## Filter results for basins inside Juruena and Tapajos
  JuruenaPercDam <- juruenabasins[!is.na(PercLoss_All),]
  TapajosPercDam <- tapajosbasins[!is.na(PercLoss_All),]
  
  message("Order")
  message(JuruenaPercDam[order(JuruenaPercDam$PercLoss_SHP), ]) 
  message(TapajosPercDam[order(TapajosPercDam$PercLoss_SHP), ]) 
  
  message("Number of basins that will decrease connectivity in more than 50%")
  message(sum(JuruenaPercDam$PercLoss_SHP < -50))
  message(sum(TapajosPercDam$PercLoss_SHP < -50))
  
  message("Get mean and range")
  message( mean(JuruenaPercDam$PercLoss_All))
  message(mean(JuruenaPercDam$PercLoss_SHP))
  message(mean(JuruenaPercDam$PercLoss_LHP))
  
  message("Overall average for Juruena")
  message( sum(JuruenaPercDam$PercLoss_SHP)/(length(JuruenaPercDam$PercLoss_SHP) + JuruenaFree))
}

mapDCI(DCIname = 'DCIp', in_bas = BasinsShape)
mapDCI(DCIname = 'DCIi', in_bas = BasinsShape)
