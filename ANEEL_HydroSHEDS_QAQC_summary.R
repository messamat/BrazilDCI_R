#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

#Import formatted data
DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT


######################## CHECK DAMS' REPORTED DRAINAGE AREA VS. HYDROSHEDS NETWORK DA #################
#AREA_DREN is the drainage area reported in ANEEL database
#UPLAND_SKM is the drainage area at the pourpoint of the reach that each dam is associated with in the RiverATLAS network
#up_area_skm_15s is the drainage area at the 15 s pixel from the HydroSHEDS accumulation area raster

DamAttributes[is.na(manualsnap), manualsnap := '-9999'] %>% #Change NA to -9999 for dams that were no checked for locational accuracy
  .[is.na(AREA_DREN), AREA_DREN := 0] #For records without AREA_DREN in ANEEL database, change it to 0 to display on graph

#Compute absolute difference in ratio between pixel-based and reach-based drainage area
#Sometimes reach-based one is not representative as it represents pourpoint of reach even if dam may be located further upstream
#Sometimes pixel-based one is not representative as point was snapped to line which may go through the corner of a pixel which does not correspond to its drainage area
#So select whichever of these two measures is the closest to the reported drainage area
#An alternative would be to route through the network, etc. but would take a while
DamAttributes[, `:=`(pxaccuratiodiff = abs((AREA_DREN/up_area_skm_15s)-1), 
                     lineaccuratiodiff = abs((AREA_DREN/UPLAND_SKM)-1))] 
DamAttributes[is.na(lineaccuratiodiff) | pxaccuratiodiff<lineaccuratiodiff, 
              finalHydroSHEDS_DA := up_area_skm_15s]
DamAttributes[pxaccuratiodiff>=lineaccuratiodiff, 
              finalHydroSHEDS_DA := UPLAND_SKM]

DamAttributes[,finalHydroSHEDS_ratio := finalHydroSHEDS_DA/AREA_DREN]

#Scatterplot of ANEEL-reported drainage area for dam and HydroSHEDS network/DEM drainage area of 
#point after semi-automatic association to HydroSHEDS network.
areacomp_plot <- ggplot(DamAttributes, aes(finalHydroSHEDS_DA+1, AREA_DREN+1)) + 
  geom_abline(slope=1) +
  # geom_abline(intercept=0.25, slope=1, alpha=1/2) +
  # geom_abline(intercept=-0.25, slope=1, alpha=1/2) +
  geom_point(aes(color=manualsnap), size=3, alpha=0.5) +
  scale_x_log10(name='Drainage area according to HydroSHEDS river network',
                breaks=c(1 , 10, 100, 1000, 10000, 100000, 1000000),
                labels=c(0 , 10, 100, 1000, 10000, 100000, 1000000)) + 
  scale_y_log10(name='Drainage area from ANEEL Hydropower plant database',
                breaks=c(1, 10, 100, 1000, 10000, 100000, 1000000),
                labels=c(0 , 10, 100, 1000, 10000, 100000, 1000000)) +
  scale_color_manual(name='Locational accuracy control status',
                     labels=c('Checked - not adjusted',
                              'Checked - manually adjusted',
                              'Not checked'),
                     values=c('#4575b4', '#d73027', 'grey')) +
  scale_size_continuous(name = 'Potential Capacity (KW)') +
  theme_classic() +
  theme(text=element_text(size=14),
        legend.position=c(0.2, 0.9))
areacomp_plot

tiff(filename = file.path(figdir, paste0("Revision_DAcomparison.tiff")),
     height = 3000, width = 3000, res = 300, compression = c("lzw"))
areacomp_plot
dev.off()

######################## CHECK PRECISION OF COORDINATES PROVIDED IN ANEEL DATABASE AND SNAPPING DISTANCE #################
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

#Remove all kinds of special characters from dam lat and lon
DamAttributes[, `:=`(ANEEL_LATdam = gsub('[â€™Â]|(009d)', '', LAT_EIXO_G),
                     ANEEL_LONdam = gsub('[â€™Â]|(009d)', '', LONG_EIXO_))] %>%
  .[, `:=`(ANEEL_LATdam = gsub('^\\s*$', NA, ANEEL_LATdam),
          ANEEL_LONdam = gsub('^\\s*$', NA, ANEEL_LONdam))] %>%
  .[, `:=`(ANEEL_LATdam = gsub('[,]', '.', ANEEL_LATdam),
           ANEEL_LONdam = gsub('[,]', '.', ANEEL_LONdam))]

#Parse lat and lon into decimal degrees
DamAttributes[, `:=`(ANEEL_LATdamparse = as.numeric(parse_lat(ANEEL_LATdam)),
                     ANEEL_LONdamparse = as.numeric(parse_lon(ANEEL_LONdam)))]

#Parse and divide lat and lon into deg, min, and sec to check precision
newcols_lat <- c('ANEEL_LATdam_deg', 'ANEEL_LATdam_min', 'ANEEL_LATdam_sec')
newcols_lon <- c('ANEEL_LONdam_deg', 'ANEEL_LONdam_min', 'ANEEL_LONdam_sec')
DamAttributes[, (newcols_lat) := tstrsplit(gsub('(".*)|S', '', ANEEL_LATdam), 
                                           "[°\\s']", fixed = FALSE)[1:3]] 
DamAttributes[, (newcols_lon) := tstrsplit(gsub('(°".*)|S', '', ANEEL_LONdam), 
                                           "[°\\s']", fixed = FALSE)[1:3]] 

DamAttributes[ANEEL_LATdam_sec %in% c('0','',' ', 'NA'), ANEEL_LATdam_sec := NA]
DamAttributes[ANEEL_LONdam_sec %in% c('0','',' ', 'NA'), ANEEL_LONdam_sec := NA]
DamAttributes[ANEEL_LATdam_min %in% c('',' ','NA'), ANEEL_LATdam_min := NA]
DamAttributes[ANEEL_LONdam_min %in% c('',' ','NA'), ANEEL_LONdam_min := NA]


DamAttributes[!(is.na(as.numeric(ANEEL_LATdam_min))) & is.na(ANEEL_LATdam_sec),
              ANEEL_LATprecision := 
                (111000/60)*10^(-sapply(as.numeric(ANEEL_LATdam_min), decimalplaces))]
DamAttributes[!(is.na(ANEEL_LATdam_min) | is.na(ANEEL_LATdam_sec)),
              ANEEL_LATprecision := 
                (111000/3600)*10^(-sapply(as.numeric(ANEEL_LATdam_sec), decimalplaces))]
DamAttributes[!(is.na(as.numeric(ANEEL_LONdam_min))) & is.na(ANEEL_LONdam_sec),
              ANEEL_LONprecision := 
                (111000/60)*10^(-sapply(as.numeric(ANEEL_LONdam_min), decimalplaces))]
DamAttributes[!(is.na(ANEEL_LONdam_min) | is.na(ANEEL_LONdam_sec)),
              ANEEL_LONprecision := 
                (111000/3600)*10^(-sapply(as.numeric(ANEEL_LONdam_sec), decimalplaces))]

#Assess distance between dam coordinates and final coordinates
DamAttributes[!(is.na(ANEEL_LATdamparse) | is.na(ANEEL_LONdamparse)),
                finaldist := geodesic_inverse(c(ANEEL_LONdamparse, ANEEL_LATdamparse),
                                 c(POINT_X, POINT_Y))[1],
              by= seq_len(nrow(DamAttributes[!(is.na(ANEEL_LATdamparse) | is.na(ANEEL_LONdamparse)),]))]

damtopoint_distplot <- 
  ggplot(DamAttributes[finaldist<49000,], aes(x=finaldist/1000)) + 
  geom_histogram(aes(fill=manualsnap)) +
  scale_x_log10(name = 'Distance between reported dam location and network-adjusted location (km)',
                breaks=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
                labels=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
                expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(name='Locational accuracy control status',
                    labels=c('Checked - not adjusted',
                             'Checked - manually adjusted',
                             'Not checked'),
                    values=c('#4575b4', '#d73027', 'grey')) +
  theme_classic() + 
  theme(legend.position=c(0.2, 0.9),
        text=element_text(size=14))

damtopoint_distplot

tiff(filename = file.path(figdir, paste0("Revision_damtopointdist.tiff")),
     height = 3000, width = 3000, res = 300, compression = c("lzw"))
damtopoint_distplot
dev.off()

DamAttributes[finaldist<49000 & !is.na(finaldist),.N]
DamAttributes[is.na(finaldist), .N]
DamAttributes[finaldist<49000 & !is.na(finaldist), mean(finaldist)]
DamAttributes[finaldist<49000 & !is.na(finaldist), median(finaldist)]

distDAratio_plot <-
  ggplot(DamAttributes[AREA_DREN > 0,], aes(x=finaldist/1000, y=finalHydroSHEDS_ratio)) + 
  geom_hline(yintercept=1) +
  geom_point(aes(color=manualsnap), size=3, alpha=0.5) + 
  scale_x_log10(name = 'Distance between reported dam location and network-adjusted location (km)',
                breaks=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 1000),
                labels=c(0.0001, 0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 1000),
                expand=c(0,0)) +
  scale_y_log10(name='HydroSHEDS drainage area/ANEEL-reported drainage area',
                breaks=c(0.001, 0.01, 0.65, 1, 1.5, 10, 100, 10000),
                labels=c(0.001, 0.01, 0.65, 1, 1.5, 10, 100, 10000),
                expand=c(0,0))+
  scale_color_manual(name='Locational accuracy control status',
                     labels=c('Checked - not adjusted',
                              'Checked - manually adjusted',
                              'Not checked'),
                     values=c('#4575b4', '#d73027', 'grey')) +
  coord_cartesian(clip='off') +
  theme_classic() + 
  theme(legend.position=c(0.85, 0.9),
        text=element_text(size=14))
distDAratio_plot

tiff(filename = file.path(figdir, paste0("Revision_distDAratio.tiff")),
     height = 3000, width = 3000, res = 300, compression = c("lzw"))
distDAratio_plot
dev.off()
