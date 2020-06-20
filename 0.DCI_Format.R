#Import packages
source('00.DCI_packages.R')
#Import directory structure and functions
source('00.DCI_functions.R')

## Import network and dams dataset
NetworkBRAZILCrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='networkattributes'))
DamAttributesCrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='damattributes'))

#------ DEAL WITH COASTAL DRAINAGES ------
#Keep only network within basins that contain dams, excluding coastal segments without dams
#(only keep segments that are either in UpSeg or DownSeg of a dam and within a basin with a dam)
Basin08WithDams <- DamAttributesCrude[!is.na(HYBAS_ID08), unique(HYBAS_ID08)] #vector of all basins that contain a dam
NetworkBRAZIL <- NetworkBRAZILCrude[which((paste(NetworkBRAZILCrude$SEGID, NetworkBRAZILCrude$HYBAS_ID08) %in% 
                                             paste(DamAttributesCrude$UpSeg, DamAttributesCrude$HYBAS_ID08)) |
                                            (paste(NetworkBRAZILCrude$SEGID, NetworkBRAZILCrude$HYBAS_ID08) %in% 
                                               paste(DamAttributesCrude$DownSeg, DamAttributesCrude$HYBAS_ID08)) &
                                            NetworkBRAZILCrude$HYBAS_ID08 %in% Basin08WithDams),]
NetworkBRAZIL$SEGIDBAS <- with(NetworkBRAZIL, paste(SEGID, HYBAS_ID08, sep='_'))


#Create a new ID for dams that combines segment IDs and basin IDs
DamAttributesCrude$DownSegBAS <- with(DamAttributesCrude, paste(DownSeg, HYBAS_ID08, sep='_'))
DamAttributesCrude$UpSegBAS <- with(DamAttributesCrude, paste(UpSeg, HYBAS_ID08, sep='_'))
DamAttributesCrude$SEGIDBAS <- with(DamAttributesCrude, paste(SEGID, HYBAS_ID08, sep='_'))

#Identify basins with multiple networks that have dams 
#(Find those basins that have multiple segments with no downstream dam i.e. multiple outlets)
multibas <- as.data.frame(table(NetworkBRAZIL[(NetworkBRAZIL$SEGID %in% DamAttributesCrude$DownSeg) & 
                                                !(NetworkBRAZIL$SEGID %in% DamAttributesCrude$UpSeg), 'HYBAS_ID08']))

#For those, assign an extra digit to basins to account for those that have multiple networks with dams
coastdamSEG <- NetworkBRAZIL[(NetworkBRAZIL$SEGID %in% DamAttributesCrude$DownSeg) & 
                               !(NetworkBRAZIL$SEGID %in% DamAttributesCrude$UpSeg) &
                               NetworkBRAZIL$HYBAS_ID08 %in% multibas[multibas$Freq > 1, 'Var1'],]
coastdamSEG <- ddply(coastdamSEG, .(HYBAS_ID08), mutate, 
                     HYBAS_ID08ext = paste0(HYBAS_ID08, seq_along(SEGID)))

#Find connected segments within that basin to identify them as part of that subbasin
seg_list <- c()
for (segnum in seq_along(coastdamSEG$HYBAS_ID08)) {
  seg_mouth <- coastdamSEG[segnum,]
  seg_iter <- seg_mouth
  while (nrow(seg_iter)>0) {
    seg_list <- c(seg_list, paste(seg_iter$SEGID, seg_mouth$HYBAS_ID08ext, sep='_'))
    seg_iter <- DamAttributesCrude[DamAttributesCrude$DownSegBAS %in% seg_iter$SEGIDBAS,] 
  }
}

#Assign new subbasin ID too all river segments (simply HYBAS_ID08 + '1' if not multiple networks)
NetworkBRAZIL <- merge(NetworkBRAZIL, 
                       data.frame(HYBAS_ID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                                  SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1)), 
                       by='SEGIDBAS', all.x=T)
NetworkBRAZIL[is.na(NetworkBRAZIL$HYBAS_ID08ext), "HYBAS_ID08ext"] <-  NetworkBRAZIL[is.na(NetworkBRAZIL$HYBAS_ID08ext), 
                                                                                     paste0(HYBAS_ID08, 1)]

#Assign new subbasin ID too all dams (simply HYBAS_ID08 + '1' if not multiple networks)
DamAttributesCrude <- merge(DamAttributesCrude, 
                            data.frame(HYBAS_ID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                                       SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1)), 
                            by='SEGIDBAS', all.x=T)
DamAttributesCrude[is.na(HYBAS_ID08ext), "HYBAS_ID08ext"] <-  DamAttributesCrude[is.na(HYBAS_ID08ext),
                                                                                 paste0(HYBAS_ID08, 1)]

## Remove a basin with problems (Probably the dam is too close to the upstream edge)
# DamAttributesCrude <- DamAttributesCrude[-which(DamAttributesCrude$HYBAS_ID08 == 6080595090),]
# NetworkBRAZIL <- NetworkBRAZIL[-which(NetworkBRAZIL$HYBAS_ID08 == 6080595090),]

#------ FORMAT DATA ------
## Final Dataframes to run DCI analyses
DamAttributes <- DamAttributesCrude
NetworkBRAZIL <- NetworkBRAZIL

# Organize the matrix based on type and stage of each dam
DamAttributes[, ESTAGIO_1 := ifelse(ESTAGIO_1 != "OperaÃ§Ã£o", 'Planned', 'Operation')]
DamAttributes[, Tipo_1 := ifelse(Tipo_1 != "UHE", "SHP", "LHP")]

## Fix mistakes in the dataset
DamAttributes$Tipo_1[which(DamAttributes$Tipo_1 == "SHP" & DamAttributes$POT_KW > 30000)] <- "LHP"    #UHE Buritizal(Das Mortes), UHE Resplendor(Doce), UHE Pouso Alto (Sucuriú)
DamAttributes$Tipo_1[which(DamAttributes$Tipo_1 == "LHP" & DamAttributes$POT_KW < 30000 &
                             DamAttributes$AREA_NA_MA < 13.0 & DamAttributes$ESTAGIO_1 == "Planned")] <- "SHP"   #Keep old dams as UHEs and new ones as SHPs

#Convert dataframes to data tables to speed up analysis
DamAttributes <- as.data.table(DamAttributes)
NetworkBRAZIL <- as.data.table(NetworkBRAZIL)

DamAttributes[, `:=`(Allcurrent = ifelse(DamAttributes$ESTAGIO_1 == 'Operation', 0.1, 1),
                     Allfuture = 1,
                     SHPcurrent = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                           DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHPcurrent = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                           DamAttributes$Tipo_1 == 'LHP', 0.1, 1),
                     SHPfuture = ifelse(DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHPfuture = ifelse(DamAttributes$Tipo_1 == 'LHP', 0.1, 1)
)]

setnames(DamAttributes, 'HYBAS_ID08ext', 'DAMBAS_ID08ext')

write.fst(DamAttributes, file.path(resdir, 'DamAttributes.fst'))
write.fst(NetworkBRAZIL, file.path(resdir, 'NetworkBRAZIL.fst'))
