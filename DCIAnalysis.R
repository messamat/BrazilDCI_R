#Developers: Thiago B. Couto; Mathis L. Messager
#Creation date: 
#Purpose: Compute riverine fragmentation by hydropower development in Brazil with the Dendritic Connectivity Index 
#Requirements: Output from GIS fragmentation analysis

library(rprojroot)
library(sf)
library(igraph) # library graph theory
library(data.table)
library(tictoc)
library(plyr)
library(magrittr)
library(microbenchmark)
library(devtools)
# devtools::install_github("hadley/lineprof")
library(lineprof)
library(Rcpp)
library(stringr)

#To get Rcpp to work - install RBuildTools
# Sys.setenv(Path = paste("C:/RBuildTools/3.4/bin", Sys.getenv("Path"), sep=";"))
# Sys.setenv(BINPREF = "C:/RBuildTools/3.4/mingw_64/bin/")

########################################### Directory structure ########################################
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
datadir <- file.path(rootdir, "data")
dcigdb <- file.path(resdir, 'dci.gdb')

########################################### Functions ########################################
DCIp <- function(d2, d3, print = NULL){
  "d2 = a data frame containing the links between patches and passability for each link
  d2$id1 = initial segment (FROM) connection; MUST BE FACTOR!
  d2$id2 = final segment (to) connection; MUST BE FACTOR!
  d2$pass = the passability from one segment to the next
  d3 = data frame containing all segments and its respective sizes
  d3$id = the ID of each segment; MUST BE FACTOR!
  d3$l = the SIZE of segments formed by barriers within a given river
  
  OBS: d2 and d3 must have the same segment names (ID)
  d = edges; vertices = nodes"
  
  require(igraph)
  
  graph <- graph.data.frame(d2, directed = F); #plot(graph)
  
  #calculating l_i/L ratio
  d3$l_L <- d3$l/sum(d3$l)
  
  patches <- as.vector(d3$id)
  resu <- outer(patches, patches, FUN = Vectorize(function(i,k){
    caminho_edge <- shortest_paths(graph, from = i, to = k, output = "both")$epath[[1]]
    pass <- E(graph)[caminho_edge]$pass
    c_ik <- prod(pass);
    
    d_subI <- d3[d3$id %in% i ,]
    d_subK <- d3[d3$id %in% k ,]
    return(c_ik * prod(d_subI$l_L, d_subK$l_L) * 100)
  } ))

  ## Plot Graphs while calculating DCI
  
  # Last correct 11/06
  # resu <- outer(patches,patches,FUN = Vectorize(function(i,k){
  #   caminho_edge <- shortest_paths(graph,from = i, to = k, output = "both")$epath[[1]];  caminho_edge
  #   pass <- E(graph)[caminho_edge]$pass
  #   c_ik <- prod(pass);
  #
  #   d_subI <-d3[d3$id %in% i ,]; d_subI
  #   d_subK <-d3[d3$id %in% k ,]; d_subK
  #   return(c_ik * prod(d_subI$l_L,d_subK$l_L) * 100)
  # } ))
  #resu;sum(resu)
  
  # #to print figure network
  # if(!is.null(print)){
  # 
  #   graph <- graph.data.frame(d2, directed = T)
  #   graph$layout <- layout_with_kk;#plot(graph)
  #   # graph$layout <- layout.reingold.tilford
  #   # graph$layout <- layout_as_tree
  #   # graph$layout <- layout_with_fr
  #   #sort data.frame according to graph
  #   d3 = d3[match(V(graph)$name, d3$id), ]
  #   #include node size on graph object
  #   V(graph)$l = (d3$l/sum(d3$l)) * 30
  #   V(graph)$frame.color = "white"
  #   # print(plot(graph,edge.width=E(graph)$pass*5,vertex.size=(V(graph)$l),
  #   #            main=paste("DCIp = ",round(sum(resu),2),sep=""),
  #   #            vertex.label.dist=1.5,vertex.color="darkgreen",
  #   #            edge.arrow.size=0.01
  #   #            ))
  
  ## Plot networks while running
  #   print(plot(graph,
  #              vertex.label = V(graph)$name, #vertex.shape="none",
  #              vertex.label.cex = 1.5,
  #              # vertex.label.cex = (V(graph)$l)/20,
  #              vertex.label.dist = 1,
  #              vertex.label.color = "black",
  #              vertex.size = (V(graph)$l),
  #              vertex.color = "yellow",
  #              vertex.label.font = 3,
  #              edge.width = E(graph)$pass * 10,
  #              edge.arrow.size = 0.05,
  #              edge.color = "gray45",
  #              main = paste("DCIp = ",round(sum(resu), 2), sep = "")
  #   ))
  # }
  
  #returning DCIp
  return(sum(resu))
}

#Function to get all combinations of 2 vectors (faster version of combn)
#From https://stackoverflow.com/questions/26828301/faster-version-of-combn/26828486
cppFunction('
            Rcpp::DataFrame combi2inds(const Rcpp::CharacterVector inputVector){
            const int len = inputVector.size();
            const int retLen = len * (len-1) / 2;
            Rcpp::IntegerVector outputVector1(retLen);
            Rcpp::IntegerVector outputVector2(retLen);
            int indexSkip;
            for (int i = 0; i < len; ++i){
            indexSkip = len * i - ((i+1) * i)/2;
            for (int j = 0; j < len-1-i; ++j){
            outputVector1(indexSkip+j) = i+1;
            outputVector2(indexSkip+j) = i+j+1+1;
            }
            }
            return(Rcpp::DataFrame::create(Rcpp::Named("xid") = outputVector1,
            Rcpp::Named("yid") = outputVector2));
            };
            ')

#DCI for potadromous species
DCIp_opti <- function(d2, d3, print = NULL){
  "d2 = a data frame containing the links between patches and passability for each link
  d2$id1 = initial segment (FROM) connection; MUST BE FACTOR!
  d2$id2 = final segment (to) connection; MUST BE FACTOR!
  d2$pass = the passability from one segment to the next
  d3 = data frame containing all segments and its respective sizes
  d3$id = the ID of each segment; MUST BE FACTOR!
  d3$l = the SIZE of segments formed by barriers within a given river
  
  OBS: d2 and d3 must have the same segment names (ID)
  d = edges; vertices = nodes"
  
  require(igraph)
  
  graph <- graph.data.frame(d2, directed = F); #plot(graph)
  
  #calculating l_i/L ratio - for each segment, the sum of that segment to that of all segments in basin
  d3$l_L <- d3$l/sum(d3$l)
  
  #Get lower triangle of matrix of all pairwise combination of segments
  indices <- combi2inds(d3$id)
  lowertricomb <- setDT(list(i=d3$id[indices$xid], k=d3$id[indices$yid])) #Does not include matrix diagonal
  
  DCImatp <- function(x) {
    #Compute connectivity for all unique pairs of segments in basin (for potadromous species)
    return(prod(shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]$pass) *
             prod(d3[d3$id==x[1],'l_L'], d3[d3$id==x[2],'l_L']) *
             100) 
  }
  
  DCImatpself <- function(x) {
    #Compute connectivity within segments
    return(1 * (d3[d3$id == x , 'l_L']^2) * 100)
  }
  
  resutri <- apply(lowertricomb, 1, DCImatp)
  resudia <- sapply(d3$id, DCImatpself)
  #Compute DCI by sum connectivity for full matrix of pairwise combination 
  return(sum(resutri)*2+sum(resudia))
}

#DCI for diadromous species
#Calculate connectivity in terms of the probability that a fish can move in both directions between 
#the mouth of the river and another section of the river network
DCId_opti <- function(d2, d3, print = NULL){
  #Main difference with DCIp is that not passability is not assessed for all combinations of segments but between
  #mouth segment and all segments (included itself) [and only one l_i/L ratio is used in the calculation]
  
  "d2 = a data frame containing the links between patches and passability for each link
  d2$id1 = initial segment (FROM) connection; MUST BE FACTOR!
  d2$id2 = final segment (to) connection; MUST BE FACTOR!
  d2$pass = the passability from one segment to the next
  d3 = data frame containing all segments and its respective sizes
  d3$id = the ID of each segment; MUST BE FACTOR!
  d3$l = the SIZE of segments formed by barriers within a given river
  
  OBS: d2 and d3 must have the same segment names (ID)
  d = edges; vertices = nodes"
  
  require(igraph)
  graph <- graph.data.frame(d2, directed = F); #plot(graph)
  
  #calculating l_i/L ratio - for each segment, the sum of that segment to that of all segments in basin
  d3$l_L <- d3$l/sum(d3$l)
  
  #get segment at mouth of network
  mouthsegID <- unique(d2$id1)[!(unique(d2$id1) %in% unique(d2$id2))]
  
  #Get lower triangle of matrix of all pairwise combination of segments
  mouthcomb <- data.table(i= rep(mouthsegID, length(d3$id)), j=d3$id)
  
  DCImatd <- function(x) {
    #Compute connectivity for all unique pairs of segments in basin (for diadromous species)
    #Note: when from= and to= are the same in shortest_paths, returns an empty numeric that, when fed to prod, returns 1
    return(prod(shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]$pass) *
             d3[d3$id==x[2],'l_L'][[1]] *
             100) 
  }
  
  return(sum(apply(mouthcomb, 1, DCImatd)))
}

########### RUN DCI - All levels dataset #################################################################

##### LEVEL 8 ######
## Loop over scenarios
level <- 8
passprob=0.1

# Identify possible scenatios
# scenarios <- c("All_current.txt", "All_future.txt", "SHP_current.txt", "SHP_future.txt", "LHP_current.txt", "LHP_future.txt")

tic()
#Read in all scenarios
netcrude <- as.data.table(sf::st_read(dsn = dcigdb, layer='networkattributes'))

# NetworkBRAZILcrude <- lapply(seq_along(scenarios), function(i) {
#   scenariotab <- fread(file.path(resdir, paste0("Brazil_L", level, "_", scenarios[i])),
#                        stringsAsFactors=T, data.table=T, integer64="numeric")
#   scenariotab[, scenario := scenarios[i]]
# }) %>%
#   rbindlist

DamAttributes <- as.data.table(sf::st_read(dsn = dcigdb, layer='damattributes'))
# 
# DamAttributes <- lapply(seq_along(scenarios), function(i) {
#   scenariotab <- fread(file.path(resdir, paste0("DamAttributes_L", level, "_", scenarios[i])),
#                        stringsAsFactors=T, data.table=T, integer64="numeric")
#   scenariotab[, scenario := scenarios[i]]
# }) %>%
#   rbindlist

# Organize the matrix based on type and stage of each dam
DamAttributes[, ESTAGIO_1 := factor(ifelse(ESTAGIO_1 == "OperaÃ§Ã£o", 'Operation', 'Planned'), 
                                    levels=c('Operation', 'Planned'))]
DamAttributes[, Tipo_1 := factor(ifelse(Tipo_1 == "UHE", 'LHP', 'SHP'), 
                                 levels=c('SHP', 'LHP'))]

DamAttributes[, `:=`(All_current = ifelse(DamAttributes$ESTAGIO_1 == 'Operation', 0.1, 1),
                     All_future = 1,
                     SHP_current = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                            DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHP_current = ifelse(DamAttributes$ESTAGIO_1 == 'Operation' & 
                                            DamAttributes$Tipo_1 == 'LHP', 0.1, 1),
                     SHP_future = ifelse(DamAttributes$Tipo_1 == 'SHP', 0.1, 1),
                     LHP_future = ifelse(DamAttributes$Tipo_1 == 'LHP', 0.1, 1)
                     )]

#Slightly rename basin id for dams to avoid conflicts in data.table processing
setnames(DamAttributes, c('HYBAS_ID08'), c('DAMBASID08'))

#Correct few glitches from exceptions that can be automatically corrected later
#(here, the dam sits at the downstream confluence of a mainstem and a first order stream)
bug1 <- DamAttributes[DAMBASID08 == 6080440420 & DownSeg==1676,]
bug1[, UpSeg := 1679]
DamAttributes <- rbind(DamAttributes, bug1)

#Keep only network within basins that contain dams, excluding coastal segments without dams
#(only keep segments that are either in UpSeg or DownSeg of a dam and within a basin with a dam)
dambas08 <- DamAttributes[!is.na(DAMBASID08), unique(DAMBASID08)] #vector of all basins that contain a dam
netdams <- netcrude[((paste(SEGID, HYBAS_ID08) %in% paste(DamAttributes$UpSeg, DamAttributes$DAMBASID08)) |
                       (paste(SEGID, HYBAS_ID08) %in% paste(DamAttributes$DownSeg, DamAttributes$DAMBASID08))) &
                      HYBAS_ID08 %in% dambas08,][
                        , SEGIDBAS := paste(SEGID, HYBAS_ID08, sep='_')] 

#Identify basins with multiple networks that have dams
#Assign an extra digit to basins to account for those that have multiple networks with dams
multibas <- netdams[(SEGID %in% DamAttributes$DownSeg) & 
                      !(SEGID %in% DamAttributes$UpSeg),][
                        , .N, by=HYBAS_ID08][
                          N>1, HYBAS_ID08]
DamAttributes[, `:=`(DownSegBAS = paste(DownSeg, DAMBASID08, sep='_'),
                     UpSegBAS = paste(UpSeg, DAMBASID08, sep='_'),
                     SEGIDDAMBAS = paste(SEGID, DAMBASID08, sep='_'))]

coastdamSEG <- netdams[(SEGID %in% DamAttributes$DownSeg) & 
                         !(SEGID %in% DamAttributes$UpSeg) &
                         HYBAS_ID08 %in% multibas,][
                             , HYBAS_ID08ext := paste0(HYBAS_ID08, seq_len(.N)), by = HYBAS_ID08]

seg_list <- c()
for (segnum in seq_along(coastdamSEG$HYBAS_ID08)) {
  seg_mouth <- coastdamSEG[segnum,]
  seg_iter <- seg_mouth
  while (nrow(seg_iter)>0) {
    seg_list <- c(seg_list, paste(seg_iter$SEGID, seg_mouth$HYBAS_ID08ext, sep='_'))
    seg_iter <- DamAttributes[DownSegBAS == seg_iter$SEGIDBAS,] 
  }
}

netdams <- merge(netdams, 
                 data.frame(HYBAS_ID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                            SEGIDBAS = substr(seg_list, 1, nchar(seg_list)-1)), 
                 by='SEGIDBAS', all.x=T)
netdams[is.na(HYBAS_ID08ext), HYBAS_ID08ext := paste0(HYBAS_ID08, 1)]


DamAttributes <- merge(DamAttributes, 
                 data.frame(DAMBASID08ext = str_split(seg_list, '_', simplify=T)[,2], 
                            SEGIDDAMBAS = substr(seg_list, 1, nchar(seg_list)-1)), 
                 by='SEGIDDAMBAS', all.x=T)
DamAttributes[is.na(DAMBASID08ext), DAMBASID08ext := paste0(DAMBASID08, 1)]


#Compute DCI
DCI_L8_current <- netdams[,
                          list(DCI = DCIp_opti(
                            DamAttributes[DAMBASID08ext == HYBAS_ID08ext,
                                          list(
                                            id1 = DownSeg,
                                            id2 = UpSeg,
                                            pass = All_current
                                          )],
                            .SD[, list(id=as.character(SEGID),
                                       l=Shape_Length)],
                            print = F)),
                          by=.(HYBAS_ID08ext)]
toc()

#To do going forward:
#Join back to basins, all others are 100 (no dams) or Null (if not in netcrude, because coastal)
#To add multiple scenarios, can melt DamAttributes by scenario
#To run in parallel, could divide scenarios into 4-8 chunks and run as a foreach loop or parallel apply

#For each region, create a list scenario with all possible permutations
Ndamsbas <- DamAttributes[ESTAGIO_1 == "Planned",.N, by=DAMBASID08]


allscenarios <- ldply(unique(DamAttributes[ESTAGIO_1 == "Planned" & !is.na(DAMBASID08), DAMBASID08]), function(reg) {
  #Create all combinations of future dams
  permut <- as.data.table(DamAttributes[ESTAGIO_1 == "Planned" & DAMBASID08 == reg, 
                                  expand.grid(rep(list(c(0.1, 1)), .N))]) 
  colnames(permut) <- as.character(DamAttributes[ESTAGIO_1 == "Planned" & DAMBASID08 == reg, DAMID])
  
  #Add all existing dams
  permut[, as.character(DamAttributes[ESTAGIO_1 == "Operation" & DAMBASID08 == reg, DAMID]) := 0.1]
  
  #Melt and assign to basin-scenario ID
  scenarios_reg <- melt(cbind(permut, 
                              DAMBASID08_scenario = paste(reg, seq_len(nrow(permut)), sep='_')), 
                        id.var="DAMBASID08_scenario")
  return(scenarios_reg)
})
setnames(allscenarios, 'variable', 'DAMID')

allscenarios_seg <- merge(allscenarios, DamAttributes[, .(DAMID, DownSeg, UpSeg)], by='DAMID')

################### IN PROGRESS #############
#Number of unique scenarios
length(unique(allscenarios$DAMBASID08_scenario))


d2 <- netdams[HYBAS_ID08 == 6080704460, DamAttributes[DAMBASID08 == 6080704460,
                                                                            list(
                                                                                id1 = DownSeg,
                                                                                id2 = UpSeg,
                                                                                pass = All_current
                                                                              )]]
#To run on subset
DCI_L8_all<- netdams[HYBAS_ID08 %in% dambas08,
                           list(DCI = DCIp_opti(
                             DamAttributes[DAMBASID08 == HYBAS_ID08,
                                           list(
                                             id1 = DownSeg,
                                             id2 = UpSeg,
                                             pass = All_current
                                           )],
                             .SD[, list(id=as.character(SEGID),
                                        l=Shape_Length)],
                             print = F)),
                           by=.(HYBAS_ID08)]


############# EXTRA STUFF
# #Vectorize join that
# ### Attributes of the basin (N of dams, cumulative MW)
# # basinNDams_L8[, x + 1] <- dim(DamX)[1]
# # basinMWDams_L8[, x + 1] <- sum(DamX$POT_KW, na.rm = T)/1000 #kw to MW
# 
# ## Create a csv file with the output of DCI analysis for all the scenarios
# write.csv(DCI_BRAZIL_L8, file = "DCI_Brazil_L8.csv")
# write.csv(basinNDams_L8, file = "basinNDams_L8.csv")
# write.csv(basinMWDams_L8, file = "basinMWDams_L8.csv")


#Get every possible future scenario
damsall<- fread(file.path(resdir,"DamAttributes_L8_All_future.txt"),
                stringsAsFactors=T, data.table=T, integer64="numeric")
length(damsall[, unique(batRegion)])

# DCIall <- NetworkBRAZIL[Numbsegments > 1,
#                         list(DCI = DCIp_opti(
#                           DamAttributes[batRegion == Region8 & scenario==scenario,
#                                         list(
#                                           id1 = paste("seg", Min_batNet, sep =""), #DownSeg
#                                           id2 = paste("seg", Max_batNet, sep =""), #UpSeg
#                                           pass = passprob
#                                         )],
#                           .SD[, list(id=paste0("seg", batNetID),
#                                      l=sum(Shape_Leng)), by=batNetID],
#                           print = F)),
#                         by=.(Region8, scenario)]


#############################################################################################################
#############################################################################################################
#############################################################################################################
#####################################     PLOTS FIGURE 1     ################################################

## Import data
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8.csv", header = T))

## Filter results of basins with both SHPs and LHPs in the future
resuDCI_Both <- resuDCI[resuDCI$SHP_fut != 100 & resuDCI$LHP_fut != 100, ]

## Calculate differences
DCILoss_All <- resuDCI$All_fut - resuDCI$All_curr
DCILoss_SHP <- resuDCI$SHP_fut - resuDCI$SHP_curr
DCILoss_LHP <- resuDCI$LHP_fut - resuDCI$LHP_curr

DCILoss_Both_All <- resuDCI_Both$All_fut - resuDCI_Both$All_curr
DCILoss_Both_SHP <- resuDCI_Both$SHP_fut - resuDCI_Both$SHP_curr
DCILoss_Both_LHP <- resuDCI_Both$LHP_fut - resuDCI_Both$LHP_curr


##### Plot

tiff(filename = "Figure1.tiff", height = 2196, width = 2900, res = 300, compression = c("lzw"))

mat <- matrix(c(rep(1, times = 2), rep(2, times = 2)), ncol = 2, nrow = 2, byrow = F)
layout(mat = mat, widths = c(4, 4), heights = c(rep(0.5, times = 2)))

# Par
par (oma = c(5.5, 6, 4, 4.0), mar = c(0.5, 0, 0, 0), bty = "n")

par(bty = "n")

### Plot data for all basins

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI)[1])

plot(x = treatment, y = DCILoss_All, ylim = c(-90, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment, y = DCILoss_All, bg = "#99999908", col = c("#42424220"), cex = 3, pch = 21)
points(x = treatment * 2, y = DCILoss_SHP, bg = "#4876FF08", col = c("#42424220"), cex = 3, pch = 21)
points(x = treatment * 3, y = DCILoss_LHP, bg = "#8B000008", col = c("#42424220"), cex = 3, pch = 21)

# Plot 95% CI
UppCI_All <- quantile(DCILoss_All, 0.975, na.rm = T)
LowCI_All <- quantile(DCILoss_All, 0.025, na.rm = T)
UppCI_SHP <- quantile(DCILoss_SHP, 0.975, na.rm = T)
LowCI_SHP <- quantile(DCILoss_SHP, 0.025, na.rm = T)
UppCI_LHP <- quantile(DCILoss_LHP, 0.975, na.rm = T)
LowCI_LHP <- quantile(DCILoss_LHP, 0.025, na.rm = T)

segments(y0 = UppCI_All, x0 = treatment, 
         y1 = LowCI_All, x1 = treatment, lwd = 4.3, col = "black")
segments(y0 = UppCI_SHP, x0 = treatment * 2, 
         y1 = LowCI_SHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_LHP, x0 = treatment * 3, 
         y1 = LowCI_LHP, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot medians
points(x = treatment[1], y = mean(DCILoss_All), pch = 18, col = "black", cex = 3.0)
points(x = treatment[1] * 2, y = mean(DCILoss_SHP), pch = 18, col = "black", cex = 3.0)
points(x = treatment[1] * 3, y = mean(DCILoss_LHP), pch = 18, col = "black", cex = 3.0)

# Plot axes and labels
axis(side = 2, at = c(0, -20, -40, -60, -80), cex.axis = 2.3, line = - 3)
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("All", "SHP", "LHP"),
     cex.axis = 2.3, line = 1)
mtext("DCI loss", side = 2, cex = 2.5, line = 1.5)



### Plot data just for basins with both SHPs and LHPs

# Position dam treatments in the x axis
treatment <- rep(0.5, times = dim(resuDCI_Both)[1])

plot(x = treatment, y = DCILoss_Both_All, ylim = c(-90, 0), xlim = c(0, 1.8), type = "n", ylab = "",
     xlab = "", xaxt = "n", yaxt = "n")

# Plot poins (each one is a basin)
points(x = treatment, y = DCILoss_Both_All, bg = "#99999908", col = c("#42424220"), cex = 3, pch = 21)
points(x = treatment * 2, y = DCILoss_Both_SHP, bg = "#4876FF08", col = c("#42424220"), cex = 3, pch = 21)
points(x = treatment * 3, y = DCILoss_Both_LHP, bg = "#8B000008", col = c("#42424220"), cex = 3, pch = 21)

# Plot 95% CI
UppCI_Both_All <- quantile(DCILoss_Both_All, 0.975, na.rm = T)
LowCI_Both_All <- quantile(DCILoss_Both_All, 0.025, na.rm = T)
UppCI_Both_SHP <- quantile(DCILoss_Both_SHP, 0.975, na.rm = T)
LowCI_Both_SHP <- quantile(DCILoss_Both_SHP, 0.025, na.rm = T)
UppCI_Both_LHP <- quantile(DCILoss_Both_LHP, 0.975, na.rm = T)
LowCI_Both_LHP <- quantile(DCILoss_Both_LHP, 0.025, na.rm = T)

segments(y0 = UppCI_Both_All, x0 = treatment, 
         y1 = LowCI_Both_All, x1 = treatment, lwd = 4.3, col = "black")
segments(y0 = UppCI_Both_SHP, x0 = treatment * 2, 
         y1 = LowCI_Both_SHP, x1 = treatment * 2, lwd = 4.3, col = "black")
segments(y0 = UppCI_Both_LHP, x0 = treatment * 3, 
         y1 = LowCI_Both_LHP, x1 = treatment * 3, lwd = 4.3, col = "black")

# Plot medians
points(x = treatment[1], y = mean(DCILoss_Both_All), pch = 18, col = "black", cex = 3.0)
points(x = treatment[1] * 2, y = mean(DCILoss_Both_SHP), pch = 18, col = "black", cex = 3.0)
points(x = treatment[1] * 3, y = mean(DCILoss_Both_LHP), pch = 18, col = "black", cex = 3.0)

# Plot axes and labels
axis(side = 1, at = c(treatment[1], treatment[1] * 2, treatment[1] * 3), labels = c("All", "SHP", "LHP"),
     cex.axis = 2.3, line = 1)

mtext("Level 8", side = 3, cex = 3.0, line = 0)


dev.off()




#############################################################################################################
#############################################################################################################
#############################################################################################################
#####################################     PLOTS FIGURE 3     ################################################

## Import data on sub-basin areas and other attributes for each level
BasinData_L10 <- read.csv("Brazil_L10_All_future.txt", header = T)
BasinData_L8 <- read.csv("Brazil_L8_All_future.txt", header = T)
BasinData_L6 <- read.csv("Brazil_L6_All_future.txt", header = T)
BasinData_L4 <- read.csv("Brazil_L4_All_future.txt", header = T)

## Get the attributes of unique sub-basins - checked if it matches with basinList order (sum(basinList != SUBAreas$Region8))
SUBAttributes_L10 <- BasinData_L10[duplicated(BasinData_L10$Region10) == F, ]
SUBAreas_L10 <- SUBAttributes_L10$SUB_AREA

SUBAttributes_L8 <- BasinData_L8[duplicated(BasinData_L8$Region8) == F, ]
SUBAreas_L8 <- SUBAttributes_L8$SUB_AREA

SUBAttributes_L6 <- BasinData_L6[duplicated(BasinData_L6$Region6) == F, ]
SUBAreas_L6 <- SUBAttributes_L6$SUB_AREA

SUBAttributes_L4 <- BasinData_L4[duplicated(BasinData_L4$Region4) == F, ]
SUBAreas_L4 <- SUBAttributes_L4$SUB_AREA


############## Number of dams in the sub-basin

## Level 10
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L10.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L10.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L10.csv", header = T))

# plot(resuNDams$N_All_curr, resuDCI$All_curr, ylim = c(0, 100), xlim = c(0, 9), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_All_fut, resuDCI$All_fut, ylim = c(0, 100), xlim = c(0, 9), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#99999930")

#plot(resuNDams$N_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), xlim = c(0, 9), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_SHP_fut, resuDCI$SHP_fut, ylim = c(0, 100), xlim = c(0, 9), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 10")
#plot(resuNDams$N_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), xlim = c(0, 9), ylab = "DCI", xlab = "Number of dams")
points(resuNDams$N_LHP_fut, resuDCI$LHP_fut, ylim = c(0, 100), xlim = c(0, 9), ylab = "DCI", xlab = "Number of dams",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 8
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8.csv", header = T))

#plot(resuNDams$N_All_curr, resuDCI$All_curr, ylim = c(0, 100), xlim = c(0, 40), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_All_fut, resuDCI$All_fut, ylim = c(0, 100), xlim = c(0, 40), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), xlim = c(0, 40), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_SHP_fut, resuDCI$SHP_fut, ylim = c(0, 100), xlim = c(0, 40), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 8")
#plot(resuNDams$N_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), xlim = c(0, 40), ylab = "DCI", xlab = "Number of dams")
points(resuNDams$N_LHP_fut, resuDCI$LHP_fut, ylim = c(0, 100), xlim = c(0, 40), ylab = "DCI", xlab = "Number of dams",
       pch = 21, col = "#42424220", bg = "#8B000030")

## Level 6
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L6.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L6.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L6.csv", header = T))

#plot(resuNDams$N_All_curr, resuDCI$All_curr, ylim = c(0, 100), xlim = c(0, 150), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_All_fut, resuDCI$All_fut, ylim = c(0, 100), xlim = c(0, 150), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), xlim = c(0, 150), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_SHP_fut, resuDCI$SHP_fut, ylim = c(0, 100), xlim = c(0, 150), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 6")
# plot(resuNDams$N_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), xlim = c(0, 150), ylab = "DCI", xlab = "Number of dams")
points(resuNDams$N_LHP_fut, resuDCI$LHP_fut, ylim = c(0, 100), xlim = c(0, 150), ylab = "DCI", xlab = "Number of dams",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 4
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L4.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L4.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L4.csv", header = T))

#plot(resuNDams$N_All_curr, resuDCI$All_curr, ylim = c(0, 100), xlim = c(0, 500), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_All_fut, resuDCI$All_fut, ylim = c(0, 100), xlim = c(0, 500), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), xlim = c(0, 500), ylab = "DCI", xlab = "Number of dams")
plot(resuNDams$N_SHP_fut, resuDCI$SHP_fut, ylim = c(0, 100), xlim = c(0, 500), ylab = "DCI", xlab = "Number of dams",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 4")
#plot(resuNDams$N_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), xlim = c(0, 500), ylab = "DCI", xlab = "Number of dams")
points(resuNDams$N_LHP_fut, resuDCI$LHP_fut, ylim = c(0, 100), xlim = c(0, 500), ylab = "DCI", xlab = "Number of dams",
       pch = 21, col = "#42424220", bg = "#8B000020")



############ Density of dams in the sub-basin

## Level 10
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L10.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L10.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L10.csv", header = T))

#plot(resuNDams$N_All_curr/SUBAreas_L10, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_All_fut/SUBAreas_L10, resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr/SUBAreas_L10, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_SHP_fut/SUBAreas_L10, resuDCI$SHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 10")
#plot(resuNDams$N_LHP_curr/SUBAreas_L10, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
points(resuNDams$N_LHP_fut/SUBAreas_L10, resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 8
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8.csv", header = T))

#plot(resuNDams$N_All_curr/SUBAreas_L8, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_All_fut/SUBAreas_L8, resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr/SUBAreas_L8, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_SHP_fut/SUBAreas_L8, resuDCI$SHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 8")
#plot(resuNDams$N_LHP_curr/SUBAreas_L8, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
points(resuNDams$N_LHP_fut/SUBAreas_L8, resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 6
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L6.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L6.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L6.csv", header = T))

#plot(resuNDams$N_All_curr/SUBAreas_L6, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_All_fut/SUBAreas_L6, resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr/SUBAreas_L6, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_SHP_fut/SUBAreas_L6, resuDCI$SHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 6")
#plot(resuNDams$N_LHP_curr/SUBAreas_L6, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
points(resuNDams$N_LHP_fut/SUBAreas_L6, resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 4
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L4.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L4.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L4.csv", header = T))

#plot(resuNDams$N_All_curr/SUBAreas_L4, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_All_fut/SUBAreas_L4, resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuNDams$N_SHP_curr/SUBAreas_L4, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
plot(resuNDams$N_SHP_fut/SUBAreas_L4, resuDCI$SHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 4")
#plot(resuNDams$N_LHP_curr/SUBAreas_L4, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)")
points(resuNDams$N_LHP_fut/SUBAreas_L4, resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Dernsity of dams (km2)",
       pch = 21, col = "#42424220", bg = "#8B000020")





############# Summed capacity (MW) in the sub-basin

## Level 10
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L10.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L10.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L10.csv", header = T))

#plot(resuMWDams$MW_All_curr, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_All_fut), resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuMWDams$MW_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_SHP_fut), resuDCI$SHP_fut, xlim = c(-7, 10), ylim = c(0, 100), ylab = "DCI", xlab = "log Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 10")
#plot(resuMWDams$MW_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
points(log(resuMWDams$MW_LHP_fut), resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 8
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L8.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L8.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L8.csv", header = T))

#plot(resuMWDams$MW_All_curr, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_All_fut), resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuMWDams$MW_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_SHP_fut), resuDCI$SHP_fut, xlim = c(-7, 10), ylim = c(0, 100), ylab = "DCI", xlab = "log Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 8")
#plot(resuMWDams$MW_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
points(log(resuMWDams$MW_LHP_fut), resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 6
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L6.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L6.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L6.csv", header = T))

#plot(resuMWDams$MW_All_curr, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_All_fut), resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuMWDams$MW_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_SHP_fut), resuDCI$SHP_fut, xlim = c(-5, 10), ylim = c(0, 100), ylab = "DCI", xlab = "log Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 6")
#plot(resuMWDams$MW_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
points(log(resuMWDams$MW_LHP_fut), resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
       pch = 21, col = "#42424220", bg = "#8B000020")

## Level 4
resuDCI <- as.data.frame(read.csv("DCI_Brazil_L4.csv", header = T))
resuNDams <- as.data.frame(read.csv("basinNDams_L4.csv", header = T))
resuMWDams <- as.data.frame(read.csv("basinMWDams_L4.csv", header = T))

#plot(resuMWDams$MW_All_curr, resuDCI$All_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_All_fut), resuDCI$All_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#99999930")
#plot(resuMWDams$MW_SHP_curr, resuDCI$SHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
plot(log(resuMWDams$MW_SHP_fut), resuDCI$SHP_fut, xlim = c(0, 10), ylim = c(0, 100), ylab = "DCI", xlab = "log Capacity (MW)",
     pch = 21, col = "#42424220", bg = "#4876FF40", main = "Level 4")
#plot(resuMWDams$MW_LHP_curr, resuDCI$LHP_curr, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)")
points(log(resuMWDams$MW_LHP_fut), resuDCI$LHP_fut, ylim = c(0, 100), ylab = "DCI", xlab = "Capacity (MW)",
       pch = 21, col = "#42424220", bg = "#8B000020")



