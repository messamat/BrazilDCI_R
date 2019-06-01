#Developers: Thiago B. Couto; Mathis L. Messager
#Creation date: 
#Purpose: Compute riverine fragmentation by hydropower development in Brazil with the Dendritic Connectivity Index 
#Requirements: Output from BAT tool

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
  
  DCImat <- function(x) {
    #Compute connectivity for all unique pairs of segments in basin
    return(prod(shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]$pass) *
             prod(d3[d3$id==x[1],'l_L'], d3[d3$id==x[2],'l_L']) *
             100) 
    }
  
  DCImatself <- function(x) {
    #Compute connectivity within segments
    return(1 * (d3[d3$id == x , 'l_L']^2) * 100)
  }
  resutri <- apply(lowertricomb, 1, DCImat)
  resudia <- sapply(d3$id, DCImatself)
  #Compute DCI by sum connectivity for full matrix of pairwise combination 
  return(sum(resutri)*2+sum(resudia))
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
setnames(DamAttributes, c('HYBAS_ID04', 'HYBAS_ID06', 'HYBAS_ID08'), c('DAMBASID04','DAMBASID06','DAMBASID08'))

#Correct few glitches from exceptions that can be automatically corrected later
#(here, the dam sits at the downstream confluence of a mainstem and a first order stream)
bug1 <- DamAttributes[DAMBASID08 == 6080440420 & DownSeg==1676,]
bug1[, UpSeg := 1679]
DamAttributes <- rbind(DamAttributes, bug1)

#Compute DCI
dambas08 = DamAttributes[!is.na(DAMBASID08), unique(DAMBASID08)]
DCI_L8 <- netcrude[HYBAS_ID08 %in% dambas08,
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
toc()


#To do going forward:
#Join back to basins, all others are 100 (no dams) or Null (if not in netcrude, because coastal)
#To add multiple scenarios, can melt DamAttributes by scenario
#To run in parallel, could divide scenarios into 4-8 chunks and run as a foreach loop or parallel apply

#Vectorize join that
### Attributes of the basin (N of dams, cumulative MW)
# basinNDams_L8[, x + 1] <- dim(DamX)[1]
# basinMWDams_L8[, x + 1] <- sum(DamX$POT_KW, na.rm = T)/1000 #kw to MW

## Create a csv file with the output of DCI analysis for all the scenarios
write.csv(DCI_BRAZIL_L8, file = "DCI_Brazil_L8.csv")
write.csv(basinNDams_L8, file = "basinNDams_L8.csv")
write.csv(basinMWDams_L8, file = "basinMWDams_L8.csv")


#Get every possible future scenario
damsall<- fread(file.path(resdir,"DamAttributes_L8_All_future.txt"),
      stringsAsFactors=T, data.table=T, integer64="numeric")
length(damsall[, unique(batRegion)])

#For each region, create a list scenario with all possible permutations
table(damsall[ESTAGIO_1 != "OperaÃ§Ã£o",.N, by=batRegion]$N)


allscenarios <- ldply(unique(damsall[ESTAGIO_1 != "OperaÃ§Ã£o", batRegion]), function(reg) {
  permut <- damsall[ESTAGIO_1 != "OperaÃ§Ã£o" & batRegion == reg, expand.grid(rep(list(c(0.1, 1)), .N))]
  names(permut) <- damsall[ESTAGIO_1 != "OperaÃ§Ã£o" & batRegion == reg, FID]
  scenarios_reg <- melt(cbind(region1, batRegion = paste(reg, rownames(region1), sep='_')), id.var="batRegion")
  return(scenarios_reg)
})

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



##### LEVEL 6 ######

# library graph theory
library(igraph)

# Identify possible scenatios
scenarios <- c("All_current.txt", "All_future.txt", "SHP_current.txt", "SHP_future.txt", "LHP_current.txt", "LHP_future.txt")

## Loop over scenarios
for (x in 1: length(scenarios)){
  tic = time()
  scenarioFileBasin <- file.path(resdir, paste("Brazil_L6_", scenarios[x], sep =""))
  scenarioFileDams <- file.path(resdir, paste("DamAttributes_L6_", scenarios[x], sep =""))
  
  ## Import BAT outputs
  # Import network data (fragments - nodes)
  NetworkBRAZILcrude <- read.csv(scenarioFileBasin, header = T)
  # Import dams data (links - edges)
  DamAttributes <- read.csv(scenarioFileDams, header = T)
  
  # Exclude reach IDs with problems (bugs in the polylines)
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61019316),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60934104),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60945892),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60901182),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60727967),]
  NetworkBRAZIL <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60700186),]
  
  # Organize the matrix based on type and stage of each dam
  levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) != "OperaÃ§Ã£o"] <- "Planned"
  levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) == "OperaÃ§Ã£o"] <- "Operation"
  levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) != "UHE"] <- "SHP"
  levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) == "UHE"] <- "LHP"
  
  ## Determine the basins present in the data       
  basinList <- unique(NetworkBRAZIL$Region6)
  
  ## Create matrixes to fill out with data just for the first scenario (not necessary for the others)
  if(x == 1){
    
    ## Create a matrix of basins and DCIs               
    DCI_BRAZIL_L6 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(DCI_BRAZIL_L6) <- c("HYBAS_ID", "All_curr", "All_fut", "SHP_curr", "SHP_fut", "LHP_curr", "LHP_fut")
    DCI_BRAZIL_L6 [ , 1] <- basinList
    
    ## Create two matrixes to fill with dam atributes of each basin (N dams, MW production)
    basinNDams_L6 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(basinNDams_L6) <- c("HYBAS_ID", "N_All_curr", "N_All_fut", "N_SHP_curr", "N_SHP_fut", "N_LHP_curr", "N_LHP_fut")
    basinNDams_L6 [ , 1] <- basinList
    
    basinMWDams_L6 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(basinMWDams_L6) <- c("HYBAS_ID", "MW_All_curr", "MW_All_fut", "MW_SHP_curr", "MW_SHP_fut", "MW_LHP_curr", "MW_LHP_fut")
    basinMWDams_L6 [ , 1] <- basinList
    
  }
  
  
  ## Loop over basins
  for (j in 1: length(basinList)){
    
    ## filter attributes of the basin j         
    BasinX <- NetworkBRAZIL[NetworkBRAZIL$Region6 == basinList[j], ]
    DamX <- DamAttributes[DamAttributes$batRegion == basinList[j], ]
    
    # Create a sequence ranging from 1 to the maximum number of segments
    Lsegments <- min(unique(BasinX$batNetID))
    Numbsegments <- length(unique(BasinX$batNetID))
    
    ## Determine that basins with no dams have DCI = 100
    if (Numbsegments == 1){
      DCI <- 100
    }
    print(Numbsegments)
    
    ## For basins with dams, DCI is calculated below
    if (Numbsegments > 1){
      
      # Paste the name seg before each segment number and make a vector
      # Create a vector of fragment lenghts
      listSeg <- rep (NA, times = Numbsegments)
      listLengths <- rep (NA, times = Numbsegments)
      
      ## Loop over segments
      # Fill these two vectors
      for (i in 1: Numbsegments){
        
        listSeg[i] <- paste("seg", Lsegments + i - 1, sep ="")
        listLengths[i] <- sum(BasinX$Shape_Leng[BasinX$batNetID == Lsegments + i - 1])
      }
      
      ## Compile the data on the edges (dam attributes)
      DowSeg <- paste("seg", DamX$Min_batNet, sep ="")
      UpSeg <- paste("seg", DamX$Max_batNet, sep ="")
      Type <- DamX$Tipo_1
      Situation <- DamX$ESTAGIO_1
      ID_number <- DamX$FID
      ID_name <- DamX$NOME
      Basin <- DamX$batRegion
      
      EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin)
      
      
      # Create a passibility vector based on the number of segments of the basin
      passVec <- rep(0.1, times = dim(EdgesData)[1])
      
      # List of edges and passability = 100
      d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[, 2], pass = passVec)
      
      # attributes of nodes; node sizes
      d3 = data.frame(id = listSeg, l = listLengths)
      
      
      # Run DCI analysis for the basin
      DCI <- DCIp (d2, d3, print = F)
      
      
    }
    
    ## Compute the DCI of the Basin in the matrix
    DCI_BRAZIL_L6[ j, x + 1] <- DCI
    
    ### Attributes of the basin (N of dams, cumulative MW)
    basinNDams_L6[j, x + 1] <- dim(DamX)[1]
    basinMWDams_L6[j, x + 1] <- sum(DamX$POT_KW, na.rm = T)/1000 #kw to MW
    
    print(x)
    print(j)
    
  }
}

## Create a csv file with the output of DCI analysis for all the scenarios
write.csv(DCI_BRAZIL_L6, file = "DCI_Brazil_L6.csv")
write.csv(basinNDams_L6, file = "basinNDams_L6.csv")
write.csv(basinMWDams_L6, file = "basinMWDams_L6.csv")





##### LEVEL 10 ######

# library graph theory
library(igraph)

# Identify possible scenatios
scenarios <- c("All_current.txt", "All_future.txt", "SHP_current.txt", "SHP_future.txt", "LHP_current.txt", "LHP_future.txt")

## Loop over scenarios
for (x in 1: length(scenarios)){
  scenarioFileBasin <- paste("Brazil_L10_", scenarios[x], sep ="")
  scenarioFileDams <- paste("DamAttributes_L10_", scenarios[x], sep ="")
  
  ## Import BAT outputs
  # Import network data (fragments - nodes)
  NetworkBRAZILcrude <- read.csv(scenarioFileBasin, header = T)
  # Import dams data (links - edges)
  DamAttributes <- read.csv(scenarioFileDams, header = T)
  
  # Exclude reach IDs with problems (bugs in the polylines)
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61436191),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61429804),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61406674),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61407347),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61401754),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61298170),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61195139),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61148586),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61007902),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60682708),]
  NetworkBRAZIL <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60352823),]
  
  # Organize the matrix based on type and stage of each dam
  levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) != "OperaÃ§Ã£o"] <- "Planned"
  levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) == "OperaÃ§Ã£o"] <- "Operation"
  levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) != "UHE"] <- "SHP"
  levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) == "UHE"] <- "LHP"
  
  ## Determine the basins present in the data       
  basinList <- unique(NetworkBRAZIL$Region10)
  
  ## Create matrixes to fill out with data just for the first scenario (not necessary for the others)
  if(x == 1){
    
    ## Create a matrix of basins and DCIs               
    DCI_BRAZIL_L10 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(DCI_BRAZIL_L10) <- c("HYBAS_ID", "All_curr", "All_fut", "SHP_curr", "SHP_fut", "LHP_curr", "LHP_fut")
    DCI_BRAZIL_L10 [ , 1] <- basinList
    
    ## Create two matrixes to fill with dam atributes of each basin (N dams, MW production)
    basinNDams_L10 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(basinNDams_L10) <- c("HYBAS_ID", "N_All_curr", "N_All_fut", "N_SHP_curr", "N_SHP_fut", "N_LHP_curr", "N_LHP_fut")
    basinNDams_L10 [ , 1] <- basinList
    
    basinMWDams_L10 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(basinMWDams_L10) <- c("HYBAS_ID", "MW_All_curr", "MW_All_fut", "MW_SHP_curr", "MW_SHP_fut", "MW_LHP_curr", "MW_LHP_fut")
    basinMWDams_L10 [ , 1] <- basinList
    
  }
  
  
  ## Loop over basins
  for (j in 1: length(basinList)){
    
    ## filter attributes of the basin j         
    BasinX <- NetworkBRAZIL[NetworkBRAZIL$Region10 == basinList[j], ]
    DamX <- DamAttributes[DamAttributes$batRegion == basinList[j], ]
    
    # Create a sequence ranging from 1 to the maximum number of segments
    Lsegments <- min(unique(BasinX$batNetID))
    Numbsegments <- length(unique(BasinX$batNetID))
    
    ## Determine that basins with no dams have DCI = 100
    if (Numbsegments == 1){
      DCI <- 100
    }
    
    
    ## For basins with dams, DCI is calculated below
    if (Numbsegments > 1){
      
      # Paste the name seg before each segment number and make a vector
      # Create a vector of fragment lenghts
      listSeg <- rep (NA, times = Numbsegments)
      listLengths <- rep (NA, times = Numbsegments)
      
      ## Loop over segments
      # Fill these two vectors
      for (i in 1: Numbsegments){
        
        listSeg[i] <- paste("seg", Lsegments + i - 1, sep ="")
        listLengths[i] <- sum(BasinX$Shape_Leng[BasinX$batNetID == Lsegments + i - 1])
      }
      
      ## Compile the data on the edges (dam attributes)
      DowSeg <- paste("seg", DamX$Min_batNet, sep ="")
      UpSeg <- paste("seg", DamX$Max_batNet, sep ="")
      Type <- DamX$Tipo_1
      Situation <- DamX$ESTAGIO_1
      ID_number <- DamX$FID
      ID_name <- DamX$NOME
      Basin <- DamX$batRegion
      
      EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin)
      
      
      # Create a passibility vector based on the number of segments of the basin
      passVec <- rep(0.1, times = dim(EdgesData)[1])
      
      # List of edges and passability = 100
      d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[, 2], pass = passVec)
      
      # attributes of nodes; node sizes
      d3 = data.frame(id = listSeg, l = listLengths)
      
      
      # Run DCI analysis for the basin
      DCI <- DCIp (d2, d3, print = F)
      
      
    }
    
    ## Compute the DCI of the Basin in the matrix
    DCI_BRAZIL_L10[ j, x + 1] <- DCI
    
    ### Attributes of the basin (N of dams, cumulative MW)
    basinNDams_L10[j, x + 1] <- dim(DamX)[1]
    basinMWDams_L10[j, x + 1] <- sum(DamX$POT_KW, na.rm = T)/1000 #kw to MW
    
    print(x)
    print(j)
    
  }
}


## Create a csv file with the output of DCI analysis for all the scenarios
write.csv(DCI_BRAZIL_L10, file = "DCI_Brazil_L10.csv")
write.csv(basinNDams_L10, file = "basinNDams_L10.csv")
write.csv(basinMWDams_L10, file = "basinMWDams_L10.csv")






##### LEVEL 4 ######


# library graph theory
library(igraph)

# Identify possible scenatios
scenarios <- c("All_current.txt", "All_future.txt", "SHP_current.txt", "SHP_future.txt", "LHP_current.txt", "LHP_future.txt")

## Loop over scenarios
for (x in 1: length(scenarios)){
  scenarioFileBasin <- paste("Brazil_L4_", scenarios[x], sep ="")
  scenarioFileDams <- paste("DamAttributes_L4_", scenarios[x], sep ="")
  
  ## Import BAT outputs
  # Import network data (fragments - nodes)
  NetworkBRAZILcrude <- read.csv(scenarioFileBasin, header = T)
  # Import dams data (links - edges)
  DamAttributes <- read.csv(scenarioFileDams, header = T)
  
  ### Exclude reach IDs with problems (bugs in the polylines)
  # To many segments (I created a loop for excluding 209 segments)
  Lixo <- read.csv("DELETE.txt", header = T)
  for (i in 1:length(Lixo$REACH_ID)){
    NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == Lixo$REACH_ID[i]),]
  }
  
  # Exclude reach IDs with problems (bugs in the polylines)
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61308381),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61308380),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309965),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310036),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310095),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310566),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310637),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310483),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310094),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310214),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310093),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61310035),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309964),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309600),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309599),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309745),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309744),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309453),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61309452),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61308695),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61308694),]
  
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61323794),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61323795),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61323796),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61323858),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61324506),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61324565),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61324938),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61324937),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325218),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325219),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325133),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325351),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325486),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325683),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325555),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325485),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325484),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325554),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325553),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325737),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325990),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325736),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325735),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325920),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61326511),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61325989),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61326451),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61326362),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61327000),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61327055),]
  
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61295884),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61295753),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61295593),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61295658),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61298462),]
  
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60802917),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60584495),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 61154949),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60727967),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60576506),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60295190),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60227082),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60367435),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60282029),]
  NetworkBRAZILcrude <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60504992),]
  NetworkBRAZIL <- NetworkBRAZILcrude[-which(NetworkBRAZILcrude$REACH_ID == 60512237),]
  
  
  # Organize the matrix based on type and stage of each dam
  levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) != "OperaÃ§Ã£o"] <- "Planned"
  levels(DamAttributes$ESTAGIO_1)[levels(DamAttributes$ESTAGIO_1) == "OperaÃ§Ã£o"] <- "Operation"
  levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) != "UHE"] <- "SHP"
  levels(DamAttributes$Tipo_1)[levels(DamAttributes$Tipo_1) == "UHE"] <- "LHP"
  
  ## Determine the basins present in the data       
  basinList <- unique(NetworkBRAZIL$Region4)
  
  ## Create matrixes to fill out with data just for the first scenario (not necessary for the others)
  if(x == 1){
    
    ## Create a matrix of basins and DCIs               
    DCI_BRAZIL_L4 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(DCI_BRAZIL_L4) <- c("HYBAS_ID", "All_curr", "All_fut", "SHP_curr", "SHP_fut", "LHP_curr", "LHP_fut")
    DCI_BRAZIL_L4 [ , 1] <- basinList
    
    ## Create two matrixes to fill with dam atributes of each basin (N dams, MW production)
    basinNDams_L4 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(basinNDams_L4) <- c("HYBAS_ID", "N_All_curr", "N_All_fut", "N_SHP_curr", "N_SHP_fut", "N_LHP_curr", "N_LHP_fut")
    basinNDams_L4 [ , 1] <- basinList
    
    basinMWDams_L4 <- matrix(NA, nrow = length(basinList), ncol = 7)
    colnames(basinMWDams_L4) <- c("HYBAS_ID", "MW_All_curr", "MW_All_fut", "MW_SHP_curr", "MW_SHP_fut", "MW_LHP_curr", "MW_LHP_fut")
    basinMWDams_L4 [ , 1] <- basinList
    
  }
  
  
  ## Loop over basins
  for (j in 1: length(basinList)){
    
    ## filter attributes of the basin j         
    BasinX <- NetworkBRAZIL[NetworkBRAZIL$Region4 == basinList[j], ]
    DamX <- DamAttributes[DamAttributes$batRegion == basinList[j], ]
    
    # Create a sequence ranging from 1 to the maximum number of segments
    Lsegments <- min(unique(BasinX$batNetID))
    Numbsegments <- length(unique(BasinX$batNetID))
    
    ## Determine that basins with no dams have DCI = 100
    if (Numbsegments == 1){
      DCI <- 100
    }
    
    
    ## For basins with dams, DCI is calculated below
    if (Numbsegments > 1){
      
      # Paste the name seg before each segment number and make a vector
      # Create a vector of fragment lenghts
      listSeg <- rep (NA, times = Numbsegments)
      listLengths <- rep (NA, times = Numbsegments)
      
      ## Loop over segments
      # Fill these two vectors
      for (i in 1: Numbsegments){
        
        listSeg[i] <- paste("seg", Lsegments + i - 1, sep ="")
        listLengths[i] <- sum(BasinX$Shape_Leng[BasinX$batNetID == Lsegments + i - 1])
      }
      
      ## Compile the data on the edges (dam attributes)
      DowSeg <- paste("seg", DamX$Min_batNet, sep ="")
      UpSeg <- paste("seg", DamX$Max_batNet, sep ="")
      Type <- DamX$Tipo_1
      Situation <- DamX$ESTAGIO_1
      ID_number <- DamX$FID
      ID_name <- DamX$NOME
      Basin <- DamX$batRegion
      
      EdgesData <- data.frame(DowSeg, UpSeg, Type, Situation, ID_number, ID_name, Basin)
      
      
      # Create a passibility vector based on the number of segments of the basin
      passVec <- rep(0.1, times = dim(EdgesData)[1])
      
      # List of edges and passability = 100
      d2 = data.frame (id1 = EdgesData[, 1], id2 = EdgesData[, 2], pass = passVec)
      
      # attributes of nodes; node sizes
      d3 = data.frame(id = listSeg, l = listLengths)
      
      
      # Run DCI analysis for the basin
      DCI <- DCIp (d2, d3, print = F)
      
      
    }
    
    ## Compute the DCI of the Basin in the matrix
    DCI_BRAZIL_L4[ j, x + 1] <- DCI
    
    ### Attributes of the basin (N of dams, cumulative MW)
    basinNDams_L4[j, x + 1] <- dim(DamX)[1]
    basinMWDams_L4[j, x + 1] <- sum(DamX$POT_KW, na.rm = T)/1000 #kw to MW
    
    print(x)
    print(j)
    
    
  }
}


## Create a csv file with the output of DCI analysis for all the scenarios
write.csv(DCI_BRAZIL_L4, file = "DCI_Brazil_L4.csv")
write.csv(basinNDams_L4, file = "basinNDams_L4.csv")
write.csv(basinMWDams_L4, file = "basinMWDams_L4.csv")




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



