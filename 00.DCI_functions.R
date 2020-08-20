#Developers: Thiago B. Couto; Mathis L. Messager
#Creation date: 
#Purpose: Functions to compute DCI for a given river basin


#Define directory structure
rootdir <- find_root(has_dir("src"))
resdir <- file.path(rootdir, "results")
figdir <- file.path(resdir, 'figures')
datadir <- file.path(rootdir, "data")
dcigdb <- file.path(resdir, 'dci.gdb')

if (!file.exists(figdir)) {
  dir.create(figdir)
}

#DCI for potadromous species
DCIp_opti5 <- function(d2, d3, print = NULL){
  "d2 = a data frame containing the links between patches and passability for each link
  d2$id1 = initial segment (FROM) connection; MUST BE FACTOR!
  d2$id2 = final segment (to) connection; MUST BE FACTOR!
  d2$pass = the passability from one segment to the next
  
  d3 = data frame containing all segments and its respective sizes
  d3$id = the ID of each segment; MUST BE FACTOR!
  d3$l = the SIZE of segments formed by barriers within a given river
  
  OBS: d2 and d3 must have the same segment names (ID)
  d = edges; vertices = nodes"
  
  graph <- igraph::graph.data.frame(d2, directed = F); #plot(graph)
  
  #calculating l_i/L ratio - for each segment, the sum of that segment to that of all segments in basin
  d3[,l_L := l/sum(l)] %>%
    setkey(id)
  
  if (max(d2$pass) == 0) {
    
    return(d3[,100*sum(l_L^2)])
    
  } else {
    
    #Get pair-wise DCI for lower triangle of matrix of all pairwise combination of segments (excluding diagonal)
    lowertricomb <- d3[, list(RcppAlgos::comboGeneral(id, 2))] %>% #Get all combinations of ids
      setDT %>% setnames(c('V1', 'V2'), c('i', 'k')) %>%
      setkey(i) %>%  .[d3] %>%
      setkey(k) %>%  .[d3] %>%
      .[!is.na(i)] %>%
      .[, 100 * l_L * i.l_L * prod(igraph::shortest_paths(graph, from = i, to = k, output = "epath")$epath[[1]]$pass**2), #Compute DCI from each segment i to all other segments
        by=1:nrow(.)]
    
    #Compute DCI by sum connectivity for full matrix of pairwise combination 
    return(2*lowertricomb[,sum(V1)]+d3[,100*sum(l_L^2)]) #Sum DCI across matrix of pairwise combinations of segments (2*lower triangle + diagonal)
  }
}

#DCI for diadromous species
#Calculate connectivity in terms of the probability that a fish can move in both directions between 
#the mouth of the river and another section of the river network
DCIi_opti <- function(d2, d3, print = NULL){
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
  
  graph <- igraph::graph.data.frame(d2, directed = F); #plot(graph)
  
  #calculating l_i/L ratio - for each segment, the sum of that segment to that of all segments in basin
  d3$l_L <- d3$l/sum(d3$l)
  
  #get segment at mouth of network
  mouthsegID <- unique(d2$id1)[!(unique(d2$id1) %in% unique(d2$id2))]
  
  #Get lower triangle of matrix of all pairwise combination of mouth by segments
  mouthcomb <- data.table(i= rep(mouthsegID, length(d3$id)), j=d3$id)
  
  DCImatd <- function(x) {
    #Compute connectivity for all unique pairs of segments in basin (for diadromous species)
    #Note: when from= and to= are the same in shortest_paths, returns an empty numeric that, when fed to prod, returns 1
    return(prod(igraph::shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]$pass**2) *
             d3[d3$id==x[2],'l_L'][[1]] *
             100) 
  }
  
  return(sum(apply(mouthcomb, 1, DCImatd)))
}


##########################
# For troubleshooting
#########################
# library(rprojroot)
# library(microbenchmark)
# rootdir <- find_root(has_dir("src"))
# resdir <- file.path(rootdir, "results")
# dcigdb <- file.path(resdir, 'dci.gdb')
# DamAttributes <- read.fst(file.path(resdir, 'DamAttributes.fst')) %>% setDT
# NetworkBRAZIL <- read.fst(file.path(resdir, 'NetworkBRAZIL.fst')) %>% setDT
# 
# j <- '60808111301'
# d3 <- NetworkBRAZIL[HYBAS_ID08ext==j, list(id=as.character(SEGID), l=Shape_Length)]
# d2 <-  DamAttributes[DAMBAS_ID08ext ==j, list(id1 = DownSeg, id2 = UpSeg, pass = Allcurrent)]
# 
# microbenchmark(
#   DCIi_opti(d2, d3, print = NULL),
#   DCIp_opti5(d2, d3, print = NULL),
#   times=50L
# )
