#Developers: Thiago B. Couto; Mathis L. Messager
#Creation date: 
#Purpose: Functions to compute DCI for a given river basin

library(igraph) # library graph theory
library(data.table)
library(plyr)
library(magrittr)
library(stringr)
library(RcppAlgos)


#DCI for potadromous species
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
  resudia <- lapply(d3$id, DCImatpself)
  #Compute DCI by sum connectivity for full matrix of pairwise combination 
  return(sum(resutri)*2+sum(unlist(resudia)))
}

DCIp_opti2 <- function(d2, d3, print = NULL){
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
  d3[,l_L := l/sum(l)]
  
  #Get lower triangle of matrix of all pairwise combination of segments
  lowertricomb <- d3[, list(RcppAlgos::comboGeneral(id, 2))] %>% setDT %>% setnames(c('V1', 'V2'), c('i', 'k'))
  
  DCImatp <- function(x) {
    #Compute connectivity for all unique pairs of segments in basin (for potadromous species)
    return(prod(shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]$pass)*
             prod(d3[d3$id==x[1],'l_L'], d3[d3$id==x[2],'l_L']) *
             100)
  }
  resutri <- apply(lowertricomb, 1, DCImatp)
  
  #Compute DCI by sum connectivity for full matrix of pairwise combination 
  return(sum(resutri)*2+d3[,100*sum(l_L^2)])
}

DCIp_opti3 <- function(d2, d3, print = NULL){
  "d2 = a data frame containing the links between patches and passability for each link
  d2$id1 = initial segment (FROM) connection; MUST BE FACTOR!
  d2$id2 = final segment (to) connection; MUST BE FACTOR!
  d2$pass = the passability from one segment to the next
  
  d3 = data frame containing all segments and its respective sizes
  d3$id = the ID of each segment; MUST BE FACTOR!
  d3$l = the SIZE of segments formed by barriers within a given river
  
  OBS: d2 and d3 must have the same segment names (ID)
  d = edges; vertices = nodes"
  
  graph <- graph.data.frame(d2, directed = F); #plot(graph)
  
  #calculating l_i/L ratio - for each segment, the sum of that segment to that of all segments in basin
  d3[,l_L := l/sum(l)]
  
  #Get pair-wise DCI for lower triangle of matrix of all pairwise combination of segments (excluding diagonal)
  lowertricomb <- d3[, list(RcppAlgos::comboGeneral(id, 2))] %>% #Get all combinations of ids
    setDT %>% setnames(c('V1', 'V2'), c('i', 'k')) %>% #Convert to data.table and rename columns
    merge(d3, by.x='i', by.y='id', all.y=F) %>% #Get size of segments (may be sped up with key-based match)
    merge(d3, by.x='k', by.y='id', all.y=F) %>%
    .[, 100 * l_L.x * l_L.y * prod(shortest_paths(graph, from = i, to = k, output = "epath")$epath[[1]]$pass), #Compute DCI from each segment i to all other segments
      by=seq_len(nrow(.))]
  
  #Compute DCI by sum connectivity for full matrix of pairwise combination 
  return(lowertricomb[, 2*sum(V1)]+d3[,100*sum(l_L^2)]) #Sum DCI across matrix of pairwise combinations of segments (2*lower triangle + diagonal)
}

DCIp_opti4 <- function(d2, d3, print = NULL){
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
  d3[,l_L := l/sum(l)]
  
  #Get lower triangle of matrix of all pairwise combination of segments
  lowertricomb <- d3[, list(RcppAlgos::comboGeneral(id, 2))] %>% setDT %>% setnames(c('V1', 'V2'), c('i', 'k')) %>%
    merge(d3, by.x='i', by.y='id', all.y=F) %>%
    merge(d3, by.x='k', by.y='id', all.y=F) 
  
  lowercomp <- unlist(lapply(lowertricomb[, unique(i)], function(c) {
    lowertricomb[i==c, 100 * l_L.x * l_L.y] *
      unlist(
        lapply(shortest_paths(graph, from = c, to = lowertricomb[i==c,k], output = "epath")$epath, function(epath) {
          prod(epath$pass)
        })
      )
  })) 
  
  #Compute DCI by sum connectivity for full matrix of pairwise combination 
  return(2*sum(lowercomp)+d3[,100*sum(l_L^2)])
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