#saveRDS(list(d2, d3), 'DCIpbenchmarking.rds')
imp <- readRDS('DCIpbenchmarking.rds')
d2 <- imp[[1]]
d3 <- imp[[2]]
library(microbenchmark)

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




DCIp_all <- function(d2, d3, print = NULL){
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
  patches <- as.character(d3$id)
  
  resu <- outer(patches, patches, FUN = Vectorize(function(i,k){
    caminho_edge <- shortest_paths(graph, from = i, to = k, output = "epath")$epath[[1]]
    pass <- E(graph)[caminho_edge]$pass
    c_ik <- prod(pass);
    
    d_subI <- d3[d3$id == i ,]
    d_subK <- d3[d3$id == k ,]
    return(c_ik * prod(d_subI$l_L, d_subK$l_L) * 100)
  } ))
  return(sum(resu))
}

DCIp_lowertri <- function(d2, d3, print = NULL){
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
  
  #Get all combinations (lower triangular matrix)
  indices <- combi2inds(d3$id)
  lowertricomb <- setDT(list(i=d3$id[indices$xid], k=d3$id[indices$yid])) #Missing matrix diagonal
  
  DCImat <- function(i, k) {
    caminho_edge <- shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]
    pass <- E(graph)[caminho_edge]$pass
    c_ik <- prod(pass)
    
    d_subI <- d3[d3$id == x[1] ,]
    d_subK <- d3[d3$id == x[2] ,]
    return(c_ik * prod(d_subI$l_L, d_subK$l_L) * 100)
  }
  
  DCImatself <- function(x) {
    return(1 * (d3[d3$id == x , 'l_L']^2) * 100)
  }
  resutri <- apply(lowertricomb, 1, DCImat)
  resudia <- sapply(d3$id, DCImatself)
  return(sum(resutri)*2+sum(resudia))
}

d3dt <- as.data.table(d3)
DCIp_lowertridt <- function(d2, d3, print = NULL){
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
  d3[, l_L := l/sum(l)]
  
  #Get all combinations (lower triangular matrix)
  indices <- combi2inds(d3$id)
  lowertricomb <- setDT(list(i=d3$id[indices$xid], k=d3$id[indices$yid])) #Missing matrix diagonal
  
  DCImat <- function(i, k) {
    caminho_edge <- shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]
    pass <- E(graph)[caminho_edge]$pass
    c_ik <- prod(pass)
    return(c_ik * prod(d3[id==x[1],l_L], d3[id==x[2],l_L]) * 100)
  }
  
  DCImatself <- function(x) {
    return(1 * (d3[id == x , 'l_L']^2) * 100)
  }
  resutri <- apply(lowertricomb, 1, DCImat)
  resudia <- sapply(d3$id, DCImatself)
  return(sum(resutri)*2+sum(resudia))
}

DCIp_lowertrishort <- function(d2, d3, print = NULL){
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
  
  #Get all combinations (lower triangular matrix)
  indices <- combi2inds(d3$id)
  lowertricomb <- setDT(list(i=d3$id[indices$xid], k=d3$id[indices$yid])) #Missing matrix diagonal
  
  DCImat <- function(x) {
    return(prod(shortest_paths(graph, from = x[1], to = x[2], output = "epath")$epath[[1]]$pass) *
             prod(d3[d3$id==x[1],'l_L'], d3[d3$id==x[2],'l_L']) *
             100) 
    }
  
  DCImatself <- function(x) {
    return(1 * (d3[d3$id == x , 'l_L']^2) * 100)
  }
  resutri <- apply(lowertricomb, 1, DCImat)
  resudia <- sapply(d3$id, DCImatself)
  return(sum(resutri)*2+sum(resudia))
}


microbenchmark(
  DCIp1 <- DCIp_all(d2, d3, print = F),
  DCIp2 <- DCIp_lowertri(d2, d3, print = F),
  DCIp3 <- DCIp_lowertridt(d2, d3dt, print = F),
  DCIp4 <- DCIp_lowertrishort(d2, d3, print = F),
  times=10
)