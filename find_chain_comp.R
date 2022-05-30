library(gRbase)
library(ggm)
library(igraph)

find.chain.comp = function(ess){
  
  # Given an EG with adjacency matrix "ess" returns its set of chain components Tau
  
  ess = as(ess,"igraph")
  
  amat = as.matrix(get.adjacency(ess))  # if the argument is ess.obs (the essential graph, a graph object)
  wmat = matrix(as.integer((amat + t(amat)) > 1), nrow = nrow(amat))
  wg = graph.adjacency(wmat, mode = "undirected")
  cc = clusters(wg)
  neworder <- sort.list(cc$membership, na.last = NA)
  a = matrix(0, nrow = length(cc$csize), ncol = length(cc$csize))
  b = cumsum(cc$csize)
  wmat = amat[neworder, neworder]
  
  for(i in 1: length(cc$csize)){
    for(j in 1: length(cc$csize)){
      if(j != i){
        a[i,j] = as.integer(sum(wmat[(max(b[i-1],0)+1):b[i],
                                     (max(b[j-1],0)+1):b[j]]) > 0)
      }
    }
  }
  
  rownames(a) = colnames(a) = as.character(1:length(b))
  
  chainorder = topOrder(a)
  vertorder = c()
  chainsize = c()
  
  for(k in 1:length(b)){
    vertorder = c(vertorder, which(cc$membership == chainorder[k]))
    chainsize = c(chainsize, cc$csize[chainorder[k]])
  }
  
  q = list(vert.order=vertorder,chain.size=chainsize)
  order = q$vert.order
  size = c(1,q$chain.size)
  
  Tau = list(0)
  for(i in 2:(length(size))){
    Tau[[i-1]] = order[sum(size[1:(i-1)]):(sum(size[2:i]))]
  }
  
  return(Tau)
  
}
