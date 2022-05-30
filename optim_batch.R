###################################
## Main Algorithm (He Geng 2008) ##
###################################

find_optim_sets = function(G_tau){
  
  library(abind)
  
  # INPUT:
  
  # G_tau : adjacency matrix of G_tau decomposable sub-graph of G for chain component tau
  
  # OUTPUT:
  
  # S.optimal : a list collecting the optimal sets
  
  h = dim(G_tau)[1]
  
  nodes = 1:h # by default node-labels are re-assigned as 1,...,h
  
  if(h == 1){stop("No edges to orient: chain component has size one")}
  
  if(h == 2){
    
    out = matrix(nodes, nrow = 1)
    
    out_list = split(out, rep(1:ncol(out), each = nrow(out)))
    
    return(list(G_tau = G_tau, S.optimal = out_list))
    
  }
  
  if(h > 2){
    
    count.duplicates <- function(X){
      
      if(is.null(dim(X))){
        
        return(table(X))
        
      }
      
      x <- do.call('paste', c(X, sep = '\r'))
      ox <- order(x)
      rl <- rle(x[ox])
      
      cbind(X[ox[cumsum(rl$lengths)],,drop=FALSE], count = rl$lengths)  
      
    }
    
    # Find all DAGs compatible with G_tau
    
    all.dags = pdag2allDags(gm = G_tau)$dags
    
    all.dags.tab = array(t(all.dags), c(h, h, nrow(all.dags)))
    
    k = 1; is.optim.k = FALSE
    
    while(is.optim.k == FALSE){
      
      # Find all sets of size k
      
      sets = as.matrix(combn(x = 1:h, m = k)) # each column denotes one possible set (sequence of manipulated nodes) S
      
      is.optim.all.k = rep(FALSE, ncol(sets))
      
      for(i in 1:ncol(sets)){
        
        S = sets[,i]
        
        out_array_tmp = abind(array(all.dags.tab[S,,], c(k, h, nrow(all.dags))),
                              array(all.dags.tab[,S,], c(k, h, nrow(all.dags))), along = 1)
        
        out_mat_tmp = t(apply(out_array_tmp, 3L, c))
        
        is.optim.all.k[i] = all(count.duplicates(data.frame(out_mat_tmp))$count == 1)
        
      }
      
      # I check at the end of the loop since for the same size k there may be more than one sufficient set,
      # but the optimal must have the smallest size; therefore I stop if any set of size k is sufficient
      
      is.optim.k = any(is.optim.all.k == TRUE)
      
      k = k + 1
      
    }
    
    out = as.matrix(sets[,is.optim.all.k])
    
    out_list = split(out, rep(1:ncol(out), each = nrow(out)))
    
    return(list(G_tau = G_tau, S.optimal = out_list))
    
  }
  
}
