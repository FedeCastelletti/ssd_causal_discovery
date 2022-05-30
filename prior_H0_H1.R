########################################################
## Evaluating prior probability for j -> i and j <- i ##
########################################################

find_prior_01 = function(G_tau, i, j){
  
  # INPUT:
  
  # G_tau : adjacency matrix of G_tau decomposable sub-graph of G for chain component tau
  
  # i,j labels of the two nodes
  
  # OUTPUT:
  
  # p(j -> i), p(j <- i); prior probabilities under H0 and H1 respectively
  
  h = dim(G_tau)[1]
  
  nodes = 1:h # by default node-labels are re-assigned as 1,...,h
  
  if(G_tau[i,j] == 0){stop("No undirected link between the two nodes")}
  
  if(h == 1){stop("No edges within the chain component: chain component has size one")}
  
  if(h == 2){
    
    out = c(0.5,0.5)
    
    names(out) = c("H0", "H1")
    
    out_list = split(out, rep(1:length(out)))
    
    return(out_list)
    
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
    
    
    # Select orientation of edge between i and j in the DAGs
    
    out_array_tmp = array(all.dags.tab[c(i,j),c(i,j),], c(2, 2, nrow(all.dags)))
    
    out_mat_tmp = t(apply(out_array_tmp, 3L, c))
    
    out_tab_tmp = count.duplicates(data.frame(out_mat_tmp))
    
    index_ji = which(out_tab_tmp[,2] == 1) # j -> i
    index_ij = which(out_tab_tmp[,2] == 0) # j <- i
    
    out = c()
    
    out[1] = out_tab_tmp$count[index_ji]/sum(out_tab_tmp$count)
    out[2] = out_tab_tmp$count[index_ij]/sum(out_tab_tmp$count)
    
    names(out) = c("H0", "H1")
    
    out_list = split(out, rep(1:length(out)))
    
    return(out_list)
    
  }
  
}