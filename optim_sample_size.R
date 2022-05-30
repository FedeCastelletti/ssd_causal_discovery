find_BF_thresholds = function(prior_01, gamma_0, gamma_1){
  
  # INPUT:
  
  # prior_01 : list with two elements representing the prior probability of H0 (j -> i) and H1 (j <- i)
  # gamma_0, gamma_1 : thresholds for posterior probability of H0 and H1, e.g. both equal to 0.9
  
  # OUTPUT:
  
  # k_0, k_1 : BF thresholds under H0 and H1 respectively
  
  # N.B. Thresholds are k_0 and 1/k_1
  
  omega = prior_01[[2]]/prior_01[[1]]; names(omega) = NULL # prior ratio
  
  k_0 = omega*gamma_0/(1 - gamma_0)
  k_1 = omega*gamma_1/(1 - gamma_1)
  
  return(list(k_0 = k_0, k_1 = k_1))
  
}

n_optim = function(nodes, k_01, G_tau, beta, H = H, N = N, S_tau, g = g){
  
  source("prior_H0_H1.r")
  source("predictive_BF.r")
  
  # INPUT :
  
  # nodes : (2,1) vector with labels of nodes i and j respectively
  # k_01  : BF thresholds (under H0 and H1 respectively)
  # beta  : threshold for probability of decisive evidence
  # G_tau : adjacency matrix of G_tau decomposable sub-graph of G for chain component tau
  # N     : number of historical data (Z)
  # S_tau : (q,q) matrix t(Z)%*%Z restricted to chain component tau
  
  # OUTPUT :
  
  # n_star : optimal sample size
  
  k0 = k_01[[1]]
  k1 = k_01[[2]]
  
  card.tau = dim(G_tau)[1]
  
  
  # Construction of all DAGs compatible with chain component can be done outside this algorithm and just one time
  # (no need to repeat it for each i,j)
  
  out_prob_comp = function(...,n){
    
    out_tmp_0 = exp(BF.predictive.ji(ij = nodes, card.tau, N = N, S_tau, a_omega = card.tau, g = g, n = n, H = H, log.BF = T))
    
    out_tmp_1 = exp(BF.predictive.ij(ij = nodes, card.tau, N = N, S_tau, a_omega = card.tau, g = g, n = n, H = H, log.BF = T))
    
    prob_decisive_BF_0 = sum(out_tmp_0 > k0)/length(out_tmp_0)
    
    prob_decisive_BF_1 = sum(out_tmp_1 < 1/k1)/length(out_tmp_1)
    
    prior_01 = find_prior_01(G_tau = G_tau, i = nodes[1], nodes[2])
    
    p0 = prior_01[[1]]; p1 = prior_01[[2]]
    
    out_prob_decisive = p0*prob_decisive_BF_0 + p1*prob_decisive_BF_1
    
    return(out_prob_decisive)
    
  }
  
  u = 100
  l = 2
  m = 2
  
  out_prob_decisive_l = out_prob_comp(n = l)
  out_prob_decisive_u = out_prob_comp(n = u)
  
  while(out_prob_decisive_u < beta){
    
    u = 2*u
    out_prob_decisive_u = out_prob_comp(n = u)
    
  }
  
  #while((out_prob_decisive_l < beta) & (out_prob_decisive_u > beta) & (l != u)){
  while((out_prob_decisive_l < beta) & (out_prob_decisive_u > beta) & (abs(l-u)>1)){
    m = round((l+u)/2)
    
    out_prob_decisive_m = out_prob_comp(n = m)
    
    if(out_prob_decisive_m > beta){
      
      u = m
      l = l
      
    }
    
    if(out_prob_decisive_m < beta){
      
      u = u
      l = m
      
    }
    
  }
  
  return(list(n_star = m, nodes = nodes))
  
}


n_optim_chain_comp = function(G_tau, gamma_0, gamma_1, beta, H = H, N = N, S_tau = S_tau, g = g){
  
  # beta       : threshold for probability of correct and decisive evidence
  
  # gamma_0, gamma_1 : thresholds for posterior probabilities of (correctly accepting) H0 and H1
  
  source("optim_batch.r")
  source("prior_H0_H1.r")
  
  optim_sets = find_optim_sets(G_tau)
  
  L = length(optim_sets$S.optimal)
  
  out_all = as.list(rep(NA, L))
  
  for(l in 1:L){
    
    K = length(optim_sets$S.optimal[[l]])
    
    out_l = list()
    
    for(h in 1:K){ # for each intervened node i in the sequence of manipulated variables:
      
      i = optim_sets$S.optimal[[l]][h]
      
      i_ne = which(G_tau[i,] > 0)
      
      out_l_i = c()
      
      for(s in 1:length(i_ne)){ # for each neighbor j of i:
        
        j = i_ne[s]
        
        prior_01 = find_prior_01(G_tau = G_tau, i = i, j = j)
        
        k_01 = find_BF_thresholds(prior_01, gamma_0 = gamma_0, gamma_1 = gamma_1)
        
        out_l_i = c(out_l_i, n_optim(nodes = c(i,j), k_01 = k_01, G_tau = G_tau, beta = beta, H = H, N = N, S_tau = S_tau, g = g)$n_star)
        
      }
      
      names(out_l_i) = i_ne; out_l_i = list(out_l_i)
      
      out_l = append(out_l, out_l_i)
      
    }
    
    names(out_l) = optim_sets$S.optimal[[l]]
    
    out_all[[l]] = out_l
    
  }
  
  names(out_all) = c(lapply(optim_sets$S.optimal, paste, collapse = " "))
  
  return(list(out_all = out_all, optim_sets = optim_sets))
  
}

