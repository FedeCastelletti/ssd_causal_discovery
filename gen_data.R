library(pcalg)
library(mvtnorm)

gen.dag.parameters = function(q, prob, inf.lim, sup.lim, var.cond){
  
  # INPUT:
  
  # q    : number of variables/nodes
  # prob : probability of edge inclusion
  
  # inf.lim  : lower bound for random generation of elements in L (regression coefficients)
  # sup.lim  : upper bound for random generation of elements in L (regression coefficients)
  # var.cond : conditional variances of the q variables (diagonal of D)
  
  # OUTPUT:
  
  # DAG   : (q,q) adjacency matrix of the generated DAG
  # L, D  : DAG parameters
  
  DAG = randomDAG(q, prob = prob)
  
  DAG = t(as(DAG, "matrix")); DAG[DAG != 0] = 1
  
  L = DAG*matrix(runif(q*q, inf.lim, sup.lim), q, q)*sample(c(-1, 1), size = q*q, replace = TRUE)
  
  diag(L) = 1
  
  D = diag(var.cond, q, q)
  
  return(list(DAG = DAG, L = L, D = D))
  
}


gen.obs.dataset = function(L, D, N){
  
  # INPUT:
  
  # L   : (q,q) matrix with regression coefficients in the observational data generating model X = LX + eps
  # D   : (q,q) diagonal matrix with conditional variances of of X_1, ..., X_q
  # N   : number of observations
  
  # OUTPUT:
  
  # Z   : a (N,q) observational dataset
  
  Sigma = solve(t(L))%*%D%*%solve(L)
  
  Z = rmvnorm(N, rep(0,q), Sigma)
  
  m = colMeans(Z)
  s = apply(X = Z, FUN = sd, MARGIN = 2)
  Z = t((t(Z) - m)/s)
  
  return(list(Z = Z, L = L, D = D))
  
}
