# We consider BF of i <- j against i -> j

## Algorithm 1 [Predictive distribution of BF under i -> j, the alternative Hp H1]

library(pcalg)
library(mvtnorm)

BF.predictive.ij = function(ij, card.tau, N, S, a_omega, g, n, H, log.BF = FALSE){
  
  # INPUT:
  
  # ij       : (2,1) vector with numerical labels of nodes i and j in the covariance matrix S_tau
  # card.tau : size of chain component i and j belong to
  # N        : number of historical data (number of rows of data matrix Z)
  # S_tau    : (q,q) matrix t(Z)%*%Z restricted to chain component tau
  # a_omega  : hyperparameter of the Wishart prior
  # g        : variance of the interventional density f.tilde
  # n        : sample size
  
  # OUTPUT:
  
  # posterior predictive of BF
  
  # check ?rWishart
  
  i = ij[1]; j = ij[2]
  
  
  sample.r.ij = function(...){
    
    Q.ij = rWishart(1, df = a_omega + N - (card.tau - 2), Sigma = solve(S_tau[ij,ij]))[,,1]
    
    Sigma.ij = solve(Q.ij)
    
    L.ij = -Sigma.ij[i,j]/(Sigma.ij[i,i])
    D.jj = Sigma.ij[j,j] - Sigma.ij[i,j]^2/Sigma.ij[i,i]
    
    xi.tilde = rnorm(n, 0, sqrt(g))
    
    xj = rnorm(n, -L.ij*xi.tilde, sqrt(D.jj))
    
    X.ij = cbind(xi.tilde, xj)
    
    P = t(X.ij)%*%X.ij
    
    r.ij.2 = P[i,j]^2/(P[i,i]*P[j,j])
    
    return(r.ij.2)
    
  }
  
  r.ij.2.sample = sapply(X = 1:H, FUN = sample.r.ij)
  
  if(log.BF == FALSE){
    
    BF.sample = pi^(-0.5)*gamma(n/2)/gamma((n + 1)/2)*(n + 1)*((1 - r.ij.2.sample)^((n - 1)/2))
    
  }
  
  if(log.BF == TRUE){
    
    BF.sample = -0.5*log(pi) + lgamma(n/2) - lgamma((n+1)/2) + log(n+1) + ((n-1)/2)*log(1 - r.ij.2.sample)
    
  }
  
  return(BF.sample)
  
}


## Algorithm 2 [Predictive distribution of BF under D: i <- j, the null Hp H0]

BF.predictive.ji = function(ij = NULL, card.tau = NULL, N = NULL, S = NULL, a_omega = NULL, g = NULL, n, H, log.BF = FALSE){
  
  # INPUT:
  
  # ij       : (2,1) vector with numerical labels of nodes i and j in the covariance matrix S
  # card.tau : size of chain component i and j belong to
  # N        : number of historical data (number of rows of data matrix Z)
  # S        : (q,q) matrix t(Z)%*%Z
  # a_omega  : hyperparameter of the Wishart prior
  # g        : variance of the interventional density f.tilde
  # n        : sample size
  
  # OUTPUT:
  
  # posterior predictive of BF
  
  i = ij[1]; j = ij[2]
  
  arg.r.ij.2.sample = rbeta(n = H, shape1 = (n - 1)/2, shape2 = 1/2)
  
  if(log.BF == FALSE){
    
    BF.sample = pi^(-0.5)*gamma(n/2)/gamma((n + 1)/2)*(n)*(arg.r.ij.2.sample^((n - 1)/2))
    
  }
  
  if(log.BF == TRUE){
    
    BF.sample = -0.5*log(pi) + lgamma(n/2) - lgamma((n+1)/2) + log(n+1) + ((n-1)/2)*log(arg.r.ij.2.sample)
    
  }
  
  return(BF.sample)
  
}
