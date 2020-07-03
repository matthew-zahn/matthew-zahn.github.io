# Morley Nelson Zivot (2003) uc_ur programs.R

# Replication of the GAUSS program "uc_ur.opt" by Morley, Nelson, and Zivot
# Used in Why Are Beveridge-Nelson and Unobserved-Component Decompositions of GDP So Different?
# R code written by Matthew Zahn
# Please cite the original paper if you use these programs.

# These programs were used to estimate Table 4
# Notes: The code is generalized so that it can be applied to any data frame the user wants. 
#        Scalers T and START are defined globally. If the functions are nested within other functions, keep these items in mind, as they may need to be moved and refined locally.

# T <- nrow(df)

library(numDeriv)

START <- 2        # Start up values for VEVD of likelihood  

# Extracts results and formats in a list
results <- function(model, prmtr){
  
  # Final parameter estimates
  prm_fnl = t(trans(model$par))
  
  # Use Hessian to get standard errors for estimates 
  hessn0 = model$hessian
  cov0 = solve(hessn0)
  grdn_fnl = jacobian(trans, model$par)
  cov = grdn_fnl%*%cov0%*%t(grdn_fnl)
  sd = sqrt(diag(cov))
  
  LL = (-1)*model$value
  
  f = matrix(c(prm_fnl[5], 1, prm_fnl[6], 0), 2, 2)
  eig = eigen(f)
  
  f_out = cbind(prm_fnl, sd)
  rownames(f_out) = c("SN", "SE", "SNE", "mu", "phi1", "phi2")
  colnames(f_out) = c("Parameter Estimates", "Standard Deviation")
  
  f_out1 = rbind(LL)
  colnames(f_out1) = c("Log likelihood value")
  
  f_out2 = cbind(trans(prmtr))
  colnames(f_out2) = c("Initial parameter values")
  rownames(f_out2) = c("SN", "SE", "SNE", "mu", "phi1", "phi2")
  
  f_out3 = cbind(t(model$par))
  colnames(f_out3) = c("Pre-transformed estimates")
  rownames(f_out3) = c("SN", "SE", "SNE", "mu", "phi1", "phi2")
  
  f_out4 = cbind(eig$values)
  colnames(f_out4) = c("Eigen values")
  
  f_out5 = cbind(eig$vectors)
  colnames(f_out5) = c("Eigen vectors", "Eigen vectors")
  
  results = list(results = f_out, loglik = f_out1, intpar = f_out2, ptpar = f_out3,
                 prm = prm_fnl, eigval = f_out4, eigvec = f_out5, sd = sd, ll = LL, convg = model$convergence)
  
  return(results)
  
}

# Parameter constraints
trans <- function(c0){ # Constraining values of reg. coeff.
  
  # Constraints imposed
  # 1. cov(n_t, e_t) is pd
  # 2. AR(2) is stationary
  
  c1 = c0
  
  c1[1:2] = exp(-c0[1:2])
  
  # Constrain cov(n_t, e_t) to be pd using Cholesky factorization
  p       = matrix(0,2,2)
  p[1,1]  = exp(-c0[1])
  p[2,1]  = c0[3]
  p[2,2]  = exp(-c0[2])
  varcov  = t(p)%*%p
  c1[1]   = sqrt(varcov[1,1])
  c1[2]   = sqrt(varcov[2,2])
  c1[3]   = varcov[2,1]
  
  aaa = c0[5] / (1+abs(c0[5]))
  ccc = (1-abs(aaa))*c0[6]/(1+abs(c0[6]))+abs(aaa)-aaa^2
  
  c1[5] = 2*aaa
  c1[6] = (-1)*(aaa^2+ccc)
  
  return(c1)
  
}

# Define likelihood function
lik_fcn <- function(par, df){
  
  # Transform hyperparameters to impose constraints
  par = trans(par)
  
  SN   = par[1] # s.e. of the random walk component
  SE   = par[2] # s.e. of the AR component
  SNE  = par[3] # cov(n,e)
  mu   = par[4]
  phi1 = par[5]
  phi2 = par[6]
  
  muvec = matrix(c(mu,0,0),3,1)               # Drift vector
  F = matrix(c(1,0,0,0,phi1,1,0,phi2,0),3,3)  # Transition matrix
  F1 = F[-1,-1]                               # Pick out I(0) part
  H = matrix(c(1,1,0),1,3)                    # Measurement equation
  
  Q = matrix(0,3,3)                           # E[V(t)t(V(t))]
  Q[1,1] = SN^2
  Q[2,2] = SE^2
  Q[1,2] = SNE
  Q[2,1] = SNE
  
  Q1 = Q[-1,-1] # Cov matrix of I(0) part
  Q1vec = Q1
  dim(Q1vec) = c(dim(Q1vec)[1]*dim(Q1vec)[2],1)
  R = 0         # Cov matrix of error in measurement equation
  
  beta_ll = matrix(c(df[1],0,0),3,1) # Starting values 
  
  # Variance matrix of initial state vector
  p_ll          = matrix(0,3,3)
  p_ll1         = solve((diag(4)-kronecker(F1,F1)))%*%Q1vec
  p_ll1         = matrix(p_ll1,2,2)
  p_ll[1,1]     = 100000000
  p_ll[2:3,2:3] = p_ll1
  
  lik_mat = matrix(0,T,1)
  for(j_iter in 1:T){
    
    beta_tl = muvec+F%*%beta_ll
    p_tl    = F%*%p_ll%*%t(F)+Q
    
    vt = df[j_iter]-H%*%beta_tl # Prediction error
    ft = H%*%p_tl%*%t(H)+R      # Variance of forecast error
    
    beta_tt = beta_tl+p_tl%*%t(H)%*%solve(ft)%*%vt
    p_tt = p_tl - p_tl%*%t(H)%*%solve(ft)%*%H%*%p_tl
    
    lik_mat[j_iter,1] = 0.5*(log((2*pi*ft))+vt^2/ft)
    
    beta_ll = beta_tt
    p_ll = p_tt
    
  }
  
  val = apply(as.matrix(lik_mat[START:T]),2,sum)
  
  return(val)
  
} 

# Define filter
filter <- function(par, df){
  
  beta_mat = matrix(0,T,2)
  
  lik_mat = matrix(0,T,1)
  
  par = trans(par)
  
  SN   = par[1] # SE of the random walk component
  SE   = par[2] # SE of the AR component
  SNE  = par[3] # Cov(n,e)
  mu   = par[4]
  phi1 = par[5]
  phi2 = par[6]
  
  muvec = matrix(c(mu,0,0),3,1)               # Drift vector
  F = matrix(c(1,0,0,0,phi1,1,0,phi2,0),3,3)  # Transition matrix
  F1 = F[-1,-1]                               # Pick out I(0) part
  H = matrix(c(1,1,0),1,3)                    # Measurement equation
  
  Q = matrix(0,3,3)                           # E[V(t)t(V(t))]
  Q[1,1] = SN^2
  Q[2,2] = SE^2
  Q[1,2] = SNE
  Q[2,1] = SNE
  
  Q1 = Q[-1,-1] # Cov matrix of I(0) part
  Q1vec = Q1
  dim(Q1vec) = c(dim(Q1vec)[1]*dim(Q1vec)[2],1)
  R = 0         # Cov matrix of error in measurement equation
  
  beta_ll = matrix(c(df[1],0,0),3,1) # Starting values 
  
  # Variance matrix of initial state vector
  p_ll          = matrix(0,3,3)
  p_ll1         = solve((diag(4)-kronecker(F1,F1)))%*%Q1vec
  p_ll1         = matrix(p_ll1,2,2)
  p_ll[1,1]     = 100000000
  p_ll[2:3,2:3] = p_ll1
  
  for(j_iter in 1:T){
    
    beta_tl = muvec+F%*%beta_ll
    p_tl = F%*%p_ll%*%t(F)+Q
    
    vt = df[j_iter,1]-H%*%beta_tl # Prediction error
    ft = H%*%p_tl%*%t(H)+R        # variance of forecast error
    
    beta_tt = beta_tl+p_tl%*%t(H)%*%solve(ft)%*%vt
    p_tt    = p_tl - p_tl%*%t(H)%*%solve(ft)%*%H%*%p_tl
    
    beta_mat[j_iter,] = matrix(c(beta_tt[1,1], beta_tt[2,1]),1,2)
    beta_ll = beta_tt
    p_ll = p_tt
    
  }
  
  return(beta_mat[START:T,])
  
}


# EOF
