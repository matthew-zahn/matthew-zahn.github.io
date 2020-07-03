# Morley Nelson Zivot (2003) arima212 programs.R

# Replication of the GAUSS program "arima212.opt" by Morley, Nelson, and Zivot
# Used in Why Are Beveridge-Nelson and Unobserved-Component Decompositions of GDP So Different?
# R code written by Matthew Zahn
# Please cite the original paper if you use these programs.

# These programs were used to estimate Table 2
# Notes: The code is generalized so that it can be applied to any data frame the user wants. 
# 		 Scalers T and START are defined globally. If the functions are nested within other functions, keep these items in mind, as they may need to be moved and refined locally.

# T <- nrow(df)

library(numDeriv)

START <- 1        # Start up values for VEVD of likelihood

# Extracts results and formats in a list
results <- function(model, prmtr){
  
  # Final parameter estimates
  prm_fnl = t(trans(model$par))
  
  # Use Hessian to get standard errors for estimates 
  hessn0   = model$hessian
  cov0     = solve(hessn0)
  grdn_fnl = jacobian(trans, model$par)
  cov      = grdn_fnl*cov0*t(grdn_fnl)
  sd       = sqrt(diag(cov))
  
  LL = (-1)*model$value
  
  f_out = cbind(prm_fnl, sd)
  rownames(f_out) = c("SE", "phi1", "phi2", "mu", "theta1", "theta2")
  colnames(f_out) = c("Parameter Estimates", "Standard Deviation")
  
  f_out1 = rbind(LL)
  colnames(f_out1) = c("Log likelihood value")
  
  f_out2 = cbind(trans(prmtr))
  colnames(f_out2) = c("Initial parameter values")
  rownames(f_out2) = c("SE", "phi1", "phi2", "mu", "theta1", "theta2")
  
  f_out3 = cbind(t(model$par))
  colnames(f_out3) = c("Pre-transformed estimates")
  rownames(f_out3) = c("SE", "phi1", "phi2", "mu", "theta1", "theta2")
  
  results = list(results = f_out, loglik = f_out1, intpar = f_out2, ptpar = f_out3,
                 prm = prm_fnl, sd = sd, ll = LL, convg = model$convergence)
  
  return(results)
  
}

# Parameter constraints
trans <- function(c0){
  
  c1 = c0
  
  c1[1] = exp(-c0[1])
  
  aaa = c0[2]/(1+abs(c0[2]))
  ccc = (1-abs(aaa))*c0[3]/(1+abs(c0[3]))+abs(aaa)-aaa^2
  
  c1[2] = 2*aaa
  c1[3] = (-1)*(aaa^2+ccc)
  
  return(c1)
  
}

# Define likelihood function
lik_fcn <- function(par, df){
  
  # Transform hyperparameters to impose constraints
  par = trans(par)
  
  SE     = par[1] # SE of the random walk component
  phi1   = par[2] # SE of the AR component
  phi2   = par[3] # Cov(n,e)
  mu     = par[4]
  theta1 = par[5]
  theta2 = par[6]
  
  F = matrix(c(phi1, 1, 0, 0,phi2,0,0,0,theta1,0,0,1,theta2,0,0,0),4,4) # Transition matrix
  H = matrix(c(1,0,0,0),1,4)                                            # Measurement equation
  Q = SE^2                                                              # E[V(t)t(V(t))]
  G = matrix(c(1,0,1,0),4,1)
  R = 0                                                                 # Cov matrix of error in measurement equation
  
  # Vector product of G*Q*G'
  GQtGvec = G%*%Q%*%t(G)
  dim(GQtGvec) = c(dim(GQtGvec)[1]*dim(GQtGvec)[2],1)
  
  beta_ll = matrix(0,4,1) # Starting values
  
  # Variance matrix of initial state vector
  p_ll = solve(diag(16)-kronecker(F,F))%*%GQtGvec
  p_ll = matrix(p_ll,4,4)
  
  lik_mat = matrix(0,T,1)
  
  for(j_iter in 1:T){
    
    beta_tl = F%*%beta_ll
    p_tl    = F%*%p_ll%*%t(F)+G%*%Q%*%t(G)
    
    vt = df[j_iter,1]-mu-H%*%beta_tl # Predicition error
    ft = H%*%p_tl%*%t(H)+R           # Variance of forecast errorr
    
    beta_tt = beta_tl+p_tl%*%t(H)%*%solve(ft)%*%vt
    p_tt    = p_tl-p_tl%*%t(H)%*%solve(ft)%*%H%*%p_tl
    
    lik_mat[j_iter,1] = 0.5*(log(2*pi*ft)+vt^2/ft)
    
    beta_ll = beta_tt
    p_ll    = p_tt
    
  }
  
  val = apply(as.matrix(lik_mat[START:T]),2,sum)
  
  return(val)
  
}

# Define filter
filter <- function(par, df){
  
  par = trans(par)
  
  beta_mat = matrix(0,T,1)
  lik_mat  = matrix(0,T,1)
  
  SE     = par[1] # SE of the random walk component
  phi1   = par[2] # SE of the AR component
  phi2   = par[3] # Cov(n,e)
  mu     = par[4]
  theta1 = par[5]
  theta2 = par[6]
  
  F = matrix(c(phi1, 1, 0, 0,phi2,0,0,0,theta1,0,0,1,theta2,0,0,0),4,4) # Transition matrix
  H = matrix(c(1,0,0,0),1,4)                                            # Measurement equation
  Q = SE^2                                                              # E[V(t)t(V(t))]
  G = matrix(c(1,0,1,0),4,1)
  R = 0                                                                 # Cov matrix of error in measurement equation
  
  # Vector product of G*Q*G'
  GQtGvec = G%*%Q%*%t(G)
  dim(GQtGvec) = c(dim(GQtGvec)[1]*dim(GQtGvec)[2],1)
  
  beta_ll = matrix(0,4,1) # Starting values
  
  # Variance matrix of initial state vector
  p_ll = solve(diag(16)-kronecker(F,F))%*%GQtGvec
  p_ll = matrix(p_ll,4,4)
  
  lik_mat = matrix(0,T,1)
  
  for(j_iter in 1:T){
    
    beta_tl = F%*%beta_ll
    p_tl    = F%*%p_ll%*%t(F)+G%*%Q%*%t(G)
    
    vt = df[j_iter,1]-mu-H%*%beta_tl # Predicition error
    ft = H%*%p_tl%*%t(H)+R           # Variance of forecast errorr
    
    beta_tt = beta_tl+p_tl%*%t(H)%*%solve(ft)%*%vt
    p_tt    = p_tl-p_tl%*%t(H)%*%solve(ft)%*%H%*%p_tl
    
    lik_mat[j_iter,1] = 0.5*(log(2*pi*ft)+vt^2/ft)
    
    beta_mat[j_iter,] = beta_tt[3,1]
    
    beta_ll = beta_tt
    p_ll    = p_tt
    
  }
  
  return(beta_mat[START:T,])
  
}

# EOF