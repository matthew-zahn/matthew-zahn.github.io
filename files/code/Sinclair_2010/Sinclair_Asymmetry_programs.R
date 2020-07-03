# Sinclair (2010) Asymmetry programs.R

# Replication of GAUSS programs "univariate asymmetric correlated with filter.opt" by Tara M. Sinclair
# Used in Asymmetry in the Business Cycle: Friedman's Plucking Model with Correlation Innovations Studies in Nonlinear Dynamics and Econometrics,  Vol 14: Iss. 1, Article 3 (2010)
# R code written by Matthew Zahn
# Please site original paper if you use these programs.

# This code was used to estimate Table 1 column 1 (and provides the estimated components and probabilities used in the figures).  
# This program is straightforward to modify to produce the other results in the paper except for the endogenous switching results which are available upon request.

# Notes: The code is generalized so that it can be applied to any data frame the user wants. 
#       Scalers T and START are defined globally. If the functions are nested within other functions, keep these items in mind, as they may need to be moved and refined locally.

# T <- nrow(df)

library(numDeriv)

results <- function(model, p_int){
  
  # Parameters from estimation
  prm_fnl = t(trans(model$par))
  
  # Use Hessian from MLE estimation to get the standard errors of estimates
  hessn0 <- model$hessian
  cov0 <- solve(hessn0)
  grdn_fnl <- jacobian(trans, model$par)
  cov <- grdn_fnl%*%cov0%*%t(grdn_fnl)
  sd <- sqrt(diag(cov))
  
  LL <- (-1)*model$value
  
  f_out <- cbind(prm_fnl, sd)
  rownames(f_out) <- c("SN", "SE", "mu", "phi1", "phi2", "gamma", "ppr", "qpr", "SNE")
  colnames(f_out) <- c("Parameter Estimates", "Standard Deviation")
  
  f_out1 <- rbind(LL)
  colnames(f_out1) <- c("Log likelihood value")
  
  f_out2 <- cbind(trans(p_int))
  colnames(f_out2) <- c("Initial parameter values")
  rownames(f_out2) <- c("SN", "SE", "mu", "phi1", "phi2", "gamma", "ppr", "qpr", "SNE")
  
  f_out3 <- cbind(t(model$par))
  colnames(f_out3) <- c("Pre-transformed estimates")
  rownames(f_out3) <- c("SN", "SE", "mu", "phi1", "phi2", "gamma", "ppr", "qpr", "SNE")
  
  results <- list(results = f_out, loglik = f_out1, intpar = f_out2, ptpar = f_out3,
                  prm = prm_fnl, sd = sd, ll = LL, convg = model$convergence)
  
  return(results)
  
}

lik_fcn <- function(par, df){
  
  # Transform hyperparamters to impose contraints
  par <- trans(par)
  
  SNY     = par[1]
  SEY     = par[2]
  mu      = par[3] # Drift 
  phi1y   = par[4]
  phi2y   = par[5]
  gammay  = par[6]
  ppr     = par[7] # Pr[St=1/St-1=1]
  qpr     = par[8] # Pr[St=0/St-1=1]
  
  # Correlation
  SNYEY   = par[9]*SNY*SEY
  
  # Drift vector for drift in mu and Markov Switching, St=0
  muvec0  = matrix(c(mu, 0, 0),3,1)
  # Drift vector for drift in mu and Markov Switching, St=1
  muvec1  = matrix(c(mu, gammay, 0),3,1)
  
  # Transition matrix
  F = matrix(c(1,0,0,0,phi1y,1,0,phi2y,0),3,3)
  
  # Measurement equation
  H = matrix(c(1,1,0),1,3)
  Q = matrix(0,3,3)
  # E[V(t)V(t)'] note no state-dependent variances in this example
  Q[1,1] = SNY^2
  Q[2,2] = SEY^2
  Q[1,2] = SNYEY
  Q[2,1] = SNYEY
  
  ### INITIAL PROB Pr[S0/Y0] ###
  prior_cf = matrix(c(df[1,1],0,0),3,1)
  prior_vr = diag(3)*1000
  prob_1 = (1-qpr) / (2-ppr-qpr)
  prob_0 = 1 - prob_1 # Pr[St-1=1/Yt-1], Stead State Prob
  
  # Initial values
  p_cf_0 = prior_cf
  p_cf_1 = prior_cf
  
  p_vr_0 = prior_vr
  p_vr_1 = prior_vr
  likv = 0.0
  
  ### PREDICTION ###
  for(j_iter in 1:T){
    pst_cf00 = muvec0 + F%*%p_cf_0
    pst_cf01 = muvec1 + F%*%p_cf_0
    pst_cf10 = muvec0 + F%*%p_cf_1
    pst_cf11 = muvec1 + F%*%p_cf_1
    
    pst_vr00 = F%*%p_vr_0%*%t(F)+Q
    pst_vr01 = F%*%p_vr_0%*%t(F)+Q
    pst_vr10 = F%*%p_vr_1%*%t(F)+Q
    pst_vr11 = F%*%p_vr_1%*%t(F)+Q
    
    f_cast01 = t(df[j_iter,])-H%*%pst_cf01
    f_cast11 = t(df[j_iter,])-H%*%pst_cf11
    f_cast10 = t(df[j_iter,])-H%*%pst_cf10
    f_cast00 = t(df[j_iter,])-H%*%pst_cf00
    
    ss00 = H%*%pst_vr00%*%t(H)
    ss01 = H%*%pst_vr01%*%t(H)
    ss10 = H%*%pst_vr10%*%t(H)
    ss11 = H%*%pst_vr11%*%t(H)
    
    pr_vl00 = v_prob(f_cast00, ss00)%*%qpr%*%prob_0
    pr_vl01 = v_prob(f_cast01, ss01)%*%(1-qpr)%*%prob_0
    pr_vl10 = v_prob(f_cast10, ss10)%*%(1-ppr)%*%prob_1
    pr_vl11 = v_prob(f_cast11, ss11)%*%ppr%*%prob_1
    
    pr_val = pr_vl00+pr_vl01+pr_vl10+pr_vl11
    
    lik =(-1)*log(pr_val)
    
    ### UPDATING ###
    k_gn00 = pst_vr00%*%t(H)%*%solve(ss00)
    k_gn01 = pst_vr01%*%t(H)%*%solve(ss01)
    k_gn10 = pst_vr10%*%t(H)%*%solve(ss10)
    k_gn11 = pst_vr11%*%t(H)%*%solve(ss11)
    
    p_cf00 = pst_cf00+k_gn00%*%f_cast00
    p_cf01 = pst_cf01+k_gn01%*%f_cast01
    p_cf10 = pst_cf10+k_gn10%*%f_cast10
    p_cf11 = pst_cf11+k_gn11%*%f_cast11
    
    p_vr00 = pst_vr00-k_gn00%*%H%*%pst_vr00
    p_vr01 = pst_vr01-k_gn01%*%H%*%pst_vr01
    p_vr10 = pst_vr10-k_gn10%*%H%*%pst_vr10
    p_vr11 = pst_vr11-k_gn11%*%H%*%pst_vr11
    
    pro_00 = pr_vl00*solve(pr_val)
    pro_01 = pr_vl01*solve(pr_val)
    pro_10 = pr_vl10*solve(pr_val)
    pro_11 = pr_vl11*solve(pr_val)
    
    pro_00 = as.numeric(pro_00)
    pro_01 = as.numeric(pro_01)
    pro_10 = as.numeric(pro_10)
    pro_11 = as.numeric(pro_11)
    
    prob_0 = pro_00+pro_10
    prob_1 = pro_01+pro_11
    
    p_cf_0 = (pro_00*p_cf00+pro_10*p_cf10)%*%solve(prob_0)
    p_cf_1 = (pro_01*p_cf01+pro_11*p_cf11)%*%solve(prob_1)
    
    p_vr_0 = (pro_00*(p_vr00+(p_cf_0-p_cf00)%*%t((p_cf_0-p_cf00))) +
                pro_10*(p_vr10+(p_cf_0-p_cf10)%*%t((p_cf_0-p_cf10)))) / prob_0
    
    p_vr_1 = (pro_01*(p_vr01+(p_cf_1-p_cf01)%*%t((p_cf_1-p_cf01))) +
                pro_11*(p_vr11+(p_cf_1-p_cf11)%*%t((p_cf_1-p_cf11)))) / prob_1
    
    likv = likv + lik 
    
  }
  
  return(likv)
  
}

lik_out <- function(par, df){
  
  # Transform hyperparamters to impose contraints
  par <- trans(par)
  
  SNY     = par[1]
  SEY     = par[2]
  mu      = par[3] # Drift 
  phi1y   = par[4]
  phi2y   = par[5]
  gammay  = par[6]
  ppr     = par[7] # Pr[St=1/St-1=1]
  qpr     = par[8] # Pr[St=0/St-1=1]
  
  # Correlation
  SNYEY   = par[9]*SNY*SEY
  
  # Drift vector for drift in mu and Markov Switching, St=0
  muvec0 = matrix(c(mu,0,0),3,1)
  
  # Drift vector for drift in mu and MArkov Switching, St=1
  muvec1 = matrix(c(mu,gammay,0),3,1)
  
  # Note that in this version both u and y have same value for St
  
  # Transition matrix
  F = matrix(c(1,0,0,0,phi1y,1,0,phi2y,0),3,3)
  
  # Measurement equation
  H = matrix(c(1,1,0),1,3)
  
  # E[V(t)V(t)'] note no state-dependent variances in this example
  Q = matrix(0,3,3)
  Q[1,1] = SNY^2
  Q[2,2] = SEY^2
  Q[1,2] = SNYEY
  Q[2,1] = SNYEY
  
  ### Initial prob Pr[S0/Y0] ###
  prior_cf = matrix(c(df[1,1],0,0),3,1)
  prior_vr = diag(3)*1000
  prob_1 = (1-qpr) / (2-ppr-qpr)
  # Pr[st-1=1/Yt-1], Steady State Prob.
  prob_0 = 1-prob_1 # Pr[St-1=1/Yt-1], Steady State Prob
  
  # Initial values
  p_cf_0 = prior_cf
  p_cf_1 = prior_cf
  
  p_vr_0 = prior_vr
  p_vr_1 = prior_vr
  
  out_mat = matrix(0,T,5)
  
  for(j_iter in 1:T){
    
    pr_0_l = qpr*prob_0+(1-ppr)*prob_1 # PR[St=0/Yt-1]
    pr_1_l = (1-qpr)*prob_0+ppr*prob_1 # PR[St=1/Yt-1]
    
    out_mat[j_iter,5] = pr_0_l
    
    ### PREDICTION ###
    pst_cf00 = muvec0 + F%*%p_cf_0
    pst_cf01 = muvec1 + F%*%p_cf_0
    pst_cf10 = muvec0 + F%*%p_cf_1
    pst_cf11 = muvec1 + F%*%p_cf_1
    
    pst_vr00 = F%*%p_vr_0%*%t(F)+Q
    pst_vr01 = F%*%p_vr_0%*%t(F)+Q
    pst_vr10 = F%*%p_vr_1%*%t(F)+Q
    pst_vr11 = F%*%p_vr_1%*%t(F)+Q
    
    f_cast01 = t(df[j_iter,])-H%*%pst_cf01
    f_cast11 = t(df[j_iter,])-H%*%pst_cf11
    f_cast10 = t(df[j_iter,])-H%*%pst_cf10
    f_cast00 = t(df[j_iter,])-H%*%pst_cf00
    
    ss00 = H%*%pst_vr00%*%t(H)
    ss01 = H%*%pst_vr01%*%t(H)
    ss10 = H%*%pst_vr10%*%t(H)
    ss11 = H%*%pst_vr11%*%t(H)
    
    pr_vl00 = v_prob(f_cast00, ss00)%*%qpr%*%prob_0
    pr_vl01 = v_prob(f_cast01, ss01)%*%(1-qpr)%*%prob_0
    pr_vl10 = v_prob(f_cast10, ss10)%*%(1-ppr)%*%prob_1
    pr_vl11 = v_prob(f_cast11, ss11)%*%ppr%*%prob_1
    
    pr_val = pr_vl00+pr_vl01+pr_vl10+pr_vl11
    
    ### UPDATING ###
    k_gn00 = pst_vr00%*%t(H)%*%solve(ss00)
    k_gn01 = pst_vr01%*%t(H)%*%solve(ss01)
    k_gn10 = pst_vr10%*%t(H)%*%solve(ss10)
    k_gn11 = pst_vr11%*%t(H)%*%solve(ss11)
    
    p_cf00 = pst_cf00+k_gn00%*%f_cast00
    p_cf01 = pst_cf01+k_gn01%*%f_cast01
    p_cf10 = pst_cf10+k_gn10%*%f_cast10
    p_cf11 = pst_cf11+k_gn11%*%f_cast11
    
    p_vr00 = pst_vr00-k_gn00%*%H%*%pst_vr00
    p_vr01 = pst_vr01-k_gn01%*%H%*%pst_vr01
    p_vr10 = pst_vr10-k_gn10%*%H%*%pst_vr10
    p_vr11 = pst_vr11-k_gn11%*%H%*%pst_vr11
    
    pro_00 = pr_vl00*solve(pr_val)
    pro_01 = pr_vl01*solve(pr_val)
    pro_10 = pr_vl10*solve(pr_val)
    pro_11 = pr_vl11*solve(pr_val)
    
    pro_00 = as.numeric(pro_00)
    pro_01 = as.numeric(pro_01)
    pro_10 = as.numeric(pro_10)
    pro_11 = as.numeric(pro_11)
    
    prob_0 = pro_00+pro_10
    prob_1 = pro_01+pro_11
    
    out_mat[j_iter,1] = prob_1
    out_mat[j_iter,5] = prob_0
    
    p_cf_0 = (pro_00*p_cf00+pro_10*p_cf10)%*%solve(prob_0)
    p_cf_1 = (pro_01*p_cf01+pro_11*p_cf11)%*%solve(prob_1)
    
    p_cftt = p_cf_0*prob_0+p_cf_1*prob_1
    
    #out_mat[j_iter,2:3] = t(p_cftt[1:2,1])
    
    out_mat[j_iter,2] = (matrix(c(1,0,0),1,3))%*%p_cftt
    out_mat[j_iter,3] = (matrix(c(0,1,0),1,3))%*%p_cftt
    
    p_vr_0 = (pro_00*(p_vr00+(p_cf_0-p_cf00)%*%t((p_cf_0-p_cf00))) +
                pro_10*(p_vr10+(p_cf_0-p_cf10)%*%t((p_cf_0-p_cf10)))) / prob_0
    
    p_vr_1 = (pro_01*(p_vr01+(p_cf_1-p_cf01)%*%t((p_cf_1-p_cf01))) +
                pro_11*(p_vr11+(p_cf_1-p_cf11)%*%t((p_cf_1-p_cf11)))) / prob_1
    
    #out_mat[j_iter,5] = p_vr_0[4,4]%*%prob_0+p_vr_1[4,4]%*%prob_1
    
  }
  
  return(out_mat[1:T,])
  
}

trans <- function(c0){ # Constraining values of reg. coeff
  
  # Constraints imposed:
  # 1. Cov matrix is pd
  # 2. AR(2) is stationary for both series
  
  c1 = c0
  
  # Constrain cov(n_t,e_t) to be pd using Cholesky factorization
  p       = matrix(0,2,2)
  p[1,1]  = exp(-c0[1])
  p[2,1]  = c0[9]
  p[2,2]  = exp(-c0[2])
  varcov  = t(p)%*%p
  c1[1]   = sqrt(varcov[1,1])
  c1[2]   = sqrt(varcov[2,2])
  
  # Correlation
  c1[9]   = varcov[2,1] / (c1[1]*c1[2])
  
  # Constraining the AR parameters
  # MNZ Version
  aaay    = c0[4]/(1+abs(c0[4]))
  cccy    = (1-abs(aaay))*c0[5]/(1+abs(c0[5]))+abs(aaay)-aaay^2
  c1[4]   = 2*aaay
  c1[5]   = (-1)*(aaay^2+cccy)
  # Constrain gamma to be negative 
  c1[6]   = -exp(-c0[6])
  # Constraining the probabilities
  c1[7:8] = exp(-c0[7:8])/(1+exp(-1*c0[7:8]))
  
  return(c1)
  
}

v_prob <- function(ev, he){
  
  # Calculates Pr[Yt/St, Yt-1]
  val = (1/sqrt(2*pi*he))*exp(-0.5*ev*ev/he)
  
  return(val)
  
}

smooth <- function(results, input_pr_tt0, input_pr_tl0){
  
  # pr_tt0 contains Pr[S_t|Y_t]
  # pr_tl0 contains Pr[S_t|Y_t-1]
  
  pr_tt0 = as.matrix(input_pr_tt0)
  pr_tl0 = as.matrix(input_pr_tl0)
  
  ppr = results$prm[7]
  qpr = results$prm[8]
  
  pr_tt1 = 1-pr_tt0
  pr_tl1 = 1-pr_tl0
  
  pr_sm0 = pr_tt0 #pr_sm0 will contain Pr[S_t|Y_t]
  pr_sm1 = pr_tt1
  
  end = T-1
  
  for(j_iter in end:0){
    
    j = j_iter + 1 
    
    # The following are P[S_t, S_t+1|Y_t]
    pr_sm00 = pr_sm0[j,1]%*%qpr%*%pr_tt0[j_iter,1]/pr_tl0[j,1]
    pr_sm01 = pr_sm1[j,1]%*%(1-qpr)%*%pr_tt0[j_iter,1]/pr_tl1[j,1]
    pr_sm10 = pr_sm0[j,1]%*%(1-ppr)%*%pr_tt1[j_iter,1]/pr_tl0[j,1]
    pr_sm11 = pr_sm1[j,1]%*%ppr%*%pr_tt1[j_iter,1]/pr_tl1[j,1]
    
    pr_sm0[j_iter,1] = pr_sm00+pr_sm01
    pr_sm1[j_iter,1] = pr_sm10+pr_sm11
    
  }
  
  #results = list(sm0 = pr_sm0, sm1 = pr_sm1) 
  
  return(pr_sm0)
  
}

# EOF