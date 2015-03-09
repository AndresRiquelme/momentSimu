# The following packages are required by some of this functions:
library("RcppArmadillo") # Matrix algebra functions in C [used by gmm]
# Generate Covariance Matrix -----------------------------------------------------------------
# gensigma: generates the covariance matrix for each setup
# ARGUMENTS:
#  q:      Total Number of Moments
#  s:      Known valid moments
#  si:     Invalid Moments
#  d:      Correlation for invalid moments (non-local)
#  czz:    Covariance between instruments (constant)
#  covuv:  Covariance between u,v
#  local:  Indicator for fixed or local-to-zero moments
# RETURN:  (q+2) covariance matrix
gensigma <- function(q, s, si, d, czz, covuv, local=FALSE){
  qs <- q-s                               # Moments to be tested
  sv <- qs-si                             # Valid moments (half of them)
  # The next 10 lines generates the sigma matrix of the structural model
  if(local) {                             # In the local-to-zero case:
    covzu <- c(rep(0, s+sv),              # cov of the s+sv valid
               rep(1/n, ceiling(si/3)),              # cov of the invalid at 1/n rate
               rep(1/sqrt(n), ceiling(si/3)),        # cov of the invalid at 1/sqrt(n) rate
               rep(1/n^{1/3}, si-2*ceiling(si/3)) )  # cov of the invalid at 1/n^(1/3) rate
  } else {
    covzu <- c(rep(0, s+sv), rep(d, si) ) # cov of the s+sv valid and si invalid set as d
  }
  covzz <- czz * diag(q)
  sigma <- cbind(rbind(covzz, covzu, 0), rbind(cbind(covzu,0) ,covuv))
  return(sigma)
}
# Generate Data -----------------------------------------------------------------
# gendata: generates the structured first step data
# ARGUMENTS:
#  n:         Number of Observations
#  q:         Total number of moments
#  s:         known valid moments
#  pcoeff:    Coeff of valid moments
#  strong:    Coeff of the strong moments
#  cholsigma: Cholesky decomposition of the covariance matrix
#  theta:     True theta parameter
#  het:       Heteroskedastic errors
# RETURN:     list of the response, covariates and instruments matrices
gendata <- function(n, q, s, pcoeff, strong, cholsigma, theta, het=FALSE){
  qs <- as.integer(q-s)                   # Moments to be tested
  sv <- as.integer(qs-si)                 # Valid moments (half of them)
  M  <- matrix(rnorm(n*{q+2}),n, q+2) # 
  M  <- M%*%cholsigma
  Z  <- M[, seq_len(q)]
  u  <- M[,{q+1}]
  v  <- M[,{q+2}]
  # Generate response variables
  pi <- c(rep(strong,s),rep(pcoeff,qs))
  if(het) u <- u * drop(apply(Z, 1, function(x) sqrt(sum(x^2)))) # Heteroskedastic errors
  Y2 <- as.vector(Z%*%pi + v)             # Reduced Form equation
  Y1 <- as.vector(Y2*theta + u)           # Structural equation
  # THE DATA HAS BEEN GENERATED!
  result <- list(Y1 = Y1, Y2=Y2, Z = Z)
}
# Indexes -----------------------------------------------------------------
# indexes: generate indexes of variables for the index-th combination of qs vars
# ARGUMENTS:
#  index:  number of the combination (from 1 to 2^qs)
#  qs:     total number of variables to be combined
#  known:  number of always included indexes (append at the beginning of the returned vector)
# RETURN: vector of indexes of length between 0 and 2^qs+known
indexes <- function(index, qs, known=0){
  vec   <- as.integer(NULL)
  l1    <- as.integer(1)
  l2    <- as.integer(2^{qs-1})
  for(i in seq_len(qs)){
    ifelse(index >= l1 & index <= l2, vec <- c(vec, as.integer(i)), l1 <- as.integer(l2+1))
    l2  <- as.integer(l1+2^{qs-i-1}-1)
  }
  return(c(seq_len(known), as.integer(known)+vec))
}
# Reverse Indexes -----------------------------------------------------------------
# rev.indexes: takes the vector of indexes and returns the index in the combinations list
# ARGUMENTS:
#  index:  vector of indexes (one of the 2^qs combinations)
#  qs:     total number of variables to be combined (the "unknown")
#  known:  number of always included indexes (append at the beginning of the returned vector)
# RETURN: scalar of the index of the vector in the combinations list
rev.indexes <- function(index, qs, known=0){
  index2  <- index[{known+1}:length(index)] - as.integer(known)
  l1      <- as.integer(1)         # lower limit
  l2      <- as.integer(2^{qs-1})  # Upper limit
  for(j in seq_len(qs)){
    if (!{j %in% index2}) l1 <- as.integer(l2+1)
    l2  <- as.integer(l1+2^{qs-j-1}-1)
  }
  return(l1)
}
# GMM ---------------------------------------------------------------------
# gmm: estimates the theta parameter of a linear gmm model using an efficient second step 
#      weighting matrix
# it makes use of the existing yy1, yy2 and n
# ARGUMENTS:
#  yy1: response
#  yy2: endogenous variable
#  zz:  Matrix of instruments
# RETURN: theta and j in a list
gmm <- function(zz){               # y1, y2, zz, n=nrow(zz)){
 # TWO STEP ESTIMATOR
 # First step
  pzz       <- zz%*%chol2inv(chol(t(zz)%*%zz))%*%t(zz)     # Projection matrix: z(z'z)^(-1)z'
  y2.tilde  <- pzz%*%y2                                    # Use the existing y22
  theta.fs  <- fastLmPure(y2.tilde,y1)$coefficients        # first step (fs) theta ("RcppArmadillo") 2SLS
  e.fs      <- y1 - theta.fs * y2                          # fs residual
  ze        <- sweep(zz, 1, e.fs, FUN="*" )                # Generates zi*ei (nxq)
  w         <- chol2inv(chol({1/n}*t(ze)%*%ze))            # Optimal Weighting matrix for second step
  w.chol    <- chol(w)                                     # Cholesky decomposition for data transformation
#Second step
  z.tilde   <- w.chol%*%t(zz)%*%y2
  y.tilde   <- as.vector(w.chol%*%t(zz)%*%y1)
  out.ss    <- fastLmPure(z.tilde,y.tilde)
  theta.ss  <- out.ss$coefficients    # second step(ss) theta
  e.ss      <- y1 - theta.ss * y2                          # ss error
  gbar      <- {t(zz)%*%e.ss}/n                            # Stacked Moments
  j2S       <- n*t(gbar)%*%w%*%gbar                        # j function Value
# CUE
  cue       <- optimize(gmm_CUE, c(-.5,1.5), instr = zz)   # Optimization using CUE
#Estimation of W based on cue optimum (to keep gmm_cue return scalar)
  e.cue     <- y1 - cue$minimum * y2                       # cue residual
  ze.cue    <- sweep(zz, 1, e.cue, FUN="*" )               # Generates zi*ei (nxq)
  w.cue     <- chol2inv(chol({1/n}*t(ze.cue)%*%ze.cue))    # Optimal Weighting matrix for second step
#Coverage 95% C:I:
  std.gmm   <- out.ss$stderr
  std.cue   <- sqrt(drop(({1/n}*t(y2)%*%zz%*%w.cue%*%t(zz)%*%y2)^{-1})) # Use inverse for scalar (single parameter)
  degf      <- n-ncol(zz)                                  # Dregrees of freedom
  moe.gmm   <- qt(1-alpha/2, degf)*std.gmm                 # Margin of error gmm
  moe.cue   <- qt(1-alpha/2, degf)*std.cue                 # Margin of error cue
  cover.gmm <- ifelse(theta > theta.ss    - moe.gmm & theta < theta.ss    + moe.gmm, TRUE, FALSE)
  cover.cue <- ifelse(theta > cue$minimum - moe.cue & theta < cue$minimum + moe.cue, TRUE, FALSE)
  result    <- c(theta.ss, j2S, cue$minimum, cue$objective, cover.gmm, cover.cue)
}
# CUE ---------------------------------------------------------------------
# gmm_CUE: estimates the theta parameter of a linear gmm model using a continuously updated estimator.
# This function is optimized by the gmm function.
# it makes use of the existing yy1, yy2 and n
# ARGUMENTS:
#  theta: Theta parameter
# RETURN: value of the j function at theta
gmm_CUE <- function(theta, instr){
  err    <- y1 - theta * y2
  gstack <- t(instr)%*%{err}/n
  ze     <- sweep(instr, 1, err, FUN="*" )  
  w      <- chol2inv(chol({1/n}*t(ze)%*%ze))                # Optimal Weighting matrix for second step
  jCUE   <- n*t(gstack)%*%w%*%gstack
  result <- jCUE
}
# Downward Testing --------------------------------------------------------
# dt: returns the index of the submodel selected by the Downward or Upward Testing Procedure
# ARGUMENTS: 
# j:      list of j-values
# k:      Number of instruments for each j value
# alpha:  significance level
# upward: option for the Upward Testing
# RETURN: index of the selected j value
dt <- function(j, k, alpha, upward=FALSE) {
  o         <- sort.list(k, decreasing = !upward)
  df2.vec   <- k - as.integer(1)                       # df = min(length(theta),ncol(Z))
  pval.vec  <- pchisq(j, df=df2.vec, lower.tail=FALSE) # Vector of p-values
  nmoment   <- k[o[which.max(pval.vec[o]<=alpha)]]     # Number of moments for the the first pval<=alpha (r_dt)
  sub       <- which(k==nmoment)                       # Elements with c = r_dt
  return(sub[which.min(j[sub])])
}
# Okui Shrinkage  -----------------------------------------------------------------
# ARGUMENT:
#  index:  number of the index in the 2^k list
# RETURN: scalar theta 2sls shrinkage
# Preliminary estimates
okui <- function(index, method="j2s", param=0){      # Args: index, y1, y2, theta .vec, zz, combinations.list
  if (method=="j2s")    theta.this <- theta.vec[index]
  if (method=="jcue")   theta.this <- thetacue.vec[index]
  if (param!=0)         theta.this <- param         # For adaptive lasso control
  nvar     <- length(combinations.list[[index]])
  err1      <- y1 - y2*theta.this 
  zzthis    <- zz[, combinations.list[[index]]]
  y2hat    <- fastLm(zzthis, y2)$fitted.values
  err2      <- y2-y2hat
  sigma2e   <- sum(err1^2)/{n-nvar}    
  sigma2eu  <- sum(err1*err2)/{n-nvar}
#print(c(sigma2e, sigma2eu))
  # Transformations
  zz_known  <- zz[, seq_len(s)]
  pxx       <- zz_known%*%chol2inv(chol(t(zz_known)%*%zz_known))%*%t(zz_known)
  if (nvar-s > 0) {   # If tested number of instruments > 0 after first step (method selects at least one)
  zz_test   <- zzthis[, {s+1}:nvar]
  zz_tr     <- {diag(n)-pxx} %*% zz_test          # Transformed Z: [I-Px]Z (Okui 2011, page 72)
  pzz       <- zz_tr%*%chol2inv(chol(t(zz_tr)%*%zz_tr))%*%t(zz_tr) 
  } else {pzz <- as.integer(0)*diag(n)}              # Only known valid are used
# s*
  sstar_num <- sigma2e/n * t(y2hat)%*%pzz%*%y2hat
  sstar_den1<- sigma2eu * ncol(pzz)/n
  sstar     <- as.numeric(sstar_num / {sstar_den1 + sstar_num})
  ps        <- pxx + sstar*pzz
  thetaokui <- solve(t(y2)%*%ps%*%y2)%*%t(y2)%*%ps%*%y1
#COVERAGE
  degf      <- n-nvar                                       # Degrees of freedom
  sigmae2.okui<- sum({y1 - y2*thetaokui}^2)/degf            # Error Variance 
  vc.okui   <- sigmae2.okui * drop(t(y2)%*%ps%*%y2)^{-1}    # Covariance matrix
  moe.okui  <- qt(1-alpha/2, degf)*sqrt(vc.okui)            # Margin of error 
  cover.okui<- ifelse(theta > thetaokui - moe.okui & theta < thetaokui + moe.okui, TRUE, FALSE)
#  print(c(moe.okui, thetaokui, cover.okui, thetaokui - moe.okui, thetaokui + moe.okui, moe.okui))
  return(list(theta=thetaokui, cover=cover.okui))
}
# Data Shrinkage --------------------------------------------------------
# shrinkdata: transform the data for linear gmm into data for shrinkage estimation
# ARGUMENTS:
#  yy1: response
#  yy2: endogenous variable
#  zz:  Matrix of instruments
#  n:   Number of observations
# RETURN: list with transformed response and covariates for shrinkage
shrinkdata <- function(y1, y2, zz, n){                       # use existing (y1, y2, zz, n)
  pzz       <- zz%*%chol2inv(chol(t(zz)%*%zz))%*%t(zz)       # Projection matrix: z(z'z)^(-1)z'
  y2.tilde  <- pzz%*%y2                                      # Use the existing y22
  theta.fs  <- fastLmPure(y2.tilde,y1)$coefficients          # first step (fs) theta ("RcppArmadillo") 2SLS
  e.fs      <- y1 - theta.fs * y2                            # fs error
  ze        <- sweep(zz, 1, e.fs, FUN="*" )                  # Generates zi*ei (nxq)
  w         <- chol2inv(chol({1/n}*t(ze)%*%ze))              # Optimal Weighting matrix for second step
  w.chol    <- chol(w)                                       # Cholesky decomposition for data transformation
  F.mat     <- rbind(matrix(as.integer(0), s, qs), diag(qs)) # Auxiliar F matrix
  y1z       <- w.chol%*%t(zz)%*%y1                           # Transform to linear OLS
  y2z       <- w.chol%*%cbind(t(zz)%*%y2 , n*F.mat)
  result <- list(y1z=y1z, y2z=y2z)
}
# Adaptive lasso --------------------------------------------------------
# alasso: performs the adaptive lasso estimation (without intercept)
# ARGUMENTS:
#  X:         covariates
#  Y:         response
#  epsilon:   tolerance level
#  max.steps: maximum number of iterations
# RETURN: list with multiple results
alasso <- function(
  X,                                 # Matrix of Covariates
  Y,                                 # Response
  epsilon     = .Machine$double.eps, # Tolerance level
  max.steps   = as.integer(2*min(dim(X))),      # Steps of Algorithm
  noshrink    = NULL
  )
{
  #  MAIN CODE HERE
  Y          <- as.matrix(Y)                            # In case Y is a vector
  n          <- as.integer(dim(Y)[1])
  m          <- as.integer(dim(X)[2])                   # Number of covariates in sample
  ones       <- rep(1, n)                               # Auxiliary Vector of ones (Constant)
  rootsumsq  <- rep(1, m)                               # Initial L2 norm
  rootsumsq  <- sqrt(ones %*% {X^2})                    # Centred Root of Sum of Squares
  X          <- scale(X, FALSE, rootsumsq)
  ###### BEGIN OF ADAPTIVE LASSO SCALING (Zou 2006) #####
  wj         <- as.matrix(rep(1,m))                     # Default weigths are 1 for LAR and LASSO
  if(abs(det(t(X)%*%X)) > epsilon) wj <- abs(solve(t(X)%*%X)%*%t(X)%*%Y)[1:ncol(X),] # ALASSO weights #
  X          <- scale(X, center=FALSE, scale=1/wj)      # ALASSO scale (weight is 1 for LASSO)
  ###### END OF ADAPTIVE LASSO SCALING ##################
  ### MODULE 3: PRE-ALLOCATIONS AND INITIALIZATIONS
  beta.std   <- matrix(0, max.steps , m)                # Matrix of Beta vectors
  c.mat      <- matrix(0, max.steps , m)                # Matrix of correlation vectors
  c.vec      <- drop(t(Y) %*% X)                        # Initial correlation Vector
  drops      <- NULL                                    # Dropping parameters indicator (for Lasso)
  s.vec      <- NULL                                    # Sign vector
  gamma.vec  <- NULL                                    # Vector of gammas (Equation 2.13)
  active     <- NULL                                    # Initial active set
  inactive   <- im <- seq(m)                            # Initial inactive set (necessary for getting right indexes)
  actions    <- as.list(seq(max.steps))                 # list of parameters in and out of the Active Set
  actions.names <- NULL                                 # Variable names in active set
  var.names  <- dimnames(X)[[2]]                        # Extracts var. names for output
  res        <- Y                                       # Initial residual vector
  ssy        <- ones %*% (Y^2)                          # Sum Y squared
  rss        <- ssy                                     # Initial Sum of Residuals
  r2         <- 1                                       # R squared
  lambda     <- double(1)                               # Lambda vector
  singularx  <- FALSE                                   # Flag for singularity error
  ### MODULE 4: MAIN ROUTINE
  k <- as.integer(0)                                                # MAIN LOOP STARTS HERE <---------
  while((k < max.steps) && (length(active) < min(m,n)) ) {
  action <- NULL                                      # Reset previous actions
  c.hat <- c.vec[inactive]                            # All inactive elements in C.vec
  k <- k + as.integer(1)                              # Updates the loop counter
  # Identify the largest nonactive gradient
  C.hat <- max(abs(c.hat))                            # Equation (2.8)
  if(C.hat<epsilon) {singularx <- TRUE; break}        # Prevents singularity error
  lambda[k]<- C.hat                                   # "Assigns" Lambda
  # Check if we are in a DROP situation (Valid for Lasso and ALASSO)
  if(!any(drops)) {
    new    <- abs(c.hat) >= C.hat                     # New index into the active set
    c.hat  <- c.hat[!new]                             # Remove the new index from C
    new    <- inactive[new]                           # Get right new index
    active <- c(active, new)                          # Updates active vector
    s.vec  <- c(s.vec, sign(c.vec[new]))              # Updates sign vector
    action <- c(action, new)                          # Update Actions sequence
  }
  else action <- - dropid                             # In Drop case
  Xa  <- X[, active]                                  # Active covariates
  Xa  <- t(s.vec*t(Xa))                               # Equation 2.4
  one <- matrix(abs(s.vec),length(s.vec),1)           # 1 (mx1) as defined in Equation 2.5
  Ga  <- t(Xa) %*% Xa                                 # Equation 2.5(a)
  if(abs(det(Ga)) < epsilon) {singularx <- TRUE; break} # Prevents singular Ga error
  A   <- drop(1/sqrt(t(one)%*%solve(Ga)%*% one))      # Equation 2.5(b)
  w   <- drop(A*solve(Ga)%*%one)                      # Equation 2.6(b)
  u   <- drop(Xa %*% w)                               # Equation 2.6(a)
  a   <- drop(t(X[,-active]) %*% u)                   # Equation 2.11
  gamma    <- c({C.hat - c.hat}/{A - a}, {C.hat + c.hat}/{A + a}) # Equation 2.13
  gamma.hat <- min(gamma[gamma > 0])         # Equation 2.13 (only +) 
  ###### BEGIN LASSO MODIFICATION 1 ####### (search for any drop) #########################
  dropid      <- NULL                                 # Number of dropped parameter (if any)
  dj          <- w * s.vec                            # Directional vector
  bj          <- beta.std[k, active]                  # Selects last active betas
  gamma.lasso <- -bj/dj                               # Equation 3.4
  gamma.tilde <- min(gamma.lasso[gamma.lasso > 0], gamma.hat)   # Equation 3.5
  if(gamma.tilde < gamma.hat) {                       # Theorem 1
    gamma.hat <- gamma.tilde                          # Equation 3.6(a)
    drops     <- gamma.lasso == gamma.tilde           # Equation 3.6(b)
  }  else {drops  <- NULL}                            # No drop
  ###### END LASSO MODIFICATION 1 ##########################################################    
  beta.std[k+1, ]       <- beta.std[k, ]              # Copy the last row 
  beta.std[k+1, active] <- beta.std[k+1, active] + gamma.hat * w * s.vec # Equation 3.3
  c.mat[k, ]            <- c.vec                      # Assign row of correlations  
  res                   <- res - gamma.hat * u        # Updates residuals
  c.vec     <- drop(t(X) %*% res)                     # Updates for next iteration
  gamma.vec <- c(gamma.vec, gamma.hat)                # Updates gamma sequence
  ###### BEGIN LASSO MODIFICATION 2 ####### (drop variable j from active set if necessary)
  if(any(drops)) {
    dropid <- active[drops]                           # Dropped variable
    beta.std[k + 1, dropid] <- 0                      # Dropped coef is zero
    active <- active[!drops]                          # Updated active set
    s.vec  <- s.vec[!drops]                           # Updated sign vector
  }
  ###### END LASSO MODIFICATION 2 ##########################################################
  if(!is.null(var.names)) actions.names <- c(actions.names , var.names[abs(action)]) # Use variable names in Actions
  actions[k]    <- action
  inactive <- im[ - active]
  }
  #   MAIN LOOP ENDS HERE <---------
  ### MODULE 5: OUTPUT PREPARATION
  beta.std  <- beta.std[seq(k + 1), ,drop=FALSE ]       # Cuts beta.std to the right dimension
  beta.orig <- scale(beta.std, FALSE, rootsumsq/wj)        # beta returns to its original scale
 # beta.orig <- scale(beta.orig,FALSE,1/wj)              # Adjusted ALASSO coefficients
  lambda    <- c(lambda[seq(k)], 0)                     # Cuts lambda vector 
  res       <- drop(Y) - X %*% t(beta.std)              # drop(Y) allows element-wise subtraction 
  rss       <- colSums(as.matrix(res^2,k,n))            # Computes the residuals in all the k steps
  r2        <- 1 - rss/rss[1]                           # Computes the R^2 for all steps
  actions   <- matrix(c(actions[1:k]),1,k)              # Actions as matrix
  netdf     <- sapply(actions,function(X)sign(X))       # Df as "k" at each stage (Efron's) eq 4.10
  df        <- cumsum(netdf)                            # Degrees of freedom
  df        <- c(Zero=0,df)                             # Adds df at "all-zero" stage
  rss.ols   <- rss[k + 1]                               # 
  df.ols    <- n - m                                    # OLS df (including constant)
  sigma2    <- rss.ols/df.ols                           # Sigma^2 OLS
  Cp        <- rss/sigma2 + 2 * df - n                  # Equation 4.5
  bic       <- rss/sigma2 + log(n)* df                  # BIC in Zou et al (2007)
  aic       <- rss/sigma2 + 2 * df                      # AIC
  ### MODULE 5: CHECk IF THERE ARE PARAMETERS NOT SHRUNK and deals with singularity error
  if (is.vector(noshrink)){
    temp.index <- matrix(as.integer(beta.orig!=0), k+1, m) * col(beta.orig) # Indexes by step
    tvec <- apply(temp.index,1, function(x) noshrink %in% x)          # tvec=TRUE when all choosen vars are in active set
    if (length(noshrink) > 1) tvec <- colSums(tvec)==length(noshrink)  
    bic <- {.Machine$double.xmax * !tvec}*2 + bic      # InfCrit is infinite if parameter is not in the active set
    aic <- {.Machine$double.xmax * !tvec}*2 + aic
    Cp  <- {.Machine$double.xmax * !tvec}*2 + Cp
  }
  if(singularx) beta.orig[1] <- NaN   # Put theta=NaN when singularity error occurs
  ### MODULE 6: RETURNED VALUES
  result    <- list(
  nobs        = n,                    # Observations
  betas.orig  = beta.orig,            # betas in original scale
  lambda      = lambda,               # lambda sequence
  sigma2      = sigma2,               # Full model variance
  cp          = Cp,                   # Mallows' Cp
  bic         = bic,                  # BIC in Zou et al (2007)
  aic         = aic                   # AIC
  )
}
# R squared --------------------------------------------------------
# r2struct: estimate the R2 for an structural regression using the variables in [index]
# ARGUMENTS:
#  index:  number of the index in the 2^k list
# RETURN: R2 of    yy2 = B zz + e
r2struct <- function(index){
 zz.this <- zz[, combinations.list[[index]]] 
 model3  <- fastLm(zz.this, y2)
 SSR     <- sum(model3$residuals^2)
 SST     <- sum(y2^2)
 r2      <- 1 - SSR/SST
 return(r2)
}
