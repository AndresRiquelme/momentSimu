try(setwd(paste("C:/Users/", Sys.info()[[7]], "/momentSimu", sep="")),silent = TRUE) # Sync all my computers on Skydrive
try(setwd(paste("C:/Users/", Sys.info()[[7]], "/momentSimu", sep="")),silent = TRUE) # Sync all my computers on Onedrive

#libraries
#My files
source("v2.R")
# USER SETINGS -----------------------------------------------------------------
# Monte Carlo Setings
reps    <- as.integer(1000)  # Number of repetitions for Monte Carlo simulation
start   <- as.integer(    1) # starting simulations (for replication with multiple computers)
seed    <- as.integer( 1116) # Initial seed.  Each repetition increases the seed in 1 unit (1115 works)
parallel<- TRUE              # Execute parallel simulations
# Data Setings
n_      <- as.integer(c(50, 100, 250)) # ***********Number of observations
strong  <- as.integer(    2) # Coeff of the strong moments
pcoeff_  <- c(.2,2)          # **********Coeff of valid moments [.2 2]
d       <-  .2               # Correlation for invalid moments (non-local)
theta   <-  .5               # True theta parameter
# Covariance Setings
czz_    <- c(0.5,1)          # *********Covariance between instruments (constant) [.5, 1]
covue_  <- c(0.5, 0.9)       # Covariance between errors (.5, .9)
sigmae  <- 1.2               # Structural variance
sigmau  <- 1                 # Reduced form variance
local   <- TRUE              # Indicator for fixed or local-to-zero invalid moments
het     <- TRUE              # Heteroskedatic error
# Misc Settings
nmoment <- NULL              # If NULL set the number of moments growth rate automatically, else c(q, s, si) c(11,3,4)  
alpha   <- .05               # Significance level for Downward Testing
# Main Code -----------------------------------------------------------------

for(covue in covue_){
for(czz in czz_) {
for(pcoeff in pcoeff_){
  for(n in n_){
    print(paste("Process started at:", Sys.time()," for Czz=", czz, ", pcoeff=", pcoeff, ", n=", n))
    ptm <- proc.time() ### Starts the clock
# Autoset. change q, s and si manually for other setups
    q  <- as.integer(ceiling(sqrt(n)))  # Total number of moments grow at the rate sqrt(n)
    s  <- as.integer(ceiling(sqrt(q)))  # known to be valid moments
    si <- as.integer(ceiling({q-s}/2))  # Invalid moments
    if(length(nmoment)==3){q <- as.integer(nmoment[1]); s <- as.integer(nmoment[2]); si <- as.integer(nmoment[3])} # Override autoset
    qs <- as.integer(q-s)               # Moments to be tested
    sv <- as.integer(qs-si)             # Number of (unknown) valid moments
    truth <- seq_len(s + sv)
    k  <- rep(as.integer(s), 2^qs)      # Dimension of Z at each iteration
    for(ii in seq_len(qs)) k <- k + rep(c(rep(as.integer(1), 2^{qs-ii} ), rep(as.integer(0), 2^{qs-ii} )), length=2^qs )
    gc()                                # k vector increases the memory usage exponentially
    c_b       <- k-as.integer(1)              # (c-b) in andrews c:#instruments b:#parameters. b=length(theta)=1
    covuv     <- matrix(c(sigmae, covue, covue, sigmau), nrow=2, ncol=2) # Covariance between errors
    cholsigma <- chol(gensigma(q=q, s=s, si=si, d=d, czz=czz, covuv=covuv, local=local))
    load(paste("combs", n, "auto.list", sep="")) # Load the indexes for combinations
    if(length(nmoment)==3) assign( "combinations.list", lapply(seq_len(2^qs), function(x) indexes(x, qs, s))) # for custom number of moments
#Define methods
    jbase.methods <- c("full","min","aic","bic","hqic","dt","ut","oracle")      # J methods
    j.methods     <- paste("j", jbase.methods, sep="")      # J methods
    jcue.methods  <- paste("jcue", jbase.methods, sep="")
    s.methods     <- c("aic", "bic", "cp")                            # adaptive lasso methods
    js.methods    <- s.methods                                        # J with alasso moments 
    jscue.methods <- js.methods
    index.lst     <- c(s.methods, j.methods, jcue.methods)            # List of indexes
    index.mat     <- matrix(as.integer(0), reps, length(index.lst))   # Matrix with Indexes
    r2.mat        <- index.mat                                        # Matrix R2 structural equation
    theta1.lst    <- c(s.methods, j.methods, js.methods, jcue.methods, jscue.methods)
    o.methods     <- theta1.lst                         # Okui
    theta.lst     <- c(theta1.lst, o.methods )
    theta.mat     <- matrix(           0 , reps, length(theta.lst))   # Matrix with theta hats
    cover.mat     <- matrix(as.integer(0), reps, length(theta.lst))# Matrix of coverage of the alpha C.I.
    if(parallel) {                                                # Begin Parallel Job
      library(parallel)                                           # Load the library only if needed
      cl   <- makePSOCKcluster(detectCores())                     # Create Parallel Sock Cluster
      setDefaultCluster(cl)                                       # Unique cluster for all parallel operations
      clusterEvalQ( NULL, c(library(RcppArmadillo)))              # Auto-Executed commands in the Clusters
      clusterExport(NULL, c('gmm', 'gmm_CUE','n','q','s','qs','combinations.list','alpha','theta')) # Variables exported to the Clusters
    }
    start0 <- start
    st <- paste(start0, "-", (start0+reps-1), sep="") # Save the starting value for output
# Iterations begin here:
    for (iter in seq_len(reps)){
    	if(iter%%100==0) print(paste("Iteration", iter))
      set.seed(seed + start0)
      start0<- start0 + 1
      data  <- gendata(n=n, q=q, s=s, pcoeff=pcoeff, strong=strong, cholsigma=cholsigma, theta=theta, het=het)
      y1    <- data$Y1            # This two do not vary across submodels
      y2    <- data$Y2            # 
      zz    <- data$Z             # This one does
      data.shrk <-shrinkdata(y1, y2, zz, n)# Transformed data for shrinkage estimation
      if(parallel) {                                                # Begin Parallel Job
  	    clusterExport(NULL, c('y1', 'y2', 'zz'))                    # Variables exported to the Clusters
        res1 <- parLapply(NULL,seq_len(2^qs),function(x) gmm(zz[, combinations.list[[x]] ])) # parLapply function
      } else res1 <- lapply(seq_len(2^qs),   function(x) gmm(zz[, combinations.list[[x]] ])) # Regular lapply command
# Prepare the results for analysis
      resmat      <- do.call(rbind, res1)                           # Turns the results into a 2^qs x 2 matrix
      theta.vec   <- resmat[ ,1]                                    # Vector of theta for each combination
      js.vec      <- resmat[ ,2]                                    # Vector of J-value for each combination
      thetacue.vec<- resmat[ ,3]                                    # Vector of thetaCUE for each combination
      jcue.vec    <- resmat[ ,4]                                    # Vector of JCUE-value for each combination
      jscover.vec <- resmat[ ,5]                                    # Vector of Coverage of GMM for each combination
      jcuecover.vec <- resmat[ ,6]                                    # Vector of Coverage of CUE for each combination
      model2      <- alasso(data.shrk$y2z, data.shrk$y1z, noshrink=1) # Shrinkage estimation (Adaptive LASSO)
# Selected Indexes  (SECOND STAGE J)                            # noshrink: do not shrink the reduced form parameter
      jfull.index <- 1
      jmin.index  <- which.min( js.vec                          )
      jaic.index  <- which.min( js.vec - 2  * c_b               )   # Andrews (2000) AIC
      jbic.index  <- which.min( js.vec -      c_b * log(n)      )   # Andrews (2000) BIC
      jhqic.index <- which.min( js.vec - 2.1* c_b * log(log(n)) )   # Andrews (2000) HQIC
      jdt.index   <- dt(js.vec, k, alpha)
      jut.index   <- dt(js.vec, k, alpha, upward=TRUE)
      joracle.index<- rev.indexes(truth, qs, s)         # Oracle estimate
# Selected Indexes  (CUE J) 
      jcuefull.index <- 1
      jcuemin.index  <- which.min( jcue.vec                          )
      jcueaic.index  <- which.min( jcue.vec - 2  * c_b               )   # Hong (2003) AIC
      jcuebic.index  <- which.min( jcue.vec -      c_b * log(n)      )   # Hong (2003) BIC
      jcuehqic.index <- which.min( jcue.vec - 2.1* c_b * log(log(n)) )   # Hong (2003) HQIC
      jcuedt.index   <- dt(jcue.vec, k, alpha)
      jcueut.index   <- dt(jcue.vec, k, alpha, upward=TRUE)
      jcueoracle.index<- rev.indexes(truth, qs, s)         # Oracle estimate
# For adaptive lasso indexes. seq_len(s) is for including the known.
      for(smethod in s.methods) assign(paste("alasso", smethod, ".index", sep=""), 
       rev.indexes(c(seq_len(s),which({model2$betas.orig[which.min(eval(parse(text=paste("model2$", smethod, sep="")))),]==0}[2:{qs+1}])+s), qs,s)) 
#equivalent for aic:
#alassoaic       <- model2$betas.orig[which.min(model2$aic) ,]
#alassoaic.index <- rev.indexes(c(seq_len(s), which({alassoaic!=0}[2:{qs+1}])+s), qs,s)
# Adaptive lasso coverage (not available)
#COVERAGE
alasso.cover  <- rep(FALSE, length(s.methods))
names(alasso.cover) <-c("alassoaic", "alassobic", "alassocp")
j.cover <- sapply(paste(j.methods, ".index", sep=""), function(x) eval(parse(text=paste("jscover.vec[",x,"]", sep="" ))))
names(j.cover) <- j.methods
jcue.cover <- sapply(paste(jcue.methods, ".index", sep=""), function(x) eval(parse(text=paste("jcuecover.vec[",x,"]", sep="" ))))
names(jcue.cover) <- jcue.methods

# JMETHODS
# Get the parameters
      for(jmethod in j.methods) assign(paste(jmethod,".theta", sep=""), theta.vec[get(paste(jmethod,".index",sep=""))])
      for(jcuemethod in jcue.methods) assign(paste(jcuemethod,".theta", sep=""), thetacue.vec[get(paste(jcuemethod,".index",sep=""))])
      for(smethod in s.methods) assign(paste("alasso",smethod,".theta", sep=""), 
                    model2$betas.orig[which.min(eval(parse(text=paste("model2$", smethod, sep="")))),][1] )
#equivalent for aic: alassoaic.theta <- alassoaic[1]
# J using alasso moments:
      for(jsmethod in js.methods) assign(paste("js",jsmethod,".theta", sep=""),
                    theta.vec[eval(parse(text=paste("alasso", jsmethod, ".index", sep="")))])
#equivalent for aic: jsaic.theta <- theta.vec[alassoaic.index]
# Coverage
     for(jsmethod in js.methods) assign(paste("js",jsmethod,".cover", sep=""),
                   jscover.vec[eval(parse(text=paste("alasso", jsmethod, ".index", sep="")))])
     js.cover <- sapply(paste("js",js.methods, ".cover", sep=""), function(x) eval(parse(text=paste(x))))

# CUE using alasso moments:
      for(jscuemethod in jscue.methods) assign(paste("jscue",jscuemethod,".theta", sep=""),
                    thetacue.vec[eval(parse(text=paste("alasso", jscuemethod, ".index", sep="")))])
#equivalent for aic: jsaic.theta <- theta.vec[alassoaic.index]
# Coverage
for(jscuemethod in jscue.methods) assign(paste("jscue",jscuemethod,".cover", sep=""),
                   jcuecover.vec[eval(parse(text=paste("alasso", jscuemethod, ".index", sep="")))])
     jscue.cover <- sapply(paste("jscue",js.methods, ".cover", sep=""), function(x) eval(parse(text=paste(x))))

# Okui's shrinkage: 
## Methods
      oalasso.vec   <- c(paste("alasso", s.methods, ".index", sep="")) # "alassoaic.index" "alassobic.index" "alassocp.index"
      oj.vec        <- c(paste(j.methods, ".index", sep=""))           # "jmin.index"      "jaic.index"      ... "joracle.index"
      ojs.vec       <- c(paste("js", s.methods, ".index", sep=""))     # "jsaic.index"     "jsbic.index" "jscp.index"
      ojcue.vec     <- c(paste(jcue.methods, ".index", sep="") )       # "jcuemin.index"   "jcueaic.index"   ... "jcueoracle.index"
      ojscue        <- c(paste("jscue", s.methods, ".index", sep=""))  # "jscueaic.index" "jscuebic.index" "jscuecp.index"
## Parameters and Coverage
### Adaptive lasso
      oalassoaic  <- okui(alassoaic.index, param=alassoaic.theta)
       oalassoaictheta <- oalassoaic$theta
       oalassoaiccover <- oalassoaic$cover
      oalassobic  <- okui(alassobic.index, param=alassobic.theta)
       oalassobictheta <- oalassobic$theta
       oalassobiccover <- oalassobic$cover
      oalassocp   <- okui(alassocp.index,  param=alassocp.theta)
       oalassocptheta  <- oalassocp$theta
       oalassocpcover  <- oalassocp$cover
      oalassotheta.vec <- c(oalassoaictheta, oalassobictheta, oalassocptheta)
      oalassocover.vec <- c(oalassoaiccover, oalassobiccover, oalassocpcover)
      names(oalassotheta.vec) <- oalasso.vec
      names(oalassocover.vec) <- oalasso.vec
### 2 Step GMM
      ojtheta.vec           <- sapply(oj.vec,      function(x) okui(eval(parse(text=x)), method="j2s"   )$theta)
      ojcover.vec           <- sapply(oj.vec,      function(x) okui(eval(parse(text=x)), method="j2s"   )$cover)
### Adaptive lasso and 2 Step GMM 
      ojalassoaic <- okui(alassoaic.index, param=jsaic.theta)
       ojalassoaictheta <- ojalassoaic$theta
       ojalassoaiccover <- ojalassoaic$cover
      ojalassobic <- okui(alassobic.index, param=jsbic.theta)
       ojalassobictheta <- ojalassobic$theta
       ojalassobiccover <- ojalassobic$cover
      ojalassocp  <- okui(alassocp.index,  param=jscp.theta)
       ojalassocptheta  <- ojalassocp$theta
       ojalassocpcover  <- ojalassocp$cover
      ojalassotheta.vec <- c(ojalassoaictheta, ojalassobictheta, ojalassocptheta)
      ojalassocover.vec <- c(ojalassoaiccover, ojalassobiccover, ojalassocpcover)
      names(ojalassotheta.vec) <- ojs.vec
      names(ojalassocover.vec) <- ojs.vec
### Continuous updating GMM
      ojcuetheta.vec        <- sapply(ojcue.vec,   function(x) okui(eval(parse(text=x)), method="jcue"  )$theta)
      ojcuecover.vec        <- sapply(ojcue.vec,   function(x) okui(eval(parse(text=x)), method="jcue"  )$cover)
### Adaptive lasso and Continuous updating GMM
      ojcuealassoaic <- okui(alassoaic.index, param=jscueaic.theta)
       ojcuealassoaictheta <- ojcuealassoaic$theta
       ojcuealassoaiccover <- ojcuealassoaic$cover
      ojcuealassobic <- okui(alassobic.index, param=jscuebic.theta)
       ojcuealassobictheta <- ojcuealassobic$theta
       ojcuealassobiccover <- ojcuealassobic$cover
      ojcuealassocp  <- okui(alassocp.index,  param=jscuecp.theta)
       ojcuealassocptheta <- ojcuealassocp$theta
       ojcuealassocpcover <- ojcuealassocp$cover
      ojcuealassotheta.vec<- c(ojcuealassoaictheta, ojcuealassobictheta, ojcuealassocptheta)
      ojcuealassocover.vec<- c(ojcuealassoaiccover, ojcuealassobiccover, ojcuealassocpcover)
      names(ojcuealassotheta.vec) <- ojscue
      names(ojcuealassocover.vec) <- ojscue
## All Okui's estimators:
      o.theta <- c(oalassotheta.vec, ojtheta.vec, ojalassotheta.vec, ojcuetheta.vec, ojcuealassotheta.vec)
      o.cover <- c(oalassocover.vec, ojcover.vec, ojalassocover.vec, ojcuecover.vec, ojcuealassocover.vec)
      okui.names.indexes <- names(o.theta)
      names(o.theta) <- unlist(strsplit(names(o.theta), split=".", fixed=TRUE))[{{1:length(names(o.theta))}*2}-1] # Delete the .index part
      names(o.theta) <- paste("o_", names(o.theta), sep="")
      names(o.cover) <- names(o.theta)
# Indexes and names
      all.theta.str <- c( paste("alasso",s.methods, sep=""),     # adaptive lasso
                        j.methods,                                 # GMM 2s
                         paste("js", js.methods, sep=""),          # Alasso -> GMM 2s
                           jcue.methods,                           # GMM CUE
                             paste("jscue", js.methods, sep="")       # Alasso -> GMM CUE
      )
      all.index.str <- c( paste("alasso",s.methods, sep=""),     # adaptive lasso
                            j.methods,                             # GMM 2s
                              jcue.methods                         # GMM CUE
      )
      all.indexes <- sapply(paste(all.index.str, ".index", sep=""), function(x) eval(parse(text=x)))
      names(all.indexes) <- all.index.str
      all.thetas <- sapply(paste(all.theta.str, ".theta", sep=""), function(x) eval(parse(text=x)))
      all.covers <- c(alasso.cover, j.cover, jcue.cover, js.cover, jscue.cover, o.cover)
#      names(all.thetas) <- all.theta.str
# R2
      r2.vec      <- sapply(all.indexes, function(x) r2struct(x))
# output
      index.mat[iter,] <- all.indexes     # Indexes
      theta.mat[iter,] <- c(all.thetas, o.theta)
      cover.mat[iter,] <- as.integer(all.covers)
      r2.mat[iter, ]   <- r2.vec
      if (n >=300) print(paste("Iteration ", iter, "of", reps))
    } # <--------------- Iterations end here
    colnames(index.mat) <- all.index.str
    colnames(theta.mat) <- c(all.theta.str, names(o.theta))
    colnames(cover.mat) <- colnames(theta.mat)
    colnames(r2.mat)    <- all.index.str
    truth <- c(rep(TRUE,s),rep(TRUE,sv), rep(FALSE,si))
    if(parallel) stopCluster(cl)                                  # Kills the Cluster / End Parallel job
# Save the results:
    if(length(nmoment==0)) setup <- 1 else setup <- 2
    if(pcoeff > 1) ws <- "strong" else ws <- "weak"
    if(local)  lc <- "local" else lc <- "constant"
    scalars <- list(n=n, corr=d, pcoeff=pcoeff, theta=theta, czz=czz, covue=covue, local=local, het=het, alpha=alpha, totalmoment=q, knownvalid=s, valid=sv, reps=reps, start=st)
    output  <- list(scalars=scalars, truth=truth, theta.mat=theta.mat, index.mat=index.mat, r2.mat=r2.mat, cover.mat=cover.mat)
    string  <- paste("setup.",setup, "&n.", n, "&q.", q, "&valid.", sv, "&known.", s, "&strength.", ws,"&czz", czz, "&covue", covue,"&type.",  lc, "&sims.", st, sep="")
    save(output, file=string)    # save output to the working directory
# Timer
    totalseconds <- round((proc.time() - ptm)[3],2)### Stops the clock
    hours        <- trunc(totalseconds/3600)
    minutes      <- trunc(totalseconds/60) - hours*60
    seconds      <- round(totalseconds - minutes*60 - hours*3600, 2)
    timer        <- paste(" Elapsed Time:", hours, "hours,", minutes, "minutes and ", seconds, "seconds")
    if(hours == 0) timer <- paste(" Elapsed Time:", minutes, "minutes and ", seconds, "seconds")
    if(hours == 0 & minutes == 0) timer <- paste(" Elapsed Time:", seconds, "seconds")
    print(timer)
  }   # <----------------- Iterations for "n" end here!
}   # <----------------- Iterations for "pcoeff" end here!
} # <----------------- Iterations for "czz" end here!
}# <----------------- Iterations for "covue" end here!
