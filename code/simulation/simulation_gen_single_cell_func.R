
############## Run Simulation for manuscript ###############
source('code/CIBERSORT.R')
library(MEAD)
library(gtools)
library(Matrix)
library(MuSiC)
#'@description Generate data: based on real Xin data, generate uncorrelated or correlated individual reference matrices.
#'Compare methods: ols, adjust for measurement error, adjust for both measurement error and correlation, adjust for both measurement error and correlation and add weights.
#'For simplicity, we only simulate data at individual level. - only generate X from U;

#'@param ref,sigma2 population gene expression matrix U, of dimension G by K; and variance matrix
#'@param p cell type proportion
#'@param n_bulk number of bulk data
#'@param dirichlet whether generate random p
#'@param n_bulk_for_cor number of bulk data generated for inferring correlations
simu_study = function(ref,
                      sigma2,
                      p,
                      n_bulk=100,
                      dirichlet = TRUE,
                      dirichlet.scale = 10,
                      d=500,
                      nreps = 100,
                      bulk_lib_size = 500,
                      n_ref = 10,
                      printevery=10,
                      n_bulk_for_cor = 100,
                      cross_plat_bias_gamma_var = 0.3,
                      n_single_cell_for_ref = 50,
                      nfold = 10,
                      bulk_nb_theta = 10,
                      single_cell_nb_theta = 5,
                      softplus_a = 10
){
  
  
  genp = function(K){
    p = runif(K)
    p/sum(p)
  }
  
  G = nrow(ref)
  K = ncol(ref)
  
  is.indep = (d==0)
  if(d!=0){
    R = matrix(0,nrow=G,ncol=G)
    for(i in 1:G){
      for(j in i:min(i+d,G)){
        R[i,j] = max(1-abs(i-j)/d,0)
      }
    }
    R = R+t(R) - diag(G)
    R = Matrix(R,sparse = TRUE)
  }
  
  
  
  
  
  
  
  
  
  ########
  # rmse #
  ########
  
  # rmse:
  ## 1. ols
  ## 2. MEAD, NO WEIGHT
  ## 3. MEAD Marker
  ## 4. MEAD
  ## 5. MuSiC
  ## 6. CIBERSORT
  
  p_hat_ols = array(dim = c(K,n_bulk,nreps))
  p_hat = array(dim = c(K,n_bulk,nreps))
  p_hat_marker = array(dim = c(K,n_bulk,nreps))
  p_hat_weight = array(dim = c(K,n_bulk,nreps))
  
  p_hat_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_marker_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_sp = array(dim = c(K,n_bulk,nreps))
  
  p_hat_music = array(dim = c(K,n_bulk,nreps))
  p_hat_ciber = array(dim = c(K,n_bulk,nreps))
  
  ############
  # Variance #
  ############
  
  # variance:
  ## 1. OLS
  ## 2. MEAD, no weight, no cor
  ## 3. MEAD, no corr
  ## 4. MEAD, known cor
  ## 5. MEAD, est cor 0.1
  ## 6. MEAD, est cor 0.3
  ## 7. MEAD, est cor 0.5
  ## 8. MEAD+cv, known cor
  ## 9. MEAD+cv, est cor 0.1
  ## 10. MEAD+cv, est cor 0.3
  ## 11. MEAD+cv, est cor 0.5
  
  p_hat_ols_se = array(dim = c(K,n_bulk,nreps))
  p_hat_se = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_knowncor = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor01 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor03 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor05 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_knowncor = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor01 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor03 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor05 = array(dim = c(K,n_bulk,nreps))
  
  p_hat_se_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_knowncor_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor01_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor03_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor05_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_knowncor_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor01_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor03_sp = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor05_sp = array(dim = c(K,n_bulk,nreps))
  
  
  
  
  gene_names = rownames(ref)
  celltypes = colnames(ref)
  #true_betas = matrix(nrow=nreps,ncol=n_bulk*K)
  ## pre calculate MLN and generate independent normal
  
  norm.ref = matrix(nrow=G,ncol=K)
  norm.Sigma.chol = list()
  norm.Sigma = matrix(nrow=G,ncol=K)
  if(!is.indep){
    chol.R = chol(R)
  }
  for(k in 1:K){
    norm.ref[,k] = log(ref[,k]^2/sqrt(ref[,k]^2+sigma2[,k]))
    norm.s = sqrt(log(1+sigma2[,k]/ref[,k]^2))
    if(!is.indep){
      norm.Sigma.chol[[k]] = t(norm.s*t(chol.R))
    }else{
      norm.Sigma[,k] = norm.s^2
    }
  }
  
  # true R01
  cor.idx = which(R!=0,arr.ind = T)
  R01_true = sparseMatrix(i = cor.idx[,1],j = cor.idx[,2],dims=c(G,G))
  
  # estimated R01
  X_array = array(dim=c(G,K,n_bulk_for_cor))
  for(k in 1:K){
    if(is.indep){
      X_array[,k,] = exp(matrix(rnorm(G*n_bulk_for_cor,norm.ref[,k],sqrt(norm.Sigma[,k])),ncol=n_bulk_for_cor))
    }else{
      X_array[,k,] = t(exp(mvnfast::rmvn(n_bulk_for_cor,mu = norm.ref[,k],sigma = norm.Sigma.chol[[k]],isChol = TRUE)))
    }
  }
  mb = apply(X_array,3,function(z){z%*%genp(K)})
  thetab = apply(mb,2,function(z){z/sum(z)})
  bulk_for_cor = matrix(rnbinom(G*n_bulk_for_cor,size = single_cell_nb_theta,prob = single_cell_nb_theta/(single_cell_nb_theta + bulk_lib_size*G*thetab)),nrow=G)
  bulk_for_cor = apply(bulk_for_cor,2,function(z){z/sum(z)*1e6})
  rownames(bulk_for_cor) = gene_names
  #cor.idx = get_cor_pairs2(bulk_for_cor,alpha=alpha.cor,method=cor_method)
  R01_01 = get_R01(bulk_for_cor,0.1)
  R01_03 = get_R01(bulk_for_cor,0.3)
  R01_05 = get_R01(bulk_for_cor,0.5)
  
  # all_fit = list()
  true_p = array(dim = c(K,n_bulk,nreps))
  
  for(reps in 1:nreps){
    
    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}
    
    ## generate group p's
    
    if(dirichlet){
      b = t(gtools::rdirichlet(n_bulk,p*dirichlet.scale))
    }else{
      b = matrix(p,ncol=nb1,nrow=K)
    }
    true_p[,,reps] = b
    
    ## generate individual reference matrices
    
    n.temp = n_ref+n_bulk
    X_array = array(dim=c(G,K,n.temp))
    
    for(k in 1:K){
      if(is.indep){
        X_array[,k,] = exp(matrix(rnorm(G*n.temp,norm.ref[,k],sqrt(norm.Sigma[,k])),ncol=n.temp))
      }else{
        X_array[,k,] = t(exp(mvnfast::rmvn(n.temp,mu = norm.ref[,k],sigma = norm.Sigma.chol[[k]],isChol = TRUE)))
      }
    }
    
    #browser()
    
    X_array_bulk = X_array[,,1:n_bulk]
    
    X_array_ref = X_array[,,(n_bulk+1):(n_bulk+n_ref)]
    
    
    # # add variance to mimic the average of cell counts
    # for(r in 1:n_ref){
    #   for(k in 1:K){
    #     X_array_ref[,k,r] = rgamma(G,shape=gamma_shape,scale = X_array_ref[,k,r]/gamma_shape)
    #   }
    # }
    
    # For each individual,  for each cell type ,generate single cells, then take average as X_array_ref. 
    for (r in 1:n_ref){
      X_ref_indi = X_array_ref[,,r]
      temp_ref = 0
      for(cc in 1:n_single_cell_for_ref){
        temp_ref = temp_ref + matrix(rnbinom(G*K,size = single_cell_nb_theta,prob = single_cell_nb_theta/(single_cell_nb_theta + X_ref_indi)),nrow=G)
      }
      X_array_ref[,,r] = temp_ref / n_single_cell_for_ref
    }
    
    print("prop of 0's in indi ref matrix")
    temp = 0
    for(r in 1:n_ref){
      temp = temp + apply(X_array_ref[,,r],2,function(xx){sum(xx==0)/G})
    }
    print(sum(temp)/n_ref/K)
    
    # select marker gene
    ref_samples = c()
    for(k in 1:K){
      temp = X_array_ref[,k,]
      colnames(temp) = rep(colnames(ref)[k],n_ref)
      ref_samples = cbind(ref_samples,temp)
    }
    rownames(ref_samples) = gene_names
    sig_mat = build_signature_matrix_CIBERSORT(ref_samples)
    sig_mat = sig_mat[,match(celltypes,colnames(sig_mat))]
    sig_genes = rownames(sig_mat)
    marker.idx = match(sig_genes,gene_names)
    
    
    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
    mb = do.call(cbind,mb)
    # mb = mb*rnorm(G,1,bulk_bias_sd)
    mb = mb*rgamma(G,shape = 1/cross_plat_bias_gamma_var, rate = 1/cross_plat_bias_gamma_var)
    thetab = apply(mb,2,function(z){z/sum(z)})
    
    
    #browser()
    y = matrix(rnbinom(G*n_bulk,size = bulk_nb_theta,prob=bulk_nb_theta / (bulk_nb_theta + bulk_lib_size*G*thetab)),nrow=G)
    rownames(y) = gene_names
    
    
    # fit model
    
    #browser()
    
    X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
    rownames(X) = gene_names
    colnames(X) = 1:K
    V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref
    cat("the rows in X that are NA:", which(is.na(rowSums(X))), "\n")
    cat("the rows in X that are all 0's:", which(rowSums(X) == 0), "\n")
    cat("the rows in V that are all 0's:", which(rowSums(V) == 0), "\n")
    cat("the rows in V that are NA:", which(is.na(rowSums(V))), "\n")
    
    
    # ols fit
    fit.ols = ols_hc3(y,X)
    p_hat_ols[,,reps] = fit.ols$p_hat
    p_hat_ols_se[,,reps] = fit.ols$p_hat_se
    
    
    # adjust for measurement error, do not adjust for correlation
    fit.err = MEAD_est(y=y,
                       X=X,
                       Vg=V,
                       w=1,
                       hc.type='hc3',
                       R01=NULL)
    p_hat[,,reps] = fit.err$p_hat
    p_hat_se[,,reps] = fit.err$p_hat_se
    
    fit.err.sp = MEAD_est(y=y,
                       X=X,
                       Vg=V,
                       w=1,
                       hc.type='hc3',
                       R01=NULL,
                       beta_tilde_transformation = list(method='softplus', a = softplus_a))
    p_hat_sp[,,reps] = fit.err.sp$p_hat
    p_hat_se_sp[,,reps] = fit.err.sp$p_hat_se
    
    
    # use marker, no weight
    fit.marker = MEAD_est(y=y[marker.idx,],
                          X=X[marker.idx,],
                          Vg = V[marker.idx,],
                          w=1,
                          calc_var = FALSE)
    p_hat_marker[,,reps] = fit.marker$p_hat
    
    # use marker, no weight
    fit.marker.sp = MEAD_est(y=y[marker.idx,],
                          X=X[marker.idx,],
                          Vg = V[marker.idx,],
                          w=1,
                          calc_var = FALSE,
                          beta_tilde_transformation = list(method='softplus', a = softplus_a))
    p_hat_marker_sp[,,reps] = fit.marker.sp$p_hat
    
    
    # add weight, no cor
    fit.vash = vashr::vash(sqrt(rowSums(V)),df=n_ref-1)
    w = 1/(fit.vash$sd.post)^2
    
    fit.mead.nocor = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=NULL)
    p_hat_weight[,,reps] = fit.mead.nocor$p_hat
    p_hat_weight_se[,,reps] = fit.mead.nocor$p_hat_se
    
    fit.mead.nocor.sp = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=NULL,
                              beta_tilde_transformation = list(method='softplus', a = softplus_a))
    p_hat_weight_sp[,,reps] = fit.mead.nocor.sp$p_hat
    p_hat_weight_se_sp[,,reps] = fit.mead.nocor.sp$p_hat_se
    
    
    fit.mead.knowncor = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='hc3',
                                 R01=R01_true)
    p_hat_weight_se_knowncor[,,reps] = fit.mead.knowncor$p_hat_se
    
    fit.mead.knowncor.sp = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='hc3',
                                 R01=R01_true,
                                 beta_tilde_transformation = list(method='softplus', a = softplus_a))
    p_hat_weight_se_knowncor_sp[,,reps] = fit.mead.knowncor.sp$p_hat_se
    
    fit.mead.cor01 = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_01)
    p_hat_weight_se_cor01[,,reps] = fit.mead.cor01$p_hat_se
    
    fit.mead.cor01.sp = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_01,
                              beta_tilde_transformation = list(method='softplus', a = softplus_a))
    p_hat_weight_se_cor01_sp[,,reps] = fit.mead.cor01.sp$p_hat_se
    
    
    fit.mead.cor03 = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_03)
    
    p_hat_weight_se_cor03[,,reps] = fit.mead.cor03$p_hat_se
    
    fit.mead.cor03.sp = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_03,
                              beta_tilde_transformation = list(method='softplus', a = softplus_a))
    
    p_hat_weight_se_cor03_sp[,,reps] = fit.mead.cor03.sp$p_hat_se
    
    fit.mead.cor05 = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_05)
    
    p_hat_weight_se_cor05[,,reps] = fit.mead.cor05$p_hat_se
    
    fit.mead.cor05.sp = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_05,
                              beta_tilde_transformation = list(method='softplus', a = softplus_a))
    
    p_hat_weight_se_cor05_sp[,,reps] = fit.mead.cor05.sp$p_hat_se
    
    folds = cut(1:G,breaks = nfold,labels = F)
    
    fit.mead.cv.knowncor = MEAD_est(y=y,
                                    X=X,
                                    Vg=V,
                                    w=w,
                                    hc.type='cv_indep',
                                    R01=R01_true,
                                    folds=folds)
    
    p_hat_weight_se_cv_knowncor[,,reps] = fit.mead.cv.knowncor$p_hat_se
    
    fit.mead.cv.knowncor.sp = MEAD_est(y=y,
                                    X=X,
                                    Vg=V,
                                    w=w,
                                    hc.type='cv_indep',
                                    R01=R01_true,
                                    folds=folds,
                                    beta_tilde_transformation = list(method='softplus', a = softplus_a))
    
    p_hat_weight_se_cv_knowncor_sp[,,reps] = fit.mead.cv.knowncor.sp$p_hat_se
    
    fit.mead.cv.cor01 = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_01,
                                 folds=folds)
    
    p_hat_weight_se_cv_cor01[,,reps] = fit.mead.cv.cor01$p_hat_se
    
    fit.mead.cv.cor01.sp = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_01,
                                 folds=folds,
                                 beta_tilde_transformation = list(method='softplus', a = softplus_a))
    
    p_hat_weight_se_cv_cor01_sp[,,reps] = fit.mead.cv.cor01.sp$p_hat_se
    
    
    
    fit.mead.cv.cor03 = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_03,
                                 folds=folds)
    
    p_hat_weight_se_cv_cor03[,,reps] = fit.mead.cv.cor03$p_hat_se
    
    fit.mead.cv.cor03.sp = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_03,
                                 folds=folds,
                                 beta_tilde_transformation = list(method='softplus', a = softplus_a))
    
    p_hat_weight_se_cv_cor03_sp[,,reps] = fit.mead.cv.cor03.sp$p_hat_se
    
    
    fit.mead.cv.cor05 = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_05,
                                 folds=folds)
    
    p_hat_weight_se_cv_cor05[,,reps] = fit.mead.cv.cor05$p_hat_se
    
    fit.mead.cv.cor05.sp = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_05,
                                 folds=folds,
                                 beta_tilde_transformation = list(method='softplus', a = softplus_a))
    
    p_hat_weight_se_cv_cor05_sp[,,reps] = fit.mead.cv.cor05.sp$p_hat_se
    
    
    # run music
    Sigma.music = t(apply(X_array_ref,c(1),function(z){diag(cov(t(z),use = 'complete.obs'))}))
    Est.prop.allgene = NULL
    for(bb in 1:n_bulk){
      fit.music = music.basic(y[,bb],X,S=1,Sigma.music,iter.max=1000,nu=1e-4,eps=0.01)
      Est.prop.allgene = cbind(Est.prop.allgene, fit.music$p.weight)
    }
    p_hat_music[,,reps] = Est.prop.allgene
    
    # run cibersort
    fit_cibersort = CIBERSORT(sig_mat,data.frame(GeneSymbol = rownames(y),y))
    c.order = match(celltypes,colnames(fit_cibersort[,1:K]))
    p_hat_ciber[,,reps] = t(fit_cibersort[,1:K])[c.order,]
    
    rmse = function(x,y){sqrt(mean((x-y)^2))}
    rmse(fit.ols$p_hat,b)
    rmse(fit.err$p_hat,b)
    rmse(fit.marker$p_hat,b)
    rmse(fit.mead.nocor$p_hat,b)
    rmse(Est.prop.allgene,b)
    rmse(t(fit_cibersort[,1:K])[c.order,],b)
    
  }
  
  return(list(p_hat_ols = p_hat_ols,
              p_hat = p_hat,
              p_hat_marker = p_hat_marker,
              p_hat_weight = p_hat_weight,
              p_hat_sp = p_hat_sp,
              p_hat_marker_sp = p_hat_marker_sp,
              p_hat_weight_sp = p_hat_weight_sp,
              p_hat_music = p_hat_music,
              p_hat_ciber = p_hat_ciber,
              
              p_hat_ols_se = p_hat_ols_se,
              p_hat_se = p_hat_se,
              p_hat_weight_se = p_hat_weight_se,
              p_hat_weight_se_knowncor = p_hat_weight_se_knowncor,
              p_hat_weight_se_cor01 = p_hat_weight_se_cor01,
              p_hat_weight_se_cor03 = p_hat_weight_se_cor03,
              p_hat_weight_se_cor05 = p_hat_weight_se_cor05,
              p_hat_weight_se_cv_knowncor = p_hat_weight_se_cv_knowncor,
              p_hat_weight_se_cv_cor01 = p_hat_weight_se_cv_cor01,
              p_hat_weight_se_cv_cor03 = p_hat_weight_se_cv_cor03,
              p_hat_weight_se_cv_cor05 = p_hat_weight_se_cv_cor05,
              
              p_hat_se_sp = p_hat_se_sp,
              p_hat_weight_se_sp = p_hat_weight_se_sp,
              p_hat_weight_se_knowncor_sp = p_hat_weight_se_knowncor_sp,
              p_hat_weight_se_cor01_sp = p_hat_weight_se_cor01_sp,
              p_hat_weight_se_cor03_sp = p_hat_weight_se_cor03_sp,
              p_hat_weight_se_cor05_sp = p_hat_weight_se_cor05_sp,
              p_hat_weight_se_cv_knowncor_sp = p_hat_weight_se_cv_knowncor_sp,
              p_hat_weight_se_cv_cor01_sp = p_hat_weight_se_cv_cor01_sp,
              p_hat_weight_se_cv_cor03_sp = p_hat_weight_se_cv_cor03_sp,
              p_hat_weight_se_cv_cor05_sp = p_hat_weight_se_cv_cor05_sp,
              
              
              true_p = true_p,
              p = p,
              parameters = list(n_bulk=n_bulk,
                                dirichlet = dirichlet,
                                dirichlet.scale = dirichlet.scale,
                                d=d,
                                nreps = nreps,
                                bulk_lib_size = bulk_lib_size,
                                n_ref = n_ref,
                                n_bulk_for_cor = n_bulk_for_cor,
                                cross_plat_bias_gamma_var = cross_plat_bias_gamma_var,
                                n_single_cell_for_ref = n_single_cell_for_ref,
                                nfold = nfold,
                                bulk_nb_theta = bulk_nb_theta,
                                single_cell_nb_theta = single_cell_nb_theta,
                                softplus_a = softplus_a)))
  
}




get_rmse = function(p_hat,b){
  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  nreps = dim(p_hat)[3]
  rmses = c()
  for(i in 1:nb){
    err = c()
    for(j in 1:nreps){
      err[j] = sum((p_hat[,i,j]-b[,i])^2)
    }
    rmses[i] = sqrt(mean(err))
  }
  names(rmses) = paste('bulk',1:nb)
  round(rmses,3)
  
}

get_rmse_array = function(x,y){
  nn = dim(x)[3]
  ses = c()
  for(i in 1:nn){
    ses[i] = mean((x[,,i]-y[,,i])^2)
  }
  sqrt(mean(ses))
}

rmse = function(x,y){
  sqrt(mean((x-y)^2))
}

get_coverage_p = function(p_hat,p_hat_se,b_array){
  
  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  z = array(dim = dim(p_hat))
  for(i in 1:dim(z)[3]){
    z[,,i] = (p_hat[,,i]-b_array[,,i])/p_hat_se[,,i]
  }
  crg = apply(z,c(1,2),function(z){round(mean(abs(z)<1.96,na.rm=T),3)})
  rownames(crg) = paste('cell',1:K)
  colnames(crg) = paste('bulk',1:nb)
  crg
}
