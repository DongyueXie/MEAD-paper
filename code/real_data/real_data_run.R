
###### run real data analysis ##########

source('real_data_func.R')
indis_ref = readRDS('~/MEAD-paper/data/real_data/indis_ref_12400by6by97.rds')


n_rep = 100
dirichlet.scales = c(5,20)
n_bulks = c(86,500,1000)
n_refs = 11
cases = c('null','all_diff')
cor_fdrs = c(0.3)
for(n_bulk in n_bulks){
  for(case in cases){
    for(aa in dirichlet.scales){
      for(n_ref in n_refs){



        ref_idx_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_ref_idx.rds',sep=''))
        bulk_p_array = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_p.rds',sep=''))
        groups_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_groups_mat.rds',sep=''))
        bulk_idx_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_idx_mat.rds',sep=''))

        for(cor_fdr in cor_fdrs){
          print(paste('Running:',case,'dirichlet.scale=',aa,"n_bulk=",n_bulk,'n_ref=',n_ref,'cor_fdr',cor_fdr))


          R01 = readRDS("data/neuron/R01_neuron_G12400_alpha03.rds")
          R01.dist = as.dist(1-R01)
          clusters01 = cluster::pam(R01.dist,10,pamonce=5)
          folds01 = clusters01$clustering
          table(folds01)
          rm(R01.dist)


          neuron_simu_study(indis_ref,
                            ref_idx_mat,
                            bulk_idx_mat,
                            bulk_p_array,
                            groups_mat,
                            case,
                            R01 = R01,
                            folds = folds01,
                            cor_fdr = cor_fdr,
                            add_bulk_bias=TRUE,
                            verbose = TRUE,
                            method_list = c('my'))

        }

      }
    }
  }
}

########## analyze the result #############

###############################################################################
### 1. Generate RMSE HEATMAP plot, left panel
###############################################################################
library(reticulate)
np <- import("numpy")
celltypes = c("DA" ,"Epen1","Sert","FPP","P_FPP","U_Neur")
celltypes_sieve = c('DA', 'Epen1', 'FPP', 'P_FPP', 'Sert', 'U_Neur')

n_ref=11
K = 6
n_rep = 100
n_bulk = 86

### change 5 to 20 for the plot in appendix
dirichlet.scales = c(5)
cases = c("null")
rmses = c()
method_list = c('MEAD','MuSiC','CIBERSORT','RNA-Sieve','OLS')


for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa,"n_bulk=",n_bulk,'n_ref=',n_ref))




    meth.order=c()

    # my fit
    if('MEAD'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'MEAD')
    }


    # my_unweighted fit
    if('my_unweighted'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_unweighted_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'my_unweighted')
    }

    if('my.QP'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_QP_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'my.QP')
    }


    # music
    if('MuSiC'%in%method_list){
      out = readRDS(paste("output/manuscript/real/music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
        p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_music,p_true_array))
      meth.order = c(meth.order,'MuSiC')
    }



    # cibersort
    if('CIBERSORT'%in%method_list){
      out = readRDS(paste("output/manuscript/real/cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
        p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_cibersort,p_true_array))
      meth.order = c(meth.order,'CIBERSORT')
    }



    # rna-sieve
    if('RNA-Sieve'%in%method_list){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_sieve,p_true_array))
      meth.order = c(meth.order,'RNA-Sieve')
    }


    # ols
    if('OLS'%in%method_list){
      out = readRDS(paste("output/manuscript/real/ols/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_ols = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_ols[,,r] = out$ols_fit[[r]]$p_hat
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_ols,p_true_array))
      meth.order = c(meth.order,'OLS')


    }


  }
}

rownames(rmses) = celltypes
colnames(rmses) = meth.order

library(ggplot2)
library(reshape2)

plot.meth.order = c('OLS','MuSiC','CIBERSORT','RNA-Sieve','MEAD')
df = melt(rmses[,match(meth.order,plot.meth.order)])
colnames(df) = c('Cell.type', 'Method','RMSE')
plot1 = ggplot(df, aes(x = Method, y = Cell.type, fill = RMSE)) +
  scale_fill_gradient2(
    low = 'steelblue', high = "red", mid = 'white', midpoint  = quantile(df$RMSE,0.4), limit = range(df$RMSE))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 11),axis.text.y = element_text(size = 11),axis.title.x = element_blank(),axis.title.y = element_blank())


###############################################################################
########## two group testing ####################
#########  two group difference figure and TABLE
###############################################################################



n_bulk = 86
n_rep = 100
dirichlet.scales = c(5)
cases = c("all_diff")

add_rnasieve = T

res=c()
res_mat = c()
res2 = c()
for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa))


    # my fit
    out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))
    p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
    p_true_array = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
      p_true_array[,,r] = out$rep_info[[r]]$b
    }
    diff_hat_my = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_my[i,] = rowMeans(p_hat_array_my[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_array_my[,out$rep_info[[i]]$groups==2,i])
    }

    # music
    out = readRDS(paste("output/manuscript/real/music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
    p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
      p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
    }
    diff_hat_music = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_music[i,] = rowMeans(p_hat_array_music[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_array_music[,out$rep_info[[i]]$groups==2,i])
    }

    # cibersort
    out = readRDS(paste("output/manuscript/real/cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
    p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
      p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
    }
    diff_hat_cibersort = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_cibersort[i,] = rowMeans(p_hat_array_cibersort[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_array_cibersort[,out$rep_info[[i]]$groups==2,i])
    }

    # rna-sieve
    if(add_rnasieve){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
      }

      diff_hat_sieve = matrix(nrow=n_rep,ncol=K)
      for(i in 1:n_rep){
        diff_hat_sieve[i,] = rowMeans(p_hat_sieve[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_sieve[,out$rep_info[[i]]$groups==2,i])
      }
    }else{
      diff_hat_sieve = NA
    }


    diff_hat_true_p = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_true_p[i,] = rowMeans(p_true_array[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_true_array[,out$rep_info[[i]]$groups==2,i])
    }


    dif = p1-p2
    cover.naive0 = matrix(nrow = n_rep,ncol=K)
    sd.naive0 = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_true_array[k,out$rep_info[[i]]$groups==1,i],p_true_array[k,out$rep_info[[i]]$groups==2,i])
        sd.naive0[i,k] = temp$stderr
        cover.naive0[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }
    mean(cover.naive0)

    cover.naive = matrix(nrow = n_rep,ncol=K)
    sd.naive = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_hat_array_my[k,out$rep_info[[i]]$groups==1,i],p_hat_array_my[k,out$rep_info[[i]]$groups==2,i])
        sd.naive[i,k] = temp$stderr
        cover.naive[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.music = matrix(nrow = n_rep,ncol=K)
    sd.naive.music = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_hat_array_music[k,out$rep_info[[i]]$groups==1,i],p_hat_array_music[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.music[i,k] = temp$stderr
        cover.naive.music[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.cibersort = matrix(nrow = n_rep,ncol=K)
    sd.naive.cibersort = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_hat_array_cibersort[k,out$rep_info[[i]]$groups==1,i],p_hat_array_cibersort[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.cibersort[i,k] = temp$stderr
        cover.naive.cibersort[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    if(add_rnasieve){
      cover.naive.rnasieve = matrix(nrow = n_rep,ncol=K)
      sd.naive.rnasieve = matrix(nrow = n_rep,ncol=K)
      for(i in 1:n_rep){
        for(k in 1:K){
          temp = t.test(p_hat_sieve[k,out$rep_info[[i]]$groups==1,i],p_hat_sieve[k,out$rep_info[[i]]$groups==2,i])
          sd.naive.rnasieve[i,k] = temp$stderr
          cover.naive.rnasieve[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
        }
      }
    }else{
      cover.naive.rnasieve = NA
    }


    diff_true = matrix(dif,nrow=n_rep,ncol=K,byrow = T)


    diff_res = cbind(sqrt(colMeans(diff_hat_true_p - diff_true)^2),
                     sqrt(colMeans(diff_hat_music- diff_true)^2),
                     sqrt(colMeans(diff_hat_cibersort - diff_true)^2),
                     sqrt(colMeans(diff_hat_sieve - diff_true)^2),
                     sqrt(colMeans(diff_hat_my - diff_true)^2))

    rownames(diff_res) = celltypes
    colnames(diff_res) = c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD')



    cc = c(mean(cover.naive0,na.rm=TRUE),
           mean(cover.naive.music,na.rm=TRUE),
           mean(cover.naive.cibersort,na.rm=TRUE),
           mean(cover.naive.rnasieve,na.rm=TRUE),
           mean(cover.naive,na.rm=TRUE))

    res_mat = cbind(colMeans(cover.naive0),
                    colMeans(cover.naive.music),
                    colMeans(cover.naive.cibersort),
                    colMeans(cover.naive.rnasieve),
                    colMeans(cover.naive))

    res = cbind(res,cc)


  }
}

rownames(res)=(c('t+truep','t+music','t+cibersort','t+rnasieve','t+mea.err+weight'))
rownames(res_mat) = celltypes
colnames(res_mat) = c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD')

datax = data.frame(Coverage = c(res_mat),Cell.type = rep(celltypes,5),Method = rep(c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD'),each=6))
datax$Method = factor(datax$Method,levels = unique(datax$Method))
ggplot(datax, aes(x=Method, y=Coverage, group=Cell.type)) +
  geom_hline(yintercept = 0.95,linetype = "dashed")+geom_point(aes(shape=Cell.type,color = Cell.type))+
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11))
ggplot(datax, aes(x = Method, y = Coverage, fill = Cell.type)) +
  geom_bar(stat = "identity", position = "dodge") +  # Create grouped bars
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +  # Add a dashed horizontal line at Coverage = 0.95
  labs(x = "Method", y = "Coverage", fill = "Cell Type") +  # Set axis labels and legend title
  theme_minimal() +  # Use a minimal theme for a cleaner look
  scale_fill_brewer(palette = "PiYG") +
  coord_cartesian(ylim = c(min(datax$Coverage),1))


###############################################################################
##### Generate RMSE HEATMAP plot about group difference, right panel
###############################################################################

datax = data.frame(rmse = c(diff_res),cell.type = rep(celltypes,5),method = rep(c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD'),each=6))

datax = data.frame(RMSE = c(diff_res),Cell.type = rep(celltypes,5),Method = rep(c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD'),each=6))
datax$Method = factor(datax$Method,levels = unique(datax$Method))


plot2 =ggplot(datax, aes(x = Method, y = Cell.type, fill = RMSE)) +scale_fill_gradient2(
  low = 'steelblue', high = "red", mid = 'white', midpoint  = quantile(datax$RMSE,0.9), limit = range(datax$RMSE))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11),axis.text.y = element_text(size = 11),axis.title.x = element_blank(),axis.title.y = element_blank())

library(gridExtra)
grid.arrange(plot1, plot2, ncol=2)


###############################################################################
###################### Generate coverage table ################################
###############################################################################


### coverage of p

dirichlet.scales = c(5)
cases = c("null")
method_list = c('my','rnasieve')

res=c()
res2 = c()
for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa))

    meth.order=c()
    if('my'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr03_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      #get_rmse_array(p_hat_array_my,p_true_array)

      res = rbind(res, rowMeans(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array),1,sd))
      meth.order = c(meth.order,'my')
    }


    if('my_unweighted'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_unweighted_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      #get_rmse_array(p_hat_array_my,p_true_array)

      res = rbind(res, rowMeans(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array),1,sd))
      meth.order = c(meth.order,'my_unweighted')
    }

    if('ols'%in%method_list){
      out = readRDS(paste("output/manuscript/real/ols/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_ols = array(dim = c(K,n_bulk,n_rep))
      p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_ols[,,r] = out$ols_fit[[r]]$p_hat
        p_hat_array_se[,,r] = out$ols_fit[[r]]$p_hat_se
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      #get_rmse_array(p_hat_array_my,p_true_array)

      res = rbind(res, rowMeans(get_coverage_p(p_hat_array_ols,p_hat_array_se,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_array_ols,p_hat_array_se,p_true_array),1,sd))
      meth.order = c(meth.order,'ols')
    }

    if('rnasieve'%in%method_list){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      p_hat_l_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat_ci_l.npy',sep=''))
      p_hat_r_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat_ci_r.npy',sep=''))

      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
        p_hat_l_sieve[,,r] = p_hat_l_sieve[match(celltypes,celltypes_sieve),,r]
        p_hat_r_sieve[,,r] = p_hat_r_sieve[match(celltypes,celltypes_sieve),,r]
      }

      res = rbind(res, rowMeans(get_coverage_p_sieve(p_hat_sieve,p_hat_l_sieve,p_hat_r_sieve,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p_sieve(p_hat_sieve,p_hat_l_sieve,p_hat_r_sieve,p_true_array),1,sd))
      meth.order = c(meth.order,'rnasieve')
    }


  }
}

colnames(res) = celltypes
colnames(res2) = celltypes
rownames(res) = meth.order
rownames(res2) = meth.order
res
res2


################################
######## plot for CI length#####

get_coverage_for_one_rep = function(p_hat,p_hat_se,true_p,alpha = 0.05){
  lower_array = p_hat - qnorm(1-alpha/2)*p_hat_se
  lower_array = ifelse(lower_array>0,lower_array,0)

  upper_array = p_hat + qnorm(1-alpha/2)*p_hat_se
  upper_array = ifelse(upper_array>1,1,upper_array)

  return(list(lower=lower_array,upper = upper_array,true_p = true_p))
}
res <- readRDS("output/manuscript/real/my/add_bulk_bias/neuron_ref11_rep100_bulk86_dirichlet5_corfdr03_null.rds")
iter=9
lu = get_coverage_for_one_rep(res$my_fit[[iter]]$p_hat,res$my_fit[[iter]]$p_hat_se,res$rep_info[[iter]]$b)

library(ggplot2)
library(gridExtra)


plot_list = list()
celltypes = c("DA" ,"Epen1","Sert","FPP","P_FPP","U_Neur")
for(k in 1:6){

  datax = data.frame(lower = lu$lower[k,],
                     upper = lu$upper[k,],
                     true_p = lu$true_p[k,])
  datax = datax[order(datax$true_p),]
  datax$indi = 1:86

  plot_list[[k]] = ggplot(datax,aes(x = indi,y=true_p))+
    geom_line()+
    geom_point(size = 1,color = 'red')+
    xlab("Sorted bulk samples") +
    ylab("p")+
    geom_ribbon(data=datax,aes(ymin=lower,ymax=upper),alpha=0.3)+
    ggtitle(paste('Cell type:', celltypes[k]))+
    ylim(0,1)
}

grid.arrange(grobs=plot_list,ncol=2)



