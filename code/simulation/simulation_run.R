
################################
################################
############ run simulation ####

source('simulation_func.R')

xin = readRDS('~/MEAD-paper/data/simulation/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
set.seed(12345)
simu_out = simu_study(ref,sigma2,p=c(0.5,0.3,0.1,0.1),n_bulk = 50,dirichlet.scale=10,printevery = 1,gamma_shape = 5)


################################
# the RMSE of six different methods
for(m in 1:6){
  print(paste(names(simu_out)[m],':',get_rmse_array(simu_out[[m]],simu_out$true_p)))
}

# coverage of different methods
print(paste(names(simu_out)[7],':',mean(get_coverage_p(simu_out$p_hat_ols,simu_out[[7]],simu_out$true_p))))
print(paste(names(simu_out)[8],':',mean(get_coverage_p(simu_out$p_hat,simu_out[[8]],simu_out$true_p))))
for(m in 9:17){
  print(paste(names(simu_out)[m],':',mean(get_coverage_p(simu_out$p_hat_weight,simu_out[[m]],simu_out$true_p))))
}

# sd of the coverage of different methods
print(paste(names(simu_out)[7],':',sd(get_coverage_p(simu_out$p_hat_ols,simu_out[[7]],simu_out$true_p))))
print(paste(names(simu_out)[8],':',sd(get_coverage_p(simu_out$p_hat,simu_out[[8]],simu_out$true_p))))
for(m in 9:17){
  print(paste(names(simu_out)[m],':',sd(get_coverage_p(simu_out$p_hat_weight,simu_out[[m]],simu_out$true_p))))
}

################################
######## plot for CI length#####

get_coverage_for_one_rep = function(p_hat,p_hat_se,true_p,alpha = 0.05){
  lower_array = p_hat - qnorm(1-alpha/2)*p_hat_se
  lower_array = ifelse(lower_array>0,lower_array,0)

  upper_array = p_hat + qnorm(1-alpha/2)*p_hat_se
  upper_array = ifelse(upper_array>1,1,upper_array)

  return(list(lower=lower_array,upper = upper_array,true_p = true_p))
}
lu = get_coverage_for_one_rep(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cv_cor01,simu_out$true_p)

library(ggplot2)
library(gridExtra)


plot_list = list()
iter=1
celltypes = c('alpha', 'beta', 'delta', 'gamma')
for(k in 1:4){

  datax = data.frame(lower = lu$lower[k,,iter],
                     upper = lu$upper[k,,iter],
                     true_p = lu$true_p[k,,iter])
  datax = datax[order(datax$true_p),]
  datax$indi = 1:50

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
