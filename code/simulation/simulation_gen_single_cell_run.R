
################################
################################
############ run simulation ####

source('code/simulation/simulation_gen_single_cell_func.R')

xin = readRDS('data/simulation/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
set.seed(12345)
simu_out = simu_study(ref,sigma2,p=c(0.5,0.3,0.1,0.1),n_bulk = 50,dirichlet.scale=10,printevery = 1, 
                      bulk_nb_theta = 10,
                      cross_plat_bias_gamma_var = 0.3)
saveRDS(simu_out,file='simu_res_gen_single_try0.rds')

################################
# the RMSE of six different methods
for(m in 1:9){
  print(paste(names(simu_out)[m],':',get_rmse_array(simu_out[[m]],simu_out$true_p)))
}

# coverage of different methods
print(paste(names(simu_out)[10],':',mean(get_coverage_p(simu_out$p_hat_ols,simu_out[[10]],simu_out$true_p))))
print(paste(names(simu_out)[11],':',mean(get_coverage_p(simu_out$p_hat,simu_out[[11]],simu_out$true_p))))
for(m in 12:20){
  print(paste(names(simu_out)[m],':',mean(get_coverage_p(simu_out$p_hat_weight,simu_out[[m]],simu_out$true_p))))
}
for(m in 21:30){
  print(paste(names(simu_out)[m],':',mean(get_coverage_p(simu_out$p_hat_weight,simu_out[[m]],simu_out$true_p))))
}

# sd of the coverage of different methods
print(paste(names(simu_out)[10],':',sd(get_coverage_p(simu_out$p_hat_ols,simu_out[[10]],simu_out$true_p))))
print(paste(names(simu_out)[11],':',sd(get_coverage_p(simu_out$p_hat,simu_out[[11]],simu_out$true_p))))
for(m in 12:20){
  print(paste(names(simu_out)[m],':',sd(get_coverage_p(simu_out$p_hat_weight,simu_out[[m]],simu_out$true_p))))
}
for(m in 21:30){
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
lu = get_coverage_for_one_rep(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cor05,simu_out$true_p)

library(ggplot2)
library(gridExtra)
library(grid)

plot_list = list()
iter = 1
celltypes = c('alpha', 'beta', 'delta', 'gamma')

for (k in 1:4) {
  datax = data.frame(lower = lu$lower[k,,iter],
                     upper = lu$upper[k,,iter],
                     true_p = lu$true_p[k,,iter])
  datax = datax[order(datax$true_p), ]
  datax$indi = 1:50
  
  plot_list[[k]] = ggplot(datax, aes(x = indi, y = true_p)) +
    geom_line() +
    geom_point(size = 1, color = 'red') +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    ylab("p") +
    ggtitle(paste('Cell type:', celltypes[k])) +
    ylim(0, 1) +
    theme(axis.title.x = element_blank())  # Remove individual x-axis labels
}

# Shared x-axis label
x.grob <- textGrob("Sorted bulk samples", gp = gpar(fontsize = 14))

# Arrange plots with shared x-axis label
grid.arrange(grobs = plot_list, ncol = 4, bottom = x.grob)

######################################

####### softplus CI ##################

get_coverage_for_one_rep = function(p_hat,p_hat_se,true_p,alpha = 0.05){
  lower_array = p_hat - qnorm(1-alpha/2)*p_hat_se
  lower_array = ifelse(lower_array>0,lower_array,0)
  
  upper_array = p_hat + qnorm(1-alpha/2)*p_hat_se
  upper_array = ifelse(upper_array>1,1,upper_array)
  
  return(list(lower=lower_array,upper = upper_array,true_p = true_p))
}
lu = get_coverage_for_one_rep(simu_out$p_hat_weight_sp,simu_out$p_hat_weight_se_cor05_sp,simu_out$true_p)

library(ggplot2)
library(gridExtra)
library(grid)

plot_list = list()
iter = 1
celltypes = c('alpha', 'beta', 'delta', 'gamma')

for (k in 1:4) {
  datax = data.frame(lower = lu$lower[k,,iter],
                     upper = lu$upper[k,,iter],
                     true_p = lu$true_p[k,,iter])
  datax = datax[order(datax$true_p), ]
  datax$indi = 1:50
  
  plot_list[[k]] = ggplot(datax, aes(x = indi, y = true_p)) +
    geom_line() +
    geom_point(size = 1, color = 'red') +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    ylab("p") +
    ggtitle(paste('Cell type:', celltypes[k])) +
    ylim(0, 1) +
    theme(axis.title.x = element_blank())  # Remove individual x-axis labels
}

# Shared x-axis label
x.grob <- textGrob("Sorted bulk samples", gp = gpar(fontsize = 14))

# Arrange plots with shared x-axis label
grid.arrange(grobs = plot_list, ncol = 4, bottom = x.grob)