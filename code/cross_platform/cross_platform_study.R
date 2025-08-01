########### Create bulk data ################
xin_raw <- readRDS("data/real_data_platform/xin_raw.rds")
# create bulk data
indis = unique(xin_raw$individual)
bulk = c()
for(ind in indis){
  print(ind)
  idx = which(xin_raw$individual == ind)
  bulk = cbind(bulk, rowSums(xin_raw@assays$data$counts[,idx]))
}
colnames(bulk) = indis

true_p = table(xin_raw$individual,xin_raw$cell_type)
true_p = true_p / rowSums(true_p)

##########Fit MuSiC###############
library(MuSiC)
XinT2D.sce = readRDS('data/real_data_platform/XinT2Dsce.rds')
dim(XinT2D.sce)
XinT2D.construct.full = bulk_construct(XinT2D.sce, clusters = 'cellType', samples = 'SubjectName')
XinT2D.construct.full$prop.real = relative.ab(XinT2D.construct.full$num.real, by.col = FALSE)
XinT2D.construct.full$prop.real = XinT2D.construct.full$prop.real[,c(2,1,3,4)]

EMTAB.sce <- readRDS("data/real_data_platform/EMTABsce_healthy.rds")
Est.prop.Xin = music_prop(bulk.mtx = XinT2D.construct.full$bulk.counts, sc.sce = EMTAB.sce,
                          clusters = 'cellType', samples = 'sampleID',
                          select.ct = c('alpha', 'beta', 'delta', 'gamma'))

Eval_multi(prop.real = data.matrix(XinT2D.construct.full$prop.real),
           prop.est = list(data.matrix(Est.prop.Xin$Est.prop.weighted),
                           data.matrix(Est.prop.Xin$Est.prop.allgene)),
           method.name = c('MuSiC', 'NNLS'))

W = Est.prop.Xin$Weight.gene

########### Create reference data ################
# only use healthy individuals
library(MEAD)
segerstolpe_raw <- readRDS("data/real_data_platform/segerstolpe_raw.rds")
seger_healthy = segerstolpe_raw[,segerstolpe_raw$disease=="normal"]
datax = MEAD_preprocessing(bulk,seger_healthy,cell_types = c("alpha","beta", "delta", "gamma"),filter.gene = F,marker_gene = rownames(W))
refs = MEAD_getX(datax$ref,datax$cell_types,datax$individuals)
# total numer of genes used
sum(rowSums(refs$X)!=0)

# get correlated genes
R01_01 = get_R01(datax$bulk,0.1)
R01_03 = get_R01(datax$bulk,0.3)
R01_05 = get_R01(datax$bulk,0.5)

# fitted_default = MEAD_est(datax$bulk,refs$X,refs$V,w=refs$w)
fitted_default = MEAD_est(datax$bulk,refs$X,refs$V,w=refs$w,R01=R01_05)

########## fit cibersort############
source('code/CIBERSORT.R')
# ref_samples = c()
# for(k in 1:K){
#   temp = X_array_ref[,k,]
#   colnames(temp) = rep(colnames(ref)[k],n_ref)
#   ref_samples = cbind(ref_samples,temp)
# }
# rownames(ref_samples) = gene_names

ref_samples = segerstolpe_raw[,segerstolpe_raw$cell_type1%in%c("alpha","beta", "delta", "gamma")]
ref_mat = counts(ref_samples)
colnames(ref_mat) = ref_samples$cell_type1
sig_mat = build_signature_matrix_CIBERSORT(ref_mat)
# sig_mat = sig_mat[,match(celltypes,colnames(sig_mat))]
fit_cibersort = CIBERSORT(sig_mat,data.frame(GeneSymbol = rownames(bulk),bulk))
temp_row_names = rownames(fit_cibersort)
for(i in 1:length(temp_row_names)){
  temp_row_names[i] = gsub("\\.", " ", temp_row_names[i])
}
rownames(fit_cibersort) = temp_row_names
fit_cibersort = fit_cibersort[,match(colnames(true_p),colnames(fit_cibersort))]
fit_cibersort = fit_cibersort[match(rownames(true_p),rownames(fit_cibersort)),]

########## heat plot of RMSE############
library(ggplot2)
library(reshape2)

indi_names = rownames(true_p)

# Find the global min and max p values
global_min_p <- 0
global_max_p <- 1

df = melt(true_p)
colnames(df) = c("Var1",  "Var2",  "p")
# Create heatmap plot with ggplot2
p1 = ggplot(df, aes(Var1, Var2, fill = p)) +
  geom_tile() +
  scale_fill_gradient2(low = 'steelblue' , high = "red", mid = "white", midpoint  = quantile(df$p,0.7),
                       limits = c(global_min_p, global_max_p)) +  # Custom color gradient
  #labs(x = "Cell types", y = "Individuals") +
  labs(x = "Individuals", y = "Cell types") +
  ggtitle("True cell type proportions") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())


df = melt(t(fitted_default$p_hat)[match(indi_names,rownames(t(fitted_default$p_hat))),])
colnames(df) = c("Var1",  "Var2",  "p")
# Create heatmap plot with ggplot2
p3 = ggplot(df, aes(Var1, Var2, fill = p)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", high = "red", mid = 'white', midpoint  = quantile(df$p,0.7),
                       limits = c(global_min_p, global_max_p)) +  # Custom color gradient
  #labs(x = "Cell types", y = "Individuals") +
  labs(x = "Individuals", y = "Cell types") +
  ggtitle("MEAD") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df = melt(Est.prop.Xin$Est.prop.weighted[match(indi_names,rownames(Est.prop.Xin$Est.prop.weighted)),])
colnames(df) = c("Var1",  "Var2",  "p")
# Create heatmap plot with ggplot2
p4 = ggplot(df, aes(Var1, Var2, fill = p)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", high = "red", mid = 'white', midpoint  = quantile(df$p,0.7),
                       limits = c(global_min_p, global_max_p)) +  # Custom color gradient
  #labs(x = "Cell types", y = "Individuals") +
  labs(x = "Individuals", y = "Cell types") +
  ggtitle("MuSiC") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

df = melt(fit_cibersort[match(indi_names,rownames(fit_cibersort)),])
colnames(df) = c("Var1",  "Var2",  "p")
# Create heatmap plot with ggplot2
p5 = ggplot(df, aes(Var1, Var2, fill = p)) +
  geom_tile() +
  scale_fill_gradient2(low = "steelblue", high = "red", mid = 'white', midpoint  = quantile(df$p,0.7),
                       limits = c(global_min_p, global_max_p)) +  # Custom color gradient
  #labs(x = "Cell types", y = "Individuals") +
  labs(x = "Individuals", y = "Cell types") +
  ggtitle("CIBERSORT") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

gridExtra::grid.arrange(p1,p5,p4,p3,nrow=4, heights = c(1, 1, 1,1.3))

#############rmse, mad, R#################
sqrt(mean((t(fitted_default$p_hat)-data.matrix(XinT2D.construct.full$prop.real))^2))
sqrt(mean((fit_cibersort-data.matrix(XinT2D.construct.full$prop.real))^2))
sqrt(mean((data.matrix(Est.prop.Xin$Est.prop.weighted)-data.matrix(XinT2D.construct.full$prop.real))^2))
sqrt(mean((data.matrix(Est.prop.Xin$Est.prop.allgene)-data.matrix(XinT2D.construct.full$prop.real))^2))

(mean(abs(t(fitted_default$p_hat)-data.matrix(XinT2D.construct.full$prop.real))))
(mean(abs(fit_cibersort-data.matrix(XinT2D.construct.full$prop.real))))
(mean(abs(data.matrix(Est.prop.Xin$Est.prop.weighted)-data.matrix(XinT2D.construct.full$prop.real))))
(mean(abs(data.matrix(Est.prop.Xin$Est.prop.allgene)-data.matrix(XinT2D.construct.full$prop.real))))


################two group difference rmse ################
true_diff = colMeans(XinT2D.construct.full$prop.real[1:12,]) - colMeans(XinT2D.construct.full$prop.real[13:18,])
mead_diff = rowMeans(fitted_default$p_hat[,1:12]) - rowMeans(fitted_default$p_hat[,13:18])
music_diff = colMeans(data.matrix(Est.prop.Xin$Est.prop.weighted)[1:12,]) - colMeans(data.matrix(Est.prop.Xin$Est.prop.weighted)[13:18,])
nnls_diff = colMeans(data.matrix(Est.prop.Xin$Est.prop.allgene)[1:12,]) - colMeans(data.matrix(Est.prop.Xin$Est.prop.allgene)[13:18,])
cibersort_diff = colMeans(fit_cibersort[1:12,]) - colMeans(fit_cibersort[13:18,])

sqrt(mean(()^2))
sqrt(mean(()^2))
sqrt(mean(()^2))
############## coverage ######################
get_coverage_for_one_rep = function(p_hat,p_hat_se,true_p,alpha = 0.05){
  lower_array = p_hat - qnorm(1-alpha/2)*p_hat_se
  lower_array = ifelse(lower_array>0,lower_array,0)

  upper_array = p_hat + qnorm(1-alpha/2)*p_hat_se
  upper_array = ifelse(upper_array>1,1,upper_array)

  return(list(lower=lower_array,upper = upper_array,true_p = true_p))
}

indi_names = rownames(true_p)

phat = t(fitted_default$p_hat)[match(indi_names,rownames(t(fitted_default$p_hat))),]
phat_se = t(fitted_default$p_hat_se)[match(indi_names,rownames(t(fitted_default$p_hat_se))),]
lu = get_coverage_for_one_rep(phat,phat_se,true_p)

covered = ((lu$lower <= lu$true_p) & (lu$upper >= lu$true_p))
covered
mean(covered)

plot_list = list()
celltypes = c('alpha', 'beta', 'delta', 'gamma')
for(k in 1:4){

  datax_temp = data.frame(lower = lu$lower[,k],
                     upper = lu$upper[,k],
                     true_p = lu$true_p[,k])
  datax_temp = datax_temp[order(datax_temp$true_p),]
  datax_temp$indi = 1:18

  plot_list[[k]] = ggplot(datax_temp,aes(x = indi,y=true_p))+
    geom_line()+
    geom_point(size = 1,color = 'red')+
    xlab("Sorted bulk samples") +
    ylab("p")+
    geom_ribbon(data=datax_temp,aes(ymin=lower,ymax=upper),alpha=0.3)+
    ggtitle(paste('Cell type:', celltypes[k]))+
    ylim(0,1)
}

gridExtra::grid.arrange(grobs=plot_list,ncol=2)


################################
###########softplus#############
################################


fitted_default_softplus = MEAD_est(datax$bulk,refs$X,refs$V,w=refs$w,R01 = R01_05,beta_tilde_transformation = list(method='softplus',a=10))
indi_names = rownames(true_p)
phat = t(fitted_default_softplus$p_hat)[match(indi_names,rownames(t(fitted_default_softplus$p_hat))),]
phat_se = t(fitted_default_softplus$p_hat_se)[match(indi_names,rownames(t(fitted_default_softplus$p_hat_se))),]
lu = get_coverage_for_one_rep(phat,phat_se,true_p)

covered = ((lu$lower <= lu$true_p) & (lu$upper >= lu$true_p))
covered
mean(covered)

plot_list = list()
celltypes = c('alpha', 'beta', 'delta', 'gamma')
for(k in 1:4){

  datax_temp = data.frame(lower = lu$lower[,k],
                     upper = lu$upper[,k],
                     true_p = lu$true_p[,k])
  datax_temp = datax_temp[order(datax_temp$true_p),]
  datax_temp$indi = 1:18

  plot_list[[k]] = ggplot(datax_temp,aes(x = indi,y=true_p))+
    geom_line()+
    geom_point(size = 1,color = 'red')+
    xlab("Sorted bulk samples") +
    ylab("p")+
    geom_ribbon(data=datax_temp,aes(ymin=lower,ymax=upper),alpha=0.3)+
    ggtitle(paste('Cell type:', celltypes[k]))+
    ylim(0,1)
}

gridExtra::grid.arrange(grobs=plot_list,ncol=2)








