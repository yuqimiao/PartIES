library(tidyverse)
library(igraph)
library(SNFtool)
source("R/PartIES.R")


# simulation for 4 clusters, data 1 separates 12/3/4, data 2 separates 1/2/34, data 3 separates 1/2/3/4 but with vague division ----
# output
# data matrix with n_sub*n_feat

get_data_4clust_dim1_numsig = function(n_sub = 200,
                                       n_feat = 1000,
                                       n_signal, # number of signals
                                       mu_vec = c(1, 2, 10, 20), # feature mean for 4 clusters
                                       signal_sd = 1,
                                       noise_sd = 1){
  data = matrix(0, n_sub, n_feat)
  signal_ind = 1:n_signal
  noise_ind = (n_signal+1):n_feat


  n_sub_perc = n_sub/4
  cluster_index = lapply(1:4, function(i) ((i-1)*n_sub_perc+1):(i*n_sub_perc))

  # noise simulation
  data[1:n_sub, noise_ind]=rnorm(n = n_sub*length(noise_ind), mean = 0, sd = noise_sd)


  # signal_simulation
  for(i in 1:4){
    cluster_index_cur = cluster_index[[i]]
    mu_cur = mu_vec[i]
    data[cluster_index_cur, signal_ind]=rnorm(n = n_sub_perc*length(signal_ind), mean = mu_cur, sd = signal_sd)
  }

  return(data)
}


# heatmap_gg(data, "sim1") # heatmap_gg in visualization_functions.R
# simulation par collection

# parameters ----
# simulation parameters: 4 data types gives different cluster structure
mu1 = c(1,1,3,3) # data 1 mean
mu2 = c(1,2,3,3)# data 2 mean
mu3 = c(3,3,1,2) # data 3 mean
feature_comb_ls = list(c(10000, 10000, 10000))
neig_single_vec = c(2,3,3)
sigma = 2 # kernel_parameter
n = 200
k = 50
c = 4

# generate data ----
set.seed(123)
n_signal = 100
n_feat_ls = list(n_feat1 = 10000,
                 n_feat2 = 10000,
                 n_feat3 = 10000)

mu_ls = list(mu1 = mu1, mu2 = mu2, mu3 = mu3)


data_list = list()
for(i in 1:3){
  if(!is.na(n_feat_ls[[i]])){
    print(mu_ls[[i]])
    data = get_data_4clust_dim1_numsig(n_feat = n_feat_ls[[i]],
                                       mu_vec = mu_ls[[i]],
                                       n_signal = n_signal)
    data_list = c(data_list, list(data))
  }
}

distance_list = lapply(data_list, function(x) dist2(x,x))
kernel_list = lapply(distance_list, function(x) kernel_calculation(distance = x, k = k, sigma = sigma ))
diff_kernel_list = lapply(kernel_list, function(x) diffusion_enhancement(kernel = x, alpha = 1, k = k , diffusion_form = "L1"))

# estimate the number of clusters for each data type using eigengap
nc_single_est = map_dbl(diff_kernel_list, function(s) estimateNumberOfClustersGivenGraph(s)[[1]])

# Perform PartIES ----
res_part_cimlr = parties(diff_kernel_list,
                            k = k,
                            neig_single = nc_single_est,
                            c = c, n_iter = 50)

igraph::compare(res_part_cimlr$cluster, rep(1:4, each = 50), "nmi")



