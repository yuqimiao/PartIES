# Dependency
library(SNFtool)
dyn.load("R/projsplx_R.so")
## parties
# Input:
# * kernel_list: kernels/affinity matrix generated from kernel_list_generation function, or any other similarity matrix, symmetric, semi-pos-def
# * k: numbers of neighbors in KNN, indicating the structure of single data
# * c: number of clusters estimated for integrated similarity matrix
# * neig_single: a s vector containing the number of eigenvec to use in each data type
# * update_neig: whether we update c based on the new weighted sum of the partition information
# * num_eig: candidate number of eigen vec to use in single data partition information F_s
# Output:
# * cluster: cluster result for each subject
# * Z_ls: similarity matrices list for data types
# * F_ls: partition information matrices list for data types
# * Y: integrated partition information
# * w: empirical weight for each data type

# Notes:
# lambda in the code is the gamma, beta in overleaf
# generate kernel list before optimization
parties = function(kernel_list,
                   k = 10,
                   c,
                   neig_single,
                   n_iter = 30,
                   thresh=1e-10
){
  # kernel distance calculation ----
  S = length(kernel_list) # number of data types
  n = nrow(kernel_list[[1]]) # number of samples
  kernel_distance_list = lapply(kernel_list, function(kernel) dist_kernels(kernel))

  # initialization ----
  # par initial
  initial_list = lapply(kernel_distance_list, function(distance) initial(distance, k = k))
  Z0 = lapply(kernel_distance_list, function(x){
    x = as.matrix(x)
    max(x)-x
  })

  F0 = lapply(1:S, function(i){
    L = diag(rowSums(Z0[[i]]))-Z0[[i]]
    F_eig1 = eig1(L, c = neig_single[i], isMax = 0)$eigvec
    # F_eig1 = dn.cimlr(F_eig1, "ave")
    F_eig1
  })


  Z_cur = Z0
  F_cur = F0
  w_cur = rep(1/S,S)
  weight_L = matrix(0,n,n)
  for(s in 1:S){
    L = diag(1,n)-F_cur[[s]]%*%t(F_cur[[s]])*2
    weight_L = weight_L+w_cur[[s]]*L
  }
  Y0 = eig1(weight_L, c=c, isMax=0)$eigvec
  Y_cur = Y0
  # optimization ----
  converge = 100
  converge_cur = converge
  rr = Reduce("+", lapply(initial_list, function(x) x$lambda))/S
  for(i in 1:n_iter){
    # Update Z ----
    Z_pre = Z_cur
    for(s in 1:S){
      # updata each data type separately
      F_eig1 = F_cur[[s]]
      distX = initial_list[[s]]$distX
      distX1 = initial_list[[s]]$distX1
      idx = initial_list[[s]]$idx
      lambda = initial_list[[s]]$lambda
      r = initial_list[[s]]$lambda
      distf = L2_distance_1(t(F_eig1),t(F_eig1))
      A = array(0,c(n,n))
      b = idx[,2:dim(idx)[2]]
      a = apply(array(0,c(n,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:n })
      inda = cbind(as.vector(a),as.vector(b)) # rank of each row aligned
      ad = (distX[inda]+rr*distf[inda])/2/r
      dim(ad) = c(n,ncol(b))

      # call the c function for the optimization
      c_input = -t(ad)
      c_output = t(ad)
      ad = t(.Call("projsplx_R",c_input,c_output))
      A[inda] = as.vector(ad)
      A[is.nan(A)] = 0
      A = (A + t(A)) / 2
      # Z_cur[[s]] = (1 - eta) * Z_cur[[s]] + eta * A
      Z_cur[[s]] = A
      # if(network_diffusion){Z_cur[[s]] = network.diffusion(Z_cur[[s]], k)}
      # Z_cur[[s]] = dn(Z_cur[[s]],"ave")
    }

    # Update F ----
    if(i>1){
      F_pre = F_cur
    }else{
      F_cur = list()
    }
    for(s in 1:S){
      # updata each data type separately
      L = normalized_GL(Z_cur[[s]])
      Lz = normalized_GL(Y_cur %*% t(Y_cur))
      F_eig1 = eig1(L+w_cur[[s]]*Lz, isMax = 0, c = neig_single[[s]] )$eigvec
      F_cur[[s]] = F_eig1
    }


    # Update Y ----
    w_pre = w_cur
    weight_L = matrix(0,n,n)
    for(s in 1:S){
      L = diag(1,n)-F_cur[[s]]%*%t(F_cur[[s]])*2
      weight_L = weight_L+w_cur[[s]]*L
    }
    # update neig_all, number of eigenvec to use in integrated partition, if not update, the same as c
    # if(update_neig){
    #   neig_all = est_nclust(S = weight_L, num_eig = num_eig, is_GL = T)
    # }else{
    #   neig_all = c
    # }
    neig_all = c
    Y_cur = eig1(weight_L, isMax = 0, c = neig_all)$eigvec

    # Update w empirically
    w_pre = w_cur
    YYt = Y_cur %*% t(Y_cur)
    for(s in 1:S){
      FsFst = F_cur[[s]]%*%t(F_cur[[s]])
      w_cur[[s]] = 1/(2*sqrt(sum((YYt-FsFst)^2)))
    }
    w_cur = w_cur/sum(w_cur)

    # convergence judging ----
    converge_pre = converge_cur
    converge_cur = 0
    for(s in 1:S){
      distX = initial_list[[s]]$distX
      beta = initial_list[[s]]$lambda
      gamma = Reduce("+", lapply(initial_list, function(x) x$lambda))/S
      term1 = sum(distX * Z_cur[[s]])
      term2 = sum(Z_cur[[s]]^2) * beta
      term3 = sum(diag( t(F_cur[[s]])%*% (diag(apply(Z_cur[[s]],1,sum))-Z_cur[[s]]) %*%F_cur[[s]])) * gamma
      term4 = gamma* w_cur[[s]] * sum((Y_cur %*% t(Y_cur)-F_cur[[s]]%*%t(F_cur[[s]]))^2)
      converge_cur = converge_cur+term1+term2+term3+term4
    }
    converge = abs(converge_cur-converge_pre)
    if(converge<thresh){
      break
    }
    print(paste0("Iteration: ",i))
  }

  ## Output ----
  cluster = kmeans(Y_cur, c, nstart = 200)$cluster
  return(list(cluster = cluster,
              Z_ls = Z_cur,
              F_ls = F_cur,
              Y = Y_cur,
              w = w_cur))
}


# Utils ----

kernel_calculation = function(distance, k, sigma){
  distance = (distance)^(1/2)
  # sort
  d_sort = t(apply(distance,2,sort))
  # calculate variance matrix
  TT = apply(d_sort[,1:k],MARGIN=1,FUN=mean) + .Machine$double.eps
  ## calc mean distance for the first k neighbors of every subjects(row), length = k
  TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
  Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
  Sig[Sig<0] = 0
  Sig = (Sig + t(Sig))/2
  # calculate kernel value
  W = dnorm(distance,0,sigma*Sig)*sigma*Sig*sqrt(2*pi)
  # symmetry

  W = (W+t(W))/2

  return(W)
}


diffusion_enhancement = function(kernel,
                                 alpha,
                                 k,
                                 diffusion_form){
  "
  This function aims to enhance the local structure of the kernel
  input:
    * kernel,
    * alpha: importance of local structure
    * k: number of neighbors
  output:
    * kernel_enh: the enhanced kernel
  "
  # get the KNN normalized transition matrix P
  P = dominateset(kernel, k)
  P = dn.cimlr(P, "ave")
  # perform transformation
  if(diffusion_form == "L1") L = kernel %*% t(P)
  if(diffusion_form == "L2") L = kernel %*% t(P) %*% t(P)
  if(diffusion_form == "L3") L = P %*% kernel %*%  t(P)
  if(diffusion_form == "L4") L = kernel %*%  t(P) %*%  t(P) %*%  t(P)

  # combine with the initial global K
  L = (L+t(L))/2
  kernel_enh = (1-alpha)* kernel + alpha * L

  return(kernel_enh)
}

dominateset <- function(xx,KK=20) {
  ### This function outputs the top KK neighbors.
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
  }

  return(A)
}
# kernel normalization
"dn.cimlr" = function( w, type ) {

  # compute the sum of any column
  w = w * dim(w)[1]
  D = apply(abs(w),MARGIN=1,FUN=sum)
  D = 1 / D
  D_temp = matrix(0,nrow=length(D),ncol=length(D))
  D_temp[cbind(1:length(D),1:length(D))] = D
  D = D_temp

  # type "ave" returns D^-1*W
  if(type=="ave") {
    wn = D %*% w
  }else if(type == "sym"){
    wn = sqrt(D) %*% w %*% sqrt(D)
  }
  else {
    stop("Invalid type!")
  }

  return(wn)

}

dist_kernels = function(kernel){
  K = kernel
  k = 1/sqrt(diag(K)+1)
  ## diagnal of the matrix is the highest similarity
  G = K * (k %*% t(k))
  ## with diag(K) all ~ 0, k is just n 1s, thus G is K
  G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
  G2 = t(G1)
  ## Use the average difference btw 2 self similarities and pairwise similarity as kenels
  D_Kernels_tmp = (G1 + G2 - 2*G)/2
  D_Kernels = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
  D_Kernels = (as.matrix(D_Kernels)+t(as.matrix(D_Kernels)))/2
  return(as.matrix(D_Kernels))
}

initial = function(distance,k){
  n = nrow(distance)
  # rank every row of the kernel distance
  res = apply(distance,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
  distance_sort = array(0,c(nrow(distance),ncol(distance)))
  idx = array(0,c(n,n))
  for(i in 1:nrow(distance)) {
    distance_sort[i,] = res[[i]]$x # ith row is the ranked ith row of kernel distance
    idx[i,] = res[[i]]$ix # index of the ranks of ith row
  }
  A = array(0,c(n,n))
  di = distance_sort[,2:(k+2)]
  rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum)) # parameter estimation
  id = idx[,2:(k+2)]
  lambda = max(mean(rr),0)
  return(list(distX = distance,
              distX1 = distance_sort,
              lambda = lambda,
              idx = idx
              ))
}

normalized_GL = function(affinity, type = 3){
  d <- rowSums(affinity)
  d[d == 0] <- .Machine$double.eps
  D <- diag(d)
  L <- D - affinity
  if (type == 1) {
    NL <- L
  }
  else if (type == 2) {
    Di <- diag(1/d)
    NL <- Di %*% L
  }
  else if (type == 3) {
    Di <- diag(1/sqrt(d))
    NL <- Di %*% L %*% Di
  }
  return(NL)
}
"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {

  # set the needed parameters
  if(is.na(c)) {
    c = dim(A)[1]
  }
  if(c>dim(A)[1]) {
    c = dim(A)[1]
  }
  if(is.na(isMax)) {
    isMax = 1
  }
  if(is.na(isSym)) {
    isSym = 1
  }

  # compute the eigenvalues and eigenvectors of A
  if(isSym==1) {
    eigen_A = eigen(A,symmetric=TRUE)
  }
  else {
    eigen_A = eigen(A)
  }
  v = eigen_A$vectors
  d = eigen_A$values

  # sort the eigenvectors
  if(isMax == 0) {
    eigen_A_sorted = sort(d,index.return=TRUE)
  }
  else {
    eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
  }
  d1 = eigen_A_sorted$x
  idx = eigen_A_sorted$ix
  idx1 = idx[1:c]

  # compute the results
  eigval = d[idx1]
  eigvec = Re(v[,idx1])
  eigval_full = d[idx]

  return(list(eigval=eigval,eigvec=eigvec,eigval_full=eigval_full))

}

"L2_distance_1" <- function( a, b ) {

  if(dim(a)[1] == 1) {
    a = rbind(a,rep(0,dim(a)[2]))
    b = rbind(b,rep(0,dim(b)[2]))
  }

  aa = apply(a*a,MARGIN=2,FUN=sum)
  bb = apply(b*b,MARGIN=2,FUN=sum)
  ab = t(a) %*% b
  d1 = apply(array(0,c(length(t(aa)),length(bb))),MARGIN=2,FUN=function(x){ x = t(aa) })
  d2 = t(apply(array(0,c(length(t(bb)),length(aa))),MARGIN=2,FUN=function(x){ x = t(bb) }))
  d = d1 + d2 - 2 * ab
  d = Re(d)
  d = matrix(mapply(d,FUN=function(x) { return(max(max(x),0)) }),nrow=nrow(d),ncol=ncol(d))
  d_eye = array(1,dim(d))
  diag(d_eye) = 0
  d = d * d_eye

  return(d)

}





