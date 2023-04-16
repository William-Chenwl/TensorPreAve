#' @title Eigenvalue Plot of a Random Sample
#' @description Function to plot the eigenvalues of the sample covariance matrix of a randomly chosen sample.
#' @details Input a tensor time series and a mode index, output the plot of eigenvalues of the sample covariance matrix of the given mode,
#' with a randomly chosen sample of the mode-\eqn{k} fibres. This helps users to choose the parameter \code{eigen_j} in function \code{pre_est}.
#' A large dip should be observed at the (\eqn{r_k+1})-th position of the plot,
#' and user can choose \code{eigen_j} to be a bit larger than the position of dip observed to avoid missing potential weak factors. If such a dip
#' is not observed, try to run the function for a few times until it can be observed.
#' @param X A 'Tensor' object defined in package \pkg{rTensor} with \eqn{K+1} modes. Mode-1 should correspond to the time mode.
#' @param k The mode to plot the eigenvalues for.
#' @export
#' @import rTensor MASS
#' @examples
#' # Example of a real data set
#' set.seed(800)
#' pre_eigenplot(value_weight_tensor, k = 2)
#'
#' \donttest{
#' # Example using generated data
#' K = 2
#' T = 100
#' d = c(40,40)
#' r = c(2,2)
#' re = c(2,2)
#' eta = list(c(0,0),c(0,0))
#' u = list(c(-2,2),c(-2,2))
#' set.seed(10)
#' Data_test = tensor_data_gen(K,T,d,r,re,eta,u)
#' X = Data_test$X
#' pre_eigenplot(X, k = 1)
#' }


################################# Function to Observe Eigenvalue Ratio for Pre-Averaging Procedure ##########################
pre_eigenplot = function(  X                                       # a 'Tensor' object defined in package 'rTensor'. Mode-1 should correspond to the time mode.
                          ,k                                       # the mode to plot the eigenvalues

)
{

  K = X@num_modes - 1                             # number of modes for the tensor time series
  d = X@modes[2:(K+1)]                                # dimensions in each mode
  T = X@modes[1]                                      # length of time series
  Z = diag(1,T) - rep(1,T) %*% t(rep(1,T))/T




  # z_k = z[k]
  d_k = d[k]
  dmk = d[(1:K)[-k]]

  mat_X_k = rTensor::unfold(X,row_idx = c(k + 1, 1), col_idx = c(1:(K+1))[-c(1, k + 1)])

  # record eigenvalue ratios for each random sample generated
  ## create random row index for each A_l with l \neq k
  random_index = list()
  # random_index_c = list()
  for (j in (1:K)[-k])
  {
    ri = sample.int(d[j],d[j]/2)  # random sample half indexes of each A_l
    random_index[[j]] = ri
    # random_index_c[[j]] = setdiff(1:d[j],ri)  # complementary sampling
  }


  # Store eigenvalue ratio for specific random sample
  index_comb = expand.grid(random_index[-k])  # combination of indexes
  total_n = dim(index_comb)[1]
  random_fibre = numeric(total_n)

  index_mul_n = append(rev(cumprod(rev(dmk))),1)[-1]   # index multiplication
  ## create random fibres based on Cartesian products of row indexes from A_l with l \neq k
  for (cc in 1:total_n)
  {
    index_cc = index_comb[cc,]
    rr = 0
    for (dd in 1:(K-1)){
      rr = rr + (index_cc[[dd]]-1) * index_mul_n[dd]
    }
    random_fibre[cc] = rr + 1
  }


  ## For K = 2, we can run the following lines (instead of the above paragraph) to create random_fibre
  # d_minus_k = prod(d)/d_k
  # n = d_minus_k / 2                   # number of mode-k fibres to sum in each random sample
  # random_fibre = sample.int(d_minus_k, n)




  tilde_X_m = matrix(rowSums(mat_X_k[,random_fibre]@data),c(d_k,T))
  tilde_Sigma_cov_m = tilde_X_m %*% Z %*% t(tilde_X_m) / T           # sample covariance matrix for each random sample

  eigen_tilde_Sigma_m = eigen(tilde_Sigma_cov_m)
  plot(eigen_tilde_Sigma_m$values, ylab = 'Eigenvalues')

}
