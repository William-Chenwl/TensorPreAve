#' @title Pre-Averaging Estimator
#' @description Function for the initial Pre-Averaging Procedure.
#' @details Input a tensor time series and return the estimated factor loading matrices using pre-averaging method.
#' @param X A 'Tensor' object defined in package \pkg{rTensor} with \eqn{K+1} modes. Mode-1 should correspond to the time mode.
#' @param z (Estimated) Rank of the core tensor, written as a vector of length \eqn{K}. For iterative projection purpose, we only need this to be 1's. Default is 1's.
#' @param M0 Number of random samples to generate, should be a positive integer. Usually set as 200, or \eqn{min(d_{-k}^2 /4, 1000)} for potential more samples. Default is 200.
#' @param M Number of chosen samples for pre-averaging, should be a positive integer. Usually can be set as constants (5 or 10) or 2.5 percents of \code{M0}. Default is 5.
#' @return A list of \eqn{K} estimated factor loading matrices.
#' @export
#' @import rTensor MASS
#' @examples
#' # Example of a real data set
#' set.seed(10)
#' Q_PRE = pre_est(value_weight_tensor)
#' Q_PRE
#'
#' set.seed(10)
#' Q_PRE_2 = pre_est(value_weight_tensor, z = c(2,2))
#' Q_PRE_2
#'
#' \donttest{
#' # Example using generated data
#' K = 2
#' T = 100
#' d = c(40,40)
#' r = c(2,2)
#' re = 10
#' eta = list(c(0,0),c(0,0))
#' u = list(c(-2,2),c(-2,2))
#' set.seed(10)
#' Data_test = tensor_data_gen(K,T,d,r,re,eta,u)
#' X = Data_test$X
#' Q_PRE = pre_est(X, z = r)
#' Q_PRE
#' }



################################# Function for Initial Pre-Averaging Procedure ##########################
pre_est = function(    X                                       # a 'Tensor' object defined in package 'rTensor'. Mode-1 should correspond to the time mode.
                     , z = rep(1,length(X@modes) - 1)          # (estimated) rank of the core tensor. For iterative projection purpose, we only need this to be 1's.
                     , M0 = 200                                # number of random samples. Usually set as 200, or min(d_minus_k^2 /4, 1000) for potential more samples
                     , M = 5                                   # number of chosen samples for pre-averaging. Usually can be set as 5 or 10 or 2.5% of M0.
)
  # output : a list of K matrices, where output[[k]] is the mode-k estimated factor loading matrix (or vector)
{

  K = length(X@modes) - 1                             # number of modes for the tensor time series
  d = X@modes[2:(K+1)]                                # dimensions in each mode
  T = X@modes[1]                                      # length of time series
  Z = diag(1,T) - rep(1,T) %*% t(rep(1,T))/T



  hat_Q = list()

  for (k in 1:K)
  {
    z_k = z[k]
    d_k = d[k]


    tilde_Sigma_cov_m_record = array(numeric(),c(M0,d_k,d_k))
    ER_record = numeric(M0)
    mat_X_k = unfold(X,row_idx = c(k + 1, 1), col_idx = c(1:(K+1))[-c(1, k + 1)])

    # record eigenvalue ratios for each random sample generated
    for (m in 1:M0)
    {


      ## create random row index for each A_l with l \neq k
      random_index = list()
      for (j in (1:K)[-k])
      {
        random_index[[j]] = sample.int(d[j],d[j]/2)  # random sample half indexes of each A_l
      }
      index_comb = expand.grid(random_index[-k])  # combination of indexes
      total_n = dim(index_comb)[1]
      random_fibre = numeric(total_n)
      dmk = d[(1:K)[-k]]
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

      # Calculate Eigenvalue Ratio
      eigen_tilde_Sigma_m = eigen(tilde_Sigma_cov_m)
      ER = eigen_tilde_Sigma_m$values[1] / eigen_tilde_Sigma_m$values[d_k/2]
      tilde_Sigma_cov_m_record[m,,] = tilde_Sigma_cov_m
      ER_record[m] = ER
    }

    ## Pre-Average Estimator of Q_k Using the largest M eigenvalue ratio samples

    tilde_Sigma_cov_pre = matrix(0, d_k, d_k)
    for (mm in order(-ER_record)[1:M])
    {
      tilde_Sigma_cov_pre = tilde_Sigma_cov_pre + tilde_Sigma_cov_m_record[mm,,]
    }
    tilde_Sigma_cov_pre = tilde_Sigma_cov_pre / M                       # the aggregrated sample covariance matrix
    eigen_tilde_Sigma_pre = eigen(tilde_Sigma_cov_pre)


    hat_Q_k_pre = array(eigen_tilde_Sigma_pre$vectors[,1:z_k],c(d_k,z_k))  # the pre-averaging estimator
    hat_Q[[k]] = hat_Q_k_pre
  }
  return(hat_Q)

}
