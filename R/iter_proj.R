#' @title Iterative Projection Estimator.
#' @description Function for Iterative Projection to re-estimate the factor loading matrices.
#' @details Input a tensor time series and initial projection direction, return the estimated factor loading matrices using iterative projection.
#' @param X A 'Tensor' object defined in package \pkg{rTensor} with \eqn{K+1} modes. Mode-1 should correspond to the time mode.
#' @param initial_direction Initial direction for projection, written in a list of \eqn{K} vectors. This can be obtained from the pre-averaging procedure.
#' @param proj_N Number of iterations, should be a positive integer. Default is 30.
#' @param z (Estimated) Rank of the core tensor, written as a vector of length \eqn{K}. Can be set as 1's when we only need to do rank estimation based on projected data. Default is 1's.
#' @return A list of \eqn{K} estimated factor loading matrices.
#' @export
#' @import rTensor MASS
#' @examples
#' \dontrun{
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
#' Q_PRE = pre_est(X)
#' Q_PROJ = iter_proj(X, initial_direction = Q_PRE, z = r)
#' Q_PROJ
#'
#' # Example of a real data set
#' set.seed(10)
#' Q_PRE = pre_est(value_weight_tensor)
#' Q_PROJ = iter_proj(value_weight_tensor, initial_direction = Q_PRE)
#' Q_PROJ
#'
#' set.seed(10)
#' Q_PRE = pre_est(value_weight_tensor)
#' Q_PROJ_2 = iter_proj(value_weight_tensor, initial_direction = Q_PRE, z = c(2,2))
#' Q_PROJ_2
#' }





############################## Function for Iterative Projection to Re-estimate Factor Loadings ###################################
iter_proj = function(    X                                         # the target 'Tensor' object, where mode-1 is the time mode
                       , initial_direction                         # initial direction for projection written in a list of K vectors, can be obtained from the pre-averaging procedure
                       , proj_N = 30                               # number of iterations
                       , z = rep(1,length(X@modes) - 1)            # (estimated) rank of the core tensor; can be set as 1's when we need to do rank estimation based on projected data.
)
  # output : a list of K matrices, where output[[k]] is the mode-k estimated factor loading matrix (or vector)
{

  K = length(X@modes) - 1
  d = X@modes[2:(K+1)]
  T = X@modes[1]

  # Demean X
  X_mean = apply(X@data, 2:(K+1), mean)
  X_mean_T = aperm(array(X_mean, c(d,T)), c(K+1,1:K))
  X_demean = as.tensor(X@data - X_mean_T)


  # Initialize the projection direction
  for (k in 1:K)
  {
    hat_Q_k_1 = initial_direction[[k]]
  }


  hat_Q_proj = list()

  ### Iterative Projection for proj_N times
  for (proj_n in 1:proj_N)
  {
    for (k in 1:K)
    {
      d_k = d[k]
      z_k = z[k]

      # Calculate q_minus_k
      q_minus_k = 1
      for (j in rev((1:K)[-k]))
      {
        q_minus_k = kronecker(q_minus_k,initial_direction[[j]])
      }

      # Calculate projected data Y_k
      mat_X_k_demean = unfold(X_demean,row_idx = c(k + 1, 1), col_idx = c(1:(K+1))[-c(1, k + 1)])
      Y_k = matrix(mat_X_k_demean@data %*% q_minus_k, c(d_k,T))

      eigen_Y_k = eigen(Y_k %*% t(Y_k))$vectors


      # Reproduce the projection direction
      hat_Q_k_1 = eigen_Y_k[,1]


      if (proj_n == proj_N)
      {
        hat_Q_k_proj = eigen_Y_k[,1:z_k]
        hat_Q_proj[[k]] = hat_Q_k_proj}
    }

  }
  return (hat_Q_proj)

}
