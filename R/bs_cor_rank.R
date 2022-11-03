#' @title Bootstrap Rank Estimation.
#' @description Function to estimate the rank of the core tensor by Bootstrapped Correlation Thresholding.
#' @details Input a tensor time series and estimated projection directions, return the estimated rank of core tensor.
#' @param X A 'Tensor' object defined in package \pkg{rTensor} with \eqn{K+1} modes. Mode-1 should correspond to the time mode.
#' @param initial_direction Initial direction for projection, written in a list of \eqn{K} vectors. This can be obtained from the iterative projection procedure.
#' @param r_range Approximate range of \eqn{r_k} (number of factors) to search from, written in a list of \eqn{K} vectors (e.g. \code{z = list(c(1,10),c(1,10))} for \eqn{K = 2}). Default range is 1 to 10 for all modes.
#' @param C_range The range of constant C for calculating threshold. Default is \code{seq(0,100,0.1)}.
#' @param B Number of bootstrap samples. Default is 50.
#' @return A vector of length \eqn{K}, indicating estimated number of factors in each mode.
#' @export
#' @import rTensor MASS
#' @examples
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
#' Q_PROJ = iter_proj(X, initial_direction = Q_PRE)
#' bs_rank = bs_cor_rank(X, Q_PROJ)
#' bs_rank
#'
#' # Example of real data set
#' set.seed(10)
#' Q_PRE = pre_est(value_weight_tensor)
#' Q_PROJ = iter_proj(value_weight_tensor, initial_direction = Q_PRE)
#' bs_rank = bs_cor_rank(value_weight_tensor, Q_PROJ)
#' bs_rank




################ Function for Bootstrapped Rank Estimation for the Core Tensor ###################################
bs_cor_rank = function(    X                                                  # the target 'Tensor' object, where mode-1 is the time mode
                         , initial_direction                                  # direction for projection written in a list of K vectors, can be obtained from the iterative projection
                         , r_range = NULL                                     # the approximate range of r_k (number of factors), written in a list of K vectors (e.g. list(c(1,5),c(1,5)) for matrix time series). The default range is c(1,10) for all modes.
                         , C_range = seq(0,100,0.1)[-1]                       # the range of C for calculating threshold
                         , B = 50                                             # number of bootstrap samples
)
  # output : a vector of K elements of estimated number of factors in each mode
{
  ### Create educated grid from the un-Bootstrapped data
  K = length(X@modes) - 1
  d = X@modes[2:(K+1)]
  T = X@modes[1]

  # Demean X
  X_mean = apply(X@data, 2:(K+1), mean)
  X_mean_T = aperm(array(X_mean, c(d,T)), c(K+1,1:K))
  X_demean = as.tensor(X@data - X_mean_T)


  r_sw_hat = numeric(K)

  if (is.null(r_range))
  {
    r_range = rep(list(c(1,10)),K)
  }

  # Eigenvalues Record
  # eigen_R_sw_record = list()


  for (k in 1:K)
  {

    d_k = d[k]
    d_minus_k = prod(d)/d_k

    # Calculate q_minus_k
    q_minus_k = 1
    for (j in rev((1:K)[-k]))
    {
      q_minus_k = kronecker(q_minus_k,initial_direction[[j]])
    }

    # Calculate projected data Y_k
    mat_X_k_demean = unfold(X_demean,row_idx = c(k + 1, 1), col_idx = c(1:(K+1))[-c(1, k + 1)])
    Y_k = matrix(mat_X_k_demean@data %*% q_minus_k, c(d_k,T))

    # Calculate Correlation Matrix
    Sigma_y_k = Y_k %*% t(Y_k) / T
    Corr_y_k = diag(diag(Sigma_y_k)^-0.5) %*% Sigma_y_k %*% diag(diag(Sigma_y_k)^-0.5)
    eigen_R_k = eigen(Corr_y_k)$values

    #### Create educated grid for Bootstrap
    C0 = C_range

    # Set the maximum and minimum possible rank for search
    r_k_min = r_range[[k]][1]
    r_k_max = r_range[[k]][2]


    grid_result_k = numeric(length(C0))
    for (c in 1:length(C0))
    {
      grid_result_k[c] = sum(eigen_R_k > 1 + C0[c]/sqrt(T))
    }


    # define educated grid of constanc C
    grid_index_min_k = which(grid_result_k <= r_k_max)[1]
    grid_index_max_k = min(which(grid_result_k < r_k_min),length(C0) + 1) - 1
    C_k = C0[grid_index_min_k:grid_index_max_k]


    mat_X_k_demean = unfold(X_demean,row_idx = c(k + 1, 1), col_idx = c(1:(K+1))[-c(1, k + 1)])



    ####### Bootstrap: Sampling with Replacement ########

    sw_result_k = matrix(0,B,length(C_k))
    eigen_R_k_sw_record = matrix(0,B,d_k)

    for (b in 1:B)
    {

      swr = sample(1:d_minus_k,d_minus_k,replace = TRUE)
      sw_weight = diag(table(data.frame(type = factor(swr,1:d_minus_k), quantity = swr)$type))


      # create bootstrapped projected data
      Y_k_sw = matrix(mat_X_k_demean@data %*% sw_weight %*% q_minus_k, c(d_k,T))

      # Calculate boostrap correlation matrix
      Sigma_y_k_sw = Y_k_sw %*% t(Y_k_sw) / T
      Corr_y_k_sw = diag(diag(Sigma_y_k_sw)^-0.5) %*% Sigma_y_k_sw %*% diag(diag(Sigma_y_k_sw)^-0.5)

      eigen_R_k_sw = eigen(Corr_y_k_sw)$values
      eigen_R_k_sw_record[b,] = eigen_R_k_sw

      for (c in 1:length(C_k))
      {
        sw_result_k[b,c] = sum(eigen_R_k_sw > 1 + C_k[c]/sqrt(T))
      }


    }


    # result for rank estimation
    c_k_index = which.min(apply(sw_result_k,2,var))         # choose constanc C with minimum sample variance
    r_k_hat = getmode(sw_result_k[,c_k_index])              # estimate r_k
    r_sw_hat[k] = r_k_hat


  }

  return(r_sw_hat)

}





############# function to get mode of a vector ####################################
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
