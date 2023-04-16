#' @title Rank and Factor Loadings Estimation
#' @description The complete procedure to estimate both rank and factor loading matrices simultaneously.
#' @details Input a tensor time series and return the estimated factor loading matrices and rank of core tensor.
#' @param X A 'Tensor' object defined in package \pkg{rTensor} with \eqn{K+1} modes. Mode-1 should correspond to the time mode.
#' @param proj_N Number of iterations for iterative projection. Default is 30.
#' @param r_range Approximate range of \eqn{r_k} (number of factors) to search from, written in a list of \eqn{K} vectors (e.g. \code{z = list(c(1,10),c(1,10))} for \eqn{K = 2}). Default range is 1 to 10 for all modes.
#' @param C_range The range of constant C for calculating threshold. Default is \code{seq(0,100,0.1)}.
#' @param M0 Number of random samples to generate in pre-averaging procedure. Usually set as 200, or \eqn{min(d_{-k}^2 /4, 1000)} for potential more samples. Default is 200.
#' @param M Number of chosen samples for pre-averaging. Usually can be set as constants (5 or 10) or 2.5 percents of \code{M0}. Default is 5.
#' @param B Number of bootstrap samples for estimating rank of core tensor by bootstrapped correlation thresholding. Default is 50.
#' @param input_r The rank of core tensor if it is already know, written as a vector of length \eqn{K}. If no input, it will be estimated. Default is \code{NULL}.
#' @return A list containing the following: \cr
#'     \code{rank}: A vector of \eqn{K} elements, indicating the estimated number of factors in each mode \cr
#'     \code{loadings}: A list of \eqn{K} estimated factor loading matrices.
#' @export
#' @import rTensor MASS
#' @importFrom stats var
#' @examples
#' # Example of real data set
#' set.seed(10)
#' results = rank_factors_est(value_weight_tensor)
#' results
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
#' results = rank_factors_est(X)
#' results
#' }





############ The complete procedure for Rank and Factor Loadings Estimation simultaneously ###################
rank_factors_est = function(    X                                                  # a 'Tensor' object defined in package 'rTensor'. Mode-1 should correspond to the time mode.
                              , proj_N = 30                                        # number of iterations for iterative projection
                              , r_range = NULL                                     # the approximate range of r_k (number of factors), written in a list of K vectors (e.g. list(c(1,5),c(1,5)) for matrix time series). The default range is c(1,10) for all modes.
                              , C_range = seq(0,100,0.1)[-1]                       # the range of constant C for calculating threshold for bootstrap rank estimator
                              , M0 = 200                                           # number of random samples. Usually set as 200, or min(d_minus_k^2 /4, 1000) for potential more samples
                              , M = 5                                              # number of chosen samples for pre-averaging. Usually can be set as 5 or 10 or 2.5% of M0.
                              , B = 50                                             # number of bootstrap samples
                              , input_r = NULL                                     # the rank of core tensor if it is already know (e.g. c(2,2)), the length of the vector should be the same as number of modes. If no input, it will be estimated.
)
  # output: a list of two elements
  #                  output$rank : a vector of K elements of estimated number of factors in each mode, where K is the number of modes for the tensor time series
  #                  output$loadings: a list of K estimated factor loading matrices
{

  K = length(X@modes) - 1                             # number of modes for the tensor time series
  d = X@modes[2:(K+1)]                                # dimensions in each mode


  ## pre-averaging estimator for strongest factor
  Q_pre = pre_est(X, M0 = M0, M = M)

  if (is.null(input_r))
  {
    ## iterative projection for estimating the rank of core tensor
    initial_direction = list()
    for (k in 1:K)
    {
      initial_direction[[k]] = (Q_pre[[k]][,1])
    }
    Q_proj = iter_proj(X,initial_direction, proj_N, z = d)

    Q_proj_1 = list()
    for (k in 1:K)
    {
      Q_proj_1[[k]] = (Q_proj[[k]][,1])
    }

    bs_rank = bs_cor_rank(X, Q_proj_1, r_range, C_range, B)


    ## Estimation of factor loading spaces
    Q_hat = list()
    for (k in 1:K)
    {
      Q_hat[[k]] = Q_proj[[k]][,1:bs_rank[k]]
    }
  }

  else
  {
    ## iterative projection for estimating the factor loading spaces
    initial_direction = list()
    for (k in 1:K)
    {
      initial_direction[[k]] = (Q_pre[[k]][,1])
    }
    Q_proj = iter_proj(X,initial_direction, proj_N, z = d)

    Q_hat = list()
    for (k in 1:K)
    {
      Q_hat[[k]] = Q_proj[[k]][,1:input_r[k]]
    }
  }

  if (is.null(input_r))
  {
    output_rank = bs_rank
  }
  else
  {
    output_rank = input_r
  }
  return(list('rank' = output_rank, 'loadings' = Q_hat))

}
