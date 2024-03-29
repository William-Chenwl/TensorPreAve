% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs_cor_rank.R
\name{bs_cor_rank}
\alias{bs_cor_rank}
\title{Bootstrap Rank Estimation.}
\usage{
bs_cor_rank(X, initial_direction, r_range = NULL, C_range = NULL, B = 50)
}
\arguments{
\item{X}{A 'Tensor' object defined in package \pkg{rTensor} with \eqn{K+1} modes. Mode-1 should correspond to the time mode.}

\item{initial_direction}{Direction corresponds to the strongest factors, written in a list of \eqn{K} vectors. This can be obtained from the iterative projection procedure by using function \code{iter_proj}.}

\item{r_range}{Approximate range of \eqn{r_k} (number of factors) to search from, written in a list of \eqn{K} vectors (e.g. \code{z = list(c(1,10),c(1,10))} for \eqn{K = 2}). Default range is 1 to 10 for all modes.}

\item{C_range}{The range of constant C for calculating threshold. Default is \code{seq(0,100,0.1)}, and set to be automatically tuned as data-driven.}

\item{B}{Number of bootstrap samples. Default is 50. Can be set as 10 to save time when dimension is large.}
}
\value{
A vector of length \eqn{K}, indicating estimated number of factors in each mode.
}
\description{
Function to estimate the rank of the core tensor by Bootstrapped Correlation Thresholding.
}
\details{
Input a tensor time series and estimated directions corresponding to the strongest factors, return the estimated rank of core tensor.
}
\examples{
# Example of real data set
set.seed(10)
Q_PRE = pre_est(value_weight_tensor)
Q_PROJ = iter_proj(value_weight_tensor, initial_direction = Q_PRE)
bs_rank = bs_cor_rank(value_weight_tensor, Q_PROJ)
bs_rank

\donttest{
# Example using generated data
K = 2
T = 100
d = c(40,40)
r = c(2,2)
re = c(2,2)
eta = list(c(0,0),c(0,0))
u = list(c(-2,2),c(-2,2))
set.seed(10)
Data_test = tensor_data_gen(K,T,d,r,re,eta,u)
X = Data_test$X
Q_PRE = pre_est(X)
Q_PROJ = iter_proj(X, initial_direction = Q_PRE)
bs_rank = bs_cor_rank(X, Q_PROJ)
bs_rank
}
}
