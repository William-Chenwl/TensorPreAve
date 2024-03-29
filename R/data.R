#' @title Value weighted Fama-French portfolio returns data.
#' @name value_weight_tensor
#' @docType data
#' @description Value weighted Fama-French portfolio returns data formed on size and operating profitability of Chen and Lam (2023).
#'
#' @details Stocks are categorized into 10 different sizes (market equity, using NYSE market equity deciles) and 10 different
#' operating profitability (OP) levels (using NYSE OP deciles. OP is annual revenues minus cost of goods
#' sold, interest expense, and selling, general, and administrative expenses divided by book equity for the last
#' fiscal year end). The stocks in each of the 10 × 10 categories form a portfolio using value
#' weighted. We use monthly data from July 1973 to June 2021, so that T = 576, and each
#' data tensor we have thus has size 10 × 10 × 576. Since the market factor is certainly pervasive in financial returns, we use the CAPM to remove
#' its effects and facilitate detection of potentially weaker factors.
#'
#'
#' @format A 576 × 10 × 10 'Tensor' object defined in package \pkg{rTensor}, where mode-1,2,3 correspond to time, OP levels and size levels, respectively.
#'
#' @references Chen, W. and Lam, C. (2023). Rank and Factor Loadings Estimation in Time Series Tensor Factor Model by Pre-averaging. Manuscript.
NULL



#' @title Equal weight Fama-French portfolio returns data.
#' @name equal_weight_tensor
#' @docType data
#' @description Equal weight Fama-French portfolio returns data formed on size and operating profitability of Chen and Lam (2023).
#'
#' @details Stocks are categorized into 10 different sizes (market equity, using NYSE market equity deciles) and 10 different
#' operating profitability (OP) levels (using NYSE OP deciles. OP is annual revenues minus cost of goods
#' sold, interest expense, and selling, general, and administrative expenses divided by book equity for the last
#' fiscal year end). The stocks in each of the 10 × 10 categories form a portfolio by equal weight. We use monthly data from July 1973 to June 2021, so that T = 576, and each
#' data tensor we have thus has size 10 × 10 × 576. Since the market factor is certainly pervasive in financial returns, we use the CAPM to remove
#' its effects and facilitate detection of potentially weaker factors.
#'
#'
#' @format A 576 × 10 × 10 'Tensor' object defined in package \pkg{rTensor}, where mode-1,2,3 correspond to time, OP levels and size levels, respectively.
#'
#' @references Chen, W. and Lam, C. (2023). Rank and Factor Loadings Estimation in Time Series Tensor Factor Model by Pre-averaging. Manuscript.
NULL
