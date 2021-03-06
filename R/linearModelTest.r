#' A test function for linear models
#'
#' This function tests screening performance on linear models with user specified dimension. The function generates synthetic data using linear model and prints the screened variable for all method to the screen.
#' @param n the sampel size
#' @param p the model dimension
#' @param beta.not.null the non-zero coefficient indexes
#' @param num.select the same as in "screening"
#' @param ebic the same as in "screening"
#'
linearModelTest <- function(n, p, beta.not.null = c(1, 2, 3), num.select = 5 * length(beta.not.null), ebic = FALSE) {
    n = 50
    p = 200
    x = matrix(rnorm(n * p), n, p)
    beta = rep(0, p)
    for (i in beta.not.null) {
        beta[i] = 2 * ((runif(1) > 0.5) - 0.5) * (abs(rnorm(1)) + 0.2)
    }
    y = x %*% beta + rnorm(n)
    family = 'gaussian'

    print('The true non-zero coefficients are:')
    print(beta.not.null)

    ##############
    ## SIS
    print('Start SIS...')
    result = screening(x, y, 'sis', num.select = num.select, family = family, ebic = ebic)
    print('SIS selects the following coefficients:')
    print(result$screen)

    ##############
    ## holp
    print('Start HOLP...')
    result = screening(x, y, 'holp', num.select = num.select, family = family, ebic = ebic)
    print('HOLP selects the following coefficients:')
    print(result$screen)

    ##############
    ## rrcs
    print('Start RRCS...')
    result = screening(x, y, 'rrcs', num.select = num.select, family = family, ebic = ebic)
    print('RRCS selects the following coefficients:')
    print(result$screen)

    ##############
    ## forward regression
    print('Start forward regression...')
    result = screening(x, y, 'forward', num.select = num.select, family = family, ebic = ebic)
    print('Forward regression selects the following coefficients:')
    print(result$screen)
}
