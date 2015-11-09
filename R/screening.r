#' An efficient variable screening method
#'
#' This function implements 4 different screening methods (SIS, HOLP, RRCS and Forward regression) for linear models and 3 (excluding RRCS) for generalized linear models.
#' @param x the predictor variables, each row corresponds to an observation. Should be a numeric matrix instead of a data.frame
#' @param y the observation.
#' @param method the screening method to use. Choices are "sis", "holp", "rrcs", "forward". Default to "holp".
#' @param num.select the number of variables to keep after screening. Default to half of the sample size. It will not be used if ebic is set to be TRUE.
#' @param family the model type choices are the same as glmnet. Default to be 'gaussian'.
#' @param ebic Indicate whether the extended BIC should be used to determine the number of variables to keep. If ebic is TRUE, then the algorithm will use ebic to terminate the screening procedure and num.select will be ignored.
#' @param ebic.gamma tunning parameter for ebic (between 0 and 1). Gamma = 0 corresponds to the usual BIC. default to be 1.
#' @return a list of two variables "screen" and "method". "screen" contains the index of the selected variables and "method" indicates the method of the screening.
#'
#' @examples There are one unit test function and two integrated test functions. Two integrated function test on linear model and logistic model. User specify the sample size, dimension and the true indexes. The two function generate simulate data and coefficients and print the screening results for all methords.
#'
#' linearModelTest(n = 50, p = 100, beta.not.null = c(1, 2, 3), num.select = 20)
#' logisticTest(n = 50, p = 100, beta.not.null = c(1, 2, 3), nums.select = 20)
#'
#' @references
#' Fan, Jianqing, and Jinchi Lv. "Sure independence screening for ultrahigh dimensional feature space." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.5 (2008): 849-911.
#'
#' Wang, Xiangyu, and Chenlei Leng. "High-dimensional ordinary least-squares projection for screening variables." arXiv preprint arXiv:1506.01782 (2015).
#'
#' Li, Gaorong, et al. "Robust rank correlation based screening." The Annals of Statistics 40.3 (2012): 1846-1877.
#'
#' Wang, Hansheng. "Forward regression for ultra-high dimensional variable screening." Journal of the American Statistical Association 104.488 (2009): 1512-1524.
#'
screening <- function(x, y, method = 'holp', num.select = floor(dim(x)[1]/2), family = 'gaussian', ebic = FALSE, ebic.gamma = 1) {

    # standardize
    x = as.matrix(x)
    X = scale(x)
    if (family == 'gaussian'){
        Y = y - mean(y)
    }
    else {
        Y = y
    }
    if (is.null(dim(X))) {
        p = 1
        n = length(X)
    } else {
        n = dim(X)[1]
        p = dim(X)[2]
    }

    # if p is smaller than the required number of variables, return all
    if (p == 1 || (p < num.select && !ebic)) {
        selectedVariable = 1 : p
    }
    else {
        # for the linear case, it is easy to compute everything.
        if (family == 'gaussian') {
            if (method == 'holp') {
                OLS = t(X) %*% solve(X %*% t(X) + diag(n) * 1, Y)
                ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'sis') {
                OLS = t(X) %*% Y
                ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'rrcs') {
                OLS = .rankScreening(X, Y)
                ranking = sort(abs(OLS), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'forward') {
                if (ebic) {
                    ranking = .forwardRegression(X, Y, n - 1, family, ebic.gamma)
                    bestModel = which.min(ranking$bic)
                    result = ranking$select[1 : bestModel]
                }
                else {
                    result = .forwardRegression(X, Y, num.select, family, ebic.gamma)$select
                }
            }
            else {
                stop('The method is unknown and not supported/')
            }
        }
        else {
            if (method == 'holp') {
                require(glmnet)
                model = glmnet(x = X, y = Y, family = family, alpha = 0, lambda = 1, intercept = FALSE)
                coefs = coef(model)[-1]
                ranking = sort(abs(coefs), index.return = TRUE, decreasing = TRUE)
                if (ebic) {
                    result = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'sis') {
                residuals = rep(0, p)
                for (i in 1 : p) {
                    model = glm(Y ~ X[, i] - 1, family = family)
                    residuals[i] = model$deviance
                }
                ranking = sort(residuals, index.return = TRUE)
                if (ebic) {
                    results = .ebicRanking(X, Y, ranking$ix, family, ebic.gamma)$select
                }
                else {
                    result = ranking$ix[1 : num.select]
                }
            }
            else if (method == 'forward') {
                if (ebic) {
                    ranking = .forwardRegression(X, Y, n - 1, family, ebic.gamma)
                    bestModel = which.min(ranking$bic)
                    result = ranking$select[1 : bestModel]
                }
                else {
                    result = .forwardRegression(X, Y, num.select, family, ebic.gamma)$select
                }
            }
            else if (method == 'rrcs') {
                stop('rrcs is not supported for GLM')
            }
            else {
                stop('The method is unknown and not supported.')
            }
        }
    }
    return (list(screen = result, method = method))
}

#############################################################################

.rankScreening <- function(X, Y) {
  n = dim(X)[1]
  p = dim(X)[2]
  w = rep(0,p)
  Ynew = sort(Y, index.return = T)
  for(j in 1 : n)
      for(k in j : n)
          w = w + (X[Ynew$ix[k], ] > X[Ynew$ix[j], ])
  w = w / n / (n-1)
  w = w - 1 / 4
  return (abs(w))
}

#############################################################################

.forwardRegression <- function(X, Y, num.select, family, ebic.gamma = 1) {

    # get the dimension
    n = dim(X)[1]
    p = dim(X)[2]

    # store the used variables, including good variable and bad variable
    usedVariables = NULL

    # store the selected variables
    selectedVariables = NULL

    # store the residual sum of squares
    rss = rep(0, p)

    # store the bic values
    bic = rep(Inf, num.select)

    iteration = 0
    while (iteration < num.select) {

        # to compute for each variable the deviance of the model if they were added
        for (i in setdiff(1 : p, usedVariables)) {
            activeVariables = c(selectedVariables, i)
            model = try(glm(Y ~ X[, activeVariables] - 1, family = family))
            if (inherits(model, "try-error")) {
                rss[i] = Inf
                usedVariables = c(usedVariables, i)
            } else {
                rss[i] = model$deviance
            }
        }
        if (min(rss) == Inf) {
            break
        }

        # select the variabel that gives the smallest deviance to the model
        toAdd = which.min(rss)
        selectedVariables = c(selectedVariables, toAdd)
        usedVariables = c(usedVariables, toAdd)

        # record the corresponding bic value
        iteration = iteration + 1
        bic[iteration] = .ebic(rss[toAdd], p, n, iteration, ebic.gamma)
        rss[toAdd] = Inf
    }

    return (list(select = selectedVariables, bic = bic))
}

##############################################################################

.ebic <- function(deviance, model.size, sample.size, num.select, ebic.gamma) {
    return (deviance + num.select * (log(sample.size) + 2 * ebic.gamma * log(model.size)))
}

##############################################################################

.ebicRanking <- function(X, Y, sortedVariables, family, ebic.gamma) {
    # get the dimension
    n = dim(X)[1]
    p = dim(X)[2]

    # store the currently selected variables
    selectedVariables = NULL

    # store the bic values
    bic = rep(Inf, n - 1)

    iteration = 0
    while (iteration < n - 1) {

        iteration = iteration + 1

        # to compute for each variable the deviance of the model if they were added
        i = sortedVariables[iteration]
        selectedVariables = c(selectedVariables, i)
        model = try(glm(Y ~ X[, selectedVariables] - 1, family = family))
        if (inherits(model, "try-error")) {
            rss = Inf
        } else {
            rss = model$deviance
        }

        if (rss == Inf) {
            break
        }

        # record the corresponding bic value
        bic[iteration] = .ebic(rss, p, n, iteration, ebic.gamma)
    }

    bestModel = which.min(bic)
    return (list(select = sortedVariables[1 : bestModel], bic = bic))
}
