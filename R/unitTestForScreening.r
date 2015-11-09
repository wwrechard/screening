p = 100
n = 100
X = matrix(rnorm(n * p), n, p)
Y = runif(n) > 0.5

testRankScreening <- function() {
    .rankScreening(X, Y)
}

testForwardRegression <- function() {
    result = .forwardRegression(X, Y, 20, 'binomial', 1)
    stopifnot(length(result$select) == 20)
}

testEbicRanking <- function() {
    .ebicRanking(X, Y, sample(1 : p, p), 'binomial', 1)
}


testScreeningForLinearModel <- function() {
    num.select = 20
    family = 'gaussian'
    result = screening(X, Y, 'holp', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'holp', family = family, ebic = TRUE)

    result = screening(X, Y, 'sis', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'sis', family = family, ebic = TRUE)

    result = screening(X, Y, 'rrcs', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'rrcs', family = family, ebic = TRUE)

    result = screening(X, Y, 'forward', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'forward', family = family, ebic = TRUE)
}

testScreeningForGLM <- function() {
    num.select = 20
    family = 'binomial'
    result = screening(X, Y, 'holp', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'holp', family = family, ebic = TRUE)

    result = screening(X, Y, 'sis', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'sis', family = family, ebic = TRUE)

    result = screening(X, Y, 'forward', num.select, family, FALSE)
    stopifnot(length(result$screen) == num.select)
    screening(X, Y, 'forward', family = family, ebic = TRUE)
}


unitTest <- function() {
    testRankScreening()
    testForwardRegression()
    testEbicRanking()
    testScreeningForLinearModel()
    testScreeningForGLM()

    print('All tests passed!')
}
