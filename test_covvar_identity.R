#' This function computes the 4 arithmetic means estimators
#' from \cite{fisherTestingIdentityCovariance2012}
#'
#' @return a1 a2 a3 a4
compute_arithmetic_means <- function(X) {
    p <- ncol(X)
    n <- nrow(X)

    X <- as.matrix(X)
    # cov_X is named S in the paper
    S <- cov(X)

    tr <- function(X) {
        sum(diag(X))
    }

    a1 <- (1 / p) * tr(S)
    a2 <- (n^2 / ((n + 1) * (n + 2) * p)) * (tr(S^2) - 1 / n * tr(S)^2)
    tau <- n^4 / ((n - 1) * (n - 2) * (n + 2) * (n + 4))
    a3 <- tau / p * (tr(S^3) - 3 / n * tr(S^2) * tr(S) + 2 / n^2 * (tr(S))^3)
    gamma <- (n^5 * (n^2 + n + 2)) / ((n + 1) * (n + 2) * (n + 4) * (n + 6)*(n - 1) * (n - 2) * (n - 3))
    a4 <- gamma / p * (tr(S^4) - 4 / n * tr(S^3) * tr(S) - (2 * n^2 + 3 * n - 6) / (n * (n^2 + n + 2)) * (tr(S^2))^2 + (2*(5 * n + 6) / (n * (n^2 + n + 2))) * tr(S^2) * (tr(S))^2 - ((5 * n + 6) / (n^2 * (n^2 + n + 2))) * (tr(S))^4)

    return(setNames(c(a1, a2, a3, a4), c("a1", "a2", "a3", "a4")))
}

#' Computes the T1 statistic from the paper
#'
#' @return T1 statistic for normality test
compute_T1_statistic <- function(X) {
    n <- nrow(X)
    p <- ncol(X)

    c <- p / n

    ai_vec <- compute_arithmetic_means(X)

    a1 <- ai_vec[1]
    a2 <- ai_vec[2]
    a3 <- ai_vec[3]
    a4 <- ai_vec[4]

    T1 <- (n / (c * sqrt(8))) * (a4 - 4 * a3 + 6 * a2 - 4 * a1 + 1)
    return(T1)
}

#' Computes the T2 statistic from the paper
#'
#' @return T2 statistic for normality test
compute_T2_statistic <- function(X) {
    n <- nrow(X)
    p <- ncol(X)

    c <- p / n

    ai_vec <- compute_arithmetic_means(X)

    a2 <- ai_vec[2]
    a4 <- ai_vec[4]

    T2 <- (n / sqrt(8 * (c^2 + 12 * c + 8))) * (a4 - 2 * a2 + 1)
    return(T2)
}
