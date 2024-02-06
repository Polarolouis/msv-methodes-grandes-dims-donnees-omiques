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

    S_svd <- svd(S)
    S_eigenval <- S_svd$d
    trS <- sum(S_eigenval)

    S2_eigenval <- (S_eigenval)^2
    trS2 <- sum(S2_eigenval)

    S3_eigenval <- (S_eigenval)^3
    trS3 <- sum(S3_eigenval)

    S4_eigenval <- (S_eigenval)^4
    trS4 <- sum(S4_eigenval)

    a1 <- (1 / p) * trS
    a2 <- (n^2 / ((n + 1) * (n + 2) * p)) * (trS2 - 1 / n * trS^2)
    tau <- n^4 / ((n - 1) * (n - 2) * (n + 2) * (n + 4))
    a3 <- tau / p * (trS3 - 3 / n * trS2 * trS + 2 / n^2 * (trS)^3)
    gamma <- (n^5 * (n^2 + n + 2)) / ((n + 1) * (n + 2) * (n + 4) * (n + 6)*(n - 1) * (n - 2) * (n - 3))
    a4 <- gamma / p * (trS4 - 4 / n * trS3 * trS - (2 * n^2 + 3 * n - 6) / (n * (n^2 + n + 2)) * (trS2)^2 + (2*(5 * n + 6) / (n * (n^2 + n + 2))) * trS2 * (trS)^2 - ((5 * n + 6) / (n^2 * (n^2 + n + 2))) * (trS)^4)

    return(setNames(c(a1, a2, a3, a4), c("a1", "a2", "a3", "a4")))
}

#' Computes the T1 statistic from the paper
#'
#' @return T1 statistic for normality test
compute_T1_statistic <- function(X) {
    n <- nrow(X)
    p <- ncol(X)

    c <- p / n

    ai_vec <- unname(compute_arithmetic_means(X))

    a1 <- ai_vec[1]
    a2 <- ai_vec[2]
    a3 <- ai_vec[3]
    a4 <- ai_vec[4]

    T1 <- (n / (c * sqrt(8))) * (a4 - 4 * a3 + 6 * a2 - 4 * a1 + 1)
    return(T1)
}

T1_pvalue <- function(X) {
    return(1 - pnorm(compute_T1_statistic(X)))
}

#' Computes the T2 statistic from the paper
#'
#' @return T2 statistic for normality test
compute_T2_statistic <- function(X) {
    n <- nrow(X)
    p <- ncol(X)

    c <- p / n

    ai_vec <- unname(compute_arithmetic_means(X))

    a2 <- ai_vec[2]
    a4 <- ai_vec[4]

    T2 <- (n / sqrt(8 * (c^2 + 12 * c + 8))) * (a4 - 2 * a2 + 1)
    return(T2)
}

T2_pvalue <- function(X) {
    return(1 - pnorm(compute_T2_statistic(X)))
}