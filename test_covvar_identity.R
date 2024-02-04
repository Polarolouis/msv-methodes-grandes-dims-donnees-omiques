#' This function computes the 4 arithmetic means estimators
#' from \cite{fisherTestingIdentityCovariance2012}
#'
#' @return a1 a2 a3 a4
compute_arithmetic_means <- function(X) {
    p <- ncol(X)
    n <- nrow(X)
    # cov_X is named S in the paper
    cov_X <- cov(X)

    tr <- function(X) {
        sum(diag(X))
    }

    a1 <- (1 / p) * tr(cov_X)
    a2 <- (n^2 / ((n + 1) * (n + 2) * p)) * (tr(S^2) - 1 / n * tr(S)^2)
    tau <- n^4/((n-1)*(n-2)*(n+2)*(n+4))
    a3 <- tau / p * (tr(S^3) - 3 / n * tr(S^2) * tr(S) + 2 / n^2 * (tr(S))^3)
    gamma <- (n^5 * (n^2 + n + 2)) / ((n + 1) * (n + 2) * (n + 4) * (n + 6)(n - 1) * (n - 2) * (n - 3))
    a4 <- gamma / p * (tr(S^4) - 4 / n * tr(S^3) * tr(S) - (2 * n^2 + 3 * n - 6) / (n * (n^2 + n + 2)) * (tr(S^2))^2 + (2(5 * n + 6) / (n * (n^2 + n + 2))) * tr(S^2) * (tr(S))^2 - ((5 * n + 6) / (n^2 * (n^2 + n + 2)))*(tr(S))^4)
    
    return(setNames(c(a1,a2,a3,a4), c("a1","a2","a3","a4")))
}