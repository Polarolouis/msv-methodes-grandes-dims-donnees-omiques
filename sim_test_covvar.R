# Require
# library(arima)
library(dplyr)
library(Matrix)
library(ggplot2)

# Import Fisher 2011 test
source("test_covvar_identity.R")

set.seed(1234)

# ## Simulation gaussienne
# data_gauss <- t(sapply(1:n, function(id) rnorm(q)))
# T1_pvalue(data_gauss)
# T2_pvalue(data_gauss)

# TODO Faire graphique de la pvalue en fonction du c (en fixant 3 valeurs de n)
# c = p/n, p = n*c
## Simulation avec AR(1)


simulate_AR1_T1_T2_pvalues <- function(N, n_values, c_values, phi1 = 0.5) {
    result <- do.call("rbind", lapply(n_values, function(n) {
        do.call("rbind", lapply(c_values, function(c) {
            all_results <- do.call("rbind", parallel::mclapply(seq(1:N), function(sim_id) {
                # Correspond au q de notre modèle :f le nb de colonnes
                p <- n * c
                cat(paste("sim", sim_id, "n", n, "c", c, "\n"))
                data <- t(sapply(1:n, function(id) {
                    arima.sim(
                        model = list(ar = phi1),
                        n = p
                    )
                }))
                T1_no_whitening <- T1_pvalue(data)
                T2_no_whitening <- T2_pvalue(data)
                result_no_whitening <- data.frame(n = n, c = c, p = p, T1_pvalue = T1_no_whitening, T2_pvalue = T2_no_whitening, type = "AR(1) no whitening")

                # Compute the sigma square root inv
                phi_hat_vect <- rep(-phi1, (p - 1))
                square_root_inv_hat_Sigma <- bandSparse(p,
                    k = c(1, 0),
                    diagonals = list(phi_hat_vect, c(
                        sqrt(1 - phi1^2),
                        rep(1, (p - 1))
                    ))
                )

                T1_white <- T1_pvalue(data %*% square_root_inv_hat_Sigma)
                T2_white <- T2_pvalue(data %*% square_root_inv_hat_Sigma)
                result_whitening <- data.frame(n = n, c = c, p = p, T1_pvalue = T1_white, T2_pvalue = T2_white, type = "AR(1) white")
                return(rbind(result_no_whitening, result_whitening))
            }, mc.cores = parallel::detectCores() - 1))

            # Aggregates results for this n and c
            aggregated_results <- all_results %>%
                dplyr::group_by(n, c, p, type) %>%
                dplyr::summarise(
                    T1_pvalue_mean = mean(T1_pvalue),
                    T1_pvalue_sd = sd(T1_pvalue),
                    T2_pvalue_mean = mean(T2_pvalue),
                    T2_pvalue_sd = sd(T2_pvalue), .groups = "keep"
                ) %>%
                dplyr::ungroup()
            return(aggregated_results)
        }))
    }))
    return(result)
}

simulate_gaussian_T1_T2_pvalues <- function(N, n_values, c_values, phi1 = 0.5) {
    result <- do.call("rbind", lapply(n_values, function(n) {
        do.call("rbind", lapply(c_values, function(c) {
            all_results <- do.call("rbind", parallel::mclapply(seq(1:N), function(sim_id) {
                # Correspond au q de notre modèle :f le nb de colonnes
                p <- n * c
                cat(paste("sim", sim_id, "n", n, "c", c, "\n"))
                data <- t(sapply(1:n, function(id) {
                    rnorm(
                        p
                    )
                }))
                T1_gauss <- T1_pvalue(data)
                T2_gauss <- T2_pvalue(data)
                result_gauss <- data.frame(n = n, c = c, p = p, T1_pvalue = T1_gauss, T2_pvalue = T2_gauss, type = "Gaussian")
                return(rbind(result_gauss))
            }, mc.cores = parallel::detectCores() - 1))

            # Aggregates results for this n and c
            aggregated_results <- all_results %>%
                dplyr::group_by(n, c, p, type) %>%
                dplyr::summarise(
                    T1_pvalue_mean = mean(T1_pvalue),
                    T1_pvalue_sd = sd(T1_pvalue),
                    T2_pvalue_mean = mean(T2_pvalue),
                    T2_pvalue_sd = sd(T2_pvalue), .groups = "keep"
                ) %>%
                dplyr::ungroup()
            return(aggregated_results)
        }))
    }))
    return(result)
}

N <- 50
phi1 <- 0.5
n_values <- seq(20, 300, by = 80)
c_values <- seq(1, 12)

filename <- paste0(
    "covvar_test_n_",
    paste0(n_values, collapse = "_"),
    "_c_",
    paste0(c_values, collapse = "_"),
    "_N_", N,
    "_AR_1_phi1_", phi1, ".Rds"
)

if (!file.exists(filename)) {
    result_AR1_T1_T2 <- simulate_AR1_T1_T2_pvalues(N = N, n_values = n_values, c_values = c_values)
    save(result_AR1_T1_T2, file = filename)
} else {
    load(filename)
}

filename_gauss <- paste0(
    "covvar_test_n_",
    paste0(n_values, collapse = "_"),
    "_c_",
    paste0(c_values, collapse = "_"),
    "_N_", N,
    "_gaussian.Rds"
)

if (!file.exists(filename_gauss)) {
    result_gauss <- simulate_gaussian_T1_T2_pvalues(N = N, n_values = n_values, c_values = c_values)
    save(result_gauss, file = filename_gauss)
} else {
    load(filename_gauss)
}

# Omit NAs
result <- result_AR1_T1_T2[complete.cases(result_AR1_T1_T2), ]

result$c <- as.numeric(result$c)
result$p <- as.numeric(result$p)
result$n <- as.numeric(result$n)
result$type <- as.factor(result$type)

result_gauss$type <- as.factor(result_gauss$type)

result <- rbind(result, result_gauss)

ggplot(result) +
    aes(x = c, y = T1_pvalue_mean, group = type) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(
        ymin = T1_pvalue_mean - T1_pvalue_sd,
        ymax = T1_pvalue_mean + T1_pvalue_sd
    ))+
    geom_hline(aes(yintercept = 0.05, color = "red")) +
    labs(color = "Seuil = 0.05") +
    facet_grid(vars(result$n), vars(result$type)) +
    ggtitle(paste0("Test sur la statistique T1 | N = ", N))

ggplot(result) +
    aes(x = c, y = T2_pvalue_mean, group = type) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(
        ymin = T2_pvalue_mean - T2_pvalue_sd,
        ymax = T2_pvalue_mean + T2_pvalue_sd
    ))+
    geom_hline(aes(yintercept = 0.05, color = "red")) +
    labs(color = "Seuil = 0.05") +
    facet_grid(vars(result$n), vars(result$type)) +
    ggtitle(paste0("Test sur la statistique T2 | N = ", N))

