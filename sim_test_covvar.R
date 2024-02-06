# Require
# library(arima)
library(Matrix)

# Import Fisher 2011 test
source("test_covvar_identity.R")

set.seed(1234)

# Simulations
n <- 30
q <- 200

phi1 <- 0.5

data <- t(sapply(1:n, function(id) arima.sim(model = list(ar = phi1), n = q)))

# Rejet de H0
T1_pvalue(data)
T2_pvalue(data)

# Calcul de la Sigma^-1/2

phi_hat_vect <- rep(-phi1, (q - 1))
square_root_inv_hat_Sigma <- bandSparse(q,
    k = c(1, 0),
    diagonals = list(phi_hat_vect, c(
        sqrt(1 - phi1^2),
        rep(1, (q - 1))
    ))
)

T1_pvalue(data %*% square_root_inv_hat_Sigma)
T2_pvalue(data %*% square_root_inv_hat_Sigma)

## Simulation gaussienne
data_gauss <- t(sapply(1:n, function(id) rnorm(q)))
T1_pvalue(data_gauss)
T2_pvalue(data_gauss)

# TODO Faire graphique de la pvalue en fonction du c (en fixant 3 valeurs de n)