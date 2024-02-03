require(MultiVarSel)
require(ggplot2)

array_to_LaTeX <- function(arr, n=nrow(arr), p=ncol(arr)){
  rows <- apply(arr, MARGIN = 1, paste, collapse = " & ")
  if(length(arr) > p) {
    rows <- apply(arr, MARGIN = 1, paste, collapse = " & ")
  }

  if(length(rows)>n){
    rows <- c(
      rows[1:ceiling(n / 2)],
      "\\vdots",
      rows[seq(length(rows) - floor(n / 2),length(rows))]
    )
  }
  matrix_string <- paste(rows, collapse = " \\\\ ")
  return(paste("\\begin{bmatrix}", matrix_string, "\\end{bmatrix}"))
}

# Import données
table <- read.table("Table_metabolome_FH_all.csv", header = TRUE, sep = ";", dec = ",")

# Retrait des 4 premières colonnes
table2 <- table[, -seq(1, 4)]

temperature <- table2$temperature
imbibition <- table2$imbibition

knitr::kable(table[1:6,5:9])

# Moyennes des colonnes
idx_null_mean <- as.integer(length(which(colMeans(table2[, 3:dim(table2)[2]]) == 0)))

# Pour retirer les metabolomes sans variations
idx_null_sd <- as.integer(length(which(apply(table2[, 3:dim(table2)[2]], 2, sd) == 0)))

# Verifier les NAs
idx_NAs <- as.integer(length(which(is.na(table2))))

knitr::asis_output(paste(
    "Il y a", idx_null_mean, "colonnes dont la moyenne est nulle,", idx_null_sd,
    "colonnes dont l'écart-type est nul et", idx_NAs, "NAs."
))



Y <- as.matrix(table2[,-c(1,2)])

# Permet de récupérer la matrice de design
X <- model.matrix(lm(Y~temperature:imbibition + 0, data=table2))

cat(paste0("$", array_to_LaTeX(X, n = 3), "$"))

p <- ncol(X)
n <- nrow(X)

q <- dim(Y)[2]

## Rescale de Y pour éviter de donner de l'importance a certine variable
Y <- scale(Y)

# Résidus
residus <- lm(as.matrix(Y) ~ X - 1)$residuals

pvalue <- whitening_test(residus)

# La pvaleur est trop grande donc les données ne sont pas blanche

## Testons des dépendances
results_whithening_choice <- whitening_choice(residus, typeDeps = c("AR1", "nonparam", "ARMA"), pAR = 1, qMA = 1)

knitr::kable(results_whithening_choice)
## Ainsi on utilise le non param

square_root_inv_hat_Sigma=whitening(residus, "nonparam", pAR = 1, qMA = 0)

nb_repli = 5000
nb.cores <- 5

filename <- paste0("Freqs_metabolomeFH_TOEPLITZ_nb_repli_", nb_repli, ".Rdata")

if (!file.exists(filename)) {
  Freqs = variable_selection(Y, X, square_root_inv_hat_Sigma, nb_repli = 5000, parallel = TRUE, nb.cores = 5)
  save(Freqs, file = filename)
} else {
  load(filename)
}

colnames(Freqs) <- c("Names_of_Y", "Names_of_X", "frequency")

plot(sort(Freqs$frequency, decreasing = TRUE), type='s')
sort(Freqs$frequency, decreasing = TRUE)[1:50]
