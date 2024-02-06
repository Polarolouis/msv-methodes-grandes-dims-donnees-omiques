require(MultiVarSel)
require(ggplot2)
require(tidyr)
require(mixOmics)

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
data(nutrimouse)
table <- nutrimouse

table2 <- cbind(
    genotype = table$genotype,
    diet = table$diet, table$gene #, table$lipid
)

# Ne s'attend pas une structure de Toeplitz mais plutot à une par block

genotype <- table2$genotype
diet <- table2$diet


knitr::kable(table2[1:6,1:5])

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
X <- model.matrix(lm(Y~genotype:diet + 0, data=table2))

cat(paste0("$", array_to_LaTeX(X, n = 3), "$"))

p <- ncol(X)
n <- nrow(X)

q <- dim(Y)[2]

## Rescale de Y pour éviter de donner de l'importance a certine variable
Y <- scale(Y)

# Résidus
residus <- lm(as.matrix(Y) ~ X - 1)$residuals

pvalue <- whitening_test(residus)

## Testons des dépendances
results_whithening_choice <- whitening_choice(residus, typeDeps = c("AR1", "nonparam", "ARMA", "no_whitening"), pAR = 1, qMA = 1)

knitr::kable(results_whithening_choice,
  caption = "Tableau de résultats des tests de Portmanteau pour les différentes méthodes"
)

nb_repli <- 5000
nb.cores <- parallel::detectCores() - 1
whitening <- "no_whitening"

square_root_inv_hat_Sigma <- whitening(residus, whitening, pAR = 1, qMA = 0)


filename <- paste0("Freqs_nutrimouse_TOEPLITZ_nb_repli_", nb_repli,"_", whitening ,".Rdata")

if (!file.exists(filename)) {
  Freqs = variable_selection(Y, X, square_root_inv_hat_Sigma, nb_repli = 5000, parallel = TRUE, nb.cores = nb.cores)
  save(Freqs, file = filename)
} else {
  load(filename)
}
