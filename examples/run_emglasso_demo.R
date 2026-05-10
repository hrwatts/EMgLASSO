library(MASS)
library(mvtnorm)
library(glasso)

repo_root <- getwd()
if (!dir.exists(file.path(repo_root, "R"))) {
  repo_root <- normalizePath("..")
}

for (path in list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(path)
}

set.seed(123)

Tau0 <- c(0.1, 0.15, 0.2, 0.25, 0.3)
K0 <- length(Tau0)
N0 <- 1000
D0 <- 10
Mu0 <- lapply(1:K0, function(k) runif(D0, -100, 100))
Sigma0 <- lapply(1:K0, function(k) rspdmatrix(D0, 0.19 * k))

Sample0 <- rmmnorm(N0, D0, Tau0, Mu0, Sigma0)
Theta0 <- list(Tau = Tau0, Mu = Mu0, Sigma = Sigma0)
Emglasso0 <- emglasso(Sample0, Tau0)

output_dir <- file.path(repo_root, "examples", "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

capture.output(
  list(
    TrueTheta0 = list(TrueTau0 = Tau0, TrueMu0 = Mu0, TrueSigma0 = Sigma0),
    Theta0Hat = Emglasso0$Theta
  ),
  file = file.path(output_dir, "emglasso_demo_output.txt")
)
