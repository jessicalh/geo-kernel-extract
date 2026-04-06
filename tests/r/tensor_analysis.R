# tensor_analysis.R
# Generated: 2026-04-03T20:42:18.007299


library(reticulate)
np <- import("numpy")

# Load all 8 calculator shielding tensors
dir <- "/tmp/nmr_regression/orca_fresh"
calcs <- c("bs", "hm", "mc", "coulomb", "hbond", "pq", "ringchi", "disp")
labels <- c("BiotSavart", "HaighMallion", "McConnell", "Coulomb", "HBond", "PiQuad", "RingChi", "Dispersion")

shielding <- list()
for (i in seq_along(calcs)) {
  f <- file.path(dir, paste0(calcs[i], "_shielding.npy"))
  shielding[[labels[i]]] <- np$load(f)
}

# DFT reference
orca <- np$load(file.path(dir, "orca_total.npy"))

cat("=== Data loaded: ", nrow(shielding[[1]]), " atoms, 8 calculators ===\n\n")

# --- T2 linear independence ---
# T2 components are columns 5:9 (0-indexed 4:8) of each (N,9) array
# Compute pairwise cosine similarity of T2 vectors (flattened across atoms)

t2_vecs <- matrix(0, nrow=8, ncol=nrow(shielding[[1]]) * 5)
for (i in 1:8) {
  t2 <- shielding[[labels[i]]][, 5:9]  # columns 5-9 (1-indexed)
  t2_vecs[i,] <- as.vector(t2)
}

cosine_sim <- matrix(0, 8, 8)
for (i in 1:8) {
  for (j in 1:8) {
    cosine_sim[i,j] <- sum(t2_vecs[i,] * t2_vecs[j,]) / 
      (sqrt(sum(t2_vecs[i,]^2)) * sqrt(sum(t2_vecs[j,]^2)) + 1e-30)
  }
}
rownames(cosine_sim) <- labels
colnames(cosine_sim) <- labels

cat("=== T2 Cosine Similarity Matrix ===\n")
print(round(cosine_sim, 3))

# --- T0 correlation with DFT ---
cat("\n=== T0 vs DFT Correlation ===\n")
dft_t0 <- orca[, 1]
for (i in 1:8) {
  calc_t0 <- shielding[[labels[i]]][, 1]
  r <- cor(calc_t0, dft_t0, use="complete.obs")
  cat(sprintf("  %-15s r = %+.4f  (mean T0 = %.4e)\n", labels[i], r, mean(calc_t0)))
}

# --- Per-ring-type signal strength ---
cat("\n=== Per-ring-type mean |T0| (Biot-Savart) ===\n")
bs_type_t0 <- np$load(file.path(dir, "bs_per_type_T0.npy"))
type_names <- c("PheBenz", "TyrPhenol", "TrpBenz", "TrpPyrr", "TrpPerim", "HisImid", "HidImid", "HieImid")
for (i in 1:8) {
  cat(sprintf("  %-12s mean|T0| = %.4e\n", type_names[i], mean(abs(bs_type_t0[,i]))))
}

# --- Tensor symmetry checks ---
cat("\n=== Tensor Symmetry Checks ===\n")
for (i in 1:8) {
  s <- shielding[[labels[i]]]
  t0 <- s[, 1]  # T0
  t1 <- s[, 2:4]  # T1[3]
  t2 <- s[, 5:9]  # T2[5]
  t1_norm <- sqrt(rowSums(t1^2))
  t2_norm <- sqrt(rowSums(t2^2))
  cat(sprintf("  %-15s mean|T1|=%.2e  mean|T2|=%.2e  max|T0|=%.2e\n",
    labels[i], mean(t1_norm), mean(t2_norm), max(abs(t0))))
}

