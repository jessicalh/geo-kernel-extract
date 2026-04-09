#!/usr/bin/env Rscript
# hidden_probe.R — Probe the MLP hidden layer dimensionality.
#
# Runs the trained MLP's scalar→hidden mapping on all atoms, captures
# the hidden activations, and computes their eigenspectrum.  Answers:
# how many of the 64 hidden dimensions is the model actually using?
#
# Also extracts the final kernel weights per atom and analyses their
# structure: how many distinct weighting patterns exist?
#
# Requires: a trained model checkpoint + the R export data.
# The heavy lifting (forward pass) is done in Python via reticulate.
#
# Usage:
#   cd learn/R
#   Rscript hidden_probe.R ../src/output/secondary/r_export ../runs/azimuthal_scalars_only_50

args <- commandArgs(trailingOnly = TRUE)
source("common.R")

if (length(args) < 2) {
  stop("Usage: Rscript hidden_probe.R <export_dir> <run_dir>")
}
export_dir <- args[1]
run_dir    <- args[2]

library(reticulate)

# ── Run the forward pass in Python to get hidden activations ─────────
cat("Running MLP forward pass via Python...\n")

py_run_string(sprintf("
import sys, json, torch, numpy as np
sys.path.insert(0, '../src')
sys.path.insert(0, '../../python')

from mutation_set.config import load_config
from mutation_set.model import ShieldingT2Model

# Load config and model
cfg = load_config('../src/calibration.toml')
header = json.load(open('%s/header.json'))
n_scalars = header['n_scalars']
n_kernels = header['n_kernels']

model = ShieldingT2Model(cfg, n_scalars, n_kernels)
state = torch.load('%s/best_model.pt', map_location='cpu', weights_only=True)
model.load_state_dict(state)
model.eval()

# Load export data
scalars_np = np.load('%s/scalars.npz')['data']
kernels_np = np.load('%s/kernels_norm.npz')['data']
N = len(scalars_np)

scalars_t = torch.tensor(scalars_np, dtype=torch.float32)
kernels_t = torch.tensor(kernels_np, dtype=torch.float32)

# Hook into the MLP hidden layers
hidden_activations = {}
def hook_fn(name):
    def fn(module, input, output):
        hidden_activations[name] = output.detach().numpy()
    return fn

# The MLP is model.mixing.mlp — FullyConnectedNet with layers
# Register hooks on each linear layer's output
for i, layer in enumerate(model.mixing.mlp):
    if hasattr(layer, 'weight'):
        layer.register_forward_hook(hook_fn(f'layer_{i}'))

# Forward pass
with torch.no_grad():
    # Get the relative weights before gating
    relative = model.mixing.mlp(scalars_t)
    magnitude = kernels_t.norm(dim=-1)
    gate = magnitude / (magnitude + model.mixing.gate_threshold)
    import math
    weights = relative * gate / math.sqrt(n_kernels)

# Save results
hidden1 = hidden_activations.get('layer_1', np.zeros((N, 1)))
hidden2 = hidden_activations.get('layer_3', np.zeros((N, 1)))
kernel_weights = weights.numpy()
relative_weights = relative.numpy()

np.savez('%s/hidden_probe.npz',
         hidden1=hidden1,
         hidden2=hidden2,
         kernel_weights=kernel_weights,
         relative_weights=relative_weights)
print(f'Hidden1: {hidden1.shape}, Hidden2: {hidden2.shape}')
print(f'Kernel weights: {kernel_weights.shape}')
", run_dir, run_dir, export_dir, export_dir, export_dir))

# ── Load results into R ──────────────────────────────────────────────
np <- import("numpy")
probe <- np$load(sprintf("%s/hidden_probe.npz", export_dir), allow_pickle = TRUE)

hidden1  <- probe$f[["hidden1"]]
hidden2  <- probe$f[["hidden2"]]
k_weights <- probe$f[["kernel_weights"]]
rel_weights <- probe$f[["relative_weights"]]

cat(sprintf("Hidden layer 1: %d atoms × %d units\n", nrow(hidden1), ncol(hidden1)))
cat(sprintf("Hidden layer 2: %d atoms × %d units\n", nrow(hidden2), ncol(hidden2)))
cat(sprintf("Kernel weights: %d atoms × %d kernels\n", nrow(k_weights), ncol(k_weights)))

# ── Eigenspectrum of hidden activations ──────────────────────────────

analyse_dim <- function(mat, label) {
  mat_c <- scale(mat, center = TRUE, scale = FALSE)
  cov_mat <- crossprod(mat_c) / nrow(mat_c)
  eig <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  eig <- pmax(eig, 0)
  cumvar <- cumsum(eig) / sum(eig)
  dim90 <- which(cumvar >= 0.90)[1]
  dim95 <- which(cumvar >= 0.95)[1]
  dim99 <- which(cumvar >= 0.99)[1]
  active <- sum(eig > max(eig) * 0.01)  # dims > 1% of max eigenvalue
  cat(sprintf("\n%s:\n", label))
  cat(sprintf("  Total dims: %d\n", ncol(mat)))
  cat(sprintf("  Active dims (>1%% of max eig): %d\n", active))
  cat(sprintf("  90%% variance: %d dims\n", dim90))
  cat(sprintf("  95%% variance: %d dims\n", dim95))
  cat(sprintf("  99%% variance: %d dims\n", dim99))

  data.frame(rank = seq_along(eig), eigenvalue = eig, cumvar = cumvar,
             label = label, stringsAsFactors = FALSE)
}

e1 <- analyse_dim(hidden1, "Hidden layer 1")
e2 <- analyse_dim(hidden2, "Hidden layer 2")
ew <- analyse_dim(k_weights, "Kernel weights (gated)")
er <- analyse_dim(rel_weights, "Relative weights (pre-gate)")

all_eigen <- rbind(e1, e2, ew, er)
all_eigen <- all_eigen[all_eigen$rank <= 30, ]

p1 <- ggplot(all_eigen, aes(x = rank, y = cumvar, colour = label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0.95, linetype = "dotted", colour = "grey50") +
  labs(title = "MLP hidden space dimensionality",
       subtitle = "How many dimensions does the model actually use?",
       x = "Eigenvalue rank", y = "Cumulative variance",
       colour = NULL) +
  coord_cartesian(ylim = c(0, 1))

save_fig(p1, "wb_hidden_probe", width = 10, height = 6)

# ── Per-kernel weight statistics ─────────────────────────────────────
cat("\n\nPer-kernel weight statistics (gated):\n")
cat(sprintf("  %-25s  %8s  %8s  %8s  %8s\n",
            "Kernel", "mean|w|", "std|w|", "frac>0", "frac_active"))

knames <- read.csv(sprintf("%s/kernel_names.csv", export_dir),
                    stringsAsFactors = FALSE)
weight_stats <- list()
for (k in seq_len(ncol(k_weights))) {
  w <- k_weights[, k]
  weight_stats[[k]] <- data.frame(
    kernel = knames$name[k],
    mean_abs = mean(abs(w)),
    std_abs = sd(abs(w)),
    frac_positive = mean(w > 0),
    frac_active = mean(abs(w) > 0.001),
    stringsAsFactors = FALSE
  )
}
ws_df <- do.call(rbind, weight_stats)
ws_df <- ws_df[order(-ws_df$mean_abs), ]
for (i in seq_len(min(20, nrow(ws_df)))) {
  r <- ws_df[i, ]
  cat(sprintf("  %-25s  %8.4f  %8.4f  %8.2f  %8.2f\n",
              r$kernel, r$mean_abs, r$std_abs, r$frac_positive, r$frac_active))
}

cat("\n══════════════════════════════════════════════════════════\n")
cat("  Done. Figure: wb_hidden_probe.pdf\n")
cat("══════════════════════════════════════════════════════════\n")
