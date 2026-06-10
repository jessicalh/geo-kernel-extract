\name{nmrstats}
\alias{nmrstats}
\alias{LADDER}
\alias{assert_disjoint_proteins}
\alias{baseline_blup}
\alias{calibration_coverage}
\alias{check_protocol}
\alias{drop_one_block_oof}
\alias{estimate_regions}
\alias{feature_blocks}
\alias{fit_element}
\alias{fit_joint}
\alias{fit_per_element}
\alias{folds_block_time}
\alias{folds_leave_cluster}
\alias{folds_leave_protein_out}
\alias{freeze_regions}
\alias{importance}
\alias{joint_views}
\alias{linear_relaimpo}
\alias{load_substrate}
\alias{make_folds}
\alias{morans_i}
\alias{n_eff_ar1}
\alias{resid_target}
\alias{run_all}
\alias{run_ladder}
\alias{run_nulls}
\alias{support_per_stratum}
\alias{within_protein_morans_i}
\title{NMR shielding row-design inferential stack}
\description{
Statistics-only analysis package for emitted \code{row_design} CSV/NPY substrates.
The public entry point is \code{run_all()}, which validates the emit, builds full-scope
folds, estimates regions, fits the ladder and per-element models, runs proxy/null
diagnostics, and writes tabular artifacts.
}
\usage{
run_all(dir, out, nfold = 3L, nthreads = 1L, load_t2 = TRUE, save_fits = FALSE)
load_substrate(dir, load_t2 = TRUE, load_embedding = FALSE)
}
\arguments{
\item{dir}{One or more emit directories.}
\item{out}{Output directory for analysis artifacts.}
\item{nfold}{Number of grouped full-scope CV folds.}
\item{nthreads}{Threads passed to \code{mgcv::bam}.}
\item{load_t2}{Whether to validate and load \code{row_design_target_T2.npy}.}
\item{load_embedding}{Whether to load the AIMNet2 embedding sidecar.}
\item{save_fits}{Whether to save large fitted \code{mgcv} objects.}
}
\value{
Functions return data tables, fitted objects, or nested result lists depending on the layer.
}
\keyword{models}
