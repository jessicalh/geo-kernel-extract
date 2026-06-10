dependency_anchor <- function() {
  # These packages are part of the declared statistical stack. Some real-pass
  # adjustments use faster shrinkage/bam paths, but keeping the namespaces
  # anchored makes DESCRIPTION honest and check-clean.
  invisible(list(
    lme4 = lme4::lmerControl,
    nlme = nlme::lmeControl
  ))
}
