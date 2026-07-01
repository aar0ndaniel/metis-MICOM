# Cross-verification against the cSEM benchmark. cSEM is a Suggests-only
# dependency used for validation, never at runtime, so these tests are skipped
# when cSEM is unavailable.

test_that("Step 2 |c| agrees with cSEM::testMICOM within 0.02", {
  skip_if_not_installed("cSEM")
  skip_if_not_installed("seminr")
  fx <- make_fixture(k = 2L)
  R <- 199; seed <- 1987; alpha <- 0.05
  s2 <- metis_micom_step2(fx$model, fx$data, "grp", permutations = R, seed = seed)

  d <- fx$data
  da <- d[d$grp == "A", setdiff(names(d), "grp")]
  db <- d[d$grp == "B", setdiff(names(d), "grp")]
  cn <- fx$constructs
  meas <- paste(vapply(cn, function(nm) {
    items <- rownames(fx$model$outer_weights)[fx$model$outer_weights[, nm] != 0]
    sprintf("%s <~ %s", nm, paste(items, collapse = " + "))
  }, character(1)), collapse = "\n")
  struct <- sprintf("%s ~ %s", cn[2], cn[1])
  model <- paste(struct, meas, sep = "\n")
  fit <- cSEM::csem(.data = list(a = da, b = db), .model = model,
                    .approach_weights = "PLS-PM", .disattenuate = FALSE,
                    .PLS_modes = stats::setNames(rep("modeA", length(cn)), cn))
  mi <- cSEM::testMICOM(fit, .R = R, .seed = seed)
  csem_c <- mi$Step2$Test_statistic[[1]]

  for (con in cn) {
    metis_c <- s2$c_value[s2$construct == con]
    expect_lt(abs(abs(metis_c) - abs(csem_c[con])), 0.02)
  }
})
