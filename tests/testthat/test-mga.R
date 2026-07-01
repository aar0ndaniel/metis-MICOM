test_that("metis_pls_mga wraps seminr and reports the two-tailed decision rule", {
  fx <- make_fixture(k = 3L)
  res <- suppressWarnings(metis_pls_mga(fx$model, fx$data, "grp", nboot = 60, alpha = 0.05))
  expect_s3_class(res, "metis_mga")
  expect_equal(res$engine, "seminr::estimate_pls_mga")
  expect_true(all(c("path", "p_value", "decision") %in% names(res$paths)))
  sig <- with(res$paths, (p_value < 0.05) | (p_value > 0.95))
  expect_equal(res$paths$decision == "significant", sig)
})

test_that("metis_pls_mga cautions when MICOM invariance is not established", {
  fx <- make_fixture(k = 3L)
  fake_micom <- list(invariance = list(partial = FALSE, full = FALSE, message = "not established"))
  res <- metis_pls_mga(fx$model, fx$data, "grp", nboot = 60, micom_result = fake_micom)
  expect_true(any(grepl("MICOM|invariance|caution", unlist(res), ignore.case = TRUE)))
})

test_that("print.metis_mga runs without error", {
  fx <- make_fixture(k = 3L)
  res <- suppressWarnings(metis_pls_mga(fx$model, fx$data, "grp", nboot = 60))
  expect_output(print(res), "MGA")
})
