test_that("metis_micom returns a metis_micom object with the documented structure", {
  fx <- make_fixture(k = 2L)
  res <- metis_micom(fx$model, fx$data, "grp", permutations = 99, seed = 1)
  expect_s3_class(res, "metis_micom")
  expect_equal(res$method, "MICOM")
  expect_true(all(c("groups", "settings", "step1", "step2", "step3",
                    "admissibility", "invariance") %in% names(res)))
  expect_true(all(c("construct", "c_value", "ci_lower", "ci_upper",
                    "p_value", "decision") %in% names(res$step2)))
  expect_setequal(res$step2$construct, fx$constructs)
})

test_that("Step 2 c-values are deterministic for a fixed seed", {
  fx <- make_fixture(k = 2L)
  a <- metis_micom_step2(fx$model, fx$data, "grp", permutations = 99, seed = 42)
  b <- metis_micom_step2(fx$model, fx$data, "grp", permutations = 99, seed = 42)
  expect_equal(a$c_value, b$c_value)
  expect_equal(a$p_value, b$p_value)
})

test_that("Step 2 c-values are bounded in [0, 1] under sign alignment", {
  fx <- make_fixture(k = 2L)
  s2 <- metis_micom_step2(fx$model, fx$data, "grp", permutations = 99, seed = 1)
  expect_true(all(s2$c_value >= 0 & s2$c_value <= 1 + 1e-8, na.rm = TRUE))
})

test_that("Step 2 decision is invariant to indicator sign reflection", {
  fx <- make_fixture(k = 2L)
  base <- metis_micom_step2(fx$model, fx$data, "grp", permutations = 199, seed = 7)
  d2 <- reflect_indicators(fx$data, seminr::multi_items("CUEX", 1:3))
  m2 <- suppressMessages(seminr::estimate_pls(d2, fx$model$measurement_model, fx$model$structural_model))
  refl <- metis_micom_step2(m2, d2, "grp", permutations = 199, seed = 7)
  expect_equal(base$decision, refl$decision)
})

test_that("admissibility report is internally consistent", {
  fx <- make_fixture(k = 2L)
  res <- metis_micom(fx$model, fx$data, "grp", permutations = 99, seed = 1)
  adm <- res$admissibility
  expect_true(all(c("construct", "requested", "admissible", "dropped", "dropped_pct") %in% names(adm)))
  expect_true(all(adm$requested == 99))
  expect_true(all(adm$dropped == adm$requested - adm$admissible))
  expect_true(all(adm$dropped >= 0))
})

test_that("metis_micom rejects an unknown grouping variable", {
  fx <- make_fixture(k = 2L)
  expect_error(metis_micom(fx$model, fx$data, "not_a_column", permutations = 10, seed = 1))
})

test_that("metis_micom requires explicit groups when more than two exist", {
  fx <- make_fixture(k = 2L)
  d3 <- fx$data
  set.seed(2); d3$grp <- sample(c("A", "B", "C"), nrow(d3), replace = TRUE)
  expect_error(metis_micom(fx$model, d3, "grp", permutations = 10, seed = 1))
  expect_s3_class(
    metis_micom(fx$model, d3, "grp", group_a = "A", group_b = "B", permutations = 10, seed = 1),
    "metis_micom")
})

test_that("print.metis_micom emits a readable summary", {
  fx <- make_fixture(k = 2L)
  res <- metis_micom(fx$model, fx$data, "grp", permutations = 30, seed = 1)
  expect_output(print(res), "MICOM")
  expect_output(print(res), "Compositional Invariance")
})
