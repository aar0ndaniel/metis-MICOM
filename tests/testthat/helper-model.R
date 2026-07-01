# Shared fixtures for tests: a small two-group seminr model on seminr::mobi.
make_fixture <- function(k = 2L, seed = 1L) {
  data <- seminr::mobi
  set.seed(seed)
  data$grp <- sample(c("A", "B"), nrow(data), replace = TRUE)
  cons <- list(
    seminr::composite("Image",       seminr::multi_items("IMAG", 1:5)),
    seminr::composite("Expectation", seminr::multi_items("CUEX", 1:3)),
    seminr::composite("Quality",     seminr::multi_items("PERQ", 1:7))
  )[seq_len(k)]
  mm <- do.call(seminr::constructs, cons)
  nm <- c("Image", "Expectation", "Quality")[seq_len(k)]
  paths <- lapply(seq_len(k - 1L), function(i) seminr::paths(from = nm[i], to = nm[i + 1L]))
  sm <- do.call(seminr::relationships, paths)
  model <- suppressMessages(seminr::estimate_pls(data, mm, sm))
  list(model = model, data = data, constructs = nm)
}

reflect_indicators <- function(data, items) {
  for (it in items) data[[it]] <- (max(data[[it]]) + min(data[[it]])) - data[[it]]
  data
}
