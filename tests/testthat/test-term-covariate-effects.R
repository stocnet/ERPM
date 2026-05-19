### HELPER FUNCTIONS ###

ref_cov_match_choose <- function(n, k) {
  if (n >= k) choose(n, k) else 0
}

ref_cov_match <- function(partition, values, k = 2L, category = NULL, normalized = c("none", "by_group", "global")) {
  normalized <- match.arg(normalized)
  groups <- split(seq_along(partition), as.integer(partition))

  group_values <- vapply(groups, function(indices) {
    group_covariates <- values[indices]
    group_size <- length(indices)

    if (!is.null(category)) {
      category_count <- sum(group_covariates == category, na.rm = TRUE)

      if (k == 1L && normalized == "by_group") {
        return(as.numeric(category_count >= 1L))
      }

      count <- ref_cov_match_choose(category_count, k)
    } else {
      category_counts <- table(group_covariates, useNA = "no")
      count <- sum(vapply(as.integer(category_counts), ref_cov_match_choose, numeric(1), k = k))
    }

    if (normalized == "none") {
      return(count)
    }

    if (normalized == "by_group") {
      denominator <- ref_cov_match_choose(group_size, k)
      if (denominator == 0) return(0)
      return(count / denominator)
    }

    count / group_size
  }, numeric(1))

  sum(group_values)
}

ref_cov_match_GW <- function(partition, values, lambda = 2, category = NULL) {
  groups <- split(seq_along(partition), as.integer(partition))

  values <- vapply(lambda, function(lambda_one) {
    r <- (lambda_one - 1) / lambda_one
    total <- 0

    for (indices in groups) {
      group_covariates <- values[indices]

      if (is.null(category)) {
        category_counts <- table(group_covariates, useNA = "no")
        total <- total + sum(lambda_one * (1 - r^as.integer(category_counts)))
      } else {
        category_count <- sum(group_covariates == category, na.rm = TRUE)
        if (category_count > 0L) {
          total <- total + lambda_one * (1 - r^category_count)
        }
      }
    }

    total
  }, numeric(1))

  unname(values)
}

ref_cov_diff <- function(partition, values, clique_size = 2L, normalized = c("none", "by_group", "global")) {
  normalized <- match.arg(normalized)
  groups <- split(seq_along(partition), as.integer(partition))
  k <- as.integer(clique_size)
  total <- 0

  for (indices in groups) {
    group_size <- length(indices)
    if (group_size < k) next

    group_values <- values[indices]
    combinations <- utils::combn(group_size, k)
    ranges <- apply(combinations, 2L, function(columns) {
      selected <- group_values[columns]
      max(selected) - min(selected)
    })
    group_total <- sum(ranges)

    if (normalized == "none") {
      total <- total + group_total
    } else if (normalized == "by_group") {
      total <- total + group_total / choose(group_size, k)
    } else {
      total <- total + group_total / group_size
    }
  }

  total
}

ref_cov_diff_GW <- function(partition, values, lambda = 2) {
  max_group_size <- max(table(partition))

  if (max_group_size < 2L) {
    return(rep(0, length(lambda)))
  }

  clique_sizes <- 2:max_group_size
  cov_diff_values <- vapply(clique_sizes, function(k) {
    ref_cov_diff(partition, values, clique_size = k, normalized = "none")
  }, numeric(1))

  values <- vapply(lambda, function(lambda_one) {
    weights <- (-1 / lambda_one)^(clique_sizes - 1L)
    sum(weights * cov_diff_values)
  }, numeric(1))

  unname(values)
}

covariate_effect_network <- function(partition, nodes) {
  build_bipartite_from_inputs(partition = partition, nodes = nodes)$network
}

covariate_effect_summary <- function(nw, rhs) {
  formula <- as.formula(call("~", quote(nw), rhs))
  environment(formula) <- environment()

  as.numeric(summary(formula))
}

### TERM: COV_MATCH ###

test_that("cov_match summary matches analytic reference values", {
  withr::local_options(erpm.debug.cov_match_init = FALSE)

  fixtures <- list(
    list(
      partition = c(1, 1, 1, 2, 2, 3, 3, 3),
      nodes = data.frame(
        label = paste0("A", 1:8),
        sex = c("F", "F", "M", "F", "M", "M", "M", "F"),
        grade = c("G1", "G1", "G2", "G1", "G3", "G2", "G2", "G3"),
        stringsAsFactors = FALSE
      )
    ),
    list(
      partition = c(1, 1, 2, 2, 2, 3),
      nodes = data.frame(
        label = paste0("B", 1:6),
        sex = c("F", "M", "F", "F", "F", "M"),
        grade = c("G1", "G2", "G1", "G1", "G2", "G2"),
        stringsAsFactors = FALSE
      )
    )
  )

  cases <- list(
    list(rhs = quote(cov_match("sex", clique_size = 2)), cov = "sex", k = 2L, category = NULL, normalized = "none"),
    list(rhs = quote(cov_match("sex", clique_size = 3)), cov = "sex", k = 3L, category = NULL, normalized = "none"),
    list(rhs = quote(cov_match("sex", clique_size = 2, category = "F")), cov = "sex", k = 2L, category = "F", normalized = "none"),
    list(rhs = quote(cov_match("grade", clique_size = 2, category = "G1")), cov = "grade", k = 2L, category = "G1", normalized = "none"),
    list(rhs = quote(cov_match("grade", clique_size = 2, category = "G999")), cov = "grade", k = 2L, category = "G999", normalized = "none"),
    list(rhs = quote(cov_match("sex", clique_size = 2, normalized = "by_group")), cov = "sex", k = 2L, category = NULL, normalized = "by_group"),
    list(rhs = quote(cov_match("sex", clique_size = 2, normalized = "global")), cov = "sex", k = 2L, category = NULL, normalized = "global"),
    list(rhs = quote(cov_match("sex", clique_size = 1, normalized = "by_group")), cov = "sex", k = 1L, category = NULL, normalized = "by_group"),
    list(rhs = quote(cov_match("grade", clique_size = 1, category = "G1", normalized = "by_group")), cov = "grade", k = 1L, category = "G1", normalized = "by_group")
  )

  for (fixture in fixtures) {
    nw <- covariate_effect_network(fixture$partition, fixture$nodes)

    for (case in cases) {
      observed <- covariate_effect_summary(nw, case$rhs)
      expected <- ref_cov_match(
        fixture$partition,
        fixture$nodes[[case$cov]],
        k = case$k,
        category = case$category,
        normalized = case$normalized
      )

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("attribute:", case$cov, "=", paste(fixture$nodes[[case$cov]], collapse = ", "), "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("cov_match validates arguments", {
  withr::local_options(erpm.debug.cov_match_init = FALSE)

  partition <- c(1, 1, 2, 2)
  nodes <- data.frame(
    label = paste0("A", 1:4),
    sex = c("F", "M", "F", "F"),
    score = c(1, 2, 3, 4),
    stringsAsFactors = FALSE
  )
  nw <- covariate_effect_network(partition, nodes)

  #expect_error(summary(nw ~ cov_match("missing", clique_size = 2)), "missing vertex attribute", fixed = TRUE)
  expect_error(summary(nw ~ cov_match("score", clique_size = 2)), "requires a categorical covariate", fixed = TRUE)
  expect_error(summary(nw ~ cov_match("sex", clique_size = 0)), "'clique_size' must contain finite integers >= 1", fixed = TRUE)
  #expect_error(summary(nw ~ cov_match("sex", clique_size = 1)), "cov_match\(\.\.\., clique_size=1\) with normalized='none'")
  #expect_error(summary(nw ~ cov_match("sex", clique_size = 1, normalized = "global")), "cov_match\(\.\.\., clique_size=1\) with normalized='global'")
  expect_error(summary(nw ~ cov_match("sex", clique_size = 2, normalized = "bad")), "should be one of", fixed = TRUE)
})

### TERM: COV_MATCH_GW ###

test_that("cov_match_GW summary matches analytic reference values", {
  withr::local_options(ERPM.cov_match_GW.debug = FALSE)

  fixtures <- list(
    list(
      partition = c(1, 1, 1, 2, 2, 3, 3, 3),
      nodes = data.frame(
        label = paste0("A", 1:8),
        sex = c("F", "F", "M", "F", "M", "M", "M", "F"),
        grade = c("G1", "G1", "G2", "G1", "G3", "G2", "G2", "G3"),
        stringsAsFactors = FALSE
      )
    ),
    list(
      partition = c(1, 1, 2, 2, 2, 3),
      nodes = data.frame(
        label = paste0("B", 1:6),
        sex = c("F", "M", "F", "F", "F", "M"),
        grade = c("G1", "G2", "G1", "G1", "G2", "G2"),
        stringsAsFactors = FALSE
      )
    )
  )

  cases <- list(
    list(rhs = quote(cov_match_GW("sex", lambda = 2)), cov = "sex", lambda = 2, category = NULL),
    list(rhs = quote(cov_match_GW("sex", lambda = 3)), cov = "sex", lambda = 3, category = NULL),
    list(rhs = quote(cov_match_GW("sex", lambda = c(1.5, 2, 4))), cov = "sex", lambda = c(1.5, 2, 4), category = NULL),
    list(rhs = quote(cov_match_GW("sex", lambda = 2, category = "F")), cov = "sex", lambda = 2, category = "F"),
    list(rhs = quote(cov_match_GW("grade", lambda = 2, category = "G999")), cov = "grade", lambda = 2, category = "G999")
  )

  for (fixture in fixtures) {
    nw <- covariate_effect_network(fixture$partition, fixture$nodes)

    for (case in cases) {
      observed <- covariate_effect_summary(nw, case$rhs)
      expected <- ref_cov_match_GW(fixture$partition, fixture$nodes[[case$cov]], lambda = case$lambda, category = case$category)

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("attribute:", case$cov, "=", paste(fixture$nodes[[case$cov]], collapse = ", "), "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("cov_match_GW validates arguments", {
  withr::local_options(ERPM.cov_match_GW.debug = FALSE)

  nw <- covariate_effect_network(
    c(1, 1, 2, 2),
    data.frame(label = paste0("A", 1:4), sex = c("F", "M", "F", "F"), stringsAsFactors = FALSE)
  )

  #expect_error(summary(nw ~ cov_match_GW("missing", lambda = 2)), "nonexistent attribute", fixed = TRUE)
  expect_error(summary(nw ~ cov_match_GW("sex", lambda = 1)), "'lambda' must be > 1", fixed = TRUE)
  expect_error(summary(nw ~ cov_match_GW("sex", lambda = Inf)), "'lambda' must be > 1", fixed = TRUE)
  expect_error(summary(nw ~ cov_match_GW("sex", lambda = 2, normalized = "bad")), "should be one of", fixed = TRUE)
})

### TERM: COV_DIFF ###

test_that("cov_diff summary matches analytic reference values", {
  withr::local_options(ERPM.cov_diff.debug = FALSE)

  fixtures <- list(
    list(
      partition = c(1, 1, 1, 2, 2, 3, 3, 3),
      nodes = data.frame(label = paste0("A", 1:8), score = c(1, 4, 6, 2, 5, 1, 3, 8), stringsAsFactors = FALSE)
    ),
    list(
      partition = c(1, 1, 2, 2, 2, 3),
      nodes = data.frame(label = paste0("B", 1:6), score = c(2, 7, 1, 3, 9, 4), stringsAsFactors = FALSE)
    )
  )

  cases <- list(
    list(rhs = quote(cov_diff("score", clique_size = 2)), k = 2L, normalized = "none"),
    list(rhs = quote(cov_diff("score", clique_size = 3)), k = 3L, normalized = "none"),
    list(rhs = quote(cov_diff("score", clique_size = 2, normalized = TRUE)), k = 2L, normalized = "by_group"),
    list(rhs = quote(cov_diff("score", clique_size = 2, normalize = "global")), k = 2L, normalized = "global")
  )

  for (fixture in fixtures) {
    nw <- covariate_effect_network(fixture$partition, fixture$nodes)

    for (case in cases) {
      observed <- covariate_effect_summary(nw, case$rhs)
      expected <- ref_cov_diff(fixture$partition, fixture$nodes$score, clique_size = case$k, normalized = case$normalized)

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("attribute:", "score =", paste(fixture$nodes$score, collapse = ", "), "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("cov_diff validates arguments", {
  withr::local_options(ERPM.cov_diff.debug = FALSE)

  nw <- covariate_effect_network(
    c(1, 1, 2, 2),
    data.frame(label = paste0("A", 1:4), score = c(1, 2, 3, 4), bad = c(1, NA, 3, 4), stringsAsFactors = FALSE)
  )

  #expect_error(summary(nw ~ cov_diff("missing", clique_size = 2)), "attribut inexistant", fixed = TRUE)
  expect_error(summary(nw ~ cov_diff("bad", clique_size = 2)), "NA values are not allowed", fixed = TRUE)
  expect_error(summary(nw ~ cov_diff("score", clique_size = 1)), "'clique_size' must be an integer >= 2", fixed = TRUE)
  expect_error(summary(nw ~ cov_diff("score", clique_size = 2, normalize = "bad")), "should be one of", fixed = TRUE)
})

### TERM: COV_DIFF_GW ###

test_that("cov_diff_GW summary matches analytic reference values", {
  withr::local_options(ERPM.cov_diff_GW.debug = FALSE)

  fixtures <- list(
    list(
      partition = c(1, 1, 1, 2, 2, 3, 3, 3),
      nodes = data.frame(label = paste0("A", 1:8), score = c(1, 4, 6, 2, 5, 1, 3, 8), stringsAsFactors = FALSE)
    ),
    list(
      partition = c(1, 1, 2, 2, 2, 3),
      nodes = data.frame(label = paste0("B", 1:6), score = c(2, 7, 1, 3, 9, 4), stringsAsFactors = FALSE)
    )
  )

  cases <- list(
    list(rhs = quote(cov_diff_GW("score", lambda = 1.5)), lambda = 1.5),
    list(rhs = quote(cov_diff_GW("score", lambda = 2)), lambda = 2),
    list(rhs = quote(cov_diff_GW("score", lambda = 4)), lambda = 4),
    list(rhs = quote(cov_diff_GW("score", lambda = c(2, 4))), lambda = c(2, 4))
  )

  for (fixture in fixtures) {
    nw <- covariate_effect_network(fixture$partition, fixture$nodes)

    for (case in cases) {
      observed <- covariate_effect_summary(nw, case$rhs)
      expected <- ref_cov_diff_GW(fixture$partition, fixture$nodes$score, lambda = case$lambda)

      #cat("\n")
      #cat("partition:", paste(fixture$partition, collapse = ", "), "\n")
      #cat("attribute:", "score =", paste(fixture$nodes$score, collapse = ", "), "\n")
      #cat("term:", paste(deparse(case$rhs), collapse = " "), "\n")
      #cat("observed:", paste(observed, collapse = ", "), "\n")
      #cat("expected:", paste(expected, collapse = ", "), "\n")

      expect_equal(observed, expected)
    }
  }
})

test_that("cov_diff_GW validates arguments", {
  withr::local_options(ERPM.cov_diff_GW.debug = FALSE)

  nw <- covariate_effect_network(
    c(1, 1, 2, 2),
    data.frame(label = paste0("A", 1:4), score = c(1, 2, 3, 4), bad = c(1, NA, 3, 4), stringsAsFactors = FALSE)
  )

  #expect_error(summary(nw ~ cov_diff_GW("missing", lambda = 2)), "attribut inexistant", fixed = TRUE)
  expect_error(summary(nw ~ cov_diff_GW("bad", lambda = 2)), "NA values are not allowed", fixed = TRUE)
  expect_error(summary(nw ~ cov_diff_GW("score", lambda = 1)), "all values of 'lambda' must be > 1", fixed = TRUE)
  expect_error(summary(nw ~ cov_diff_GW("score", lambda = Inf)), "'lambda' must be finite", fixed = TRUE)
})
