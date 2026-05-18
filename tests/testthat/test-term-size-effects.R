### HELPER FUNCTIONS ###

ref_size_effect_group_sizes <- function(partition) {
  as.integer(table(partition))
}

ref_squared_sizes <- function(partition, sizes = NULL, pow = 2) {
  group_sizes <- ref_size_effect_group_sizes(partition)
  effective_sizes <- if (is.null(sizes) || length(sizes) == 0L) seq_len(length(partition)) else as.integer(sizes)
  selected <- group_sizes[group_sizes %in% effective_sizes]
  
  if (!length(selected)) {
    return(0)
  }
  
  sum(selected^pow)
}

ref_log_factorial_sizes <- function(partition) {
  sum(lgamma(ref_size_effect_group_sizes(partition)))
}

ref_groups <- function(partition, from = 1, to = Inf) {
  group_sizes <- ref_size_effect_group_sizes(partition)
  sum(group_sizes >= from & group_sizes < to)
}

size_effect_network <- function(partition) {
  build_bipartite_from_inputs(partition = partition)$network
}

erpm_groups_summary <- function(partition, rhs) {
  formula <- as.formula(call("~", quote(partition), rhs))
  environment(formula) <- environment()
  translated <- erpm(formula, eval.call = FALSE, verbose = FALSE)
  
  as.numeric(summary(translated[[2]], constraints = ~ b1part))
}

size_effect_summary <- function(nw, rhs) {
  formula <- as.formula(call("~", quote(nw), rhs))
  environment(formula) <- environment()
  
  as.numeric(summary(formula))
}

### TERM: GROUPS ###

test_that("groups summary matches analytic reference values", {
  withr::local_options(ERPM.groups.debug = FALSE)
  
  partitions <- list(
    c(1, 2, 2, 3, 3, 3),
    c(1, 1, 2, 3, 3, 4, 4, 4),
    c(1, 1, 1, 2, 2, 3),
    c(1, 2, 3, 4, 5),
    rep(1, 6)
  )
  
  cases <- list(
    list(rhs = quote(groups), from = 1, to = Inf),
    list(rhs = quote(groups(1)), from = 1, to = 2),
    list(rhs = quote(groups(2)), from = 2, to = 3),
    list(rhs = quote(groups(size = 2)), from = 2, to = 3),
    list(rhs = quote(groups(3)), from = 3, to = 4),
    list(rhs = quote(groups(from = 2, to = 4)), from = 2, to = 4),
    list(rhs = quote(groups(from = 1, to = Inf)), from = 1, to = Inf)
  )
  
  for (partition in partitions) {
    for (case in cases) {
      observed <- erpm_groups_summary(partition, case$rhs)
      expected <- ref_groups(partition, from = case$from, to = case$to)
      
      expect_equal(observed, expected)
    }
  }
})

test_that("groups validates arguments", {
  withr::local_options(ERPM.groups.debug = FALSE)
  
  partition <- c(1, 1, 2, 3)
  
  expect_error(erpm_groups_summary(partition, quote(groups(from = -1, to = 2))), "groups\\(from\\): must be >= 0")
  expect_error(erpm_groups_summary(partition, quote(groups(from = 1.5, to = 3))), "groups\\(from\\): integer required")
  expect_error(erpm_groups_summary(partition, quote(groups(from = 2, to = 1))), "groups\\(from,to\\): requires 'from' < 'to'")
  expect_error(erpm_groups_summary(partition, quote(groups(from = Inf, to = Inf))), "groups\\(from\\): must be a finite integer >= .")
})

### TERM: SQUARED_SIZES ###

test_that("squared_sizes summary matches analytic reference values", {
  withr::local_options(ERPM.squared_sizes.debug = FALSE)
  
  partitions <- list(
    c(1, 2, 2, 3, 3, 3),
    c(1, 1, 2, 3, 3, 4, 4, 4),
    c(1, 1, 1, 2, 2, 3),
    c(1, 2, 3, 4, 5),
    rep(1, 6)
  )
  
  cases <- list(
    list(call = quote(squared_sizes), sizes = NULL, pow = 2),
    list(call = quote(squared_sizes(sizes = c(2, 3, 4))), sizes = c(2, 3, 4), pow = 2),
    list(call = quote(squared_sizes(pow = 3)), sizes = NULL, pow = 3),
    list(call = quote(squared_sizes(sizes = 1:2)), sizes = 1:2, pow = 2),
    list(call = quote(squared_sizes(sizes = 3)), sizes = 3, pow = 2)
  )
  
  for (partition in partitions) {
    nw <- size_effect_network(partition)
    
    for (case in cases) {
      observed <- size_effect_summary(nw, case$call)
      expected <- ref_squared_sizes(partition, sizes = case$sizes, pow = case$pow)
      
      expect_equal(observed, expected)
    }
  }
})

test_that("squared_sizes validates arguments", {
  withr::local_options(ERPM.squared_sizes.debug = FALSE)
  
  nw <- size_effect_network(c(1, 1, 2, 3))
  
  expect_error(summary(nw ~ squared_sizes(size = 2)), "argument 'size' is not supported", fixed = TRUE)
  expect_error(summary(nw ~ squared_sizes(sizes = c(1, NA))), "'sizes' must not contain NA", fixed = TRUE)
  expect_error(summary(nw ~ squared_sizes(sizes = 0)), "'sizes' must be integer >= 1", fixed = TRUE)
  expect_error(summary(nw ~ squared_sizes(sizes = 1.5)), "'sizes' must be integer >= 1", fixed = TRUE)
  expect_error(summary(nw ~ squared_sizes(pow = c(2, 3))), "'pow' must be of length 1", fixed = TRUE)
  expect_error(summary(nw ~ squared_sizes(pow = 0)), "'pow' must be an integer >= 1", fixed = TRUE)
  expect_error(summary(nw ~ squared_sizes(pow = 1.5)), "'pow' must be an integer >= 1", fixed = TRUE)
})


### TERM: LOG_FACTORIAL_SIZES ###

test_that("log_factorial_sizes summary matches analytic reference values", {
  withr::local_options(ERPM.log_factorial_sizes.debug = FALSE)
  
  partitions <- list(
    c(1, 1, 2, 2, 2, 3),
    c(1, 1, 1, 2, 3, 3, 3, 3),
    c(1, 2, 2, 3, 3, 4, 4, 4),
    c(1, 2, 3, 4, 5),
    rep(1, 6),
    c(1, 1, 1, 1, 2, 2, 3, 3, 3, 4),
    c(1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 4)
  )
  
  for (partition in partitions) {
    nw <- size_effect_network(partition)
    
    expect_equal(as.numeric(summary(nw ~ log_factorial_sizes)), ref_log_factorial_sizes(partition))
    expect_equal(as.numeric(summary(nw ~ log_factorial_sizes())), ref_log_factorial_sizes(partition))
  }
})

test_that("log_factorial_sizes validates arguments", {
  withr::local_options(ERPM.log_factorial_sizes.debug = FALSE)
  
  nw <- size_effect_network(c(1, 1, 2, 3))
  
  expect_error(summary(nw ~ log_factorial_sizes(sizes = 2)))
})
