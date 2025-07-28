erpm <- function(formula) {
  if (!inherits(formula, "formula")) stop("The input should be a formula.")

  lhs <- formula[[2]] # nw
  rhs <- formula[[3]] # formula

  # --- 1. Renaming table ---
  rename_map <- c(
    groups        = "b2degrange",
    cov_match     = "nodematch",
    cov_diff      = "absdiff",
    dyadcov       = "edgecov",
    squared_sizes = "squared_sizes"
  )

  # --- 2. Proj1(~ ...) encapsulation list ---
  needs_proj1 <- c("nodematch", "absdiff", "edgecov")

  # --- 3. B(~ ..., form = "nonzero") encapsulation list ---
  needs_B <- c("edgecov")  # Exemple

  # --- Find individual terms ---
  extract_terms <- function(expr) {
    if (is.call(expr) && expr[[1]] == as.name("+")) {
      c(extract_terms(expr[[2]]), extract_terms(expr[[3]]))
    } else {
      list(expr)
    }
  }

  terms <- extract_terms(rhs)
  transformed_terms <- list()

  for (term in terms) {
    if (is.call(term)) {
      original_fun <- as.character(term[[1]])
      args <- as.list(term)[-1]

      # 1. Renaming
      renamed_fun <- if (original_fun %in% names(rename_map)) rename_map[[original_fun]] else original_fun
      new_term <- as.call(c(as.name(renamed_fun), args))

      # 2. B encapsulation
      if (renamed_fun %in% needs_B) {
        new_term <- call("B", call("~", new_term), form = "nonzero")
      }

      # 3. Proj1 encapsulation
      if (renamed_fun %in% needs_proj1) {
        new_term <- call("Proj1", call("~", new_term))
      }

      transformed_terms[[length(transformed_terms) + 1]] <- new_term
    } else {
      message(sprintf(" Unknown ERPM effect : %s", deparse(term)))
      transformed_terms[[length(transformed_terms) + 1]] <- term
    }
  }

  # Reconstruc the formula
  new_rhs <- if (length(transformed_terms) == 1) {
    transformed_terms[[1]]
  } else {
    Reduce(function(x, y) call("+", x, y), transformed_terms)
  }
  form <- list(
    formula = call("~", new_rhs),
    constraints = call("~", as.name("b1part")))
  
  full_formula <- as.call(list(as.name("~"), as.name("nw"), form$formula[[2]]))

  # build the final call
  return(call("ergm", full_formula, constraints = form$constraints))

  # Make the call if wanted
  # eval(call("ergm", full_formula, constraints = form$constraints)) 
}