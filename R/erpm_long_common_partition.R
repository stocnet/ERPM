# ################################################################################
# # FILE: R/erpm_long_common_partition.R
# ################################################################################
# #' ERPM longitudinal common partition helpers (internal)
# #'
# #' @name erpm_long_common_partition
# #' @note erpm_long_common_partition.R
# #'
# #' @description
# #' Shared helpers to derive group objects and stable signatures from partitions.
# #' Kept because they are useful in both PLS/PLE designs.
# #'
# #' @keywords ERPM ERGM longitudinal internal helpers
# NULL

# #' Partition -> list of groups as integer actor index vectors
# #' @noRd
# .erpm_long_groups_from_partition <- function(p) {
#   p <- as.integer(p)
#   split(seq_along(p), p)
# }

# #' Groups -> stable signatures (character), e.g. "1,3,5"
# #' @noRd
# .erpm_long_group_signatures <- function(groups) {
#   vapply(groups, function(v) paste(sort(as.integer(v)), collapse = ","), character(1))
# }

# #' Build a co-membership matrix (pairs) from a partition.
# #'
# #' Returns an n x n integer matrix with 1 if same group, 0 otherwise, diag=0.
# #'
# #' @noRd
# .erpm_long_comembership_matrix <- function(p) {
#   p <- as.integer(p)
#   n <- length(p)
#   M <- outer(p, p, "==")
#   diag(M) <- FALSE
#   storage.mode(M) <- "integer"
#   M
# }