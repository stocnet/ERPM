ref_order_groupids <- function(partition) {
  ids <- unique(partition)
  match(partition, ids)
}

ref_group_sizes <- function(partition) {
  as.integer(table(partition))
}
