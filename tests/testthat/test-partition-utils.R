test_that("order_groupids relabels groups by first appearance", {
  expect_equal(order_groupids(c(2, 1, 1, 4, 2)), c(1, 2, 2, 3, 1))
  expect_equal(order_groupids(partition_unordered_labels), ref_order_groupids(partition_unordered_labels))
  expect_equal(order_groupids(partition_isolates), partition_isolates)
  expect_equal(order_groupids(partition_one_group), partition_one_group)
})

test_that("check_sizes validates allowed group sizes and number of groups", {
  expect_true(check_sizes(partition_mixed, sizes.allowed = c(1, 2), numgroups.allowed = 4))
  expect_false(check_sizes(partition_mixed, sizes.allowed = 2, numgroups.allowed = 4))
  expect_false(check_sizes(partition_mixed, sizes.allowed = c(1, 2), numgroups.allowed = 3))

  expect_true(check_sizes(partition_one_group, sizes.allowed = 6, numgroups.allowed = 1))
  expect_false(check_sizes(partition_one_group, sizes.allowed = 1:5, numgroups.allowed = 1))
})

test_that("find_all_partitions returns Bell-number counts for small n", {
  expect_equal(length(find_all_partitions(0)), 0)
  expect_equal(nrow(find_all_partitions(1)), 1)
  expect_equal(ncol(find_all_partitions(1)), 1)
  expect_equal(nrow(find_all_partitions(2)), 2)
  expect_equal(nrow(find_all_partitions(3)), 5)
})

test_that("count_classes counts group-size profiles for n = 3 partitions", {
  all_partitions <- find_all_partitions(3)
  counted <- count_classes(all_partitions)

  expect_equal(sum(counted$counts), nrow(all_partitions))
  expect_equal(nrow(counted$classes), 3)

  profiles <- apply(counted$classes, 1, paste, collapse = ",")
  expected_profiles <- c("3,0,0", "1,1,0", "0,0,1")
  expect_setequal(profiles, expected_profiles)
})

test_that("constrained partition counts match small known values", {
  expect_equal(Stirling2_constraints(4, 2, 2, 2), 3)
  expect_equal(Bell_constraints(4, 2, 2), 3)
  expect_equal(Bell_constraints(3, 2, 2), 0)
})
