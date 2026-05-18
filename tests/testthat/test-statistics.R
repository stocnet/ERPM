test_that("group_size and proportion_isolate compute basic partition summaries", {
  expect_equal(group_size(partition_isolates, "avg"), 1)
  expect_equal(group_size(partition_one_group, "avg"), 6)
  expect_equal(group_size(partition_mixed, "avg"), 1.5)
  
  expect_equal(group_size(partition_isolates, "sd"), 0)
  expect_true(is.na(group_size(partition_one_group, "sd")))
  expect_equal(group_size(partition_mixed, "sd"), stats::sd(ref_group_sizes(partition_mixed)))
  
  expect_equal(proportion_isolate(partition_isolates), 1)
  expect_equal(proportion_isolate(partition_one_group), 0)
  expect_equal(proportion_isolate(partition_mixed), 2 / 6)
})

test_that("icc handles constant and varying attributes", {
  expect_true(is.na(icc(partition_isolates, rep(1, 6))))
  expect_true(is.na(icc(partition_one_group, rep(1, 6))))
  expect_true(is.na(icc(partition_mixed, rep(1, 6))))

  mixed_attribute <- c(2, 3, 2, 1, 3, 4)
  expect_true(is.na(icc(partition_isolates, mixed_attribute)))
  expect_true(is.na(icc(partition_one_group, mixed_attribute)))
  expect_equal(icc(partition_mixed, mixed_attribute), 11 / 26)

  increasing_attribute <- 1:6
  expect_true(is.na(icc(partition_isolates, increasing_attribute)))
  expect_true(is.na(icc(partition_one_group, increasing_attribute)))
  expect_equal(icc(partition_mixed, increasing_attribute), 29 / 32)
})

test_that("range_attribute computes sum and average range per group", {
  partition <- c(1, 2, 2, 3, 3, 4, 4, 4, 5)
  attribute <- c(3, 5, 23, 2, 1, 0, 3, 9, 2)
  
  expect_equal(range_attribute(partition, attribute, "sum_pergroup"), 28)
  expect_equal(range_attribute(partition, attribute, "avg_pergroup"), 28 / 5)
  
  expect_equal(range_attribute(partition_isolates, attribute_numeric, "sum_pergroup"), 0)
  expect_equal(range_attribute(partition_isolates, attribute_numeric, "avg_pergroup"), 0)
  
  expect_equal(range_attribute(partition_one_group, attribute_numeric, "sum_pergroup"), 3)
  expect_equal(range_attribute(partition_one_group, attribute_numeric, "avg_pergroup"), 3)
})

test_that("number_categories counts selected and all categories", {
  partition <- c(1, 2, 2, 3, 3, 4, 4, 4, 5)
  attribute <- c(1, 0, 2, 0, 1, 1, 2, 0, 1)
  
  expect_equal(number_categories(partition, attribute, "sum", 1), 4)
  expect_equal(as.numeric(number_categories(partition, attribute, "sum", "all")[1]), 3)
  expect_equal(number_categories(partition, attribute, "avg", 1), 1 + 0 + 1 / 2 + 1 / 3 + 1)
  expect_equal(as.numeric(number_categories(partition, attribute, "avg", "all")[1]), 1 / 2 + 1 / 2 + 1 / 3)
})

test_that("same_pairs counts categorical matches within groups", {
  expect_equal(same_pairs(partition_isolates, attribute_binary, "sum_pergroup"), 0)
  expect_equal(same_pairs(partition_isolates, attribute_binary, "avg_pergroup"), 0)
  expect_equal(same_pairs(partition_isolates, attribute_binary, "sum_perind"), 0)
  expect_equal(same_pairs(partition_isolates, attribute_binary, "avg_perind"), 0)
  
  expect_equal(same_pairs(partition_one_group, attribute_binary, "sum_pergroup"), 6)
  expect_equal(same_pairs(partition_one_group, attribute_binary, "avg_pergroup"), 6)
  expect_equal(same_pairs(partition_one_group, attribute_binary, "sum_perind"), 6)
  expect_equal(same_pairs(partition_one_group, attribute_binary, "avg_perind"), 1)
  
  expect_equal(same_pairs(partition_mixed, attribute_binary, "sum_pergroup"), 1)
  expect_equal(same_pairs(partition_mixed, attribute_binary, "avg_pergroup"), 1 / 4)
  expect_equal(same_pairs(partition_mixed, attribute_binary, "sum_perind"), 2)
  expect_equal(same_pairs(partition_mixed, attribute_binary, "avg_perind"), 2 / 6)
})

test_that("similar_pairs counts numeric near-matches within groups", {
  expect_equal(similar_pairs(partition_isolates, attribute_numeric, "sum_pergroup", 1), 0)
  expect_equal(similar_pairs(partition_isolates, attribute_numeric, "avg_pergroup", 1), 0)
  expect_equal(similar_pairs(partition_isolates, attribute_numeric, "sum_perind", 1), 0)
  expect_equal(similar_pairs(partition_isolates, attribute_numeric, "avg_perind", 1), 0)
  
  expect_equal(similar_pairs(partition_one_group, attribute_numeric, "sum_pergroup", 0), 1)
  expect_equal(similar_pairs(partition_one_group, attribute_numeric, "sum_pergroup", 1), 7)
  expect_equal(similar_pairs(partition_one_group, attribute_numeric, "sum_pergroup", 8), 15)
  
  expect_equal(similar_pairs(partition_mixed, attribute_numeric, "sum_pergroup", 0), 0)
  expect_equal(similar_pairs(partition_mixed, attribute_numeric, "avg_pergroup", 0), 0)
  expect_equal(similar_pairs(partition_mixed, attribute_numeric, "sum_pergroup", 1), 1)
  expect_equal(similar_pairs(partition_mixed, attribute_numeric, "avg_pergroup", 1), 1 / 4)
  expect_equal(similar_pairs(partition_mixed, attribute_numeric, "sum_perind", 1), 2)
  expect_equal(similar_pairs(partition_mixed, attribute_numeric, "avg_perind", 1), 2 / 6)
})

test_that("number_ties counts dyadic ties within groups", {
  expect_equal(number_ties(partition_isolates, dyad_binary_6, "sum_pergroup"), 0)
  expect_equal(number_ties(partition_isolates, dyad_binary_6, "avg_pergroup"), 0)
  expect_equal(number_ties(partition_isolates, dyad_binary_6, "sum_perind"), 0)
  expect_equal(number_ties(partition_isolates, dyad_binary_6, "avg_perind"), 0)
  
  expect_equal(number_ties(partition_one_group, dyad_binary_6, "sum_pergroup"), 2)
  expect_equal(number_ties(partition_one_group, dyad_binary_6, "avg_pergroup"), 2)
  expect_equal(number_ties(partition_one_group, dyad_binary_6, "sum_perind"), 3)
  expect_equal(number_ties(partition_one_group, dyad_binary_6, "avg_perind"), 3 / 6)
})
