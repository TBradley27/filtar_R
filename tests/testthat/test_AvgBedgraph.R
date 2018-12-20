library(filtar)

context('test_bedgraph_averaging')

merged_bedgraph = tibble::tibble(chromosome=c(1,1,1,2,2,2,2),
			start=c(50,80,110,30,60,90,150),
			stop=c(60,90,120,40,70,100,160),
			cov1=c(1,2,5,6,3,4,6),
			cov2=c(1,10,3,6,4,5,6)
)

avg_bedgraph = AvgBedgraph(merged_bedgraph)

test_that("avg is a true mean average of 2 coverage columns", {
  expect_equal(avg_bedgraph$avg[1], 1)
  expect_equal(avg_bedgraph$avg[2], 6)
  expect_equal(avg_bedgraph$avg[3], 4)
  expect_equal(avg_bedgraph$avg[4], 6)
  expect_equal(avg_bedgraph$avg[5], 3.5)
  expect_equal(avg_bedgraph$avg[6], 4.5)
  expect_equal(avg_bedgraph$avg[7], 6)
})
