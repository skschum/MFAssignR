test_that("Isotope filtering works", {
  load("../../data/Raw_Neg_ML.rda")
  expected <- readRDS("test-data/isotopes.rda")

  actual <- IsoFiltR(Raw_Neg_ML)

  expect_equal(actual, expected)
})
