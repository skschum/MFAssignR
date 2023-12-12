test_that("KMDNoise works", {
  load("../../data/Raw_Neg_ML.rda")
  names(Raw_Neg_ML) <- c("mass", "intensity")
  actual <- KMDNoise(Raw_Neg_ML)
  expect_equal(actual$Noise, 346.0706, tolerance = 1e-4)
})
