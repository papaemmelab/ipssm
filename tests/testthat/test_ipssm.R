library(testthat)
library(ipssm)

print(sessionInfo())

test_that("ipssm examples are computed as expected", {

  # 1) Read and Validate File
  path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm")
  dd <- IPSSMread(path.file)

  # 2) Process User Input Variables into Model Variables
  dd.process <- IPSSMprocess(dd)

  # 3) Calculate IPSS-M
  dd.res <- IPSSMmain(dd.process)

  # 4) Annotate Results
  dd.annot <- IPSSMannotate(dd.res)

  # 5) Check expected results against actual results
  pp347 <- dd.annot[dd.annot$ID == "pp347", ]

  expect_equal(pp347$IPSSMscore, 0.31)
  expect_equal(pp347$IPSSMcat, "Moderate High")
  expect_equal(pp347$IPSSMscore_mean, 0.31)
  expect_equal(pp347$IPSSMcat_mean, "Moderate High")
  expect_equal(pp347$IPSSMscore_best, 0.3)
  expect_equal(pp347$IPSSMcat_best, "Moderate High")
  expect_equal(pp347$IPSSMscore_worst, 0.63)
  expect_equal(pp347$IPSSMcat_worst, "High")

  pp564 <- dd.annot[dd.annot$ID == "pp564", ]

  expect_equal(pp564$IPSSMscore, 4.55)
  expect_equal(pp564$IPSSMcat, "Very High")
  expect_equal(pp564$IPSSMscore_mean, 4.55)
  expect_equal(pp564$IPSSMcat_mean, "Very High")
  expect_equal(pp564$IPSSMscore_best, 4.55)
  expect_equal(pp564$IPSSMcat_best, "Very High")
  expect_equal(pp564$IPSSMscore_worst, 4.55)
  expect_equal(pp564$IPSSMcat_worst, "Very High")
})
