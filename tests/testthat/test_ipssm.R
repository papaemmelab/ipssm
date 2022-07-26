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
  expected_results = list(
    pp347 = list(
      IPSSMscore = 0.31,
      IPSSMcat = "Moderate High",
      IPSSMscore_mean = 0.31,
      IPSSMcat_mean = "Moderate High",
      IPSSMscore_best = 0.3,
      IPSSMcat_best = "Moderate High",
      IPSSMscore_worst = 0.63,
      IPSSMcat_worst = "High"
    ),
    pp564 = list(
      IPSSMscore = 4.55,
      IPSSMcat = "Very High",
      IPSSMscore_mean = 4.55,
      IPSSMcat_mean = "Very High",
      IPSSMscore_best = 4.55,
      IPSSMcat_best = "Very High",
      IPSSMscore_worst = 4.55,
      IPSSMcat_worst = "Very High"
    )
  )

  for (patient in names(expected_results)) {
    for (field in names(expected_results[[patient]])) {
      actual_result <- dd.annot[dd.annot$ID == patient, field]
      expected_result <- expected_results[[patient]][[field]]

      cat(paste("Field: ", field, "\n"))
      cat(paste("\t", expected_result, " (Expected)\n"))
      cat(paste("\t", actual_result, " (Actual)\n"))
      expect_equal(actual_result, expected_result)
    }
  }
})
