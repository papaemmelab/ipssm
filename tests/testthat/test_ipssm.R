library(ipssm)
print(sessionInfo())


assert_expected_results_were_obtained <- function(annot) {
  expected_results <- list(
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

  cat("\nRunning Tests:\n")
  for (patient in names(expected_results)) {
    for (field in names(expected_results[[patient]])) {
      actual_result <- annot[annot$ID == patient, field]
      expected_result <- expected_results[[patient]][[field]]

      cat(paste("Field: ", field, "\n"))
      cat(paste("\t", expected_result, " (Expected)\n"))
      cat(paste("\t", actual_result, " (Actual)\n"))
      expect_equal(actual_result, expected_result)
    }
  }
}


test_that("Run ipssm workflow by steps", {
  path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm")

  # 1) Read and Validate File
  dd <- IPSSMread(path.file)

  # 2) Process User Input Variables into Model Variables
  dd.process <- IPSSMprocess(dd)

  # 3) Calculate IPSS-M
  dd.res <- IPSSMmain(dd.process)

  # 4) Annotate Results
  dd.annot <- IPSSMannotate(dd.res)

  # Check expected results against actual results
  assert_expected_results_were_obtained(dd.annot)
})


test_that("Run ipssm workflow using wrapper function", {
  path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm")
  annot <- IPSSMwrapper(path.file)
  assert_expected_results_were_obtained(annot)
})
