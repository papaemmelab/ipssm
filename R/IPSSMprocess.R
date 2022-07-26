#' Processing of user-defined clinical and molecular variables
#'
#' Process clinical and molecular user-defined variables to model-based variables.
#'
#' @param patientInput a patient \code{data.frame}, one patient per row, and variables as columns.
#' @param genesRes a \code{vector} containing the names of the residual genes.
#' @param maxvafloh maximum variant allele frequency (vaf) threshold to determine theres Loss of Heterozigosity (LOH).
#' @param Nref the average reference value for min(Nres,2) where Nres is the number of mutated residual genes.
#'
#' @return A processed patient \code{data.frame}, same number of rows/patients as in \code{patientInput}, and with the processed variables as additional columns.
#'
#' @export
#'
#' @examples
#' dd <- read.csv(system.file("extdata", "IPSSMexample.csv", package = "ipssm"), header = TRUE)
#' dd.process <- IPSSMprocess(patientInput = dd)
#' print(dd.process)
#'
IPSSMprocess <- function(patientInput,
                         genesRes = c(
                           "BCOR", "BCORL1", "CEBPA", "ETNK1",
                           "GATA2", "GNB1", "IDH1", "NF1",
                           "PHF6", "PPM1D", "PRPF8", "PTPN11",
                           "SETBP1", "STAG2", "WT1"
                         ),
                         maxvafloh = 0.55,
                         Nref = 0.3880) {
  patientProcess <- patientInput

  cat("Pre-processing your input data...\n")

  # Construction of SF3B1 features i.e SF3B1_5q | SF3B1_alpha
  patientProcess$SF3B1_5q <- NA
  patientProcess$SF3B1_5q[which(patientProcess$SF3B1 == 0)] <- 0 # SF3B1 should be mutated for SF3B1_5q=1
  patientProcess$SF3B1_5q[which(patientProcess$del5q == 0 | patientProcess$del7_7q == 1 | patientProcess$complex == 1)] <- 0 # if not isolated 5q then SF3B1_5q = 0 in any case
  patientProcess$SF3B1_5q[which(patientProcess$SF3B1 == 1 & patientProcess$del5q == 1 & patientProcess$del7_7q == 0 & patientProcess$complex == 0)] <- 1

  patientProcess$SF3B1_alpha <- NA
  patientProcess$SF3B1_alpha[which(patientProcess$SF3B1 == 0)] <- 0
  patientProcess$SF3B1_alpha[which(patientProcess$SF3B1_5q == 1 | # if one of those variables = 1 then SF3B1_alpha = 0 in any case
    patientProcess$SRSF2 == 1 | patientProcess$STAG2 == 1 |
    patientProcess$BCOR == 1 | patientProcess$BCORL1 == 1 |
    patientProcess$RUNX1 == 1 | patientProcess$NRAS == 1)] <- 0
  patientProcess$SF3B1_alpha[which(patientProcess$SF3B1 == 1 & patientProcess$SF3B1_5q == 0 &
    patientProcess$SRSF2 == 0 & patientProcess$STAG2 == 0 &
    patientProcess$BCOR == 0 & patientProcess$BCORL1 == 0 &
    patientProcess$RUNX1 == 0 & patientProcess$NRAS == 0)] <- 1

  # Construction of TP53multi feature
  patientProcess$TP53loh[which(patientProcess$TP53maxvaf > maxvafloh | patientProcess$del17_17p == 1)] <- 1 # inferred LOH
  patientProcess$TP53multi <- NA
  patientProcess$TP53multi[which((patientProcess$TP53mut %in% c("0")) |
    (patientProcess$TP53mut %in% c("1") & patientProcess$TP53loh == 0))] <- 0 # mono-allelic
  patientProcess$TP53multi[which((patientProcess$TP53mut %in% c("2 or more")) |
    (patientProcess$TP53mut %in% c("1", "2 or more") & patientProcess$TP53loh == 1))] <- 1 # multi-hit

  # Transformation of clinical variables
  patientProcess$HB1 <- pmin(pmax(4, patientProcess$HB), 20) # within range 4-20 (very permissive)
  patientProcess$BLAST5 <- pmin(patientProcess$BM_BLAST, 20) / 5 # within range 0-20 i.e. MDS
  patientProcess$TRANSF_PLT100 <- pmin(patientProcess$PLT, 250) / 100 # ceiling at 250

  # Cytogenetics as a numerical vector
  patientProcess$CYTOVEC <- car::recode(patientProcess$CYTO_IPSSR,
    "'Very Good'=0; 'Good'=1; 'Intermediate'=2; 'Poor'=3; 'Very Poor'=4",
    as.factor = FALSE
  )

  # Calculation of number of residual mutations Nres2
  # with generalization to allow some missing genes in the list of residual genes
  res2 <- as.data.frame(t(apply(
    patientProcess[, genesRes, drop = F], 1,
    function(x) calculateNres2(patientRes = x, genesRes = genesRes, Nref = Nref)
  )))
  patientProcess <- cbind(patientProcess, res2)

  cat("Success\n")

  return(patientProcess)
}
