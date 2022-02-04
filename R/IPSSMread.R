#' Read User Input Data and Validate Format
#'
#' Function that reads the user input data and performs some validation for the presence, format, and values of the input data.
#'
#' @param path.file path to the input data file. Allowed formats are .csv, .tsv, .xls, and .xlsx.
#' @param sheet sheet number for excel format input file. Default is 1.
#'
#' @return A patient input \code{data.frame}.
#'
#' @export
#'
#' @examples
#' path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm")
#' dd <- IPSSMread(path.file)
#' print(dd)
#'
#' path.file2 <- system.file("extdata", "IPSSMexample.xlsx", package = "ipssm")
#' dd2 <- IPSSMread(path.file2, sheet=2) # equivalent 


IPSSMread <- function(path.file, sheet=1) {

   #~~~~~~~~~ Extract File Format
   f <- tail(strsplit(path.file,split=".",fixed=T)[[1]],1)
   if (f %in% c("","","")) {
      stop("Allowed file extensions are .csv, .tsv, .xls, and .xlsx")
   }
   if (f=="csv") {
      d <- read.csv(path.file, stringsAsFactors=F, na.strings=c(""," ","NA","N/A","n/a","na"))
   }
   if (f=="tsv") {
      d <- read.table(path.file, stringsAsFactors=F, header=T, sep="\t", na.strings=c(""," ","NA","N/A","n/a","na"))
   }
   if (f%in%c("xlsx",".xls")) {
      d <- readxl::read_excel(path.file, na = c(""," ","NA","N/A","n/a","na"), sheet=sheet)
      d <- as.data.frame(d)
   }

   #print("Did you make sure that the unit for Hemoglobin (HB) is g per dL. We expect values between 4 and 20 g per dL.")
   print("Did you make sure that the unit for Hemoglobin (HB) is g per dL.")
   print("Did you make sure that the unit for Platelets (PLT) is Giga per L.")
   print("Did you make sure that the unit for Bone Marrow Blast (BM_BLAST) is percentage.")

   #~~~~~~~~~ Variables Validation: | Presence | Format (numerical binary etc) | Range
   check.failed.numerical <- function(x) {
      return( any(!((is.numeric(x)|is.na(x)))) )
   }
   check.failed.binary <- function(x) {
      return( any(!((x%in%c(0,1)|is.na(x)))) )
   }
   check.failed.char.given <- function(x,mychar) {
      return( any(!((x%in%mychar|is.na(x)))) )
   }

   # ~~~ CLINICAL
   clinical.var <- c("HB","PLT","BM_BLAST") # clinical continuous (#"ANC","AGE" # optional continuous)
   if (any(!clinical.var %in% colnames(d))) {
      stop(paste(paste(clinical.var,collapse=" "),"should be columns of your input data"))
   }
   for (cc in clinical.var) {
      if (check.failed.numerical(d[,cc])) {
	 stop(paste(cc,"should have numerical values"))
      }
   }
   for (cc in clinical.var) {
      if (any(d[,cc]<0,na.rm=T)) {
	 stop(paste(cc,"contains negative values. Please re-check your input data!"))
      }
   }
   if (any(d[,"HB"]<4, na.rm=T)) {
      warning("You have HB values smaller than 4 g/dL. This is suspicious.")
   }
   if (any(d[,"HB"]>20, na.rm=T)) {
      warning("You have HB values larger than 20 g/dL. This is suspicious.")
   }

   # ~~~ CYTO
   cyto.var <- c(
		 "del5q","del7_7q","del17_17p","complex", # cytogenetics: 0/1
		 "CYTO_IPSSR" # cytogenetics categorical: Very Good, Good, Intermediate, Poor, Very Poor
   )
   if (any(!cyto.var %in% colnames(d))) {
      stop(paste(paste(cyto.var,collapse=" "),"should be columns of your input data"))
   }
   for (cc in cyto.var[1:4]) {
      if (check.failed.binary(d[,cc])) {
	 stop(paste(cc,"should have binary 0/1 values"))
      }
   }
   if (check.failed.char.given(as.vector(d[,"CYTO_IPSSR"]),c("Very Good", "Good", "Intermediate", "Poor", "Very Poor"))) {
      stop(paste("CYTO_IPSSR","should have values from Very Good, Good, Intermediate, Poor, Very Poor"))
   }

   # TODO: ADD VALIDATION TEST For del5q/CYTO etc


   # ~~~ TP53
   tp53.var <- c(
		 "TP53mut", # TP53 number of mutations categorical: 0, 1, 2 or more 
		 "TP53maxvaf", # continuous %: default 0
		 "TP53loh" # TP53 LOH: 0/1/NA
   )
   if (any(!tp53.var %in% colnames(d))) {
      stop(paste(paste(tp53.var,collapse=" "),"should be columns of your input data"))
   }
   d$TP53mut <- as.character(d$TP53mut)
   if (check.failed.char.given(as.vector(d[,"TP53mut"]),c("0", "1", "2 or more"))) {
      stop(paste("TP53mut","should have values from 0, 1, 2 or more"))
   }
   if (check.failed.numerical(d[,"TP53maxvaf"])) {
      stop(paste("TP53maxvaf","should have numerical values"))
   }
   if (any(d[,"TP53maxvaf"]<0,na.rm=T)) {
      stop(paste("TP53maxvaf","contains negative values. Please re-check your input data!"))
   }
   if (any(d[,"TP53maxvaf"]>1,na.rm=T)) {
      stop(paste("TP53maxvaf","contains values>1. We expect values between 0 and 1."))
   }
   if (check.failed.binary(d[,"TP53loh"])) {
      stop(paste("TP53loh","should have binary 0/1 values"))
   }

   # ~~~ GENES
   gene.var <- c(
		 "MLL_PTD", # all main effect genes below: 0/1/NA
		 "FLT3",
		 "ASXL1", 
		 "CBL",
		 "DNMT3A",
		 "ETV6",
		 "EZH2",
		 "IDH2",
		 "KRAS",
		 "NPM1",
		 "NRAS",
		 "RUNX1",
		 "SF3B1",
		 "SRSF2",
		 "U2AF1",
		 "BCOR", # all residual genes below: 0/1/NA
		 "BCORL1",
		 "CEBPA",
		 "ETNK1",
		 "GATA2",
		 "GNB1",
		 "IDH1",
		 "NF1",
		 "PHF6",
		 "PPM1D",
		 "PRPF8",
		 "PTPN11",
		 "SETBP1",
		 "STAG2",
		 "WT1"
   )
   if (any(!gene.var %in% colnames(d))) {
      stop(paste(paste(gene.var,collapse=" "),"should be columns of your input data"))
   }
   for (cc in gene.var) {
      if (check.failed.binary(d[,cc])) {
	 stop(paste(cc,"should have binary 0/1 values"))
      }
   }

   print("Data successfully imported.")

   return(d)

}
