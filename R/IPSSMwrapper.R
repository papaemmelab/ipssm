#' IPSS-M Wrapper
#'
#' Function that reads, processes, calculates and annotates to output IPSS-M results
#'
#' @param path.file path to the input data file. Allowed formats are .csv, .tsv, .xls, and .xlsx.
#' @param sheet sheet number for excel format input file. Default is 1.
#' @param genesRes a \code{vector} containing the names of the residual genes. Default is from the IPSS-M publised model.
#' @param Nref the average reference value for min(Nres,2) where Nres is the number of mutated residual genes.
#' @param betaValues a covariate vector of the model weights. Should have name attributes.
#' @param meanValues vector of average values for each covariate. Should have the same names attributes as in \code{betaValues}.
#' @param bestValues vector of best values (leading to minimal risk) for each covariate (nRes2 not needed as already taken care of in \code{IPSSMprocess}).
#' @param worstValues vector of worst values (leading to maximal risk) for each covariate (nRes2 not needed as already taken care of in \code{IPSSMprocess}).
#' @param rounding should the raw IPSS-M risk score be rounded. Default is TRUE.
#' @param rounding.digits number of digits for the rounding. Default is 2.
#' @param risk.cutpoints cutpoints to be applied to the IPSS-M risk score to create risk categories.
#' @param risk.cat names of the IPSS-M risk categories.
#' @param range.max threshold for the allowed maximal range in IPSS-M risk (worst - best) to confidently report the mean scenario as the main result. Default 1.
#'
#' @return A patient annoated results \code{data.frame}.
#'
#' @export
#'
#' @examples
#' path.file <- system.file("extdata", "IPSSMexample.csv", package = "ipssm")
#' ddres <- IPSSMwrapper(path.file)
#' print(ddres)


IPSSMwrapper <- function(path.file, 
			 sheet=1,
			 genesRes=c("BCOR","BCORL1","CEBPA","ETNK1",
				    "GATA2","GNB1","IDH1","NF1",
				    "PHF6","PPM1D","PRPF8","PTPN11",
				    "SETBP1","STAG2","WT1"),
			 Nref=0.3880,
			 betaValues = c("HB1"=-0.171, "TRANSF_PLT100"=-0.222, "BLAST5"=0.352, "CYTOVEC"=0.287, 
					"TP53multi"=1.1800, "FLT3"=0.7980, "MLL_PTD"=0.7980, "SF3B1_5q"=0.5040, "NPM1"=0.4300, 
					"RUNX1"=0.4230, "NRAS"=0.4170, "ETV6"=0.3910, 
					"IDH2"=0.3790, "CBL"=0.2950, "EZH2"=0.2700, "U2AF1"=0.2470, "SRSF2"=0.2390, "DNMT3A"=0.2210, 
					"ASXL1"=0.2130, "KRAS"=0.2020, "SF3B1_alpha"=-0.0794,
					"nRes2"=0.2310),
			 meanValues = c("HB1"=9.87, "TRANSF_PLT100"=1.41, "BLAST5"=0.922, "CYTOVEC"=1.39, 
					"TP53multi"=0.0710, "FLT3"=0.0108, "MLL_PTD"=0.0247, 
					"SF3B1_5q"=0.0166, "NPM1"=0.0112, "RUNX1"=0.1260, "NRAS"=0.0362, "ETV6"=0.0216, 
					"IDH2"=0.0429, "CBL"=0.0473, "EZH2"=0.0588, "U2AF1"=0.0866, "SRSF2"=0.1580, "DNMT3A"=0.1610, 
					"ASXL1"=0.2520, "KRAS"=0.0271, "SF3B1_alpha"=0.1860,
					"nRes2"=0.3880),
			 bestValues = c("HB1"=20, "TRANSF_PLT100"=2.5, "BLAST5"=0, "CYTOVEC"=0, 
					"TP53multi"=0, "FLT3"=0, "MLL_PTD"=0, 
					"SF3B1_5q"=0, "NPM1"=0, "RUNX1"=0, "NRAS"=0, "ETV6"=0, 
					"IDH2"=0, "CBL"=0, "EZH2"=0, "U2AF1"=0, "SRSF2"=0, "DNMT3A"=0, 
					"ASXL1"=0, "KRAS"=0, "SF3B1_alpha"=1),
			 worstValues = c("HB1"=4, "TRANSF_PLT100"=0, "BLAST5"=4, "CYTOVEC"=4,
					 "TP53multi"=1, "FLT3"=1, "MLL_PTD"=1,
					 "SF3B1_5q"=1, "NPM1"=1, "RUNX1"=1, "NRAS"=1, "ETV6"=1,
					 "IDH2"=1, "CBL"=1, "EZH2"=1, "U2AF1"=1, "SRSF2"=1, "DNMT3A"=1,
					 "ASXL1"=1, "KRAS"=1, "SF3B1_alpha"=0),
			 rounding=TRUE,
			 rounding.digits=2,
			 risk.cutpoints=c(-1.5,-0.5,0,0.5,1.5),
			 risk.cat=c("Very Low","Low","Moderate Low","Moderate High","High","Very High"),
			 range.max=1
			 ) {


   dd <- IPSSMread(path.file)   
   dd.process <- IPSSMprocess(patientInput=dd)
   dd.res <- IPSSMmain(patientProcess=dd.process)
   dd.annot <- IPSSMannotate(patientResult=dd.res)

   return(dd)

}
