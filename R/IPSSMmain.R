#' Molecular International Prognostic Scoring System IPSS-M
#'
#' Main function of the IPSS-M model. It applies the IPSS-M on the patient processed clinical and molecular variables. It calculates the IPSS-M risk score and risk categories under the best, mean, and worst scenarios if some input data are missing. If no missing data, all scenarios are equal.
#'
#' @param patientProcess a patient processed \code{data.frame}, as the result of \code{IPSSMprocess} function.
#' @param betaValues a covariate vector of the model weights. Should have name attributes.
#' @param meanValues vector of average values for each covariate. Should have the same names attributes as in \code{betaValues}.
#' @param bestValues vector of best values (leading to minimal risk) for each covariate (nRes2 not needed as already taken care of in \code{IPSSMprocess}).
#' @param worstValues vector of worst values (leading to maximal risk) for each covariate (nRes2 not needed as already taken care of in \code{IPSSMprocess}).
#' @param rounding should the raw IPSS-M risk score be rounded. Default is TRUE.
#' @param rounding.digits number of digits for the rounding. Default is 2.
#' @param risk.cutpoints cutpoints to be applied to the IPSS-M risk score to create risk categories.
#' @param risk.cat names of the IPSS-M risk categories
#'
#' @return A result patient \code{data.frame}, same number of rows/patients as in \code{patientProcess}, with additonal columns labelled \code{IPSSMscore_best}, \code{IPSSMscore_mean}, \code{IPSSMscore_worst}, \code{IPSSMcat_best}, \code{IPSSMcat_mean}, \code{IPSSMcat_worst}, corresponding to the IPSS-M risk score and category in the best, mean, and worst scenarios.

#' @export IPSSMmain
#'
#' @examples
#' dd <- read.csv(system.file("extdata", "IPSSMexample.csv", package = "ipssm"),header=T)
#' dd.process <- IPSSMprocess(patientInput=dd)
#' dd.res <- IPSSMmain(patientProcess=dd.process)
#' print(dd.res[,c(1,grep("IPSSM",colnames(dd.res)))])

IPSSMmain <- function(patientProcess,
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
		      risk.cat=c("Very Low","Low","Moderate Low","Moderate High","High","Very High")
		      ) {

   cat("Calculating IPSS-M ...\n")
   meanValuesShort <- meanValues[names(meanValues)!="nRes2"]
   lres <- lapply(1:nrow(patientProcess), function(i) {
		     x <- patientProcess[i,]
		     # Create the patientScenarioValues under the best, mean, and worst scenarios
		     # Then Calculate the IPSS-M risk score and category for those scenarios 
		     # Best Scenario
		     xbest <- FillScenario(x, imputeValues=bestValues, scenario="best")
		     rbest <- IPSSMmodel(patientValues=xbest, betaValues=betaValues,meanValues=meanValues,
					 rounding=rounding,rounding.digits=rounding.digits,
					 risk.cutpoints=risk.cutpoints,risk.cat=risk.cat
		     )

		     # Mean Scenario
		     xmean <- FillScenario(x, imputeValues=meanValuesShort, scenario="mean")
		     rmean <- IPSSMmodel(patientValues=xmean, betaValues=betaValues,meanValues=meanValues,
					 rounding=rounding,rounding.digits=rounding.digits,
					 risk.cutpoints=risk.cutpoints,risk.cat=risk.cat
		     )
		     # Worst Scenario
		     xworst <- FillScenario(x, imputeValues=worstValues, scenario="worst")
		     rworst <- IPSSMmodel(patientValues=xworst, betaValues=betaValues,meanValues=meanValues,
					  rounding=rounding,rounding.digits=rounding.digits,
					  risk.cutpoints=risk.cutpoints,risk.cat=risk.cat
		     )
		     resscore <- c(rbest$IPSSMscore, rmean$IPSSMscore, rworst$IPSSMscore)
		     rescat <- c(rbest$IPSSMcat, rmean$IPSSMcat, rworst$IPSSMcat)
		     return(list(resscore,rescat))
		      })

   dscore <-  as.data.frame(do.call("rbind",lapply(lres,function(x) x[[1]])))
   names(dscore) <- c("IPSSMscore_best","IPSSMscore_mean","IPSSMscore_worst")
   dcat <-  as.data.frame(do.call("rbind",lapply(lres,function(x) x[[2]])))
   names(dcat) <- c("IPSSMcat_best","IPSSMcat_mean","IPSSMcat_worst")

   patientResult <- cbind(patientProcess, dscore, dcat)

   cat("Sucess\n")

   return(patientResult)
}



IPSSMmodel <- function(patientValues,
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
		       rounding=TRUE,
		       rounding.digits=2,
		       risk.cutpoints=c(-1.5,-0.5,0,0.5,1.5),
		       risk.cat=c("Very Low","Low","Moderate Low","Moderate High","High","Very High")
		       ) {

   if (!"names"%in%names(attributes(betaValues)) |
       !"names"%in%names(attributes(meanValues))) { 
      stop("vectors betaValues and meanValues should have names attributes")
   }

   if (!identical(names(betaValues),names(meanValues))) {
      stop("vectors betaValues and meanValues should have the same length and names attributes")
   }

   if (any(!names(betaValues)%in%names(patientValues))) {
      stop(paste("patient dataframe patientValues should contains all attributes from betaValues, i.e.",paste(names(betaValues),collapse=" ")))
   }

   # Be On The Safe Side, columns reordering
   targetvar <- names(betaValues)
   meanValues <- meanValues[match(targetvar, names(meanValues))]
   patientValues <- patientValues[,targetvar,drop=F]

   # relative risk contribution of each variable, with scaling for interpretability
   # --> Hazard Ratio from the Average patient
   # --> Log2 HR scale
   contributions <- ((patientValues - meanValues) * betaValues)/log(2) # log2 is just a scaling factor

   # IPSS-M risk score
   IPSSMscore <- sum(contributions)

   if (rounding) {
      IPSSMscore <- round(IPSSMscore, rounding.digits)
   }

   # IPSS-M risk categories
   if (length(risk.cat)!=(length(risk.cutpoints)+1)) {
      stop("you should have 1 more risk category (risk.cat) that the number of cutpoints (risk.cutpoints)")
   }
   IPSSMcat <- risk.cat[cut(IPSSMscore, breaks=c(-Inf, risk.cutpoints, +Inf), labels=FALSE)]

   return(list(IPSSMscore=IPSSMscore, IPSSMcat=IPSSMcat, contributions=contributions))

}


FillScenario <- function(patientProcessValues,
			 imputeValues=c("HB1"=9.87, "TRANSF_PLT100"=1.41, "BLAST5"=0.922, "CYTOVEC"=1.39,
					"TP53multi"=0.0710, "FLT3"=0.0108, "MLL_PTD"=0.0247,
					"SF3B1_5q"=0.0166, "NPM1"=0.0112, "RUNX1"=0.1260, "NRAS"=0.0362, "ETV6"=0.0216,
					"IDH2"=0.0429, "CBL"=0.0473, "EZH2"=0.0588, "U2AF1"=0.0866, "SRSF2"=0.1580, "DNMT3A"=0.1610,
					"ASXL1"=0.2520, "KRAS"=0.0271, "SF3B1_alpha"=0.1860),
			 scenario="mean") {


   if (any(!names(imputeValues)%in%names(patientProcessValues))) {
      stop(paste("patient processed dataframe should contains",paste(names(imputeValues),collapse=" ")))
   }

   # Impute Missing Values based in the provided imputeValues vector
   patientScenarioValues <- patientProcessValues[names(imputeValues)]
   ina <- which(is.na(patientScenarioValues))

   if (length(ina)>0) {
      patientScenarioValues[ina] <- imputeValues[ina]
   }

   # Add the nRes2 that corresponds to the corresponding scenario
   nres2 <- patientProcessValues[paste0("nRes2",scenario)]
   names(nres2) <- "nRes2"
   patientScenarioValues <- cbind(patientScenarioValues, nres2)

   return(patientScenarioValues)

}
