#' Number of mutated residual genes
#'
#' Performs the generalized calculation for the number of mutated residual genes per patient, with a cap at 2 in accordance to observed saturation, and adjusting for missing information.
#'
#' @param patientRes a patient \code{data.frame} (with one row and residual genes in columns), or a \code{vector} with residual genes as names attributes. Expected entries are \code{0} (not mutated), \code{1} (mutated), or \code{NA} (missing).
#' @param genesRes a \code{vector} containing the names of the residual genes.
#' @param Nref the average reference value.
#'
#' @return A \code{vector} with the calculated number under the mean, best, and worst scenarios. If no missing data, all scenarios are equal.
#'
#' @export
#'
#' @examples
#' vec <- rep(0,15) ; vec[1] <-1 # 1 mutated genes and no missing data
#' names(vec) <- c("BCOR","BCORL1","CEBPA","ETNK1","GATA2","GNB1","IDH1","NF1","PHF6","PPM1D","PRPF8","PTPN11","SETBP1","STAG2","WT1")
#' calculateNres2(patientRes=vec)
#'
#' vec1 <- vec ; vec1[2:3] <- NA # 1 mutated genes and 2 missing
#' calculateNres2(patientRes=vec1) # scenarios are now different
#'
#' vec2 <- vec1 ; vec2[4] <- 1 # 2 mutated genes and 2 missing
#' calculateNres2(patientRes=vec2) # scenarios are equal, because the number is capped at 2 mutated genes

calculateNres2 <- function(patientRes,
			   genesRes=c("BCOR","BCORL1","CEBPA","ETNK1",
				      "GATA2","GNB1","IDH1","NF1",
				      "PHF6","PPM1D","PRPF8","PTPN11",
				      "SETBP1","STAG2","WT1"),
			   Nref=0.3880) {
   
   # Formatting
   if(any(!genesRes %in% names(patientRes))) {
      stop(paste("Input data does not contain all required residual genes, or is missing their attribute names. The residual genes are:",
		 paste(genesRes,collapse=" ")))
   } else {
      patientRes <- patientRes[names(patientRes) %in% genesRes]
   }

   # Number of missing genes i.e. with NA   
   M <- sum(is.na(patientRes))

   # number of sequenced genes
   S <- sum(!is.na(patientRes))

   # Number of mutated genes within S
   Ns <- sum(as.numeric(patientRes), na.rm=T)

   # Worst scenario: all missing are mutated
   nres2.worst <- min(Ns+M,2)

   # Best scenario: none missing are mutated
   nres2.best <- min(Ns,2)

   # Average scenario: generalized min(Nres,2)
   # i.e. adjustment from min(Ns,2) proportional the average Nref value and to proportion of missing genes M/(S+M)
   # if no genes missing: M=0 --> min(Ns,2) i.e. observed value
   # if all genes missing: S=0,Ns=0 --> Nref i.e. mean imputation
   nres2.mean <- min(Ns,2) + max(((2-Ns)/2)*(M/(S+M))*Nref,0)

   # Return
   res <- c(nres2.mean,nres2.best,nres2.worst)
   names(res) <- c("nRes2mean","nRes2best","nRes2worst")

   return(res)

}
