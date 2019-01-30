

#' @title rehh liftover for vcfR
#'
#' @description Converts a vcfR object to a "transposed" haplotype object for use in other packages (i.e. \{rehh})
#'
#'

vcfR2thap <- function(vcfRobj){
  # -----------------------------------------------------
  # Error catch
  #------------------------------------------------------
  if(!identical(vcfRobj, vcfRobj[vcfR::is.biallelic(vcfRobj)])){
    stop("VCF must be biallelic for this to work.")
  }

  refalt <- vcfR::getFIX(vcfRobj, getINFO = F)[,4:5] # based on vcf spec

  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\|")){ # grab first site
    stop("This tool does not support phased vcfs. See vcfR extract.haps for this functionality")
  } else if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\/")) {
    if(length( stringr::str_split(vcfR::extract.gt(vcfRobj)[1,2], "\\/", simplify = T)) > 2){
      stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by this tool")
    } else{
      gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
      gtmatrix[gtmatrix == "0/0"] <- 0
      gtmatrix[gtmatrix == "0/1"] <- NA
      gtmatrix[gtmatrix == "1/1"] <- 1
      gtmatrix[is.na(gtmatrix)] <- NA
    }
  } else {

    gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=T)
    gtmatrix[gtmatrix == 0] <- 0
    gtmatrix[gtmatrix == 1] <- 1
    gtmatrix[is.na(gtmatrix)] <- NA
  }

  thap <- matrix(NA, nrow(gtmatrix), ncol(gtmatrix))
  for(i in 1:nrow(gtmatrix)){
    for(j in 1:ncol(gtmatrix)){
      if(is.na(gtmatrix[i,j])){
        thap[i,j] <- NA
      } else if(gtmatrix[i,j] == 0){
        thap[i,j] <- refalt[i,1]
      } else if(gtmatrix[i,j] == 1){
        thap[i,j] <- refalt[i,2]
      }
    }
  }

  return(thap)

}
