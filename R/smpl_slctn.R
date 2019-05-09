#' @title select_samples
#' @description Select samples from a \code{vcfR}
#' @param smplvctr vector; A character vector that corresponds to the samples that should be extracted
#' @export

#-----------------------------------------------------
# Filter VCF based on given positions among the chromosomes
#------------------------------------------------------
select_samples <- function(vcfRobject = NULL,
                           smplvctr = NULL ){

  meta <- append(vcfRobject@meta, paste("##Specific samples were selected using vcfRmanip"))
  fix <- vcfRobject@fix
  gt <- vcfRobject@gt[ , c("FORMAT", smplvctr)]

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = vcfRobject@fix, gt = gt)

  return(newvcfR)
}



#' @title calc_loci_missingness_by_smpl
#' @description Calculate the proportion of the loci that are missing by samples from a \code{vcfR} object
#' @export

#-----------------------------------------------------
# Filter VCF based on given positions among the chromosomes
#------------------------------------------------------
calc_loci_missingness_by_smpl <- function(vcfRobject = NULL){

  gt <- vcfR::extract.gt(vcfRobject, element = "GT")
  ret <- gt %>% as.data.frame(.) %>%
    tidyr::gather(., key = "sample", value = "gt") %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n = n(),
      misscount = sum(is.na(gt)),
      missprop = sum(is.na(gt))/n
      )

  return(ret)

}

























