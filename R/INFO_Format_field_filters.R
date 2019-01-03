


#' @title vcfR2vcffilter_FORMAT
#'
#' @param prop.loci.missing Given a loci, how many samples can have missing information before that loci is dropped
#' @export

#-----------------------------------------------------
# Filter VCF based on Format Fields
#------------------------------------------------------
vcffilter_format <- function(vcfRobject = NULL,
                             formatGQ=NULL,
                             formatDP = NULL,
                             formatSP = NULL,
                             prop.loci.missing = NULL, # this is given a loci, how many samples can have missing information before it is dropped
                             biallelic = TRUE,
                             SNPs = TRUE){


  vcf <- vcfRobject # legacy
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  if(biallelic==TRUE){
    vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  }
  if(SNPs==TRUE){
    vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  }
  # store loci objects on info fields
  formatlist <- list()
  if(!is.null(formatGQ)){
    formatlist <- append(formatlist, "gt_GQ")
  }
  if(!is.null(formatDP)){
    formatlist <- append(formatlist, "gt_DP")
  }
  if(!is.null(formatSP)){
    formatlist <- append(formatlist, "gt_SP")
  }

  # extract format list information
  tidyvcf <- vcfR::vcfR2tidy(vcf)
  formatdf <- tidyvcf$gt


  #--------------------------------------------------------
  # filter loci
  #--------------------------------------------------------
  if(!is.null(formatGQ)){
    formatdf <- formatdf %>%
      dplyr::mutate(gt_GQ = ifelse(gt_GQ < formatGQ, "DROP", gt_GQ))
  }
  if(!is.null(formatDP)){
    formatdf <- formatdf %>%
      dplyr::mutate(gt_DP = ifelse(gt_DP < formatDP, "DROP", gt_DP))
  }
  if(!is.null(formatSP)){
    formatdf <- formatdf %>%
      dplyr::mutate(gt_SP = ifelse(gt_SP > formatSP, "DROP", gt_SP)) # this is a fisher-score p-value for likelihood of strand bias. Higher worse
  }

  #--------------------------------------------------------
  # apply filter
  #--------------------------------------------------------

  formatdf <- formatdf[ , colnames(formatdf) %in% c("ChromKey", "POS", "Indiv", formatlist) ]

  # NAs can arise when vcfs are merged and don't have the same INFO field parameters (i.e. cortex versus gatk)
  if(any(is.na(formatdf))){
    warning("Your VCF had NAs that were produced in the tidy2vcf call, which means your Format field parameters are inconsistent. \n You should investigate why this occuring. ")
  }


  formatdf <- formatdf %>%
    dplyr::mutate(excl = apply(., 1, function(x){any(x == "DROP")})) %>%
    dplyr::select(ChromKey, POS, Indiv, excl) %>%
    tidyr::spread(., key="Indiv", value="excl") %>%
    dplyr::select(-c(ChromKey, POS)) %>%
    cbind(FORMAT=rep(FALSE, nrow(.)), .)

  vcf@gt[as.matrix(formatdf)] <- NA

  #--------------------------------------------------------
  # Subset by loci missingness
  #--------------------------------------------------------
  if(!is.null(prop.loci.missing)){
    locimissingness <- apply(vcf@gt, 1, function(x){sum(is.na(x))})
    locimissingnessprop <- locimissingness/(ncol(vcf@gt)-1) # -1 for the format column


    # filter loci with too much missing information
    vcf@gt <- vcf@gt[ locimissingnessprop < prop.loci.missing , ]
  }

  fix <- as.matrix(vcfR::getFIX(vcf, getINFO = T)[ locimissingnessprop < prop.loci.missing ,])
  gt <- as.matrix(vcf@gt)
  meta <- append(vcf@meta, "##Additional Filters provided by NFB filter tools")

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}
