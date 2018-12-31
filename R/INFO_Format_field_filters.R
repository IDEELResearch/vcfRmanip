
#' @title vcfR2vcffilter_INFO
#'
#' @param infoDP options is based on a percentile cutoff. As a result, if you specify 0.1, the top and bottom 10th read depth percentiles will be excluded
#'
#' @export

## note if there is total agreement in the alignments, i.e. 1/1 where the reads are 0,100 0,100 etc then the
## readposranksum and mqranksum won't be calculated since there is no p-value (there is no read diff)
## this can happen when the reference genome is different from samples and there are snps where
## samples are in total agreement but diff from ref

#-----------------------------------------------------
# Filter VCF based on INFO Fields
#------------------------------------------------------

vcffilter_info <- function(vcfRobject = NULL,
                           infoMQ=NULL,
                           infoQD=NULL,
                           infoSOR=NULL,
                           infoAF = NULL,
                           infoDP = NULL, # this is a percentile cutoff
                           infoFS = NULL,
                           infoMQRankSum = NULL,
                           infoReadPosRankSum = NULL,
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
  infolist <- list()
  if(!is.null(infoMQ)){
    infolist <- append(infolist, "MQ")
  }
  if(!is.null(infoQD)){
    infolist <- append(infolist, "QD")
  }
  if(!is.null(infoSOR)){
    infolist <- append(infolist, "SOR")
  }
  if(!is.null(infoAF)){
    infolist <- append(infolist, "AF")
  }
  if(!is.null(infoDP)){
    infolist <- append(infolist, "DP")
  }
  if(!is.null(infoFS)){
    infolist <- append(infolist, "FS")
  }

  if(!is.null(infoMQRankSum)){
    infolist <- append(infolist, "MQRankSum")
  }
  if(!is.null(infoReadPosRankSum)){
    infolist <- append(infolist, "ReadPosRankSum")
  }

  # extract format list information
  tidyvcf <- vcfR::vcfR2tidy(vcf)
  infodf <- tidyvcf$fix


  #--------------------------------------------------------
  # filter loci
  #--------------------------------------------------------
  if(!is.null(infoMQ)){
    infodf <- infodf %>%
      dplyr::mutate(MQ = ifelse(MQ < infoMQ, "DROP", MQ))
  }
  if(!is.null(infoQD)){
    infodf <- infodf %>%
      dplyr::mutate(QD = ifelse(QD < infoQD, "DROP", QD))
  }
  if(!is.null(infoSOR)){
    infodf <- infodf %>%
      dplyr::mutate(SOR = ifelse(SOR > infoSOR, "DROP", SOR))
  }
  if(!is.null(infoAF)){
    infodf <- infodf %>%
      dplyr::mutate(AF = ifelse(AF< infoAF, "DROP", AF))
  }
  if(!is.null(infoDP)){
    DPpercentile <- quantile(infodf$DP, c(infoDP, 1-infoDP))
    infodf <- infodf %>%
      dplyr::mutate(DP = ifelse(DPpercentile[1] < DP & DP < DPpercentile[2],
                                DP, "DROP"))
  }
  if(!is.null(infoFS)){
    infodf <- infodf %>%
      dplyr::mutate(FS = ifelse(FS > infoFS, "DROP", FS))
  }
  if(!is.null(infoMQRankSum)){
    infodf <- infodf %>%
      dplyr::mutate(MQRankSum = ifelse(MQRankSum < infoMQRankSum, "DROP", MQRankSum))
  }
  if(!is.null(infoReadPosRankSum)){
    infodf <- infodf %>%
      dplyr::mutate(ReadPosRankSum = ifelse(ReadPosRankSum < infoReadPosRankSum, "DROP", ReadPosRankSum))
  }

  #--------------------------------------------------------
  # apply filter
  #--------------------------------------------------------

  infodf <- infodf[ , colnames(infodf) %in% c("CHROM", "POS", infolist) ]
  passedloci <- infodf %>%
    dplyr::mutate(CHROMPOS = paste0(CHROM, POS)) %>%
    dplyr::mutate(incl = apply(., 1, function(x){! any(x == "DROP") })) %>%
    dplyr::select(incl)

  passedloci <- passedloci$incl

  # NAs can arise when vcfs are merged and don't have the same INFO field parameters (i.e. cortex versus gatk)
  if(any(is.na(passedloci))){
    warning("Your VCF had NAs that were produced in the tidy2vcf call, which means your INFO field parameters are inconsistent. \n You should investigate why this occuring. ")
  }

  fix <- as.matrix(vcfR::getFIX(vcf, getINFO = T)[ passedloci ,])
  gt <- as.matrix(vcf@gt[ passedloci , ])
  meta <- append(vcf@meta, "##Additional Filters provided by NFB filter tools")

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}




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
