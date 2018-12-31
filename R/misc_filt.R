#' @title vcfR2SubsetChrom
#'
#' @description Produces a subsetted \code{vcfR} based on chromosome
#'
#' @export

vcfR2SubsetChrom <- function(vcfRobject = NULL,
                             chromvect = NULL){

  vcftidy <- vcfR2tidy(vcfRobject)

  passchrom <- vcftidy$fix %>%
    dplyr::select(CHROM) %>%
    dplyr::mutate( keep = ifelse(CHROM %in% chromvect, TRUE, FALSE) ) %>%
    dplyr::select(keep) %>%
    as.matrix(.)



  meta <- append(vcfRobject@meta, paste("##Chromsomes were subsetted"))
  fix <- vcfRobject@fix[ passchrom, ]
  gt <- vcfRobject@gt[ passchrom, ]


  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  newvcfR

}




#' @title vcfR2vcffilter_CHROMPOS
#'
#' @description Produces a subsetted \code{vcfR} with specific loci excluded as determined \code{chromposdf}.
#'
#' @param chromposdf an dataframe that contains loci to be excluded and has been produced by \code{GFF2VariantAnnotation_Short}
#'
#'
#' @export

#-----------------------------------------------------
# Filter VCF based on given positions among the chromosomes
#------------------------------------------------------
vcffilter_CHROMPOS <- function(vcfRobject = NULL,
                               chromposdf = NULL # this is expected to be from GFF2VariantAnnotation_Short
){
  chromposdflist = split(chromposdf, f=factor(chromposdf$GeneID))
  chromposlong <- do.call("rbind", parallel::mclapply(chromposdflist, function(dat){

    i <- length(seq(from=dat$start, to=dat$end))

    temp <- tibble::tibble(CHROM = rep(dat$seqname, i),
                           POS = seq(from=dat$start, to=dat$end),
                           remove = rep("Y", i)
    )
    return(temp)

  }))


  vcftidy <- vcfR2tidy(vcfRobject)

  passloci <- vcftidy$fix %>%
    dplyr::select(CHROM, POS) %>%
    dplyr::left_join(x=., y=chromposlong, by=c("CHROM", "POS")) %>%
    dplyr::mutate(keep = ifelse(is.na(remove), TRUE, FALSE)) %>%
    dplyr::select(keep) %>%
    as.matrix(.)

  meta <- append(vcfRobject@meta, paste("##Additional Filters for excluding hypervariable regions determined by user"))
  fix <- vcfRobject@fix[ passloci, ]
  gt <- vcfRobject@gt[ passloci, ]


  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  newvcfR

}



#' @title vcfR2vcffilter_SEGSITES
#'
#'
#' @export
#-----------------------------------------------------
# Read vcfR object to segregated sites
#------------------------------------------------------

vcfR2segsites <- function(vcfRobject = NULL, err = 0.025){
  vcf <- vcfRobject # legacy
  if(!identical(vcf, vcf[vcfR::is.biallelic(vcf)])){
    stop("VCF must be biallelic for this to work.")
  }

  ad <- vcfR::extract.gt(vcf, element = "AD")
  refad <- masplit(ad, record=1, sort=0, decreasing = 0)
  altad <- masplit(ad, record=2, sort=0, decreasing = 0)

  NRAF <- altad/(altad + refad)
  NRAF[is.nan(NRAF)] <- NA # the 0,0 are returning Nans

  segsites <- apply(NRAF, 1, function(x){
    if(all(is.na(x))){
      return(FALSE)
    } else {
      return( (max(x, na.rm=T) - min(x, na.rm=T)) > err )

    }

  })

  vcfRobject@gt <- vcfRobject@gt[segsites,]

  fix <- as.matrix(vcfR::getFIX(vcfRobject, getINFO = T)[segsites , ])
  gt <- as.matrix(vcfRobject@gt)
  meta <- append(vcfRobject@meta, paste("##Additional Filters for segregating sites, such that within-sample AF must vary more than:", err))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}



