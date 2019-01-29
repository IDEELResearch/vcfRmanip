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



#' @title vcfR2segsites_wsaf
#'
#' @description Read vcfR object to segregated sites based on wsaf
#' @export


vcfR2segsites_wsaf <- function(vcfRobject = NULL, err = 0.025){
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



#' @title vcfR2segsites_gt
#'
#' @description Read vcfR object to segregated sites based on GT
#' @export

vcfR2segsites_gt <- function(vcfRobj){

  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\|")){ # grab first site
    stop("This tool does not support phased vcfs")
  } else if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\/")) {
    if(length( stringr::str_split(vcfR::extract.gt(vcfRobj)[1,2], "\\/", simplify = T)) > 2){
      stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by this tool")
    } else{
      gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
      gtmatrix[gtmatrix == "0/0"] <- 0
      gtmatrix[gtmatrix == "0/1"] <- 1
      gtmatrix[gtmatrix == "1/1"] <- 2
      gtmatrix[is.na(gtmatrix)] <- NA
      gtmatrix <- apply(gtmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
    }
  } else {

    gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=T)
    gtmatrix[gtmatrix == 0] <- 0
    gtmatrix[gtmatrix == 1] <- 2
    gtmatrix[is.na(gtmatrix)] <- NA
  }


  segsites <- apply(gtmatrix, 1, function(x){
    homoref <- all(x[!is.na(x)] == 0)
    homoalt <- all(x[!is.na(x)] == 2)
    het <- all(x[!is.na(x)] == 1)

    if(sum(c(homoref, homoalt, het)) == 1){
      return(F)
    } else{
      return(T)
    }
  }
  )

  vcfRobj@gt <- vcfRobj@gt[segsites,]

  fix <- as.matrix(vcfR::getFIX(vcfRobj, getINFO = T)[segsites,])
  gt <- as.matrix(vcfRobj@gt)
  meta <- append(vcfRobj@meta, paste("##Additional Filters for segregating sites, such that GT call must be segregating within samples"))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)
}


#' @title vcfR2removesingletons_gt
#'
#' @description Read vcfR object and remove all loci that have are singletons
#' @export
#'

vcfR2removesingletons_gt <- function(vcfRobj){

  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\|")){ # grab first site
    stop("This tool does not support phased vcfs")
  } else if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\/")) {
    if(length( stringr::str_split(vcfR::extract.gt(vcfRobj)[1,2], "\\/", simplify = T)) > 2){
      stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by this tool")
    } else{
      gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
      gtmatrix[gtmatrix == "0/0"] <- 0
      gtmatrix[gtmatrix == "0/1"] <- 1
      gtmatrix[gtmatrix == "1/1"] <- 2
      gtmatrix[is.na(gtmatrix)] <- NA
      gtmatrix <- apply(gtmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
    }
  } else {

    gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=T)
    gtmatrix[gtmatrix == 0] <- 0
    gtmatrix[gtmatrix == 1] <- 2
    gtmatrix[is.na(gtmatrix)] <- NA
  }



  singleton <- apply(gtmatrix, 1, function(x){

    x <- factor(x, levels = c(0,1,2))
    return(
      all(as.data.frame(table(x))$Freq %in% c(0, 1, (length(x)-1))) # singleton must have 1, 0, and then rest #
    )
    }

  )

  vcfRobj@gt <- vcfRobj@gt[!singleton,]

  fix <- as.matrix(vcfR::getFIX(vcfRobj, getINFO = T)[!singleton,])
  gt <- as.matrix(vcfRobj@gt)
  meta <- append(vcfRobj@meta, paste("##Additional Filters for segregating sites, such that GT call must be segregating within samples"))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)
}



#' @title vcfR2removesingletons_gt
#'
#' @description Read vcfR object and convert genotype calls to 0,1,2
#' @export
#'


gtmat012 <- function(vcfRobj){

  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\|")){ # grab first site
    stop("This tool does not support phased vcfs")
  } else if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\/")) {
    if(length( stringr::str_split(vcfR::extract.gt(vcfRobj)[1,2], "\\/", simplify = T)) > 2){
      stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by this tool")
    } else{
      gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
      gtmatrix[gtmatrix == "0/0"] <- 0
      gtmatrix[gtmatrix == "0/1"] <- 1
      gtmatrix[gtmatrix == "1/1"] <- 2
      gtmatrix[is.na(gtmatrix)] <- NA
      gtmatrix <- apply(gtmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
    }
  } else {

    gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=T)
    gtmatrix[gtmatrix == 0] <- 0
    gtmatrix[gtmatrix == 1] <- 2
    gtmatrix[is.na(gtmatrix)] <- NA
  }

  return(gtmatrix)

}









#' @title Remove heterozygote calls and make monoclonal
#'
#' @description Read vcfR object and code all het sites in GT as missing
#' @export

vcfR2nohetsmonoclonal_gt <- function(vcfRobj){

  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\|")){ # grab first site
    stop("This tool does not support phased vcfs")
  } else if(stringr::str_detect(vcfR::extract.gt(vcfRobj)[1,2], "\\/")) {
    if(length( stringr::str_split(vcfR::extract.gt(vcfRobj)[1,2], "\\/", simplify = T)) > 2){
      stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by this tool")
    } else{
      gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
      hetmatrix <- as.matrix(gtmatrix == "0/1",
                            nrow = nrow(gtmatrix), ncol=ncol(gtmatrix))

      hetmatrix <- cbind(F, hetmatrix)
    }
  } else {

    stop("You passed a monoclonal file.")
  }




  vcfRobj@gt[hetmatrix] <- NA

  for(i in 2:ncol(vcfRobj@gt)){
    for(j in 1:nrow(vcfRobj@gt)){
      vcfRobj@gt[j,i] <- stringr::str_replace(vcfRobj@gt[j,i], "0/0", "0")
      vcfRobj@gt[j,i] <- stringr::str_replace(vcfRobj@gt[j,i], "1/1", "1")
    }
  }


  fix <- as.matrix(vcfR::getFIX(vcfRobj, getINFO = T))
  gt <- as.matrix(vcfRobj@gt)
  meta <- append(vcfRobj@meta, paste("##Removed all het calls and made haploid"))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)
}



