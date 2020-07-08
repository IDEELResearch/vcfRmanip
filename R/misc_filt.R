#' @title vcfR2removepoorcovloci_dp
#'
#' @description Read vcfR object to remove loci where too few samples have coverage
#' @param vcfRobject vcfR; a vcfR object
#' @param minsampleswithcoverage numeric; the minimum number of samples with coverage \(i.e. DP > 0) needed for a loci to be kept
#' @export


vcfR2removepoorcovloci_dp <- function(vcfRobject = NULL, minsampleswithcoverage = 1){


  dp <- vcfR::extract.gt(vcfRobject, element = "DP", as.numeric = T)
  dpnot0counts <- apply(dp, 1, function(x){sum(x != 0)})
  keeploci <- dpnot0counts >= minsampleswithcoverage

  vcfRobject@gt <- vcfRobject@gt[keeploci,]

  fix <- as.matrix(vcfR::getFIX(vcfRobject, getINFO = T)[keeploci , ])
  gt <- as.matrix(vcfRobject@gt)
  meta <- append(vcfRobject@meta, "##Additional Filters for missingness by loci based on DP")

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}

#' @title Filter VCF based on given loci missingness by sample
#' @description Calculate the proportion of the loci that are missing by samples from a \code{vcfR} object
#' @param vcfRobject S4 object; A vcf that has been converted to a \code{vcfR} object
#' @export
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




#' @title vcffilter_ChromPos
#'
#' @description Produces a subsetted \code{vcfR} with specific loci excluded as determined \code{chromposdf}.
#'
#' @param chromposdf an dataframe that contains loci to be excluded and is in the form of a bedfile or has
#' been produced by \code{GFF2VariantAnnotation_Short}
#' @details chromposdf must be a tibble with columns: seqname (where chromomsome information is stored), start, end, geneid. The logic
#' behding the "seqname" is allow users to carry around chromosome names that may be different in the bed file and the vcf
#'
#' @export

#-----------------------------------------------------
# Filter VCF based on given positions among the chromosomes
#------------------------------------------------------
vcffilter_ChromPos <- function(vcfRobject = NULL,
                               chromposbed = NULL
){

  colnames(chromposbed) <- tolower(colnames(chromposbed))
  chromposbed <- split(chromposbed, f=factor(chromposbed$geneid))
  chromposlong <- do.call("rbind", parallel::mclapply(chromposbed, function(dat){

    i <- length(seq(from=dat$start, to=dat$end))

    temp <- tibble::tibble(CHROM = rep(dat$seqname, i),
                           POS = seq(from=dat$start, to=dat$end),
                           remove = rep("Y", i)
    )
    return(temp)

  }))

  # in case overlapping genes or geneids
  chromposlong <- chromposlong[!duplicated(chromposlong), ]
  chromposlong$POS <- as.character(chromposlong$POS)

  passloci <- tibble::as_tibble(vcfR::getFIX(vcfRobject)) %>%
    dplyr::select(CHROM, POS) %>%
    dplyr::left_join(x=., y=chromposlong, by=c("CHROM", "POS")) %>%
    dplyr::mutate(keep = ifelse(is.na(remove), TRUE, FALSE)) %>%
    dplyr::select(keep) %>%
    as.matrix(.)


  if(nrow(passloci) == 0){
    return(NA)
    warning("There were no variants found in the requested region(s)")
  } else{

    meta <- append(vcfRobject@meta, paste("##Additional Filters for excluding hypervariable regions determined by user"))
    fix <- vcfRobject@fix[ passloci, ]
    gt <- vcfRobject@gt[ passloci, ]


    # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
    newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

    return(newvcfR)
  }

}



#' @title bivcfR2segsites_wsaf
#'
#' @description Read biallelic vcfR object to segregated sites based on wsaf
#' @param vcfRobject vcfR; a biallelic vcfR object
#' @param err numeric; the difference between the max WSAF and min WSAF among samples being compared that is tolerated as not being a segregating site.
#' @export


bivcfR2segsites_wsaf <- function(vcfRobject = NULL, err = 0.025){

  if(!identical(vcfRobject, vcfRobject[vcfR::is.biallelic(vcfRobject)])){
    stop("VCF must be biallelic for this to work.")
  }

  ad <- vcfR::extract.gt(vcfRobject, element = "AD")
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

  if(sum(segsites) == 1){
    vcfRobj@gt <- t(vcfRobj@gt[segsites,])
    fix <- t(vcfR::getFIX(vcfRobj, getINFO = T)[segsites,])
    meta <- append(vcfRobj@meta, paste("##Additional Filters for segregating sites, such that GT call must be segregating within samples"))

  } else {
    vcfRobj@gt <-  vcfRobj@gt[segsites,]
    fix <- as.matrix(vcfR::getFIX(vcfRobj, getINFO = T)[segsites,])
    meta <- append(vcfRobj@meta, paste("##Additional Filters for segregating sites, such that GT call must be segregating within samples"))

  }



  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = vcfRobj@gt)

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
  ploidy <- unique(stringr::str_extract(vcfRobj@meta, "ploidy=[0-9]+"))
  ploidy <- ploidy[!is.na(ploidy)]
  ploidy <- as.numeric(stringr::str_extract(ploidy, "[0-9]+"))

  if(ploidy == 2) {
    gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
    gtmatrix[gtmatrix == "0/0"] <- 0
    gtmatrix[gtmatrix == "0/1"] <- 1
    gtmatrix[gtmatrix == "1/1"] <- 2

    # phased (discounting this information)
    gtmatrix[gtmatrix == "0|0"] <- 0
    gtmatrix[gtmatrix == "0|1"] <- 1
    gtmatrix[gtmatrix == "1|1"] <- 2

    gtmatrix[is.na(gtmatrix)] <- NA
    gtmatrix <- apply(gtmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
  } else if(ploidy == 1) {

    gtmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=T)
    gtmatrix[gtmatrix == 0] <- 0
    gtmatrix[gtmatrix == 1] <- 2
    gtmatrix[is.na(gtmatrix)] <- NA
  } else {
    stop("ploidy not supported by this tool")
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



