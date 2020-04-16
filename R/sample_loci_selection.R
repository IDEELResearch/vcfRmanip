#' @title Select Samples from VCF
#' @description Select samples from a \code{vcfR}
#' @param smplvctr vector; A character vector that corresponds to the samples that should be extracted
#' @export

select_samples <- function(vcfRobject = NULL,
                           smplvctr = NULL ){

  meta <- append(vcfRobject@meta, paste("##Specific samples were selected using vcfRmanip"))
  fix <- vcfRobject@fix
  keep <- c(TRUE, colnames( vcfRobject@gt[,2:ncol(vcfRobject@gt)] ) %in% smplvctr)
  gt <- vcfRobject@gt[ , keep]
  # first true for FORMAT column

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = vcfRobject@fix, gt = gt)

  return(newvcfR)
}



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




#' @title vcfR2SubsetChromPos
#'
#' @description Produces a subsetted \code{vcfR} based on chromosome. Note the chromosome name must be contained in the column \code{seqname}
#'
#' @export

vcfR2SubsetChromPos <- function(vcfRobject = NULL,
                                chromposbed = NULL){

  colnames(chromposbed) <- tolower(colnames(chromposbed))
  chromposbed <- split(chromposbed, f=factor(chromposbed$geneid))
  chromposlong <- do.call("rbind", parallel::mclapply(chromposbed, function(dat){

    i <- length(seq(from=dat$start, to=dat$end))

    temp <- tibble::tibble(CHROM = rep(dat$seqname, i),
                           POS = seq(from=dat$start, to=dat$end)
    )
    return(temp)

  }))

  chromposlong$keep <- "Y"
  chromposlong$POS <- as.character(chromposlong$POS)
  passloci <- tibble::as_tibble(vcfR::getFIX(vcfRobject)) %>%
    dplyr::select(CHROM, POS) %>%
    dplyr::left_join(x=., y=chromposlong, by=c("CHROM", "POS")) %>%
    dplyr::mutate(keep = !is.na(keep)) %>%
    dplyr::select(keep) %>%
    as.matrix(.)


  meta <- append(vcfRobject@meta, paste("##Chromsome and Positiion were subsetted by user using vcfR2SubsetChromPos"))
  #......
  # catch instance of 1
  #......
  if(sum(passloci) == 1){
    # going to do some extra work so attributes are all correct for matrix that has to go into vcfR class
    fix <- matrix(vcfRobject@fix[ passloci, ], nrow=1)
    fix <- as.data.frame(fix)
    colnames(fix) <- colnames(vcfRobject@fix)
    fix <- as.matrix(fix)


    gt <- matrix(vcfRobject@gt[ passloci, ], nrow=1)
    gt <- as.data.frame(gt)
    colnames(gt) <- colnames(vcfRobject@gt)
    gt <- as.matrix(gt)

  } else {
    fix <- vcfRobject@fix[ passloci, ]
    gt <- vcfRobject@gt[ passloci, ]
  }

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)

  return(newvcfR)

}





#' @title Split a VCF to Individual Loci
#' @param vcfRobj S4object; vcfR object from the vcfR package
#' @param write_new_meta boolean; Flag for whether or not the original vcf metadata should be included.
#' @details If write_new_meta is set to False, a simplified vcf header will be returned (just file format)
#' @return List of vcfRobj

vcf_loci_splitter <- function(vcfRobj, write_new_meta = F){

  #..............................................................
  # simple track
  #..............................................................
  chrompos <- paste0(vcfR::getCHROM(vcfRobj), "_", vcfR::getPOS(vcfRobj))

  #..............................................................
  # output
  #..............................................................
  if (write_new_meta) {
    meta.split <- vcfRobj@meta[[1]] # just return file format version
  } else {
    meta.split <- vcfRobj@meta
  }

  # extract
  gt.split <- split(vcfRobj@gt, 1:nrow(vcfRobj@gt), drop = F)
  gt_col_names <- colnames(vcfRobj@gt)

  fix.split <-  split(vcfRobj@fix, 1:nrow(vcfRobj@fix), drop = F)
  fix_col_names <- colnames(vcfRobj@fix) # these can't change but for transp.

  # out function
  make_new_vcfR <- function(meta, fix, gt){
    fix <- matrix(fix, nrow = 1)
    colnames(fix) <- fix_col_names

    gt <- matrix(gt, nrow = 1)
    colnames(gt) <- gt_col_names

    ret <- new("vcfR", meta = meta, fix = fix, gt = gt)
    return(ret)
  }


  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR.list <- mapply(make_new_vcfR,
                         meta = meta.split, gt = gt.split, fix = fix.split)
  # fix names
  newvcfR.list <- lapply(newvcfR.list, function(newvcfR, fix_names, gt_names){
    colnames(newvcfR@fix) <- fix_names
    colnames(newvcfR@gt) <- gt_names
    return(newvcfR)
  }, fix_names = fix_col_names, gt_names = gt_col_names)

  # add in easier identifier information
  names(newvcfR.list) <- chrompos

  return(newvcfR.list)
}




















