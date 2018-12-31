# dab
# not exported

dab <- function(afmat){

  ret <-  afmat[,1]*(1-afmat[,2]) +  afmat[,2]*(1-afmat[,1])  #  dAB <- fA*(1-fB) + fB*(1-fA)
  return(ret)

}

#' @title DAB weight calculator
#'
#' @author Nick Brazeau
#' @description This is based on the within sample Fw gen dist function described in PMID: 26943619
#' @details Inspired by Bob Verity's corMat function

wicalc <- function(NRAFdf, window=50){
  # error handle
  if(window > nrow(NRAFdf)){
    stop("The smoothing window you have specified is larger than the number of SNPs you have on Chromosome ", paste(NRAFdf$CHROM[1]), ". Remember a window is centered a SNP and needs to be considered for i +/- L/2")
  }
  # start
  cij <- NULL # init, pos i and j -- compare pop AF for i and j
  for(i in 1:nrow(NRAFdf)){
    LW <- i - window
    UW <- i + window
    if(LW <= 0 ){
      UW <- UW + abs(LW) + 1
      LW <- 1
    } else if( nrow(NRAFdf) <= UW  ){
      LW <- LW - (UW -  nrow(NRAFdf))
      UW <- nrow(NRAFdf)
    }

    cij <- rbind(cij,
                 cbind(NRAFdf$POS[i], NRAFdf$POS[LW:UW])) # this should be in a while loop...will be faster, but the problem is the other side...so need an ifelse condtional?
  }
  cij <- cij[ cij[,1] != cij[,2] , ] # drop self comparisons
  cij.list <- split(cij, seq(nrow(cij)))

  rij.list <-  parallel::mclapply(cij.list, function(v){
    x1 <- unlist( NRAFdf[NRAFdf$POS == v[1], 3:ncol(NRAFdf)] )
    mu1 <- mean(x1,na.rm=TRUE)
    x2 <- unlist( NRAFdf[NRAFdf$POS == v[2], 3:ncol(NRAFdf)] )
    mu2 <- mean(x2,na.rm=TRUE)
    r <- sum((x1-mu1)*(x2-mu2), na.rm=TRUE)/sqrt( sum((x1-mu1)^2, na.rm=TRUE)*sum((x2-mu2)^2, na.rm=TRUE) )
    return(r)
  }
  )

  cij <- do.call("rbind", cij.list)
  rij <- do.call("rbind", rij.list)

  rij <- cbind(cij,rij)

  widf <- rij %>%
    tibble::as.tibble(.) %>%
    magrittr::set_colnames(c("POS", "WindowPos", "R")) %>%
    dplyr::mutate(Rsq = R^2 ) %>%
    dplyr::group_by(POS) %>%
    dplyr::summarise(Rsqsum = sum(Rsq)) %>%
    dplyr::mutate(wi = 1 / (1 + Rsqsum) ) # TODO speed this up

  return(widf)
}



#' @title vcfR2Fw_pairwisegendist
#'
#' @author Nick Brazeau
#'
#' @export



# TODO fix sloppy NRAF

vcfR2Fw_pairwisegendist <- function(vcfRobject = NULL, biallelicsnps=TRUE, segsites=0.001, window=50){

  if(!class(vcfRobject) %in% "vcfR"){
    stop("vcfRobject must contain  class vcfR")
  }

  if(biallelicsnps == FALSE){
    stop("Must take in biallelic SNPs")
  }


  vcf <-vcfR::extract.indels(vcfRobject, return.indels = F) # subset to SNPs
  vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  if(!identical(vcfRobject, vcf)){
    stop("Your file does not appear to be a vcf with only biallelic SNPs")
  }
  # Run segsites
  if(!is.null(segsites)){
    vcf <- NFBtools::vcfR2segsites(vcfRobject = vcf, err = segsites)
  }


  #--------------------------------------------------------
  # read in and convert vcfR object
  #--------------------------------------------------------

  # extract the within sample allele frequencies
  ad <- vcfR::extract.gt(vcfRobject, element = "AD")
  refad <- masplit(ad, record=1, sort=0, decreasing = 0)
  altad <- masplit(ad, record=2, sort=0, decreasing = 0)

  # now divide alt by ref to get NRAF
  NRAF <- altad/(altad + refad)
  NRAF[is.nan(NRAF)] <- NA #TODO fix that the 0,0 are returning Nans, the NA right now is a bit sloppy

  # get chrom & pos for coordinate information and make LD calculations not possible/independent between chromosomes
  NRAFdf <- vcfR::getFIX(vcfRobject) %>%
    as.tibble(.) %>%
    dplyr::select(CHROM, POS) %>%
    dplyr::mutate(CHROM = factor(CHROM)) %>%
    dplyr::mutate(POS = as.numeric(POS)) %>%
    cbind.data.frame(., NRAF) %>%
    as.tibble(.)

  NRAFlist <- split(NRAFdf, NRAFdf$CHROM)

  # get all pairs of samples
  pairs <- t(combn(
    colnames(NRAFlist[[1]])[3:ncol(NRAFlist[[1]])], # just colnames of samples not chrom,pos
    m=2)) # choose 2

  # get the correlation matrix
  wi <- parallel::mclapply(NRAFlist, wicalc, window=window)
  #--------------------------------------------------------
  # calculations
  #--------------------------------------------------------

  # pairings
  pairs.list <- split(pairs, seq(nrow(pairs)))

  # just need a vector of the wi results
  wiret <- do.call("rbind", wi)$wi

  dablist <- parallel::mclapply(pairs.list, function(pair, nrafdf=NRAFdf){

    dabret <- dab( nrafdf[ , colnames(nrafdf) %in% pair] )
    dabret <- 1/sum(!is.na(dabret)) * sum( (dabret * wiret), na.rm = T )

    ret <- cbind.data.frame(pair[1], pair[2], dabret)
    return(ret)

  }) # end internal loop for DAB calculations

  ### FINAL RETURN
  ret <- do.call("rbind", dablist)
  return(ret)
}





#' @title NRAFlist2Fw_pairwisegendist
#'
#' @author Nick Brazeau
#'
#' @export



# TODO eventually make this compatible with the end MIP output

NRAFdf2Fw_pairwisegendist <- function(NRAFdf = NULL, window=50){

  # setup
  colnames(NRAFdf)[1:2] <- toupper(colnames(NRAFdf)[1:2])
  NRAFlist <- split(NRAFdf, NRAFdf$CHROM)

  # get the correlation matrix
  wi <- parallel::mclapply(NRAFlist, wicalc, window=window)
  #--------------------------------------------------------
  # calculations
  #--------------------------------------------------------

  # get all pairs of samples
  pairs <- t(combn(
    colnames(NRAFlist[[1]])[3:ncol(NRAFlist[[1]])], # just colnames of samples not chrom,pos
    m=2)) # choose 2

  # pairings
  pairs.list <- split(pairs, seq(nrow(pairs)))

  # just need a vector of the wi results
  wiret <- do.call("rbind", wi)$wi

  dablist <- parallel::mclapply(pairs.list, function(pair, nrafdf=NRAFdf){

    dabret <- dab( nrafdf[ , colnames(nrafdf) %in% pair] )
    dabret <- 1/sum(!is.na(dabret)) * sum( (dabret * wiret), na.rm = T )

    ret <- cbind.data.frame(pair[1], pair[2], dabret)
    return(ret)

  }) # end internal loop for DAB calculations

  ### FINAL RETURN
  ret <- do.call("rbind", dablist)
  return(ret)
}

