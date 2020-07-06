#' @title Perform Empirical Bayes Shrinkage, Loci-by-Loci
#'
#' @param vcfRobject vcfR; a vcfR object
#' @param imputemissing boolean; Should imputation for missing loci be performed
#' @return matrix of wsaf (non-referent alelle) shrunk by a beta-binomial empirical bayes model
#' @details Loci must all be biallelic. Imputation for missing loci is performed by drawing a random number from the beta-binomial model
#' @details Note, optim will fail
#' @export

EB_for_wsaf <- function(vcfRobj, imputemissing = T, seed = 48){

  #..............................................................
  # catches
  #..............................................................
  if(!identical(vcfRobj, vcfRobj[vcfR::is.biallelic(vcfRobj),])){
    stop("VCF must be biallelic for this to work.")
  }

  #..............................................................
  # wsaf
  #..............................................................
  ad <- vcfR::extract.gt(vcfRobj, element = "AD")
  dp <- vcfR::extract.gt(vcfRobj, element = "DP", as.numeric = T)
  altad <- vcfR::masplit(ad, record=2, sort=0, decreasing = 0)
  NRAF <- altad/dp
  NRAF[is.nan(NRAF)] <- NA # the 0,0 are returning Nans

  # liftover nas to ref and alt
  altad[is.na(NRAF)] <- NA
  dp[is.na(NRAF)] <- NA

  # split out
  altad.list <- split(altad, 1:nrow(altad))
  dp.list <- split(dp, 1:nrow(dp))

  #..............................................................
  # EB
  #..............................................................
  #..................
  # get prior
  #..................
  # beta is conjugate of binomial
  # log-likelihood function
  get_EB_betabinomial <- function(ad, dp){

    if(length(unique(ad)) == 1){
      return("Site has no variation")
    } else {
      dat <- ad/dp
      dat <- dat[!is.na(dat)]
      ret <- fitdistrplus::fitdist(data =  dat,
                                   distr = "beta",
                                   method = "mme", # use method of moments, which may be a bit more robust
                                   start = list(shape1 = 0.5, shape2 = 0.5))


      alpha <- ret$estimate[["shape1"]]
      beta <- ret$estimate[["shape2"]]
      ret <- data.frame(alpha = alpha,
                        beta = beta)

      return(ret)
    }

  }

  EB.beta.params <- mapply(get_EB_betabinomial,
                           dp = dp.list,
                           ad = altad.list,
                           SIMPLIFY = F)

  #..................
  # get sample lvl "posterior" (shrinkage)
  #..................
  shrink_nrwsaf <- function(ad, dp, EB.beta.param){
    shrunk_nrwsaf <-  (ad + EB.beta.param[["alpha"]]) /  (dp + EB.beta.param[["alpha"]] + EB.beta.param[["beta"]] )
    return(shrunk_nrwsaf)
  }

  shrunken_nrwsaf <- mapply(shrink_nrwsaf,
                            ad = altad.list,
                            dp = dp.list,
                            EB.beta.param = EB.beta.params)

  # out the shrunk nrwsaf
  return(shrunken_nrwsaf)

}


#TODO extend to biallelic framework
#' @title Perform GT Imputation, Loci-by-Loci
#'
#' @param vcfRobject vcfR; a vcfR object
#' @param seed numeric; Reproducible Seed
#' @return matrix of GT (0,1,2) with missing loci imputed based on
#' @details Loci must all be biallelic. Imputation for missing loci is performed by drawing a 0 or 2 weighted by the PLAF
#' @export

binomial_draw_gt_imputation <- function(vcfRobj, seed = 48){
  set.seed(seed)
  #..............................................................
  # catches
  #..............................................................
  if(!identical(vcfRobj, vcfRobj[vcfR::is.biallelic(vcfRobj),])){
    stop("VCF must be biallelic for this to work.")
  }

  smplnames <- colnames(vcfRobj@gt)[2:ncol(vcfRobj@gt)]

  # get gt
  gt <- vcfR::extract.gt(vcfRobj, element = "GT")

  # convert to matrix
  gt012 <- gt012.refallelecounts <- vcfRmanip::gtmat012(vcfRobj)

  # drop het calls as uninformative
  gt012.refallelecounts[gt012 == "1"] <- NA

  #..............................................................
  # corner case catch -- if single loci vector, make it a matrix
  #..............................................................
  if (is.vector(gt012.refallelecounts)) {
    gt012.refallelecounts <- matrix(gt012.refallelecounts, nrow = 1,
                                    dimnames = list(NULL,
                                                    names(gt012.refallelecounts)))
    gt012 <- matrix(gt012, nrow = 1,
                    dimnames = list(NULL,
                                    names(gt012)))
  }

  #..............................................................
  # find allele counts and impute
  #..............................................................
  # find missing
  loci.missing <- apply(gt012, 1, function(x){return(sum(is.na(x)))})
  if (loci.missing > 0) {
    gt012.refallelecounts <- apply(gt012.refallelecounts, 1, function(x){
      return( mean(x == 0, na.rm = T) )
    })


    # draw from binaomial distribution if missing
    imputation.loci <- mapply(function(n.draws, success){
      ret <- rbinom(n = n.draws, size = 1, prob = success)
      ret <- gsub("0", "2", ret) # if fail, alt allele
      ret <- gsub("1", "0", ret) # if success, ref allele
      ret <- as.numeric(ret)
      return(ret)
    }, loci.missing, gt012.refallelecounts)

    #..............................................................
    # now put back in loci draws
    #..............................................................
    # not, we went by rows but R will go by columns to fill in
    # so can't do simple is.na drop in
    #
    # speed this up if case of one loci
    if (is.vector(gt012.refallelecounts)) {

      gt012.imp <- gt012
      gt012.imp[is.na(gt012.imp)] <- imputation.loci

    } else {

      gt012.imp <- split(gt012, f = 1:nrow(gt012))
      gt012.imp <- mapply(function(gt, imp){
        gt[is.na(gt)] <- unlist(imp)
        return(gt)
      }, gt012.imp, imputation.loci, SIMPLIFY = F) %>%
        do.call("rbind.data.frame", .)

    }

    # out
    colnames(gt012.imp) <- smplnames
  } else {
    # out
    gt012.imp <- gt012
  }

  return(gt012.imp)
}

