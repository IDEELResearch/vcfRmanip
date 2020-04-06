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
  if(!identical(vcfRobj, vcfRobj[vcfR::is.biallelic(vcfRobj)])){
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

