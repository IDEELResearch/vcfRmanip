
#------------------------------------------------
#' @title Plot vcfR Object GT Matrix
#'
#' @description Takes a \code{vcfRobject} and plots the GT levels
#'
#' @importFrom
#' @export

plot_vcfRobj_GT <- function(vcfRobj = NULL){
  #.....................
  # Read and check input
  #.....................

  if(!inherits(vcfRobj, "vcfR")){
    stop("vcfR object must be of class vcfR")
  }


  # check biallelic
  if (!identical(vcfRobj[vcfR::is.biallelic(vcfRobj)], vcfRobj)) {
    warning("This function only supports biallelic SNPs, subsetting now")
    vcfRobj <- vcfRobj[vcfR::is.biallelic(vcfRobj)]
  }


  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  ploidy <- stringr::str_extract(vcfRobj@meta[grepl("ploidy", vcfRobj@meta)], "ploidy=[0-9]")
  ploidy <- stringr::str_split(t(ploidy), "=", simplify = T)
  ploidy <- as.numeric(ploidy[1,2])

  if(ploidy == 2){

    snpmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
    snpmatrix[snpmatrix == "0/0"] <- 0
    snpmatrix[snpmatrix == "0/1"] <- 1
    snpmatrix[snpmatrix == "1/1"] <- 2
    snpmatrix <- apply(snpmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
  } else if(ploidy == 1){
    snpmatrix <- vcfR::extract.gt(vcfRobj, element='GT', as.numeric=T)
    snpmatrix[snpmatrix == 0] <- 0
    snpmatrix[snpmatrix == 1] <- 2
  } else {
    stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by polyIBD")
  }

  CHROM <- vcfR::getCHROM(vcfRobj)
  POS <- vcfR::getPOS(vcfRobj)
  snpdf <- cbind.data.frame(CHROM, POS, snpmatrix)
  snpdf <- snpdf %>%
    tidyr::pivot_longer(cols = -c("CHROM", "POS"), names_to="sample", values_to="GT") %>%
    dplyr::mutate(GTfact = factor(GT, levels=c(0,1,2), labels=c("Ref", "Het", "Alt")))

  # get lagged position
  snpdf <- snpdf %>%
    dplyr::group_by(CHROM, sample) %>%
    dplyr::mutate(start = lag(POS),
                  end = POS,
                  start = ifelse(is.na(start), 0, start)) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(sample = factor(sample)) %>%
    dplyr::mutate(sampleindex = as.numeric(factor(sample)))  # Will pull out levels

  # plot out
  snpdf %>%
    ggplot() +
    geom_rect(mapping=aes(xmin=start, xmax=end, ymin=sampleindex-0.49, ymax=sampleindex+0.49, fill=GTfact)) +
    scale_y_continuous(breaks = min(snpdf$sampleindex):max(snpdf$sampleindex),
                       labels = levels(factor(snpdf$sample))) +
    rplasmodium::scale_x_genome(scale = "Mb") +
    scale_fill_manual("Genotype Call", values=c("#AA0A3C", "#313695", "#fee090", "#cccccc")) +
    facet_grid(CHROM~.) +
    xlab("Genomic Coords. (Mb)") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size=9, family = "Arial", angle = 45),
          axis.title.y = element_text(size=14, face="bold", family = "Arial"),
          axis.text.y = element_text(size=12, face="bold", family = "Arial"))

}



