#' @title Plot Genetic Autocorrelation result
#'
#' @description Plots the autocorrelation between loci (calculated from population allele frequencies)
#'
#' @param genautocorrobj an object of class \code{genautocorr}, as produced by the function \code{polyIBD::genautocorr}
#'
#' @export

# TODO this was copied over from polyIBD on Aug 13, 2018. Need to determine if needed and fix up LD processing

plotgenautocorr <- function(genautocorrobj){

  plotdf <- data.frame(CHROM=rep(vcfdf_fromlist$CHROM[1], length(gendist)),
                       correlation = c(c),
                       gendist = c(gendist))

  corrplot <- plotdf %>%
    dplyr::mutate(gendisttol = gendist+1) %>%
    ggplot(aes(x=gendisttol, y=correlation)) +
    geom_point() +
    geom_hline(yintercept = 0, colour="#de2d26") +
    scale_x_log10() +
    xlab("log(base-pair distance + 1)") + ylab("correlation") +
    ggtitle(paste(vcfdf_fromlist$CHROM[1])) +
    theme_minimal()

}

