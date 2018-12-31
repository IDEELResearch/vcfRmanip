#' @title gff2geneannotation
#'
#'
#' @export


GFF2VariantAnnotation_Short <- function(gff){
  ## -----------------------------------
  #  READ IN GFF GENE INFORMATION
  # -----------------------------------
  geneid.end <- grep('##FASTA',readLines(gff))
  geneid.gff <- read.delim(file = gff,
                           nrows = geneid.end-1, comment= "#", header=F)
  colnames(geneid.gff) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "info")


  geneid.gff <- subset(geneid.gff, geneid.gff$feature == "gene") # subset to only genes, don't want the other mRNA, etc data

  # -----------------------------------
  #  EXTRACT GENEID from GFF INFO
  # -----------------------------------
  geneid.gff[,10:13] <- stringr::str_split(geneid.gff$info, ";", n=4, simplify = T) # give it to columns to parse on
  geneid.gff <- geneid.gff[,c(1:10, 12)] # drop that second column -- we don't need it
  colnames(geneid.gff)[10:11] <- c("GeneID", "Description")
  geneid.gff$GeneID <- gsub("ID=", "", geneid.gff$GeneID, fixed=T)
  # These are the gene identifiers that we care about
  geneid.gff$Description <- gsub("description=", "", geneid.gff$Description, fixed=T)
  geneid.gff$Description <- gsub("+", "_", geneid.gff$Description, fixed=T)
  geneid.gff$Description <- gsub("%", "-", geneid.gff$Description, fixed=T)
  return(geneid.gff)

}





