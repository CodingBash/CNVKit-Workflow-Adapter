#
# Retrieve a sample's FACETS segmented data from specified directory
#
retrieveCNVkitSegments <- function(sample, dir = "resources/testCNVkit/", genes = TRUE){
  cnvkit_segment_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.cns"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  if(genes == FALSE){
    cnvkit_segment_data <- cnvkit_segment_data[,-c(4)]
  }
  return(cnvkit_segment_data)
}

retrieveCNVkitBins <- function(sample, dir = "resources/testCNVkit/"){
  cnvkit_bins_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.cnr"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(cnvkit_bins_data)
}

retrieveCNVkitTargetCoverage <- function(sample, dir = "resources/testCNVkit/"){
  cnvkit_bins_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.targetcoverage.cnn"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(cnvkit_bins_data)
}

retrieveCNVkitAntiTargetCoverage <- function(sample, dir = "resources/testCNVkit/"){
  cnvkit_bins_data <- as.data.frame(read.table(paste0(dir, sample, ".ucsc.hg38.bwa.BQSR.antitargetcoverage.cnn"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote=""))
  return(cnvkit_bins_data)
}

sourceTestCode <- function(){
  setwd("~/Git-Projects/Git-Research-Projects/CNVKit-Workflow-Adapter")
  segments_table <- retrieveCNVkitSegments("hF2", dir = "resources/testCNVkit/", genes = FALSE)
  bins_table <- retrieveCNVkitBins("hF2", dir = "resources/testCNVkit/")
  scheme_table <- retrieveCNVkitScheme("hF2", dir = "resources/testCNVkit/") 
  scheme_table <- retrieveCNVkitAntiTargetCoverage("hF2", dir = "resources/testCNVkit/") 
}
