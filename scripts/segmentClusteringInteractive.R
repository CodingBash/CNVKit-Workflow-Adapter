#
# Install if needed
#
#source("https://bioconductor.org/biocLite.R")
#biocLite("Rsamtools")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")

#
# Load genome
#
library(BSgenome.Hsapiens.UCSC.hg19)

#
# Load source libraries
# TODO: Organize dependencies
#
setwd("~/Git-Projects/Git-Research-Projects/CNVKit-Workflow-Adapter")
source("./scripts/helperFunctions.R")
source("./scripts/segmentClusteringLibrary.R")
source("./scripts/cnvkitAdapterFunctions.R")

#
# Load input
#
normal_samples <- load_samples(classes = c("N"), sampleList = "./resources/sampleList.csv")

cytobands <- retrieveCytobands(dir = "./resources/cytoBand.txt")
chromosomeSizes <- generateChromosomeSizes(genome = BSgenome.Hsapiens.UCSC.hg19)

normalSegments <- do.call(rbind, lapply(normal_samples, function(sample){
  return(retrieveCNVkitSegments(sample, dir="./resources/testCNVkitNorm/", genes = FALSE))
}))

normalSegments <- cnvKitSegmentsToBedFormat(normalSegments)

target_samples <- load_samples(classes = c("T", "F", "M"), sampleList = "./resources/sampleList.csv")

# Generate norminput argument
norminput <- retrieveNormInput(normalSegments)
norminput <- filterNormInput(norminput, length_threshold=10000000)

for(target_samples in seq(1, length(target_samples))){
  sample <- target_samples[target_samples]
  
  print(paste("Analyzing sample", sample))
  
  #
  # Retrieve sample data
  #
  setwd("~/Git-Projects/Git-Research-Projects/CNVKit-Workflow-Adapter")
  cnvkit_segment_data <- retrieveCNVkitSegments(sample, dir="./resources/testCNVkit/", genes = FALSE)
  cnvkit_bins_data <- retrieveCNVkitBins(sample, dir="./resources/testCNVkit/")
  
  segment_probes_assignment <- assignProbeNumbers(cnvkit_segment_data, cnvkit_bins_data)

  # Generate seginput argument
  seginput <- retrieveSegInput(cnvkit_segment_data = cnvkit_segment_data, segment_probes_assignment = segment_probes_assignment, sample=sample, chromosomeSizes = chromosomeSizes, cytobands = cytobands)
  print(paste("Retrieved segment input for sample", sample))
  
  # Generate ratinput argument
  ratinput <- retrieveRatInput(cnvkit_bins_data, sample)
  print(paste("Retrieved ratio input for sample", sample))  

  # Run CNprep:CNpreprocessing
  try({
    #TODO: Fix error
    segtable <- runCNpreprocessing(seginput = seginput, distrib="vanilla", ratinput = ratinput, norminput = norminput, modelNames = "E")
    print(paste("Produced segtable for sample", sample))
    
    #
    # Write results out
    #
    cd_local()
    #write.table(segtable, paste("segClusteringResultsPar/", sample, "_segtable.tsv", sep = ""), row.names = F, sep = "\t", quote = FALSE)
    print(head(segtable))
    print(paste("Wrote output for sample", sample))
  }, silent=TRUE)
  print("test")
}