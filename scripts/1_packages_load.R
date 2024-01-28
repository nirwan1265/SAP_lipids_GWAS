### Packages needed:
packages <- c("tidyverse","ggplot2","gplots", "Rsamtools","GenomicAlignments",
              "rtracklayer","GenomicRanges","AnnotationHub","knitr","gtools",
              "data.table","stringi","GBJ","metap","multtest","Hmisc","devtools",
              "SNPRelate","gdsfmt","dplyr","vcfR","tidyr","AssocTests","SKAT",
              "NCmisc","ACAT","PANTHER.db","UniProt.ws","ape","sp","rgdal",
              "rworldmap","janitor","countrycode","tibble","vroom","gtools",
              "tictoc","gridExtra","VennDiagram", "foreach","doParallel",
              "ranger","tidyverse","kableExtra","parallel","stringr","purrr"
)

### Install packages that are not already installed
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if (length(new_packages) > 0) {
#   install.packages(new_packages)
# } else {
#   cat("All packages are already installed.\n")
# }
# # Install packages from Bioconductor
# BiocManager::install(c("GenomicAlignments", "rtracklayer", "multtest", "PANTHER.db", "UniProt.ws"))
# remotes::install_github("yaowuliu/ACAT")


# Load Packages
suppressMessages(invisible(lapply(packages, library, character.only = TRUE)))
