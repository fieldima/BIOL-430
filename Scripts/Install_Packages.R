# Installing necessary libraries

install.packages(c(
  "BiocManager"
))

BiocManager::install(c(
  "MetaVolcanoR",
  "edgeR",
  "limma",
  "GEOquery",
  "EnhancedVolcano"
))
