#GSE78888, GSE77193, GSE68380, GSE77192, and GSE67473

#GSE68380 and GSE67473 are in hours, rest are measured in days

#Load required packages

library(limma)
library(GEOquery)
library(pheatmap)
library(tidyverse)
library(MetaVolcanoR)
library(ggvenn)
library(patchwork)
library(EnhancedVolcano)
library(biomaRt)

#Load the data
GSE_to_load <- list("GSE78888", "GSE77193", "GSE68380", "GSE77192", "GSE67473")

dataList <- map(GSE_to_load, getGEO) %>% flatten()

GSE78888 <- dataList[1] 
dataList <- dataList[2:5]

#Compare only at time point of one day/24 hours

only1d <- function(GSE){
  res <- GSE[,GSE$"time:ch1" %in% c("1", "24", "1d")]
}

fixedList <- map(dataList, only1d)

fixedList[[2]] <- fixedList[[2]][,fixedList[[2]]$`time:ch1` == "24"]
fixedList[[4]] <- fixedList[[4]][,fixedList[[4]]$`time:ch1` == "24"]


## for GSE78888
GSE78888_fixed <- GSE78888 %>% map(only1d)

#Correlation Plot

## First get all data into one data frame
all_data <- map_df(fixedList, function(x)x@assayData$exprs) %>% bind_cols(map_df(GSE78888_fixed, function(x)x@assayData$exprs)) %>% na.omit()

## Get the correlations
all_cors <- all_data %>% cor(method = "pearson") 

## Retrieve sample meta data
all_meta <- map_df(fixedList, function(x)x@phenoData@data) %>% bind_rows(map_df(GSE78888_fixed, function(x)x@phenoData@data))

heatmap_meta <- all_meta %>% dplyr::select(`virus:ch1`, `infection:ch1`, source_name_ch1) %>%
  unite(col = "Treatment", `virus:ch1`:`infection:ch1`) %>%
  mutate(Treatment = str_remove(Treatment, "NA"), Treatment = str_remove(Treatment, "_"), Treatment = str_replace(Treatment, "mockulum", "Mock")) %>% 
  separate(col = source_name_ch1, c("Tissue"), sep = ",")

colnames(all_cors) <- rownames(heatmap_meta)
rownames(all_cors) <- rownames(heatmap_meta)

## Plot the heatmap
all_cors %>% pheatmap(cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = heatmap_meta) 

#GSM1670634-36 and GSM1670645-47 are outliers, need to remove them
fixedList[[2]] <- fixedList[[2]][,!fixedList[[2]]$geo_accession %in% c("GSEM1670634", "GSEM1670635", "GSEM1670636", "GSEM1670645", "GSEM1670646", "GSEM1670647")]

#split WNVWT from WNVE218A

pickVirus <- function(GSE, genotype){
  res <- GSE[,GSE$"virus:ch1" %in% c(genotype, "mockulum", "Mock")]
}

WNVWT <- map(fixedList, ~ pickVirus(.x, genotype = "WNVWT")) 
WNVE218A <- map(fixedList, ~ pickVirus(.x, genotype = "WNVE218A"))


## Split virus for GSE78888
pickVirus2 <- function(GSE, genotype){
  res <- GSE[,GSE$"infection:ch1" %in% c(genotype, "mockulum", "Mock")]
}

GSE78888WT <- GSE78888_fixed %>% map(~ pickVirus2(.x, genotype = "WNVWT"))
GSE78888E218A <- GSE78888_fixed %>% map(~ pickVirus2(.x, genotype = "WNVE218A"))

#Split Esets into Data and Meta Data

splitEset <- function(eset){
  Dat <- assayData(eset)$exprs
  Meta <- pData(eset) %>% dplyr::rename(time = "time:ch1", virus = "virus:ch1") %>% mutate(time = as.numeric(time))
  res <- list(eset, Dat, Meta)
  res
}

WNVWT_split <- map(WNVWT, splitEset)
WNVE218A_split <- map(WNVE218A, splitEset)

WTDat <- assayData(GSE78888WT$GSE78888_series_matrix.txt.gz)$exprs
WTMeta <- pData(GSE78888WT$GSE78888_series_matrix.txt.gz) %>% dplyr::rename(virus = "infection:ch1")

E218ADat <- assayData(GSE78888E218A$GSE78888_series_matrix.txt.gz)$exprs
E218AMeta <- pData(GSE78888E218A$GSE78888_series_matrix.txt.gz) %>% dplyr::rename(virus = "infection:ch1")

GSE_WT <- list(list(GSE78888WT$GSE78888_series_matrix.txt.gz, WTDat, WTMeta))
GSE_E218A <- list(list(GSE78888E218A$GSE78888_series_matrix.txt.gz, E218ADat, E218AMeta))


#Fit models to identify differentially expressed genes between WT and mutant virus

runLimma <- function(obj, strain){
  Fit <- lmFit(obj[[2]], model.matrix( ~ virus, data = obj[[3]])) %>% eBayes()
  try(Top1000 <- topTable(Fit, coef = strain, sort.by = "p", p.value = 0.01, lfc = 1, confint = TRUE, n = 1000) %>%
    rownames_to_column("ID") %>% inner_join(obj[[1]]@featureData@data) %>% dplyr::select(ENSEMBL_ID, logFC, adj.P.Val, CI.L, CI.R) %>%
    dplyr::rename(Symbol = ENSEMBL_ID, Log2FC = logFC, pvalue = adj.P.Val))
}

WNVWT_list <- map(WNVWT_split, ~ runLimma(.x, strain = "virusWNVWT"))
WNVE218A_list <- map(WNVE218A_split, ~ runLimma(.x, strain = "virusWNVE218A"))

#No DEGs for GSE77192 for WNVE218A
WNVE218A_list <- WNVE218A_list[c(1,2,4)]

# Limma with GSE78888

GSEWT_list <- GSE_WT %>% map(~ runLimma(.x, strain = "virusWNVWT"))
GSEE218A_list <- GSE_E218A %>% map(~ runLimma(.x, strain = "virusWNVE218A"))

rm(GSE78888, GSE78888_fixed, dataList)

GSEs <- c("GSE78888", "GSE77193", "GSE68380", "GSE77192", "GSE67473")

WT_res <- append(GSEWT_list, WNVWT_list)
names(WT_res) <- GSEs

E218A_res <- append(GSEE218A_list, WNVE218A_list)
names(E218A_res) <- GSEs[c(1,2,3,5)]

WT_symbols <- map(WT_res, function(x) x$Symbol) %>% map(~ str_subset(., pattern = ".+"))
E218A_symbols <- map(E218A_res, function(x) x$Symbol) %>% map(~ str_subset(., pattern = ".+"))

#Find DEGs common to all experiments
common_WT <- intersect(WT_symbols[[1]], WT_symbols[[2]]) %>% intersect(WT_symbols[[3]]) %>% intersect(WT_symbols[[4]]) %>% intersect(WT_symbols[[5]])
common_E218A <- intersect(E218A_symbols[[1]], E218A_symbols[[2]]) %>% intersect(E218A_symbols[[3]]) %>% intersect(E218A_symbols[[4]])

WTvenn <- WT_symbols  %>% ggvenn() + ggtitle("Wild Type")
E218Avenn <- E218A_symbols %>% ggvenn() + ggtitle("E218A") 

WTvenn + E218Avenn

## Finding DEGs consistent across studies using p-value combining via Fisher's method. 

Comb_WT <- combining_mv(diffexp=WT_res,
             pcriteria='pvalue', 
             foldchangecol='Log2FC',
             genenamecol='Symbol',
             geneidcol=NULL,
             metafc='Mean',
             metathr=0.1, 
             collaps=TRUE,
             jobname="MetaVolcano",
             outputfolder=".",
             draw='HTML')


Comb_E218A <- combining_mv(diffexp=E218A_res,
                           pcriteria='pvalue', 
                           foldchangecol='Log2FC',
                           genenamecol='Symbol',
                           geneidcol=NULL,
                           metafc='Mean',
                           metathr=0.1, 
                           collaps=TRUE,
                           jobname="MetaVolcano",
                           outputfolder=".",
                           draw='HTML')

Comb_E218A@MetaVolcano + Comb_WT@MetaVolcano

Comb_vol_WT <- Comb_WT@metaresult
Comb_vol_E218A <- Comb_E218A@metaresult

#Volcano plot using combined P-values (Fischer's Test)
EnhancedVolcano(Comb_vol_E218A, lab = Comb_vol_E218A$Symbol, x = "metafc", y = "metap", pCutoff = 0.01, FCcutoff = 1, title = "E218A", subtitle = NULL) +
  EnhancedVolcano(Comb_vol_WT, lab = Comb_vol_WT$Symbol, x = "metafc", y = "metap", pCutoff = 0.01, FCcutoff = 1, title = "Wild Type", subtitle = NULL) 

#Filter for cutoffs
filtered_WT <- Comb_vol_WT %>% filter(abs(metafc) >= 1) %>% filter(metap <= 0.01)
filtered_E218A <- Comb_vol_E218A %>% filter(abs(metafc) >= 1) %>% filter(metap <= 0.01)

#Venn diagram
list(WT = Comb_WT@metaresult %>% pull(Symbol), E218A = Comb_E218A@metaresult %>% pull(Symbol)) %>%
  ggvenn()

#Get all GO Terms
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
GO_terms <- getBM(attributes = c("ensembl_transcript_id", "go_id", "name_1006", "namespace_1003"),
      filters = "ensembl_transcript_id",
      values = c(Comb_vol_WT$Symbol, Comb_vol_E218A$Symbol),
      mart = ensembl) %>% rename(Symbol = ensembl_transcript_id) %>% dplyr::filter(namespace_1003 == "biological_process", name_1006 != "biological_process")

#Common DEGs before meta-analysis
WT_premeta_GO <- common_WT %>% data.frame(Symbol = .) %>% left_join(GO_terms)
E218A_premeta_GO <- common_E218A %>% data.frame(Symbol = .) %>% left_join(GO_terms)

#DEGs shared between two
shared_DEGs <- inner_join(filtered_WT, filtered_E218A, by = "Symbol")
shared_GO <- shared_DEGs %>% left_join(GO_terms) %>%
  mutate(upregulated = metafc.x > 0) %>%
  dplyr::select(Symbol, go_id, name_1006, upregulated) %>%
  group_by(name_1006, upregulated) %>%
  summarise(n = n()) %>%
  drop_na() %>%
  arrange(desc(n))

#DEGs unique to WT
WT_DEGs <- anti_join(filtered_WT, filtered_E218A, by = "Symbol")
WT_GO <- WT_DEGs %>% left_join(GO_terms) %>%
  mutate(upregulated = metafc > 0) %>%
  dplyr::select(Symbol, go_id, name_1006, upregulated) %>%
  group_by(name_1006, upregulated) %>%
  summarise(n = n()) %>%
  drop_na() %>%
  arrange(desc(n))

#DEGs unique to E218A
E218A_DEGs <- anti_join(filtered_E218A, filtered_WT, by = "Symbol")
E218A_GO <- E218A_DEGs %>% left_join(GO_terms) %>%
  mutate(upregulated = metafc > 0) %>%
  dplyr::select(Symbol, go_id, name_1006, upregulated) %>%
  group_by(name_1006, upregulated) %>%
  summarise(n = n()) %>%
  drop_na() %>%
  arrange(desc(n))

#Turn tables into CSVs
write_csv(shared_GO, "Outputs/Shared_GO_Terms.csv")
write_csv(WT_GO, "Outputs/WT_GO_Terms.csv")
write_csv(E218A_GO, "Outputs/E218A_GO_Terms.csv")