source("../general_code/sim_utils.R")

## https://www.sciencedirect.com/science/article/pii/S2001037023002684#se0200
## https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010295
## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061262
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10080318/
## https://proteomic.datacommons.cancer.gov/pdc/analysis/dbe94609-1fb3-11e9-b7f8-0a80fada099c?StudyName=CPTAC%20CCRCC%20Discovery%20Study%20-%20Proteome

## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10296504/ VEGFA NNMT
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1913536/ ppia tbp
## https://celldiv.biomedcentral.com/articles/10.1186/s13008-023-00103-9 GDF-15
## https://pubmed.ncbi.nlm.nih.gov/35678045/ NNMT
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6584416/ PPARA

set.seed(54535)

cancer_metadata <- read.table("GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt")
normal_metadata <- read.table("GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt")

cancer_data <- read.table("GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt")
normal_data <- read.table("GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt")
subset_col <- gsub("-", "\\.", sample(cancer_metadata[cancer_metadata[,2]=="KIRC",1], 30))
cancer_subs <- cancer_data[,subset_col]
subset_col <- gsub("-", "\\.", sample(normal_metadata[normal_metadata[,2]=="KIRC",1], 30))
normal_subs <- normal_data[,subset_col]
dim(normal_subs)
dim(cancer_subs)
sum(row.names(cancer_subs)==row.names(normal_subs))

## Filter require 3/4 non-zero in one condition.
Y <- cbind(normal_subs, cancer_subs)
other <- unlist(colSums(Y[(rowSums(cancer_subs!=0)<20)|(rowSums(normal_subs!=0)<20),]))
Y <- Y[(rowSums(cancer_subs!=0)>19)&(rowSums(normal_subs!=0)>19),]
Y <- rbind(Y, other)
Y <- as.matrix(Y)
sort(unlist(Y["GAPDH",]))
dim(Y)

data <- list()
data$Y <- Y
X <- c(rep(0, ncol(normal_subs)), rep(1, ncol(cancer_subs)))
data$X <- cbind(1, X)

indexa_res <- run_indexa(data, 250, 0, 2.5, denom=which(row.names(data$Y)=="GAPDH"))
saveRDS(indexa_res, "../output/indexa_res.RDS")

aldex2_res <- run_aldex2(data, 250, "BH", denom=which(row.names(data$Y)=="GAPDH"))
saveRDS(aldex2_res, "../output/aldex2_res.RDS")

aldex2_res <- readRDS("../output/aldex2_res.RDS")
indexa_res <- readRDS("../output/indexa_res.RDS")

indexa_res <- indexa_res[!is.na(indexa_res[,"gtt.pval.BH.adj"]),]
aldex2_res <- aldex2_res[!is.na(aldex2_res[,"we.eBH"]),]

sum(indexa_res[,"gtt.pval.BH.adj"]<0.05)
sum(aldex2_res[,"we.eBH"]<0.05)
cbind(indexa_res[c("PPIA", "TBP", "VEGFA", "NNMT" ,"GDF15", "PPARA"),"gtt.pval.BH.adj"],
      aldex2_res[c("PPIA", "TBP", "VEGFA", "NNMT", "GDF15", "PPARA"),"we.eBH"])

Yp <- apply(Y, 2, function(col) log2((col+0.5)/sum(col+0.5)))

mean(Yp["GAPDH",X==0])-mean(Yp["GAPDH",X==1])+2.5

write.csv(data.frame(indexa_res), "../../supplementary_files/aldex2_gtt_ccrcc_results.csv",
          row.names=T)
write.csv(data.frame(aldex2_res), "../../supplementary_files/aldex2_ccrcc_results.csv",
          row.names=T)
