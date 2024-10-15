source("../general_code/sim_utils.R")
library("ExtDist")

args = commandArgs(trailingOnly=TRUE)

p.method <- args[3]

benchmark_func <- function(n, iters, mc.samples, min_dep, max_dep, keep_taxa,
                           subset_species, lfcs, prevalences, epsilon_l,
                           epsilon_u) {
    all_padj <- c()
    all_lfcs <- c()
    for(i in 1:iters) {
        seq_depths <- runif(n, min_dep, max_dep)
        data <- sim_dossa2(n, "Stool", seq_depths, subset_species,
                           lfcs, prevalences)
        Y <- data$Y[keep_taxa,]
        other <- colSums(data$Y[other_taxa,])
        Y <- rbind(Y, other)
        colnames(Y) <- paste0("n", 1:ncol(Y))
        row.names(Y) <- paste0("Taxa", 1:nrow(Y))
        data$Y <- Y
        data$X <- cbind(1, data$X[,1])

        ## ILDEx2
        indexa_1_res <- run_indexa(data, mc.samples, epsilon_l[1], epsilon_u[1], denom="none")
        indexa_2_res <- run_indexa(data, mc.samples, epsilon_l[2], epsilon_u[2], denom="none")
        indexa_ctt1_pval <- indexa_1_res[,"ctt.pval"]
        indexa_ctt1_padj <- indexa_1_res[,"ctt.pval.BH.adj"]
        indexa_ctt1_lfc <- indexa_1_res[,"median.effect"]
        indexa_gtt1_pval <- indexa_1_res[,"gtt.pval"]
        indexa_gtt1_padj <- indexa_1_res[,"gtt.pval.BH.adj"]
        indexa_gtt1_lfc <- indexa_1_res[,"median.effect"]
        indexa_ctt2_pval <- indexa_2_res[,"ctt.pval"]
        indexa_ctt2_padj <- indexa_2_res[,"ctt.pval.BH.adj"]
        indexa_ctt2_lfc <- indexa_2_res[,"median.effect"]
        indexa_gtt2_pval <- indexa_2_res[,"gtt.pval"]
        indexa_gtt2_padj <- indexa_2_res[,"gtt.pval.BH.adj"]
        indexa_gtt2_lfc <- indexa_2_res[,"median.effect"]
        ## DESeq2
        deseq2_res <- run_deseq2(data, "poscounts")
        deseq2_pval <- deseq2_res$pvalue
        deseq2_padj <- p.adjust(deseq2_pval, method=p.method)
        deseq2_lfc <- deseq2_res$log2FoldChange
        ## Limma
        limma_res <- run_limma(data)
        limma_pval <- limma_res[,"P.Value"]
        limma_padj <- p.adjust(limma_pval, method=p.method)
        limma_lfc <- limma_res[,"logFC"]
        ## edgeR
        edger_res <- run_edger(data)
        edger_pval <- edger_res[,"PValue"]
        edger_padj <- p.adjust(edger_pval, method=p.method)
        edger_lfc <- edger_res[,"logFC"]
        ## ALDEx2 
        aldex2_res <- run_aldex2(data, mc.samples, p.method)
        aldex2_pval <- aldex2_res[,"we.ep"]
        aldex2_padj <- aldex2_res[,"we.eBH"]
        aldex2_lfc <- aldex2_res[,"effect"]
        
        res_padj <- cbind(deseq2_padj, limma_padj, edger_padj,
                          aldex2_padj, indexa_ctt1_padj,
                          indexa_ctt2_padj, indexa_gtt1_padj,
                          indexa_gtt2_padj)
        res_lfcs <- cbind(deseq2_lfc, limma_lfc, edger_lfc,
                          aldex2_lfc, indexa_ctt1_lfc,
                          indexa_ctt2_lfc, indexa_gtt1_lfc,
                          indexa_gtt2_lfc)
        
        all_padj <- abind(all_padj, res_padj, along=3)
        all_lfcs <- abind(all_lfcs, res_lfcs, along=3)
    }
    list(all_padj=all_padj, all_lfcs=all_lfcs)
}

##set.seed(5314463)
set.seed(3207485)
# Identify which taxa to keep, require 20% non-zeros
no <- 10000
seq_depths <- runif(no, 50000, 200000)
data_other <- sim_dossa2(no, "Stool", seq_depths)
species_names <- row.names(data_other$Y)
keep_taxa <- names(which((rowSums(data_other$Y!=0)/no)>=0.25))
other_taxa <- names(which((rowSums(data_other$Y!=0)/no)<0.25))
print(length(keep_taxa))

## LFCs
lfcs <- abs(rt(65, 4, 1))

subset_species <- sample(keep_taxa, length(lfcs))
prevalences <- lfcs/2
tpi <- match(subset_species, keep_taxa)
fpi <- (1:length(keep_taxa))[-tpi]
true_lfcs <- rep(0, length(keep_taxa)+1)
true_lfcs[tpi] <- lfcs

## Get Scale
seq_depths <- runif(20000, 50000, 200000)
data <- sim_dossa2(20000, "Stool", seq_depths,
                   subset_species, lfcs, prevalences)
mean(log(colSums(data$A)[10001:20000]))-mean(log(colSums(data$A)[1:10000]))

iters <- as.numeric(args[1])
n <- as.numeric(args[2])
res_l <- benchmark_func(n, iters, 300, 20000, 200000, keep_taxa,
                        subset_species, lfcs, prevalences,
                        epsilon_l=c(0,0), epsilon_u=c(0.5,1.75))

saveRDS(res_l, paste0("../output/sparse_dossa_out_", iters,
                      "_", n, "_", p.method, ".RDS"))

