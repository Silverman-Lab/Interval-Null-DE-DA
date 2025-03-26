library(ggplot2)
source("../general_code/sim_utils.R")
library("ExtDist")
library(gridExtra)
library(ggpubr)

set.seed(3207485)
no <- 10000
seq_depths <- runif(no, 50000, 200000)
data_other <- sim_dossa2(no, "Stool", seq_depths)
species_names <- row.names(data_other$Y)
keep_taxa <- names(which((rowSums(data_other$Y!=0)/no)>=0.25))
other_taxa <- names(which((rowSums(data_other$Y!=0)/no)<0.25))
Y <- data_other$Y[keep_taxa,]
other <- colSums(data_other$Y[other_taxa,])
Y <- rbind(Y, other)
colnames(Y) <- paste0("n", 1:ncol(Y))
row.names(Y) <- paste0("Taxa", 1:nrow(Y))


## LFCs
lfcs <- abs(rt(65, 4, 1))
subset_species <- sample(keep_taxa, length(lfcs))
prevalences <- lfcs/2
tpi <- match(subset_species, keep_taxa)
fpi <- (1:length(keep_taxa))[-tpi]
true_lfcs <- rep(0, length(keep_taxa)+1)
true_lfcs[tpi] <- lfcs

pos <- length(tpi)
neg <- length(fpi)+1 ## other category for +1

rename_vec <- c("LinDA", "ANCOM-BC", "ALDEx2", "limma", "DESeq2", "edgeR",
                "ALDEx2-GTT [0,0.5]", "ALDEx2-GTT [0,1.75]",
                "ALDEx2-CTT [0,0.5]", "ALDEx2-CTT [0,1.75]",
                "ALDEx2-CTT [1,1.75]", "ALDEx2-GTT [1,1.75]")
names(rename_vec) <- c("linda_padj", "ancom_bc_padj", "aldex2_padj", "limma_padj",
                       "deseq2_padj", "edger_padj", "indexa_gtt1_padj",
                       "indexa_gtt2_padj", "indexa_ctt1_padj", "indexa_ctt2_padj",
                       "indexa_ctt3_padj", "indexa_gtt3_padj")

lfc_sorted <- sort(true_lfcs, decreasing=T)
df <- data.frame(lfcs=lfc_sorted)
write.csv(df, "../output/fig_2_lfcs.csv")

fpr_df <- c()
power_df <- c()
reps <- 100
for(n in c(14, 30, 50, 70, 100, 120, 150, 170, 200, 250, 300)) {
    res_l <- readRDS(paste0("../output/sparse_dossa_out_", reps, "_", n, "_BH.RDS"))
    for(m in colnames(res_l$all_padj[,,1])) {
        if(m=="aldex2_0_padj") {next}
        power_avg <- 0
        fdr_avg <- 0
        fp_avg <- 0
        tp_avg <- 0
        for(r in 1:(dim(res_l$all_padj)[3])) {
            col <- res_l$all_padj[,m,r]
            col[is.na(col)] <- 1
            tp <- sum(col[tpi] < 0.05)
            fp <- sum(col[fpi] < 0.05)
            tn <- sum(col[fpi] >= 0.05)
            fn <- sum(col[tpi] >= 0.05)
            if(is.na(tp)&is.na(fp)) {
                ;
            } else if((tp==0)&(fp==0)) {
                ;
            } else {
                tp_avg <- tp_avg + tp
                fp_avg <- fp_avg + fp
                power_avg <- power_avg + (tp/pos)
                fdr_avg <- fdr_avg + (fp/(fp+tp))
            }
        }
        print(m)
        print(c(m, n, fp_avg/reps, tp_avg/reps))
        fpr_df <- rbind(fpr_df, c(n, rename_vec[m], fdr_avg/reps))
        power_df <- rbind(power_df, c(n, rename_vec[m], power_avg/reps))
    }
}

colnames(fpr_df) <- c("n", "method", "fpr")
fpr_df <- data.frame(fpr_df)
fpr_df[,1] <- as.numeric(fpr_df[,1])
#fpr_df[,2] <- as.factor(fpr_df[,2])
fpr_df[,3] <- as.numeric(fpr_df[,3])

g1 <- ggplot(fpr_df, aes(x=n, y=fpr, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 1)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000", "#666666", "#9E9E9E", "#CCCCCC",
                                       "#1f77b4", "#1777b4", "#d62728", "#d62728",
                                        "red", "red", "red", "red"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "dotted", "dotted", "dotted",
                                          "dashed", "solid", "dashed", "solid",
                                          "solid", "solid", "solid", "solid"))
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("False Positive Rate") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

write.csv(ggplot_build(g1)$data[[1]], "../output/geom_smooth_plot_2_fpr_data.csv")

colnames(power_df) <- c("n", "method", "power")
power_df <- data.frame(power_df)
power_df[,1] <- as.numeric(power_df[,1])
#power_df[,2] <- as.factor(power_df[,2])
power_df[,3] <- as.numeric(power_df[,3])

g1 <- ggplot(power_df, aes(x=n, y=power, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 1)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000", "#666666", "#9E9E9E", "#CCCCCC",
                                       "#1f77b4", "#1777b4", "#d62728", "#d62728",
                                        "red", "red", "red", "red"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "dotted", "dotted", "dotted",
                                          "dashed", "solid", "dashed", "solid",
                                          "solid", "solid", "solid", "solid"))
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("False Positive Rate") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

write.csv(ggplot_build(g1)$data[[1]], "../output/geom_smooth_plot_2_power_data.csv")



