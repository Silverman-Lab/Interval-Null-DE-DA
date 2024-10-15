library(ggplot2)
source("../general_code/sim_utils.R")
library("ExtDist")
library(gridExtra)
library(ggpubr)

set.seed(1264649)
mc.samples <- 300
p.method = "BH"
iters <- 100

LFC <- cbind(rnorm(100), c(runif(40, 0, 3.5), rep(0, 40), runif(20, -2, 0)))

lfcs <- LFC[,2]
fpi <- 41:80
tpi <- c(1:40, 81:100)

pos <- length(tpi)
neg <- length(fpi)

rename_vec <- c("ALDEx2-SM Unif", "ALDEx2-SM left-skew",
                "ALDEx2-SM right-skew", "ALDEx2", "ALDEx2-CTT",
                "ALDEx2-GTT")
names(rename_vec) <- c("aldex2_un_padj", "aldex2_bl_padj", "aldex2_br_padj",
                       "aldex2_padj", "indexa_ctt_padj", "indexa_gtt_padj")

lfc_sorted <- sort(lfcs, decreasing=T)
df <- data.frame(lfcs=lfc_sorted)
write.csv(df, "../output/sfig_1_lfcs.csv")

fpr_df <- c()
power_df <- c()
f1_df <- c()
f05_df <- c()
for(n in c(10,20,40,60,80,100,130,160,200)) {
    res_l <- readRDS(paste0("../output/supp_fig_1_100_", n, "_BH.RDS"))
    for(m in colnames(res_l$all_padj[,,1])) {
        if(m=="aldex2_0_padj") {next}
        power_avg <- 0
        fdr_avg <- 0
        f1_avg <- 0
        f05_avg <- 0
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
                power_avg <- power_avg + (tp/pos)
                ##fdr_avg <- fdr_avg + (fp/(fp+tn))
                fdr_avg <- fdr_avg + (fp/(fp+tp))
                f1_avg <- f1_avg + (tp/(tp+(0.5*(fp+fn))))
                f05_avg <- f05_avg + ( (1.25*(tp)) / ((1.25)*tp + 0.25*fn+fp) )
            }
        }
        fpr_df <- rbind(fpr_df, c(n, rename_vec[m], fdr_avg/100))
        power_df <- rbind(power_df, c(n, rename_vec[m], power_avg/100))
        f1_df <- rbind(f1_df, c(n, rename_vec[m], f1_avg/100))
        f05_df <- rbind(f05_df, c(n, rename_vec[m], f05_avg/100))
    }
}


colnames(f05_df) <- c("n", "method", "f05")
f05_df <- data.frame(f05_df)
f05_df[,1] <- as.numeric(f05_df[,1])
f05_df[,3] <- as.numeric(f05_df[,3])
write.csv(f05_df, "../output/fig_s1_f05_25_2.csv")

g1 <- ggplot(f05_df, aes(x=n, y=f05, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 1)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000", "#666666", "#9E9E9E", "#CCCCCC",
                                       "#1f77b4", "#1777b4", "#d62728", "#d62728"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "dotted", "dotted", "dotted",
                                          "dashed", "solid", "dashed", "solid"))
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("False Positive Rate") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

write.csv(ggplot_build(g1)$data[[1]], "../output/fig_s1_geom_smooth_plot_2_data.csv")

colnames(fpr_df) <- c("n", "method", "fpr")
fpr_df <- data.frame(fpr_df)
fpr_df[,1] <- as.numeric(fpr_df[,1])
fpr_df[,3] <- as.numeric(fpr_df[,3])
write.csv(fpr_df, "../output/fig_s1_fpr_25_2.csv")

g1 <- ggplot(fpr_df, aes(x=n, y=fpr, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 1)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000", "#666666", "#9E9E9E", "#CCCCCC",
                                       "#1f77b4", "#1777b4", "#d62728", "#d62728"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "dotted", "dotted", "dotted",
                                          "dashed", "solid", "dashed", "solid"))
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("False Positive Rate") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

write.csv(ggplot_build(g1)$data[[1]], "../output/fig_s1_geom_smooth_plot_2_fpr_data.csv")


colnames(power_df) <- c("n", "method", "power")
power_df <- data.frame(power_df)
power_df[,1] <- as.numeric(power_df[,1])
power_df[,2] <- factor(power_df[,2])
power_df[,3] <- as.numeric(power_df[,3])
write.csv(power_df, "../output/fig_s1_power_25_2.csv")

g1 <- ggplot(power_df, aes(x=n, y=power, color=method, linetype=method))
g1 <- g1 + geom_smooth(alpha=0, size=1.5)
g1 <- g1 + theme_bw() + ylim(0, 1)
g1 <- g1 + geom_hline(yintercept=0.05, linetype="dotted", alpha=0.5)
g1 <- g1 + scale_color_manual(values=c("#000000", "#666666", "#9E9E9E", "#CCCCCC",
                                       "#1f77b4", "#1777b4", "#d62728", "#d62728"))
g1 <- g1 + scale_linetype_manual(values=c("dotted", "dotted", "dotted", "dotted",
                                          "dashed", "solid", "dashed", "solid"))
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + ylab("False Positive Rate") + xlab("Sample Size")
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + theme(legend.key.size = unit(1, "cm"))

write.csv(ggplot_build(g1)$data[[1]], "../output/fig_s1_geom_smooth_plot_2_power_data.csv")

