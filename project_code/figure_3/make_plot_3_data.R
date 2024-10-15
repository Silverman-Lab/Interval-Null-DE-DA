source("../general_code/sim_utils.R")
library(biomformat)
library(abind)
library(stringr)
library(lmerTest)

## Load sequence count data
biom_data <- read_hdf5_biom("oral_trimmed_deblur.biom")
Y <- do.call(rbind, biom_data$data)
row.names(Y) <- sapply(biom_data$rows, function(item){item$id})

## Get metadata, add average flowcount if needed
metadata <- read.csv("oral_trimmed_metadata.csv", sep="\t")
subset_metadata <- metadata[,c("HostSubject", "Timepoint.", "X.SampleID", "brushing_event", "flow.cell.5min.1", "flow.cell.5min.2")]
subset_metadata[,"flowcount"] <- round(rowMeans(subset_metadata[,c("flow.cell.5min.1", "flow.cell.5min.2")]))

## Add taxonomy info
tax <- read.csv("taxonomy.tsv", sep="\t")
inds <- match(row.names(Y), tax[,"Feature.ID"])
row.names(Y) <- tax[inds,"Taxon"]

## Collapse to Genus Level
Y_genus <- c()
genus <- unname(sapply(row.names(Y), function(x) {strsplit(x, split=";")[[1]][6]}))
unique_genus <- unique(genus[!is.na(genus)])
all_inds <- c()
for(genus_n in unique_genus) {
    genus_inds <- which(genus%in%genus_n)
    all_inds <- c(all_inds, genus_inds)
    if(length(genus_inds)>1) {
        new_row <- colSums(Y[genus_inds,])
    } else {
        new_row <- Y[genus_inds,]
    }
    Y_genus <- rbind(Y_genus, new_row)
}
row.names(Y_genus) <- unique_genus

## Collapse rows with less than a third of samples having counts
other <- colSums(Y[which(is.na(genus)),])
other <- other + colSums(Y_genus[rowSums(Y_genus>0)<12,])
Y_genus <- Y_genus[rowSums(Y_genus>0)>=12,]
other <- other + Y_genus[" g__",]
Y_genus <- Y_genus[row.names(Y_genus)!=" g__",]
Y_genus <- rbind(Y_genus, other)

## Build Metadata Matrix for Mixed Effects Model
meta_mat <- cbind(c(log2(subset_metadata[,"flow.cell.5min.1"]),
                    log2(subset_metadata[,"flow.cell.5min.2"])),
                  c(subset_metadata[,"Timepoint."],
                    subset_metadata[,"Timepoint."]),
                  c(subset_metadata[,"HostSubject"],
                    subset_metadata[,"HostSubject"]),
                  c(subset_metadata[,"brushing_event"],
                    subset_metadata[,"brushing_event"]))
colnames(meta_mat) <- c("flowcounts", "time_point", "Host", "brushing_event")
meta_mat <- data.frame(meta_mat)
meta_mat$brushing_event <- factor(meta_mat$brushing_event, levels=c("before", "after"))
meta_mat$flowcounts <- as.numeric(meta_mat$flowcounts)
meta_mat$time_point_me <- apply(meta_mat, 1, function(row) {
    if(row[2]%in%c(1,2)) {
        return("morning")
    } else {
        return("evening")
    }
})

## Mixed Effects Model with 95% Confidence Interval
mixed_effect_model <- lmer(flowcounts~brushing_event+(1|Host)+(1|time_point_me), meta_mat)
z <- lm(flowcounts~brushing_event, data=meta_mat)
1-(2^(coef(z)[2]))
scale_ci <- confint(mixed_effect_model, level=0.95)[5,]
summary(mixed_effect_model)
1-(2^scale_ci)
log2(1-0.2572)

pt((-0.73--1.1583)/0.3695, df=51.45)-pt((-2.19--1.1583)/0.3695, df=51.45)
pt((0.87--1.1583)/0.3695, df=51.45)-pt((-0.59--1.1583)/0.3695, df=51.45)
pt((1.57--1.1583)/0.3695, df=51.45)-pt((0.11--1.1583)/0.3695, df=51.45)

pnorm((-0.73--1.1583)/0.3695)-pnorm((-2.19--1.1583)/0.3695)
pnorm((0.87--1.1583)/0.3695)-pnorm((-0.59--1.1583)/0.3695)
pnorm((1.57--1.1583)/0.3695)-pnorm((0.11--1.1583)/0.3695)

pt((-1.1583--0.43)/0.3695, df=51)-pt((-1.1583--1.89)/0.3695, df=51, lower.tail=T)
2*pt((-1.1583--2.19)/0.3695, df=51, lower.tail=F)

set.seed(54623)
## Run Analyses
data <- list()
data$Y <- Y_genus
data$X <- cbind(1, as.numeric(subset_metadata[,"brushing_event"]=="after"))

pt((-1.1583-log2(1-0.34))/0.3695, 51.45)
pt((-1.1583-log2(1+0.1))/0.3695, 51.45)

## Gold Standard Using CI & Aldex2 Sens Interval Test
aldex2_gold <- run_indexa(data, 1000, scale_ci[1], scale_ci[2], denom="none")
aldex2_gold[(aldex2_gold[,c("ctt.pval.BH.adj")]<0.05)|(aldex2_gold[,c("gtt.pval.BH.adj")]<0.05),]
aldex2_gold[c(" g__Rothia", " g__Selenomonas", " g__Schwartzia"),]

## Existing Methods
## DESeq2
deseq2_res <- run_deseq2(data, "poscounts")
## Limma
limma_res <- run_limma(data)
## edgeR
edger_res <- run_edger(data)
## ALDEx2
aldex2_res <- run_aldex2(data, 1000, "BH")

aldex2_lib <- run_indexa(data, 1000, log2(1-0.99), log2(1-0.15), denom="none")
aldex2_con <- run_indexa(data, 1000, log2(1-0.99), log2(1-0.25), denom="none")

## Sensitivity Analysis
sens_res <- c()
for(j in seq(-0.5, 2, 0.01)) {
    print(j)
    lower <- scale_ci[1]+j
    upper <- scale_ci[2]+j
    one_res <- run_indexa(data, 1000, lower, upper, denom="none")
    sens_res <- rbind(sens_res, cbind(j, one_res))
}

write.csv(sens_res, "../output/figure_3b_data.csv")

gm_scale <- apply((data$Y+0.5), 2, function(col) {
    p <- col/sum(col)
    -log2(gm_mean(p))
})
mean(gm_scale[data$X[,2]==1])-mean(gm_scale[data$X[,2]==0])

## False/True Positive/Negative
## Helper Function
get_tf_np <- function(obs_lfc, obs_pvals, gold_standard_lfc, tp_rows, tn_rows) {
    res <- rep("", length(obs_lfc))
    tp <- (sign(obs_lfc[tp_rows])==sign(gold_standard_lfc[tp_rows]))
    tp <- tp&(obs_pvals[tp_rows]<0.05)
    res[tp_rows][tp] <- "tp"
    res[tp_rows][which(obs_pvals[tp_rows]>=0.05)] <- "fn"
    res[tn_rows][which(obs_pvals[tn_rows]<0.05)] <- "fp"
    res[tn_rows][which(obs_pvals[tn_rows]>=0.05)] <- "tn"
    return(res)
}
tp_tn_fp_fn <- c()

tp_rows <- which(row.names(aldex2_gold)%in%row.names(aldex2_gold[aldex2_gold[,"gtt.pval.BH.adj"] < 0.05,]))
tn_rows <- which(row.names(aldex2_gold)%in%row.names(aldex2_gold[aldex2_gold[,"gtt.pval.BH.adj"] >= 0.05,]))

tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_gold[,"median.effect"],
                                            aldex2_gold[,"gtt.pval.BH.adj"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(deseq2_res[,"log2FoldChange"],
                                            deseq2_res[,"padj"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(limma_res[,"logFC"],
                                            limma_res[,"adj.P.Val"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(edger_res[,"logFC"],
                                            edger_res[,"FDR"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_res[,"effect"],
                                            aldex2_res[,"we.eBH"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_lib[,"median.effect"],
                                            aldex2_lib[,"gtt.pval.BH.adj"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))
tp_tn_fp_fn <- cbind(tp_tn_fp_fn, get_tf_np(aldex2_con[,"median.effect"],
                                            aldex2_con[,"gtt.pval.BH.adj"],
                                            aldex2_gold[,"median.effect"], tp_rows, tn_rows))

colnames(tp_tn_fp_fn) <- c("gold_standard", "deseq2", "limma", "edger", "aldex2",
                           "aldex2_lib", "aldex2_con")
row.names(tp_tn_fp_fn) <- row.names(aldex2_gold)
write.csv(tp_tn_fp_fn, "../output/figure_3_data.csv")

