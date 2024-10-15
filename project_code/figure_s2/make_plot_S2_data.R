source("../general_code/sim_utils.R")
library(MASS)
library(LaplacesDemon)
library(ExtDist)

rand_str <- function(n) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

set.seed(1264649)
mc.samples <- 300
p.method = "BH"
iters <- 100

LFC <- cbind(rnorm(100), c(runif(40, 0, 3.5), rep(0, 40), runif(20, -2, 0)))
Sigma <- rinvwishart(103, diag(2.2,100))

##hist(rBeta_ab(100000, 8, 1.25, -1.25, 1.25), breaks=100)

n <- 10000
X <- rbind(1, c(rep(0, n/2), rep(1, n/2)))
A <- matrix(0, nrow=100, ncol=n)
for (i in 1:ncol(X)) {
    mu <- (LFC%*%X)[,i]
    A[,i] <- mvrnorm(1, mu, Sigma)
}
## Scale
print(mean(log2(colSums(2^A)[((n/2)+1):n]))-mean(log2(colSums(2^A)[1:(n/2)])))

for(n in c(10, 20, 40, 60, 80, 100, 130, 160, 200)) {
    all_padj <- c()
    all_lfcs <- c()
    for(iter in 1:iters) {
        X <- rbind(1, c(rep(0, n/2), rep(1, n/2)))
        A <- matrix(0, nrow=100, ncol=n)
        for (i in 1:ncol(X)) {
            mu <- (LFC%*%X)[,i]
            A[,i] <- mvrnorm(1, mu, Sigma)
        }
        P <- apply(2^A, 2, function(col) col/sum(col))
        Y <- apply(P, 2, function(col) rmultinom(1, 1000000, col))
        data <- list()
        row.names(Y) <- rand_str(100)
        data$Y <- Y
        data$X <- t(X)

        m <- matrix(rep(rep(c(rep(0, n/2), rep(1, n/2))), mc.samples), ncol=mc.samples)
        gamma_a <- 2^sweep(m, 2, runif(mc.samples, -1.25, 1.25), "*")
        aldex2_un_res <- run_aldex2(data, mc.samples, p.method, gamma=gamma_a)
        aldex2_un_pval <- aldex2_un_res[,"we.ep"]
        aldex2_un_padj <- aldex2_un_res[,"we.eBH"]
        aldex2_un_lfc <- aldex2_un_res[,"effect"]
        
        m <- matrix(rep(rep(c(rep(0, n/2), rep(1, n/2))), mc.samples), ncol=mc.samples)
        gamma_b <- 2^sweep(m, 2, rBeta_ab(mc.samples, 1.25, 8, -1.25, 1.25), "*")
        aldex2_bl_res <- run_aldex2(data, mc.samples, p.method, gamma=gamma_b)
        aldex2_bl_pval <- aldex2_bl_res[,"we.ep"]
        aldex2_bl_padj <- aldex2_bl_res[,"we.eBH"]
        aldex2_bl_lfc <- aldex2_bl_res[,"effect"]       

        m <- matrix(rep(rep(c(rep(0, n/2), rep(1, n/2))), mc.samples), ncol=mc.samples)
        gamma_c <- 2^sweep(m, 2, rBeta_ab(mc.samples, 8, 1.25, -1.25, 1.25), "*")
        aldex2_br_res <- run_aldex2(data, mc.samples, p.method, gamma=gamma_c)
        aldex2_br_pval <- aldex2_br_res[,"we.ep"]
        aldex2_br_padj <- aldex2_br_res[,"we.eBH"]
        aldex2_br_lfc <- aldex2_br_res[,"effect"]
        
        aldex2_res <- run_aldex2(data, mc.samples, p.method)
        aldex2_pval <- aldex2_res[,"we.ep"]
        aldex2_padj <- aldex2_res[,"we.eBH"]
        aldex2_lfc <- aldex2_res[,"effect"]

        ## INDExA
        indexa_res <- run_indexa(data, mc.samples, -1.25, 1.25, denom="none")
        indexa_ctt_pval <- indexa_res[,"ctt.pval"]
        indexa_ctt_padj <- indexa_res[,"ctt.pval.BH.adj"]
        indexa_ctt_lfc <- indexa_res[,"median.effect"]
        indexa_gtt_pval <- indexa_res[,"gtt.pval"]
        indexa_gtt_padj <- indexa_res[,"gtt.pval.BH.adj"]
        indexa_gtt_lfc <- indexa_res[,"median.effect"]

        res_padj <- cbind(aldex2_un_padj, aldex2_bl_padj,
                          aldex2_br_padj, aldex2_padj,
                          indexa_gtt_padj, indexa_ctt_padj)
        res_lfcs <- cbind(aldex2_un_lfc, aldex2_bl_lfc,
                          aldex2_br_lfc, aldex2_lfc,
                          indexa_gtt_lfc, indexa_ctt_lfc)

        all_padj <- abind(all_padj, res_padj, along=3)
        all_lfcs <- abind(all_lfcs, res_lfcs, along=3)
    }
    l <- list(all_padj=all_padj, all_lfcs=all_lfcs)
    saveRDS(l, paste0("../output/supp_fig_1_", iters,
                      "_", n, "_", p.method, ".RDS"))
}

