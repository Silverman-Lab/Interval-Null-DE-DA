library(MASS)
library(INDExA)
library(ggplot2)
library(tidyr)

# lower -epsilon_u which is -epsilon_l
ci_ctt <- function(Z, conds, epsilon_l, epsilon_u, alpha=0.05) {
    n1 <- sum(conds==0)
    n2 <- sum(conds==1)
    s1 <- sd(Z[which(conds==0)])
    s2 <- sd(Z[which(conds==1)])

    m <- mean(Z[which(conds==1)]) - mean(Z[which(conds==0)])
    sp <- sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
    se <- sp * sqrt(1/n1 + 1/n2)
    df <- n1+n2-2
    return(list(lower=((m-epsilon_u) - abs(qt(alpha/2, df=df))*se),
                upper=((m-epsilon_l) + abs(qt(alpha/2, df=df))*se)))
}

get_critical_value <- function(Z, conds, epsilon_l, epsilon_u, alpha=0.05) {
    n1 <- sum(conds==0)
    n2 <- sum(conds==1)
    s1 <- sd(Z[which(conds==0)])
    s2 <- sd(Z[which(conds==1)])
    em <- mean(c(epsilon_l, epsilon_u))
    epsilon <- epsilon_u - em

    m <- mean(Z[which(conds==1)]) - mean(Z[which(conds==0)])
    sp <- sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2))
    se <- sp * sqrt(1/n1 + 1/n2)
    df <- n1+n2-2

    best_score_i <- 1000
    best_cv_i <- NULL
    for(cv in seq(0,10,0.1)) {
        lower <- (-cv-epsilon) / se
        upper <- (cv-epsilon) / se
        alpha_est <- 1+pt(lower, df)-pt(upper, df)
        score <- abs(alpha_est-alpha)
        if(score<best_score_i) {
            best_score_i <- score
            best_cv_i <- cv
        }
    }
    
    best_score <- 1000
    best_cv <- NULL
    for(cv in seq(max(0, best_cv_i-0.1), best_cv_i+0.1, 0.0001)) {
        lower <- (-cv-epsilon) / se
        upper <- (cv-epsilon) / se
        alpha_est <- 1+pt(lower, df)-pt(upper, df)
        score <- abs(alpha_est-alpha)
        if(score<best_score) {
            best_score <- score
            best_cv <- cv
        }
   }
    return(list(m=m, rej=(m>best_cv) | (m<(-best_cv)), cv=best_cv,
                s=best_score, se=se, df=df))
}

ci_gtt <- function(Z, conds, epsilon_l, epsilon_u, alpha=0.5) {
    em <- mean(c(epsilon_l, epsilon_u))
    g_res <- get_critical_value(Z, conds, epsilon_l, epsilon_u, alpha=0.05)
    cv <- g_res$cv
    m <- g_res$m
    se <- g_res$se
    df <- g_res$df
    return(list(lower=m - cv,
                upper=m + cv,
                lower_adj= m - em - cv,
                upper_adj= m - em + cv))
}

sim.scale <- function(X, eff, avar=0.25) {
    N <- nrow(X)
    D <- nrow(eff)
    Al <- matrix(0, ncol=N, nrow=D) 
    for(n in 1:N) {
        for(d in 1:D) {
            taxa_mean <- 0
            for(col in 1:ncol(X)) {
                taxa_mean <- taxa_mean + X[n,col]*eff[d,col]
            }
            Al[d,n] <- rnorm(1, taxa_mean, avar)
        }
    }
    A <- 2^Al
    P <- apply(A, 2, function(col) col/sum(col))
    S <- colSums(A)
    return(list(S=S, A=A, P=P))
}

set.seed(54563)

## Simulation
D <- 50
Intercept <- runif(D, 4, 14)
treatmentLFC <- c(rep(0, 10), runif(40, -1, 2))

N <- 10
treatment <- c(rep(0, N/2), rep(1, N/2))
other <- rep(c(0, 1), N/2)
X <- cbind(1, treatment)
eff <- cbind(Intercept, treatmentLFC)

no_err <- sim.scale(X, eff, 0)
true_scale <- mean(log2(no_err$S[X[,2]==1]))-mean(log2(no_err$S[X[,2]==0]))
true_p <- mean(log2((no_err$P[25,])[X[,2]==1]))-mean(log2((no_err$P[25,])[X[,2]==0]))

## Coverage
el_f <- -0.5
eu_f <- -0.1
el_t <- 0.1
eu_t <- true_scale

m <- c()
m <- cbind(seq(-4, 4, 0.01), 0, 0, 0, 0)
colnames(m) <- c("e", "fgs", "fcs", "tgs", "tcs")

iters <- 10000
for(i in 1:iters) {
    P <- sim.scale(X, eff, 1)$P
    false_ctt <- ci_ctt(log2(P[25,]), X[,2], -eu_f, -el_f)
    true_ctt <- ci_ctt(log2(P[25,]), X[,2], -eu_t, -el_t)
    false_gtt <- ci_gtt(log2(P[25,]), X[,2], -eu_f, -el_f)
    true_gtt <-ci_gtt(log2(P[25,]), X[,2], -eu_t, -el_t)

    sig_fgs <- c()
    sig_fcs <- c()
    sig_tgs <- c()
    sig_tcs <- c()
    for(e in seq(-4, 4, 0.01)) {
        val <- e
        sig_fgs <- c(sig_fgs, ((false_gtt$lower_adj < val)&(false_gtt$upper_adj)>val))
        sig_fcs <- c(sig_fcs, ((false_ctt$lower < val)&(false_ctt$upper)>val))
        sig_tgs <- c(sig_tgs, ((true_gtt$lower_adj < val)&(true_gtt$upper_adj)>val))
        sig_tcs <- c(sig_tcs, ((true_ctt$lower < val)&(true_ctt$upper)>val))
    }
    m[,2] <- m[,2] + as.numeric(sig_fgs)
    m[,3] <- m[,3] + as.numeric(sig_fcs)
    m[,4] <- m[,4] + as.numeric(sig_tgs)
    m[,5] <- m[,5] + as.numeric(sig_tcs)
}

m[,2] <- m[,2]/iters
m[,3] <- m[,3]/iters
m[,4] <- m[,4]/iters
m[,5] <- m[,5]/iters
m

mdf <- data.frame(m)
plot_data <- mdf %>%
  pivot_longer(cols = c(fgs, fcs, tgs, tcs), 
               names_to = "group", 
               values_to = "val")

g <- ggplot(data=plot_data, aes(x=e, y=val, color=group)) + geom_line()
g <- g + geom_vline(xintercept=true_p+true_scale)
g <- g + geom_vline(xintercept=true_p+true_scale)
g <- g + geom_hline(yintercept=0.95, color="red")
g <- g + geom_vline(xintercept=true_p--eu_t, color="blue")
g <- g + geom_vline(xintercept=true_p--el_t, color="blue")
g <- g + geom_vline(xintercept=true_p--eu_f, color="green")
g <- g + geom_vline(xintercept=true_p--el_f, color="green")

g

