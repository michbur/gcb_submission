library(biogram)
library(signalHsmm)
library(seqinr)
library(hmeasure)
library(dplyr)
library(foreach)
#library(doMC)
other_soft <- read.csv2("benchmark_other.csv")

group_best <- list(`1` = c("r", "n", "d", "q", "e", "h", "k"), 
     `2` = c("g", "p", "s", "t", "y"), `3` = c("i", "l", "m", "f", "w", "v"), 
     `4` = c("a", "c"))

signalHsmm2010 <- train_hsmm(read_uniprot("sp1950_2010.txt", euk = TRUE, what = "signal"), 
                             aa_group = group_best)
signalHsmm1989 <- train_hsmm(read_uniprot("sp1950_1989.txt", euk = TRUE, what = "signal"), 
                             aa_group = group_best)
benchmark_full <- read_uniprot("sp2010_2015.txt", euk = TRUE, what = "signal")
benchmark_data <- read.fasta("benchmark_data.fasta", seqtype = "AA")
benchmark_data2 <- lapply(c(benchmark_full,benchmark_data[-c(1:214)]), function(i) i[1L:ifelse(length(i) > 80, 80, length(i))])

# usuniecie 4 obserwacji negatywnych
negs <- c(432, 215, 280, 417)
benchmark_data2 <- benchmark_data2[-negs]
other_soft <- other_soft[-negs,]
real_labels <- c(rep(1, 214), rep(0, 214))

all_preds <- data.frame(other_soft,
                       signalHsmm2010 = pred2df(predict(signalHsmm2010, benchmark_data2))[["sp.probability"]],
                       signalHsmm1989 = pred2df(predict(signalHsmm1989, benchmark_data2))[["sp.probability"]])

p4 <- 0.1
p3 <- 0.25/(1-p4)
p2 <- 0.35/(1-p4)/(1-p3)
p1 <- 0.1/(1-p4)/(1-p3)/(1-p2)

sum(b[,1]%in%c(3) & b[,2]%in%c(1) & b[,3] %in% c(2) & b[,4]%in%c(1) & b[,4]%in%c(1,3))/nrow(b)
kMer1 <- list(3, 1, 2, 1, c(1,3))
pState1 <- 3
nState1 <- 4
pTrans1 <- p1

# -3|222
sum(b[,1]%in%c(2,3) & b[,2]%in%c(1,2,3) & b[,3] %in% c(2,4))/nrow(b)
kMer2 <- list(c(2,3),c(1,2), c(2,4))
pState2 <- 3
nState2 <- 4
pTrans2 <- p2

#-3|4(1,2,3,4)(2,4)
sum(b[,1]%in%c(4) & b[,2]%in%c(1,2,3) & b[,3] %in% c(2,4))/nrow(b)
kMer3 <- list(4, c(1,2,3), c(2,4))
pState3 <- 3
nState3 <- 4
pTrans3 <- p3

# kMer4
sum(b[,1]%in%c(2) & b[,2]%in%c(1) & b[,3] %in% c(1,4) & b[,4]%in%c(1,2,3,4))/nrow(b)
kMer4 <- list(2, 1, c(1,4))
pState4 <- 3
nState4 <- 4
pTrans4 <- p4

lastState1 = 4+1+length(kMer1)
lastState2 = lastState1+1+length(kMer2)
lastState3 = lastState2+1+length(kMer3)
lastState4 = lastState3+1+length(kMer4)

parametersSet2 <- add_k_mer_state(kMer1, signalHsmm2010$pipar, signalHsmm2010$tpmpar, signalHsmm2010$od, 
                                   signalHsmm2010$params, pState1, nState1, pTrans1, d=3)
parametersSet3 <- add_k_mer_state(kMer2, pipar = parametersSet2$pipar, tpmpar = parametersSet2$tpmpar, 
                                   od = parametersSet2$od, params = parametersSet2$params,
                                   pState2, nState2, pTrans2, d=3)
parametersSet4 <- add_k_mer_state(kMer3, pipar = parametersSet3$pipar, tpmpar = parametersSet3$tpmpar, 
                                   od = parametersSet3$od, params = parametersSet3$params, 
                                   pState3, nState3, pTrans3, d=3)
parametersSet5 <- add_k_mer_state(kMer4, pipar = parametersSet4$pipar, tpmpar = parametersSet4$tpmpar, 
                                   od = parametersSet4$od, params = parametersSet4$params, 
                                   pState4, nState4, pTrans4, d=3)
overall_probs_log <- signalHsmm2010$overall_probs_log
maxSignal <- 45

#registerDoMC(cores=4)
n <- length(benchmark_data2)
results <- foreach(i=1:n) %do% {
  prot <- benchmark_data2[[i]]
  deg_sample <- na.omit(as.numeric(degenerate(toupper(prot)[1L:min(maxSignal, length(prot))], group_best)))
  
  viterbi_res <- duration_viterbi(deg_sample-1, signalHsmm2010$pipar, signalHsmm2010$tpmpar, 
                                   signalHsmm2010$od, signalHsmm2010$params)
  viterbi_path <- viterbi_res[["path"]]+1
  c_site <- ifelse(any(viterbi_path == 4), max(which(viterbi_path == 3)), length(deg_sample))
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  
  viterbi_res5 <- duration_viterbi(deg_sample-1, parametersSet5$pipar, parametersSet5$tpmpar, 
                                    parametersSet5$od, parametersSet5$params)
  viterbi_path5 <- viterbi_res5[["path"]]+1
  c_site5 <- ifelse(any(viterbi_path5 == 4), min(which(viterbi_path5 == 4))-1, length(deg_sample))
  c_site5 <- ifelse(any(viterbi_path5 == lastState1), min(which(viterbi_path5 == lastState1)), c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState2), min(which(viterbi_path5 == lastState2)), c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState3), min(which(viterbi_path5 == lastState3)), c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState4), min(which(viterbi_path5 == lastState4)), c_site5)
  prob.signal5 <- viterbi_res5[["viterbi"]][c_site5, viterbi_path5[c_site5]]
  prob.non5 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site5], 0)
  prob.total5 <- exp(prob.signal5 - prob.non5)
  
  c(1 - 1/(1 + prob.total), 1 - 1/(1 + prob.total5), c_site, c_site5)
}
probs <- t(sapply(results, function (x) head(x, 2)))
cut <- t(sapply(results, function (x) tail(x, 2)))

all_preds$signalKmer = probs[,2]


calc_mcc <- function(metrics) {
  TP <- as.numeric(metrics[["TP"]])
  TN <- metrics[["TN"]]
  FP <- metrics[["FP"]]
  FN <- metrics[["FN"]]
  data.frame(metrics, mcc = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
}

bench_metrics <- calc_mcc(HMeasure(real_labels, all_preds,
                                   threshold = c(rep(0.5, 5), 0.05, 0.05, 0.5))[["metrics"]])
bench_metrics <- data.frame(rownames(bench_metrics), bench_metrics[, c("AUC", "Sens", "Spec", "mcc")])
bench_metrics
colnames(bench_metrics) <- c("Software name", "AUC", "Sensitivity", "Specificity", "MCC")

levels(bench_metrics[["Software name"]]) <- c("Philius \\citep{2008reynoldstransmembrane}", "Phobius \\citep{2004klla}", 
                                              "PrediSi \\citep{2004hillerpredisi}", 
                                              "signalHsmm-1989", 
                                              "signalHsmm-2010",
                                              "signalHsmm-2010 with k-mers",
                                              "signalP 4.1 (no tm) \\citep{2011petersensignalp}", 
                                              "signalP 4.1 (tm) \\citep{2011petersensignalp}")

library(xtable)

rws <- seq(1, nrow(bench_metrics) - 1, by = 2)
col <- rep("\\rowcolor[gray]{0.85}", length(rws))
print(xtable(bench_metrics, caption = "Comparison of Area Under the Curve, Sensitivity, Specificity and Matthews Correlation Coefficient for different classifiers.", 
             label = "tab:bench2010", digits = 4), include.rownames = FALSE, booktabs = TRUE,
      add.to.row = list(pos = as.list(rws), command = col), sanitize.text.function = identity)


# plasmodium ------------------------------------------


other_soft_plas <- read.csv2("benchmark_plas_other.csv")

benchmark_plas_data <- read.fasta("benchmark_plas_data.fasta", seqtype = "AA")

real_labels_plas <- c(rep(1, 102), rep(0, 358))

all_preds_plas <- data.frame(other_soft_plas,
                             signalHsmm2010 = pred2df(predict(signalHsmm2010, benchmark_plas_data))[["sp.probability"]],
                             signalHsmm1989 = pred2df(predict(signalHsmm1989, benchmark_plas_data))[["sp.probability"]])

#registerDoMC(cores=4)
n <- length(benchmark_plas_data)
results <- foreach(i=1:n) %do% {
  prot <- benchmark_plas_data[[i]]
  deg_sample <- na.omit(as.numeric(degenerate(toupper(prot)[1L:min(maxSignal, length(prot))], group_best)))
  
  viterbi_res <- duration_viterbi(deg_sample-1, signalHsmm2010$pipar, signalHsmm2010$tpmpar, 
                                   signalHsmm2010$od, signalHsmm2010$params)
  viterbi_path <- viterbi_res[["path"]]+1
  c_site <- ifelse(any(viterbi_path == 4), max(which(viterbi_path == 3)), length(deg_sample))
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  
  viterbi_res5 <- duration_viterbi(deg_sample-1, parametersSet5$pipar, parametersSet5$tpmpar, 
                                    parametersSet5$od, parametersSet5$params)
  viterbi_path5 <- viterbi_res5[["path"]]+1
  c_site5 <- ifelse(any(viterbi_path5 == 4), min(which(viterbi_path5 == 4))-1, length(deg_sample))
  c_site5 <- ifelse(any(viterbi_path5 == lastState1), min(which(viterbi_path5 == lastState1)), c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState2), min(which(viterbi_path5 == lastState2)), c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState3), min(which(viterbi_path5 == lastState3)), c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState4), min(which(viterbi_path5 == lastState4)), c_site5)
  prob.signal5 <- viterbi_res5[["viterbi"]][c_site5, viterbi_path5[c_site5]]
  prob.non5 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site5], 0)
  prob.total5 <- exp(prob.signal5 - prob.non5)
  
  c(1 - 1/(1 + prob.total), 1 - 1/(1 + prob.total5), c_site, c_site5)
}
probs <- t(sapply(results, function (x) head(x, 2)))
cut <- t(sapply(results, function (x) tail(x, 2)))

all_preds_plas$signalKmer = probs[,2]


bench_metrics <- calc_mcc(HMeasure(real_labels_plas, all_preds_plas,
                                   threshold = c(rep(0.5, 5), 0.05, 0.05, 0.05))[["metrics"]])

bench_metrics <- data.frame(rownames(bench_metrics), bench_metrics[, c("AUC", "Sens", "Spec", "mcc")])
colnames(bench_metrics) <- c("Software name", "AUC", "Sensitivity", "Specificity", "MCC")

levels(bench_metrics[["Software name"]]) <- c("Philius \\citep{2008reynoldstransmembrane}", "Phobius \\citep{2004klla}", 
                                              "PrediSi \\citep{2004hillerpredisi}", 
                                              "signalHsmm-1989", 
                                              "signalHsmm-2010",
                                              "signalHsmm-2010 with k-mers",
                                              "signalP 4.1 (no tm) \\citep{2011petersensignalp}", 
                                              "signalP 4.1 (tm) \\citep{2011petersensignalp}")


rws <- seq(1, nrow(bench_metrics) - 1, by = 2)
col <- rep("\\rowcolor[gray]{0.85}", length(rws))
print(xtable(bench_metrics, caption = "Comparison of Area Under the Curve, H-measure and Matthews Correlation Coefficient for different classifiers considering only proteins belonging to Plasmodiidae.", 
             label = "tab:bench2010plas", digits = 4), include.rownames = FALSE, booktabs = TRUE,
      add.to.row = list(pos = as.list(rws), command = col), sanitize.text.function = identity)
