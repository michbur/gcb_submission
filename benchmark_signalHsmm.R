library(biogram)
library(signalHsmm)
library(seqinr)
library(hmeasure)
library(dplyr)
other_soft <- read.csv2("benchmark_other.csv")

group_best <- list(`1` = c("r", "n", "d", "q", "e", "h", "k"), 
     `2` = c("g", "p", "s", "t", "y"), `3` = c("i", "l", "m", "f", "w", "v"), 
     `4` = c("a", "c"))

signalHsmm2010 <- train_hsmm(read_uniprot("sp1950_2010.txt", euk = TRUE, what = "signal"), 
                             aa_group = group_best)
signalHsmm1989 <- train_hsmm(read_uniprot("sp1950_1989.txt", euk = TRUE, what = "signal"), 
                             aa_group = group_best)
benchmark_data <- read.fasta("benchmark_data.fasta", seqtype = "AA")
benchmark_data2 <- lapply(benchmark_data, function(i) i[1L:ifelse(length(i) > 80, 80, length(i))])

real_labels <- c(rep(1, 214), rep(0, 218))

all_preds <- data.frame(other_soft,
                       signalHsmm2010 = pred2df(predict(signalHsmm2010, benchmark_data2))[["sp.probability"]],
                       signalHsmm1989 = pred2df(predict(signalHsmm1989, benchmark_data2))[["sp.probability"]])

calc_mcc <- function(metrics) {
  TP <- as.numeric(metrics[["TP"]])
  TN <- metrics[["TN"]]
  FP <- metrics[["FP"]]
  FN <- metrics[["FN"]]
  data.frame(metrics, mcc = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
}

bench_metrics <- calc_mcc(HMeasure(real_labels, all_preds,
                                   threshold = c(rep(0.5, 5), 0.05, 0.05))[["metrics"]])

bench_metrics <- data.frame(rownames(bench_metrics), bench_metrics[, c("AUC", "Sens", "Spec", "mcc")])
colnames(bench_metrics) <- c("Software name", "AUC", "Sensitivity", "Specificity", "MCC")

levels(bench_metrics[["Software name"]]) <- c("Philius \\citep{2008reynoldstransmembrane}", "Phobius \\citep{2004klla}", 
                                              "PrediSi \\citep{2004hillerpredisi}", 
                                              "signalHsmm-1989", 
                                              "signalHsmm-2010", 
                                              "signalP 4.1 (no tm) \\citep{2011petersensignalp}", 
                                              "signalP 4.1 (tm) \\citep{2011petersensignalp}")

library(xtable)

rws <- seq(1, nrow(bench_metrics) - 1, by = 2)
col <- rep("\\rowcolor[gray]{0.85}", length(rws))
print(xtable(bench_metrics, caption = "Performance measures for the best encoding. 60 repetitions of cross-validation.", 
             label = "tab:bench2010", digits = 4), include.rownames = FALSE, booktabs = TRUE,
      add.to.row = list(pos = as.list(rws), command = col), sanitize.text.function = identity)


# plasmodium ------------------------------------------

other_soft_plas <- read.csv2("benchmark_plas_other.csv")[, -1]

benchmark_plas_data <- read.fasta("benchmark_plas_data.fasta", seqtype = "AA")

real_labels_plas <- c(rep(1, 102), rep(0, 358))

all_preds_plas <- data.frame(other_soft_plas,
                             signalHsmm2010 = pred2df(predict(signalHsmm2010, benchmark_plas_data))[["sp.probability"]],
                             signalHsmm1989 = pred2df(predict(signalHsmm1989, benchmark_plas_data))[["sp.probability"]])

bench_metrics <- calc_mcc(HMeasure(real_labels_plas, all_preds_plas,
                                   threshold = c(rep(0.5, 4), 0.05, 0.05))[["metrics"]])