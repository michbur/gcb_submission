library(dplyr)
library(biogram)

other_soft <- read.csv2("benchmark_other_new.csv")

group_best <- list(`1` = c("r", "n", "d", "q", "e", "h", "k"), 
     `2` = c("g", "p", "s", "t", "y"), `3` = c("i", "l", "m", "f", "w", "v"), 
     `4` = c("a", "c"))

signalHsmm2010 <- train_hsmm(read_uniprot("sp1950_2010.txt", euk = TRUE, what = "signal"), 
                             aa_group = group_best)
signalHsmm1989 <- train_hsmm(read_uniprot("sp1950_1989.txt", euk = TRUE, what = "signal"), 
                             aa_group = group_best)
benchmark_data <- read.fasta("benchmark_data.fasta", seqtype = "AA")
benchmark_data2 <- lapply(benchmark_data, function(i) i[1L:ifelse(length(i) > 80, 80, length(i))])

real_labels <- c(rep(1, 218), rep(0, 218))

all_preds <- data.frame(other_soft,
                       signalHsmm2010 = pred2df(predict(signalHsmm2010, benchmark_data2))[["sp.probability"]],
                       signalHsmm1989 = pred2df(predict(signalHsmm1989, benchmark_data2))[["sp.probability"]])

bench_metrics <- HMeasure(real_labels, all_preds,
                          threshold = c(rep(0.5, 5), 0.05, 0.05))[["metrics"]]

all_preds[["signalHsmm2010"]] <- all_preds[["signalHsmm2010"]] > 0.05
all_preds[["signalHsmm1989"]] <- all_preds[["signalHsmm1989"]] > 0.05

FN_preds <- all_preds[1L:218, ]
rownames(FN_preds[!FN_preds[["signalHsmm2010"]] | !FN_preds[["signalHsmm1989"]], c("signalHsmm2010", "signalHsmm1989")])