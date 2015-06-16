library(seqinr)
data(aaindex)
library(signalHsmm)
library(seqinr)
library(cvTools)
library(hmeasure)
library(pbapply)

#normalized values of amino acids ----------------------------------------
aa_nvals <- t(sapply(aaindex, function(i) {
  vals <- i[["I"]]
  vals <- vals - min(vals, na.rm = TRUE)
  vals/max(vals, na.rm = TRUE)
}))

colnames(aa_nvals) <- tolower(a(colnames(aa_nvals)))


aa_props <- sapply(aaindex, function(i) i[["D"]])

write.csv2(data.frame(property = unname(aa_props)), file = "aa_props.csv")


traits <- list(size = c(63, 72, 109, 399),
               hydroph = c(54, 68, 151, 244),
               polarity = c(111, 321),
               alpha = c(3, 41, 253))

# clustering of amino acids -------------------------------------------------

all_traits_combn <- expand.grid(traits)



all_groups <- apply(all_traits_combn, 1, function(single_trait_combn) {
  #use ward, because Piotr does it
  cl <- hclust(dist(t(aa_nvals[single_trait_combn, ])), method = "ward.D2")
  gr <- cutree(cl, k = 4)
  names(gr) <- tolower(names(gr))
  agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
  names(agg_gr) <- 1L:length(agg_gr)
  agg_gr
})

# read data ---------------------------------

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "tobit" )
  pathway <- "/home/piotr//Dropbox/signal-peptide2_data (1)//"

if(Sys.info()["nodename"] == "sobczyk-pc"  )
  pathway <- "/home/sobczyk/Dropbox//signal-peptide2_data (1)//"


pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), euk = TRUE, what = "signal")
neg_seqs <- read.fasta(paste0(pathway, "nonsignal_peptides.fasta"), seqtype = "AA")
#remove sequences with atypical aminoacids
atyp_aa <- which(sapply(neg_seqs, function(i) any(i %in% c("X", "J", "Z", "B", "U"))))
too_short <- which(sapply(neg_seqs, length) < 80)
neg_seqs <- neg_seqs[-unique(c(atyp_aa, too_short))]

too_short <- which(sapply(pos_seqs, length) < 80)
pos_seqs <- pos_seqs[-c(too_short)]


# cross-validation ---------------------------------

fold_res <- pblapply(1L:20, function(dummy) {
  #assure the same proteins in each fold for each repetition
  pos_ids <- cvFolds(length(pos_seqs), K = 5)
  cv_neg <- neg_seqs[sample(1L:length(neg_seqs), length(pos_seqs))]
  lapply(all_groups, function(agg_group) {
    lapply(1L:5, function(fold) {
      model_cv <- train_hsmm(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] != fold]], agg_group)
      test_dat <- c(pos_seqs[pos_ids[[4]][,][pos_ids[[5]] == fold]],
                    cv_neg[pos_ids[[4]][,][pos_ids[[5]] == fold]])
      preds <- cbind(t(sapply(predict(model_cv, test_dat), function(single_pred)
        c(prob = single_pred[["sp_probability"]], cs_pred = single_pred[["sp_end"]]))),
        cs_real = sapply(test_dat, function(i) 
          ifelse(is.null(attr(i, "sig")[2]), NA, attr(i, "sig")[2])))
      preds
    })
  })
})


if(Sys.info()["nodename"] == "MICHALKOMP" )
  output <- "fold_res_MB1.RData"

if(Sys.info()["nodename"] == "phobos" )
  output <- "fold_res_MB1.RData"

if(Sys.info()["nodename"] == "tobit" )
  output <- "fold_res_PS1.RData"

if(Sys.info()["nodename"] == "sobczyk-pc"  )
  output <- "fold_res_PS1.RData"

save(fold_res, all_groups, file = paste0(pathway, output))
