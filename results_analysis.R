library(hmeasure)
library(reshape2)
library(seqinr)
library(dplyr)
library(magrittr)
library(pbapply)

load(paste0(pathway, "fold_res_MB1.RData"))
fold_res1 <- fold_res

load(paste0(pathway, "fold_res_PS1.RData"))
fold_res2 <- fold_res

load(paste0(pathway, "fold_res_PS4.RData"))
fold_res3 <- fold_res

fold_res <- c(fold_res1,
              fold_res2,
              fold_res3)

source("plot_tools.R")


perf_rep <- function(folds, threshold = 0.5)
  do.call(rbind, lapply(1L:length(folds), function(repetition_id) {
    res <- t(sapply(folds[[repetition_id]], function(single_fold)
      rowMeans(sapply(single_fold, function(single_group) {
        metrics <- unlist(HMeasure(as.numeric(!is.na(single_group[, "cs_real"])),
                                   single_group[, "prob"], threshold = threshold)[["metrics"]])
        TP <- metrics["TP"]
        TN <- metrics["TN"]
        FP <- metrics["FP"]
        FN <- metrics["FN"]
        
        c(metrics, 
          mcc = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
      }))))
    
    res <- melt(res)
    
    colnames(res) <- c("encoding", "measure", "value")
    
    res[["encoding"]] <- factor(res[["encoding"]])
    cbind(repetition = factor(rep(repetition_id, nrow(res))), res)
  }))


# optimize cut-off

cutoff_opt <- pblapply(c(0.5, 0.2, 0.1, 0.05, 0.01), function(i)
  perf_rep(fold_res, threshold = i) %>% filter(measure %in% c("AUC", "Sens", "Spec", "mcc")) %>%
    group_by(encoding, measure) %>% summarise(mean_value = mean(value)) %>%
    group_by(measure) %>% summarise(max(mean_value))
)

#0.05 the best specificity/sensitivity


rep_res <- perf_rep(fold_res, 0.05)

# rep_res %>% filter(measure %in% c("Sens")) %>% 
#   group_by(encoding) %>% 
#   summarise(mAUC = mean(value)) %>% 
#   ggplot(aes(x = mAUC)) + geom_density()
# 
# best_encodings <- rep_res %>% filter(measure %in% c("AUC")) %>% 
#   group_by(encoding) %>% 
#   summarise(mAUC = mean(value)) %>% 
#   filter(mAUC > quantile(mAUC, 0.9)) %>%
#   droplevels %>%
#   arrange(desc(mAUC)) %>%
#   select(encoding) %>%
#   unlist %>%
#   as.character %>%
#   as.numeric
#   
# 
# 
# p1 <- rep_res %>% filter(encoding %in% best_encodings, measure %in% c("AUC")) %>%
#   droplevels %>% 
#   mutate(encoding = factor(as.character(encoding), levels = best_encodings)) %>% 
#   ggplot(aes(x = encoding, y = value)) +
#   geom_point(size = 4, position = position_jitter(width = 0.1, height = 0), colour = "blue") +
#   geom_boxplot(fill = adjustcolor("blue", 0.25), colour = "blue") + 
#   scale_x_discrete("Encoding ID\n") + 
#   scale_y_continuous("AUC") + 
#   my_theme

mean_res <- rep_res %>% filter(measure %in% c("AUC", "Sens", "Spec")) %>%
  group_by(encoding, measure) %>% summarise(mean_value = mean(value)) %>% ungroup

best_enc <- lapply(c("AUC", "Sens", "Spec"), function(i)
  mean_res %>% filter(measure == i) %>%
    filter(mean_value > quantile(mean_value, 0.9)) %>%
    droplevels %>%
    #arrange(desc(mean_value)) %>%
    select(encoding) %>%
    unlist %>%
    as.character %>%
    as.numeric)


#levels(p1_dat[["measure"]]) <- c("AUC", "Sensitivity", "Specificity")

# p1 <- ggplot(p1_dat, aes(x = encoding, y = value)) +
#   #geom_boxplot(fill = adjustcolor("blue", 0.25), colour = "blue") + 
#   geom_point(size = 2, position = position_jitter(width = 0.1, height = 0), colour = "blue") +
#   facet_wrap(~ measure) +
#   scale_x_discrete("Encoding ID\n") +
#   scale_y_continuous("Value") + 
#   my_theme

p1_dat <- rep_res %>% filter(#encoding %in% unique(unlist(best_enc)), 
  measure %in% c("AUC", "Spec", "Sens", "mcc")) %>% droplevels %>%
  group_by(encoding, measure) %>% summarize(mean = mean(value)) %>%
  dcast(encoding ~ measure, value.var = "mean")

p1_dat <- p1_dat[!duplicated(p1_dat[, -1]), ]
p1_dat[, "encoding"] <- rep("", nrow(p1_dat))
p1_dat[p1_dat[, "Spec"] > 0.955, "encoding"] <- "2"
#p1_dat[p1_dat[, "Sens"] > 0.85 & p1_dat[, "Spec"] > 0.94, "encoding"] <- "3"
p1_dat[p1_dat[, "Sens"] > 0.93, "encoding"] <- "1"
p1_dat[, "encoding"] <- as.factor(p1_dat[, "encoding"])

p1 <- ggplot(p1_dat, aes(x = Sens, y = Spec, label = encoding, colour = encoding == "", fill = encoding == "")) +
  geom_point(size = 5, shape = 21) +
  geom_text(size = 9, hjust = -0.5, vjust = 0) +
  scale_colour_manual(values = c("red","blue")) + 
  scale_fill_manual(values = c(adjustcolor("red", 0.25), adjustcolor("blue", 0.25))) + 
  scale_x_continuous("Sensitivity\n") +
  scale_y_continuous("Specificity") + 
  my_theme +
  guides(colour = FALSE, fill = FALSE)
  


png("./figures/cvres.png", width = 2257*0.5, height = 1201*0.5)
print(p1)
dev.off()

#caption for cvres
paste0("Results of cross-validation. 1. The encoding providing the best sensitivity (AUC = ", 
       round(mean(p1_dat[p1_dat[, "encoding"] == "1", "AUC"]), 4),
       ", MCC = ", 
       round(mean(p1_dat[p1_dat[, "encoding"] == "1", "mcc"]), 4),
       "). 2. The encoding providing the best specificity (AUC = ", 
       round(mean(p1_dat[p1_dat[, "encoding"] == "2", "AUC"]), 4),
       ", MCC = ", 
       round(mean(p1_dat[p1_dat[, "encoding"] == "2", "mcc"]), 4),
       ").")

#interesting encodings -------------------------------------
int_enc <- as.numeric(rownames(p1_dat[p1_dat[, "encoding"] != "", ]))
group_worst <- all_groups[int_enc][[1]]
#re-arrange group for better comparision
group_worst <- group_worst[c(2, 3, 4, 1)]




