library(hmeasure)
library(reshape2)
library(seqinr)
library(dplyr)
library(magrittr)
library(pbapply)

load(paste0(pathway, "fold_res_MB1.RData"))

source("plot_tools.R")


perf_rep <- function(folds, threshold = 0.5)
  do.call(rbind, lapply(1L:length(folds), function(repetition_id) {
    res <- t(sapply(folds[[repetition_id]], function(single_fold)
      rowMeans(sapply(single_fold, function(single_group) {
        unlist(HMeasure(as.numeric(!is.na(single_group[, "cs_real"]) - 1),
                        single_group[, "prob"], threshold = threshold)[["metrics"]])
        }))))
    
    res <- melt(res)
    
    colnames(res) <- c("encoding", "measure", "value")
    
    res[["encoding"]] <- factor(res[["encoding"]])
    cbind(repetition = factor(rep(repetition_id, nrow(res))), res)
  }))




cutoff_opt <- pblapply(c(0.5, 0.8, 0.9, 0.95, 0.99), function(i)
  perf_rep(fold_res, threshold = i) %>% filter(measure %in% c("AUC", "Sens", "Spec")) %>%
    group_by(encoding, measure) %>% summarise(mean_value = mean(value)) %>%
    group_by(measure) %>% summarise(max(mean_value))
)

#0.95 the best specificity/sensitivity


rep_res <- perf_rep(fold_res, 0.95)

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

p1_dat <- rep_res %>% filter(encoding %in% unique(unlist(best_enc)), measure %in% c("AUC", "Spec", "Sens")) %>% droplevels %>%
  group_by(encoding, measure) %>% summarize(mean = mean(value)) %>%
  dcast(encoding ~ measure, value.var = "mean")

p1 <- ggplot(p1_dat, aes(x = Sens, y = Spec, label = encoding)) +
  geom_point(size = 5) +
  #geom_text(size = 9, position = "jitter") +
  scale_x_continuous("Sensitivity") +
  scale_y_continuous("Specificity\n") + 
  my_theme +
  coord_flip()


png("./figures/cvres.png", width = 2257*0.5, height = 1201*0.5)
print(p1)
dev.off()



