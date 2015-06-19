
signalHSMM_proteins <- read_uniprot("sp1950_2010.txt", euk = TRUE, what = "signal")
b <- lapply(signalHSMM_proteins, function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:50], group_best))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b <- matrix(unlist(b), ncol=5, byrow = TRUE)
# print.xtable(xtable(b, digits = 0), type = "html", include.rownames = FALSE)
# kable(b)
(table(apply(b, 1, function(x) paste0(x, collapse=""))))
(table(apply(b, 1, function(x) paste0(x[1:4], collapse=""))))
(table(apply(b, 1, function(x) paste0(x[4:1], collapse=""))))
(table(apply(b, 1, function(x) paste0(x[1:3], collapse=""))))
(table(apply(b, 1, function(x) paste0(x[3:1], collapse=""))))

sum(b[,1]==2 & b[,3] %in% c(2,4))/length(signalHSMM_proteins)
sum(b[,1]==4 & b[,3] %in% c(2,4))/length(signalHSMM_proteins)
sum(b[,1]%in%c(4) & b[,2] %in% c(1,2,3,4) & b[,3] %in% c(2,4))/length(signalHSMM_proteins)
sum(b[,1]%in%c(2,4) & b[,3] %in% c(4))/length(signalHSMM_proteins)
sum(b[,1]%in%c(3) & b[,2] %in% c(1,2,3) & b[,3]%in%c(2))/length(signalHSMM_proteins)

sum(b[,2]%in%c(1,2) & b[,3] %in% c(2,4) & b[,4]==3)/length(signalHSMM_proteins)
sum(b[,1]%in%c(2,3,4) & b[,2]%in%c(1,2,3) & b[,3] %in% c(2,4) & b[,4]%in%c(1))/length(signalHSMM_proteins)

sum(b[,1]%in%c(2) & b[,2]%in%c(1,2) & b[,3] %in% c(2,3))/length(signalHSMM_proteins)


sum(b[,4]==3, na.rm=T)/length(signalHSMM_proteins)

p4 <- 0
p3 <- 0.3/(1-p4)
p2 <- 0.2/(1-p4)/(1-p3)
p1 <- 0.15/(1-p4)/(1-p3)/(1-p2)

# -3|(3,4)(1,2,3)(2,4)(1,3)
kMer1 <- list(3, c(1,2,3), c(2,4), c(1))
pState1 <- 3
nState1 <- 4
pTrans1 <- p1

# -3|222
kMer2 <- list(2,c(1,2,3), 2)
pState2 <- 3
nState2 <- 4
pTrans2 <- p2

#-3|4(1,2,3,4)(2,4)
kMer3 <- list(4, c(1,2,3,4), c(2,4))
pState3 <- 3
nState3 <- 4
pTrans3 <- p3

# # -2|1_44
# kMer4 <- list(1, c(1,2,3,4), 4,4)
# pState4 <- 3
# nState4 <- 4
# pTrans4 <- p4


parametersSet2 <- add_k_mer_state2(kMer1, signalHsmm2010$pipar, signalHsmm2010$tpmpar, signalHsmm2010$od, 
                                   signalHsmm2010$params, pState1, nState1, pTrans1, d=4)
parametersSet3 <- add_k_mer_state2(kMer2, pipar = parametersSet2$pipar, tpmpar = parametersSet2$tpmpar, 
                                   od = parametersSet2$od, params = parametersSet2$params,
                                   pState2, nState2, pTrans2, d=3)
# parametersSet5 <- parametersSet3
parametersSet4 <- add_k_mer_state2(kMer3, pipar = parametersSet3$pipar, tpmpar = parametersSet3$tpmpar, 
                                   od = parametersSet3$od, params = parametersSet3$params, pState3, nState3, pTrans3, d=3)#
parametersSet5 <- parametersSet4

overall_probs_log <- signalHsmm2010$overall_probs_log

n <- length(benchmark_data2)
a <- matrix(nrow=n, ncol=3)
c <- matrix(nrow=n, ncol=3)
for(i in 1:n){
  prot <- benchmark_data2[[i]]
  deg_sample <- as.numeric(degenerate(toupper(prot)[1L:50], group_best))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  viterbi_res <- duration_viterbi2(deg_sample-1, signalHsmm2010$pipar, signalHsmm2010$tpmpar, 
                                   signalHsmm2010$od, signalHsmm2010$params)
  viterbi_path <- viterbi_res[["path"]]+1
  c_site <- ifelse(any(viterbi_path == 4), 
                   max(which(viterbi_path == 3)), 
                   length(deg_sample))
  #get probabilities of signal peptide model
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  #get probabilities of no signal peptide model
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  
  viterbi_res5 <- duration_viterbi2(deg_sample-1, parametersSet5$pipar, parametersSet5$tpmpar, 
                                    parametersSet5$od, parametersSet5$params)
  viterbi_path5 <- viterbi_res5[["path"]]+1
  c_site5 <- ifelse(any(viterbi_path5 == 4), 
                    min(which(viterbi_path5 == 4))-1, 
                    length(deg_sample))
  c_site5 <- ifelse(any(viterbi_path5 == 8), 
                    min(which(viterbi_path5 == 8)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == 13), 
                    min(which(viterbi_path5 == 13)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == 17), 
                    min(which(viterbi_path5 == 17)), 
                    c_site5)
  #get probabilities of signal peptide model
  prob.signal5 <- viterbi_res5[["viterbi"]][c_site5, viterbi_path5[c_site5]]
  #get probabilities of no signal peptide model
  prob.non5 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site5], 0)
  prob.total5 <- exp(prob.signal5 - prob.non5)
  
  c_site6 <- min(c_site5+1,length(deg_sample))
  prob.signal6 <- viterbi_res5[["viterbi"]][c_site6, viterbi_path5[c_site6]]
  #get probabilities of no signal peptide model
  prob.non6 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site6], 0)
  prob.total6 <- exp(prob.signal6 - prob.non6)
  
  a[i,1] <- 1 - 1/(1 + prob.total)
  a[i,2] <- 1 - 1/(1 + prob.total5)
  a[i,3] <- 1 - 1/(1 + prob.total6)
  c[i,1] <- c_site
  c[i,2] <- c_site5
  c[i,3] <- c_site6
  
}

a <- data.frame(a)
a$sig <- factor(real_labels)

library(ggplot2)
ggplot(a, aes(x=X1,y=X2, color=sig)) + geom_point() + xlab("No k-mer") + ylab("with k-mer") + 
  scale_color_discrete("Signal peptide") +geom_abline(slope=1)

all_preds2 = all_preds
all_preds2$signalKmer = signalKmer=a$X2

bench_metrics2 <- HMeasure(real_labels, all_preds2,
                          threshold = c(rep(0.5, 5), 0.05, 0.05, 0.5))[["metrics"]]
bench_metrics2

library(ROCR)
ROCRpredTest = prediction(a[,1], a[,4])
auc = as.numeric(performance(ROCRpredTest, "auc")@y.values)
auc 

ROCRpredTest = prediction(a[,2], a[,4])
auc = as.numeric(performance(ROCRpredTest, "auc")@y.values)
auc 
  
zle <- which(real_labels==1 & a[,1]>a[,3])
b2 <- lapply(benchmark_full[zle], function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:50], group_best))
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b2 <- matrix(unlist(b2), ncol=5, byrow = TRUE)
table(apply(b2, 1, function(x) paste0(x[3:1], collapse="")))

zle <- which(real_labels==0 & a[,1]<a[,3])
b2 <- lapply(signalHSMM_proteins[zle], function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:50], group_best))
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b2 <- matrix(unlist(b2), ncol=5, byrow = TRUE)
table(apply(b2, 1, function(x) paste0(x[3:1], collapse="")))

