
  p4 <- 0
  p3 <- 0.34/(1-p4)
  p2 <- 0.26/(1-p4)/(1-p3)
  p1 <- 0/(1-p4)/(1-p3)/(1-p2)
  
  # -3|(2,3,4)(1,2,3)(2,4)(1,3)
  sum(b[,1]%in%c(1,3) & b[,2]%in%c(1,2,3) & b[,3] %in% c(1,2,3,4) & b[,4]%in%c(1,3,4))/nrow(b)
  kMer1 <- list(c(1,3), c(1,2,3), c(1,2,3,4), c(1,3,4))
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
  
#   # -2|1_44
#   sum(b[,1]%in%c(2) & b[,2]%in%c(2,3) & b[,3] %in% c(2,4))/length(signalHSMM_proteins)
#   kMer4 <- list(2, c(2,3), c(2,4))
#   pState4 <- 3
#   nState4 <- 4
#   pTrans4 <- p4
  
  
  parametersSet2 <- add_k_mer_state2(kMer1, signalHsmm2010$pipar, signalHsmm2010$tpmpar, signalHsmm2010$od, 
                                     signalHsmm2010$params, pState1, nState1, pTrans1, d=4)
  parametersSet3 <- add_k_mer_state2(kMer2, pipar = parametersSet2$pipar, tpmpar = parametersSet2$tpmpar, 
                                     od = parametersSet2$od, params = parametersSet2$params,
                                     pState2, nState2, pTrans2, d=3)
  # parametersSet5 <- parametersSet3
  parametersSet4 <- add_k_mer_state2(kMer3, pipar = parametersSet3$pipar, tpmpar = parametersSet3$tpmpar, 
                                     od = parametersSet3$od, params = parametersSet3$params, 
                                     pState3, nState3, pTrans3, d=3)
  parametersSet5 <- parametersSet4
#   parametersSet5 <- add_k_mer_state2(kMer4, pipar = parametersSet4$pipar, tpmpar = parametersSet4$tpmpar, 
#                                      od = parametersSet4$od, params = parametersSet4$params, 
#                                      pState4, nState4, pTrans4, d=3)
  
  overall_probs_log <- signalHsmm2010$overall_probs_log
  
  n <- length(benchmark_data2)
  a <- matrix(nrow=n, ncol=3)
  c <- matrix(nrow=n, ncol=3)
  for(i in 1:n){
    prot <- benchmark_data2[[i]]
    deg_sample <- as.numeric(degenerate(toupper(prot)[1L:min(43, length(prot))], group_best))
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
    c_site5 <- ifelse(any(viterbi_path5 == 21), 
                      min(which(viterbi_path5 == 21)), 
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
  
  ggplot(a, aes(x=X1,y=X2, color=sig)) + geom_point() + xlab("No k-mer") + ylab("with k-mer") + 
    scale_color_discrete("Signal peptide") +geom_abline(slope=1)
  
  all_preds2 = all_preds
  all_preds2$signalKmer = signalKmer=a$X2
  
  bench_metrics2 <- HMeasure(real_labels, all_preds2,
                             threshold = c(rep(0.5, 5), 0.1, 0.1, 0.5))[["metrics"]]
  bench_metrics2

library(ROCR)
ROCRpredTest = prediction(a[,1], a[,4])
auc = as.numeric(performance(ROCRpredTest, "auc")@y.values)
auc 

ROCRpredTest = prediction(a[,2], a[,4])
auc = as.numeric(performance(ROCRpredTest, "auc")@y.values)
auc 

zle <- which(real_labels==1 & a[,2]<0.75)
zle <- which(real_labels==1 & a[,2]<a[,1] & a[,2]<0.99)
b2 <- lapply(benchmark_full[zle], function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:min(50, length(x))], group_best))
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b2 <- matrix(unlist(b2), ncol=5, byrow = TRUE)
cbind(a[zle,], seq=apply(b2, 1, function(x) paste0(x, collapse="")), 
      csite=sapply(benchmark_full[zle], function(x) cut <- attr(x, "sig")[2]))

which(real_labels==0 & a[,2]>0.5)
a[which(real_labels==0 & a[,2]>0.5),1]
  a[which(real_labels==0 & a[,2]>0.5),2]
