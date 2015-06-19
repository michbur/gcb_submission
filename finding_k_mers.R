#' 
#' tools for heuristic recognition of important k-mers
#' 

library(ggplot2)

signalHSMM_proteins <- read_uniprot("sp1950_2010.txt", euk = TRUE, what = "signal")
b <- lapply(benchmark_full, function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:min(70,length(x))], group_best))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b <- matrix(unlist(b), ncol=5, byrow = TRUE)
table(apply(b, 1, function(x) paste0(x, collapse="")))
table(apply(b, 1, function(x) paste0(x[1:4], collapse="")))
table(apply(b, 1, function(x) paste0(x[4:1], collapse="")))
ggplot(data.frame(table(apply(b, 1, function(x) paste0(x[4:2], collapse="")))),
       aes(x=Var1, y=Freq)) + geom_point()

# 2 (1,2,3,4) (2,3,4)
# (2,3,4) (4) (2,4)
# (2,3,4) (1,2,3) (2,4)
# _ 4 (2,4) (2,3,4) 
table(apply(b, 1, function(x) paste0(x[3:1], collapse="")))
ggplot(data.frame(table(apply(b, 1, function(x) paste0(x[1:3], collapse="")))),
       aes(x=Var1, y=Freq)) + geom_point()

sum(b[,1]%in%c(2,3,4) & b[,2] %in% c(4) & b[,3]%in%c(2,4))/length(signalHSMM_proteins)
sum(b[,1]%in%c(2,3,4) & b[,2] %in% c(1,2,3) & b[,3]%in%c(2,4))/length(signalHSMM_proteins)
sum(b[,2]%in%c(4) & b[,3] %in% c(2,4) & b[,4]%in%c(2,3,4))/length(signalHSMM_proteins)


benchmark_full <- read_uniprot("sp2010_2015.txt", euk = TRUE, what = "signal")
b2 <- lapply(benchmark_full, function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:50], group_best))
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b2 <- matrix(unlist(b2), ncol=5, byrow = TRUE)
table(apply(b2, 1, function(x) paste0(x[3:1], collapse="")))

which(real_labels==0 & a$X2>0.5)
