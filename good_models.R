p4 <- 0
p3 <- 0.2/(1-p4)
p2 <- 0.4/(1-p4)/(1-p3)
p1 <- 0/(1-p4)/(1-p3)/(1-p2)

# -3|(2,3,4)(1,2,3)(2,4)(1,3)
sum(b[,1]%in%c(1) & b[,2]%in%c(1,3) & b[,3] %in% c(1,3,4) & b[,4]%in%c(1,3,4))/nrow(b)
kMer1 <- list(1, c(1,3), c(1,3,4), c(1,3,4))
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



p4 <- 0
p3 <- 0.25/(1-p4)
p2 <- 0.35/(1-p4)/(1-p3)
p1 <- 0/(1-p4)/(1-p3)/(1-p2)


##

p4 <- 0.2
p3 <- 0.1/(1-p4)
p2 <- 0.05/(1-p4)/(1-p3)
p1 <- 0.55/(1-p4)/(1-p3)/(1-p2)

# kMer1
sum(b[,1]%in%c(2,3,4) & b[,2]%in%c(1) & b[,3] %in% c(2,4) & b[,4]%in%c(1,2,3,4))/nrow(b)
kMer1 <- list(c(2,3,4), c(1,2,3), c(2,4))
pState1 <- 3
nState1 <- 4
pTrans1 <- p1

# kMer2
sum(b[,1]%in%c(3) & b[,2]%in%c(1) & b[,3] %in% c(2) & b[,4]%in%c(1) & b[,4]%in%c(1,3))/nrow(b)
kMer2 <- list(3, 1, 2, 1, c(1,3))
pState2 <- 3
nState2 <- 4
pTrans2 <- p2

# kMer3
sum(b[,1]%in%c(2) & b[,2]%in%c(1) & b[,3] %in% c(1,4) & b[,4]%in%c(1,2,3,4))/nrow(b)
kMer3 <- list(2, 1, c(1,4))
pState3 <- 3
nState3 <- 4
pTrans3 <- p3

# kMer4
sum(b[,1]%in%c(2,3,4) & b[,2]%in%c(2,3) & b[,3] %in% c(2,4) & b[,4]%in%c(1,2,3,4))/nrow(b)
kMer4 <- list(c(2,3,4), c(2,3), c(2,4))
pState4 <- 3
nState4 <- 4
pTrans4 <- p4