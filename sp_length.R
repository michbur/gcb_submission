library(reshape2)

#distributions of lengths of dignal peptides regions

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"


pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), euk = TRUE, what = "signal")

nhc_borders <- t(sapply(pos_seqs, find_nhc))

nhc_len <- sapply(2L:ncol(nhc_borders), function(i)
  nhc_borders[, i] - nhc_borders[, i - 1])

colnames(nhc_len) <- c("n", "h", "c")
mlen <- melt(nhc_len)

levels(mlen[["Var2"]]) <- paste0(levels(mlen[["Var2"]]), "-region")


source("plot_tools.R")

p1 <- ggplot(mlen, aes(x = value)) +
  geom_density(col = "blue", fill = adjustcolor("blue", 0.25)) +
  facet_wrap(~ Var2) +
  scale_x_continuous("Length \n ") + 
  scale_y_continuous("Value") +
  my_theme


png("reglen.png", width = 2257*0.5, height = 801*0.5)
print(p1)
dev.off()
  

