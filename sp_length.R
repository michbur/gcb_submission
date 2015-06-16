library(reshape2)
library(signalHsmm)

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

whole_sp <- data.frame(value = nhc_borders[, "cs"] - 1)
whole_sp <- cbind(var = rep("Signal peptide", nrow(whole_sp)), whole_sp)

mlen <- cbind(mlen, ann = ifelse(mlen[["Var2"]] == "n-region", "B", ""))

source("plot_tools.R")


p2 <- ggplot(whole_sp, aes(x = value)) + 
  geom_density(col = "blue", fill = adjustcolor("blue", 0.25)) +
  scale_x_continuous("Length \n ") + 
  facet_wrap(~ var) +
  scale_y_continuous("Value") + 
  #annotate("text", x = 5, y = 0.09, label = "A", fontface = "bold") +
  my_theme

p1 <- ggplot(mlen, aes(x = value)) +
  geom_density(col = "blue", fill = adjustcolor("blue", 0.25)) +
  facet_wrap(~ Var2) +
  scale_x_continuous("Length \n ") + 
  scale_y_continuous("Value") +
  #geom_text(aes(x = 5, y = 0.3, label = ann), fontface = "bold") +
  my_theme


png("./figures/reglen.png", width = 2257*0.5, height = 1201*0.5)
print(arrangeGrob(textGrob("A", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), p2, 
                  textGrob("B", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), p1, 
                  nrow = 2, ncol = 2, widths = c(0.05, 0.95)))
dev.off()
  

