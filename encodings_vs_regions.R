library(reshape2)
library(signalHsmm)
library(biogram)

#distributions of lengths of dignal peptides regions

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"


pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), euk = TRUE, what = "signal")

#interesting encodings -------------------------------------
#p1_dat comes from results_analysis.R
int_enc <- as.numeric(rownames(p1_dat[p1_dat[, "encoding"] != "", ]))
group_best <- all_groups[int_enc][[1]]
#re-arrange group for better comparision
group_best <- group_best[c(2, 3, 4, 1)]
names(group_best) <- 1L:4

group_worst <- all_groups[int_enc][[2]]

nhc_borders <- t(sapply(pos_seqs, find_nhc))

deg_pos <- lapply(list(group_best, group_worst), function(single_encoding)
  lapply(1L:length(pos_seqs), function(seq_id) {
    single_seq <- pos_seqs[[seq_id]]
    aa_n <- single_seq[1L:(nhc_borders[seq_id, "start_h"] - 1)]
    aa_h <- single_seq[nhc_borders[seq_id, "start_h"]:(nhc_borders[seq_id, "start_c"] - 1)]
    aa_c <- single_seq[nhc_borders[seq_id, "start_c"]:(nhc_borders[seq_id, "cs"])]
    aa_mature <- single_seq[(nhc_borders[seq_id, "cs"] + 1):length(single_seq)]
    res <- list(aa_n = aa_n,
                aa_h = aa_h,
                aa_c = aa_c,
                aa_m = aa_mature)
    lapply(res, function(i) 
      degenerate(i, single_encoding))
  }))

deg_region <- do.call(rbind, lapply(deg_pos, function(single_encoding)
  cbind(region = unlist(lapply(c("n", "h", "c", "m"), function(i) rep(i, 4))), 
        do.call(rbind, lapply(1L:4, function(single_region) {
          res <- factor(unlist(lapply(single_encoding, function(single_seq)
            single_seq[[single_region]])), levels = as.character(1L:4))
          res_tab <- data.frame(table(res))
          cbind(res_tab, prop = res_tab[["Freq"]]/length(res))
        })))))

deg_region <- cbind(enc = unlist(lapply(c("best", "worst"), function(i) rep(i, 16))),
                    deg_region)

colnames(deg_region) <- c("enc", "region", "group", "count", "freq")

deg_region[["region"]] <- factor(deg_region[["region"]], levels = c("n", "h", "c", "m"))

source("plot_tools.R")

levels(deg_region[["region"]]) <- c("n-region", "h-region", "c-region", "Mature\nprotein")
levels(deg_region[["group"]]) <- paste0("Group ", levels(deg_region[["group"]]))

p3 <- ggplot(deg_region, aes(x = region, y = freq, fill = enc, colour = enc)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ group, nrow = 1) + 
  scale_y_continuous("Frequency") +
  scale_x_discrete("Region\n") +
  scale_colour_manual("Encoding: ", values = c("red", "blue")) +
  scale_fill_manual("Encoding: ", values = c(adjustcolor("red", 0.25), adjustcolor("blue", 0.25))) + 
  guides(colour = FALSE) +
  my_theme


#figure comparing encodings ---------------------------
data(aaindex)
group_properties <- function(group) {
  res <- do.call(rbind, lapply(1L:length(group), function(subgroup_id)
    melt(data.frame(group = paste0("Group ", as.character(rep(subgroup_id, 4))), critertion = c("size", "hydroph", "polarity", "alpha"),
                    aa_nvals[unlist(all_traits_combn[int_enc[2], ]), group[[subgroup_id]]]))))
  levels(res[["variable"]]) <- toupper(levels(res[["variable"]]))
  res
}

dat_best <- group_properties(group_best)
dat_worst <- group_properties(group_worst)

dat_bestworst <- cbind(enc = unlist(lapply(c("best", "worst"), function(i) rep(i, nrow(dat_best)))),
                       rbind(dat_best, dat_worst))

# p1 <- ggplot(dat_best, aes(x = critertion, y = value, label = variable)) +
#   geom_point(size = 5, shape = 21, colour = "blue", fill = adjustcolor("blue", 0.25)) +
#   #geom_text(hjust = -1) +
#   facet_wrap(~group, ncol = 4) +
#   scale_x_discrete("", labels = c("size" = "Size","hydroph" = "Hydroph.",
#                                   "polarity" = "Polarity","alpha" = expression(paste(alpha, "-helix")))) +
#   scale_y_continuous("Value\n") + 
#   my_theme
# 
# p2 <- ggplot(dat_worst, aes(x = critertion, y = value, label = variable)) +
#   geom_point(size = 5, shape = 21, colour = "blue", fill = adjustcolor("blue", 0.25)) +
#   #geom_text(hjust = -1) +
#   facet_wrap(~group, ncol = 4) +
#   scale_x_discrete("", labels = c("size" = "Size","hydroph" = "Hydroph.",
#                                   "polarity" = "Polarity","alpha" = expression(paste(alpha, "-helix")))) +
#   scale_y_continuous("Value") + 
#   my_theme

p1 <- ggplot(dat_bestworst, aes(x = critertion, y = value, col = enc, fill = enc)) +
  geom_point(size = 5, shape = 21, position = position_dodge(width=0.5)) +
  #geom_text(hjust = -1) +
  facet_grid(~group) +
  scale_x_discrete("Criterion\n", labels = c("size" = "Size","hydroph" = "Hydroph.",
                                  "polarity" = "Polarity","alpha" = expression(paste(alpha, "-helix")))) +
  scale_y_continuous("Value") + 
  scale_colour_manual("Encoding: ", values = c("red", "blue")) +
  scale_fill_manual("Encoding: ", values = c(adjustcolor("red", 0.25), adjustcolor("blue", 0.25))) + 
  my_theme

png("./figures/enccomp.png", width = 2257*0.5, height = 1201*0.85)
print(arrangeGrob(textGrob("A", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), p1, 
                  textGrob("B", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), p3,
                  nrow = 2, ncol = 2, widths = c(0.05, 0.95)))
dev.off()


# tables of groups for interesting encodings --------------------------------
library(xtable)

group2df <- function(group_list, caption = NULL, label = NULL) {
  tab <- data.frame(Groups = sapply(group_list, function(i)
    paste0(toupper(sort(i)), collapse = ", ")))
  rws <- seq(1, nrow(tab) - 1, by = 2)
  col <- rep("\\rowcolor[gray]{0.85}", length(rws))
  print(xtable(tab, caption = caption, label = label), include.rownames = FALSE, booktabs = TRUE,
        add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE)
}

cat(group2df(group_best,
             "Best-sensitivity (final) encoding",
             "tab:best"))

cat(group2df(group_worst,
             "Worst-sensitivity encoding",
             "tab:worst"))