library(rmarkdown)
render(input = "./literature/literature.Rmd", output_format = 'html_document')

tmp <- read.bibtex("./literature/kmer_filo.bib")
