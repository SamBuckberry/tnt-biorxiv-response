# install.packages(c("magrittr", "ggplot2", "stringr", "reshape2",
#                    "gridExtra", "cowplot", "plotly", "dplyr",
#                    "data.table", "R.utils", "kableExtra"))
#
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install(c("limma", "edgeR", "GenomicFeatures", "Rsamtools", "Gviz",
#                        "rtracklayer", "bsseq",
#                        "BSgenome.Hsapiens.UCSC.hg38",
#                        "BSgenome.Hsapiens.UCSC.hg19", "GenomicRanges",
#                        "TxDb.Hsapiens.UCSC.hg38.knownGene"))


library(limma)
library(edgeR)
library(magrittr)
library(ggplot2)
library(GenomicFeatures)
library(stringr)
library(reshape2)
library(gridExtra)
library(cowplot)
library(Rsamtools)
library(Gviz)
library(rtracklayer)
library(plotly)
library(dplyr)
library(data.table)

sams_pub_theme <- function(legend_pos = "none", x.text.angle = 45, hjust = 1,
                           y.text.angle = 0, line_point = 0.5) {
    line_mm <- line_point / 2.835
    custom_theme <- theme_bw() +
        theme(plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(size = line_mm),
              panel.border = element_rect(color = "black", size = line_mm, fill = NA), # Add border around each facet
              axis.text.x = element_text(angle = x.text.angle, hjust = hjust, size = 6, colour = 'black'),
              axis.text.y = element_text(size = 6, colour = 'black', angle = y.text.angle),
              strip.text.y = element_text(size = 6),
              text = element_text(size = 6),
              strip.background = element_blank(),
              legend.position = legend_pos,
              axis.line.x = element_line(color = 'black', size = line_mm),
              axis.line.y = element_line(color = 'black', size = line_mm),
              axis.ticks = element_line(color = 'black', size = line_mm))
    return(custom_theme)
}

#
# sams_pub_theme <- function(legend_pos = "NA", x.text.angle = 45, hjust=1,
#                            y.text.angle=0, line_point=0.5){
#     line_mm <- line_point / 2.835
#     custom_theme <- theme_bw() +
#         theme(plot.background = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.grid.major = element_line(size = line_mm),
#               panel.border = element_blank(),
#               axis.text.x = element_text(angle = x.text.angle, hjust = hjust, size = 6,
#                                          colour = 'black'),
#               axis.text.y = element_text(size=6, colour='black', angle = y.text.angle),
#               strip.text.y = element_text(size = 6),
#               text = element_text(size=6),
#               strip.background = element_blank(),
#               legend.position = legend_pos,
#               axis.line.x = element_line(color = 'black', size = line_mm),
#               axis.line.y = element_line(color = 'black', size = line_mm),
#               axis.ticks = element_line(color = 'black', size = line_mm))
#     return(custom_theme)
# }
