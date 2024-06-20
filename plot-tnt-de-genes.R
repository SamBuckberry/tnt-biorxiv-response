## Sample metadata
sample_dat <- data.table::fread("polyA_RNAseq_sample_table.tsv")

## RNA-seq counts data
count_dat <- read.table("mel1_rmdup.featureCounts.gz", header = TRUE)

colnames(count_dat) <- colnames(count_dat) %>%
    str_remove("X.data.tki_bodl.sammyb.") %>%
    str_remove("_Aligned.markdup.bam")

## Subset just to counts
dat <- count_dat[ ,colnames(count_dat) %in% sample_dat$Library]
rownames(dat) <- count_dat$Geneid

# Check library ID's match and order accordingly
stopifnot(all(colnames(dat) %in% sample_dat$Library))
sample_dat <- sample_dat[match(colnames(dat), sample_dat$Library), ]
stopifnot(all(sample_dat$Library == colnames(dat)))

# Calculate TPM
countToTpm <- function(counts, effLen){
    rate <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
}

tpm <- countToTpm(counts = dat, effLen = count_dat$Length)

gene_gr <- import("gencode.v27.primary_assembly.annotation.gtf.gz", format = "gtf")

gene_gr$ensembl_gene_id <- sub("\\..*$", "", gene_gr$gene_id)

## Clean up labels
sample_dat$Group[sample_dat$Group == "N2P"] <- "NtP"

reprog_pal <- c(Primed="#009E73", Naive="#0072B2",
                TNT="#EEBC4C", NtP="#55B3EA", ESC="darkgrey")

plot_gene <- function(gene_id){

    ensg_select <- gene_gr$gene_id[gene_gr$ensembl_gene_id %in% gene_id] %>% unique()

    tpm_select <- tpm[rownames(tpm) %in% ensg_select, ]
    tpm_select$ENSG_id <- rownames(tpm_select)

    tpm_select <- reshape2::melt(tpm_select)
    ind <- match(tpm_select$variable, sample_dat$Library)
    tpm_select$group <- sample_dat$Group[ind]

    ind2 <- match(tpm_select$ENSG_id, gene_gr$gene_id)
    tpm_select$Gene <- gene_gr$ensembl_gene_id[ind2]

    tpm_select$group <- factor(tpm_select$group, levels=c("ESC", "Fibroblast",
                                                          "Naive", "TNT", "NtP",
                                                          "Primed"))

    gg_gene_select <- ggplot(tpm_select, aes(x = group, y = value, colour=group)) +
        #stat_summary(fun = mean, geom = "bar", width = 0.2, alpha=0.5) +
        geom_point(size=2) +
        scale_colour_manual(values = reprog_pal) +
        #geom_col(width = 0.1) +
        ggtitle(str_c("MEL1: ", gene_id)) +
        ylab("TPM") +
        sams_pub_theme()

    # gg_gene_select_log <- ggplot(tpm_select, aes(x = group, y = log2(value+1))) +
    #     geom_point() +
    #     #facet_grid(Gene~., scales = "free") +
    #     ggtitle(str_c("MEL1: ", gene_id)) +
    #     ylab("Log2 TPM+1") +
    #     sams_pub_theme()

    #cowplot::plot_grid(gg_gene_select, gg_gene_select_log)
    gg_gene_select

}

ensg_select <- gene_gr$ensembl_gene_id[gene_gr$gene_name %in%
                                   c("CAT", "OXCT1", "TSPYL5")] %>% unique()

pdf("tnt-not-corrected-3-gene-plot.pdf", height = 2, width = 5)
cowplot::plot_grid(plotlist = lapply(ensg_select, plot_gene),
                   nrow = 1)
dev.off()
