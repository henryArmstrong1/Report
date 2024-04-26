library("tidyverse")
library("DESeq2")
library("ggplot2")
library("dplyr")
# install.packages(remotes)
library("remotes")
# remotes::install_github("kevinblighe/EnhancedVolcano")
library("EnhancedVolcano")
# install.packages("ggrepel")
library("ggrepel")
# install.packages("BiocManager")
library("BiocManager")
install.packages('BioCircos')
library(BioCircos)
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
install.packages("sjmisc")
library(sjmisc)


  counts <- read.table("raw_data/Hi_PRJNA293882_counts.tsv", row.names = 1, header = TRUE)
featname <- read.table("raw_data/Hi_feature_names.tsv", row.names = 1, header = TRUE)
featlocs <- read.table("raw_data/Hi_feature_locations.bed", header = FALSE, col.names = c("chr","start","end","feat_ID","biotype","strand"))
counts |> ggplot(aes(x = kw20.BHI1.F)) + geom_histogram()
# histogram shows majority of counts on the left not much use

counts |> ggplot(aes(x = log10(kw20.BHI1.F + 1))) + geom_histogram()
# Much clearer! The count data has a fairly normal distribution with the means just above a thousand counts (103 on the x axis)
# not much use for analysis tho
head(featname)
featname[grep("muA|muB|gam", featname$symbol),]
mu_feats <- rownames(featname[grep("muA|muB|gam", featname$symbol),])
mu_feats_grep <- paste(mu_feats, collapse = "|")
featlocs[grep(mu_feats_grep, featlocs$feat_ID),]
hi_prophage_region <- featlocs |> filter(between(start, 1558774, 1597183))
gam_counts <- counts |> filter(row.names(counts) %in% c("gene-HI_1483")) |> 
  t() |> 
  as.data.frame()
colnames(gam_counts) <- "gamcounts"
ggplot(gam_counts, aes(x = rownames(gam_counts), y = gamcounts)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "samples", y = "Hi gam total reads")
ggsave("plots/higam_total_reads.pdf")
# good shows genes of interest
count_tots <- data.frame(colSums(counts))
colnames(count_tots) <- "sum"
max_reads <- max(count_tots$sum)
min_reads <- min(count_tots$sum)
max_reads / min_reads

ggplot(count_tots, aes(x = rownames(count_tots), y = sum)) + 
  geom_bar(stat = "identity") + 
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "samples", y = "total reads")
ggsave("plots/totalreads.pdf")

# The TPM of any feature can be calculated by:
# dividing the feature counts by the feature length (i.e. counts per base)
# dividing counts/base by the total counts for the sample
# multiplying by 1 million to get numbers we can work with (otherwise everything is very small)
# this is done to normalise data for differential expresswion analysis
core_counts <- counts |> dplyr::select(starts_with("kw20.MIV"))
featlens <- data.frame(dplyr::mutate(featlocs, feat_len = (featlocs$end+1)-featlocs$start))
rownames(featlens) <- featlens$feat_ID
withlens <- merge(core_counts, featlens, by=0)
rownames(withlens) <- withlens$Row.names
counts_per_base <- subset(withlens, select = colnames(core_counts)) / withlens$feat_len
tpms <- data.frame(apply(counts_per_base, 2, function(x){(x/sum(x))*1000000}))
colSums(tpms)
tpms <- round(tpms, 2)


gam_tpms <- tpms |> filter(row.names(tpms) %in% c("gene-HI_1483")) |> t() |> as.data.frame()
colnames(gam_tpms) <- "gamtpms"
gamcor <- merge(gam_tpms, gam_counts, by=0)
rownames(gamcor) <- gamcor$Row.names
ggplot(gamcor, aes(x=gamtpms, y=gamcounts)) + 
  geom_point() + 
  geom_text_repel(label=rownames(gamcor)) + 
  labs(x = "gam TPMs", y = "gam read counts") +
  ylim(0, max(gamcor$gamcounts))
ggsave("plots/gamreadcounrs.pdf")
# correlation pretty good

res_pca_log10 <- princomp(log10(tpms+1))
pca_var <- round((data.frame(res_pca$sdev^2/sum(res_pca$sdev^2)))*100,2)
colnames(pca_var) <- "var_perc"
pca_var
pca_comps <- data.frame(res_pca$loadings[, c("Comp.1", "Comp.2")])
pca_x <- paste(c("PC1 ("), pca_var["Comp.1",], c("%)"))
pca_y <- paste(c("PC2 ("), pca_var["Comp.2",], c("%)"))
ggplot(pca_comps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(res_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)
ggsave("plots/Pca1.pdf")
#  good graph shows there have been changes in gene expression in response to stress

ggplot(pca_comps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(res_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)
ggsave("plots/PCA2.pdf")

ggplot(pca_comps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(res_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)
ggsave("plots/log_pca.pdf")
# log reduces the weighting caused by expression/variance
 #---------------------------------------------------------------------------- 
comp_counts <- counts[, grep("kw20.MIV0|kw20.MIV2", colnames(counts))]
sample_info <- data.frame(colnames(comp_counts))
colnames(sample_info) <- c("sample")
rownames(sample_info) <- sample_info$sample
sample_info
sample_info <- sample_info |> tidyr::separate(sample, c("genotype", "condition", "replicate"))
all(rownames(sample_info) %in% colnames(comp_counts))
all(rownames(sample_info) == colnames(comp_counts))
dds <- DESeqDataSetFromMatrix(countData = comp_counts, colData = sample_info, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "MIV0")
dds <- DESeq(dds)
dds_results <- results(dds)
summary(results(dds, alpha=0.05, lfcThreshold = 1))
dds_results

dds_results <- merge(as.data.frame(dds_results), featname, by=0)
rownames(dds_results) <- dds_results$Row.names
dds_results <- dds_results[,-1]

comp_tpms <- tpms[, grep("kw20.MIV0|kw20.MIV2", colnames(tpms))]
comp_tpms <- mutate(comp_tpms, MIV0_avg = rowMeans(select(comp_tpms, contains("MIV0"))), 
                    MIV2_avg = rowMeans(select(comp_tpms, contains("MIV2"))))
comp_tpms <- mutate(comp_tpms, log2FC = log2((MIV2_avg+1)/(MIV0_avg+1)))
comp_red_tpms <- comp_tpms |> select(-starts_with("kw20.MIV"))

dds_tpm <- merge(dds_results, comp_red_tpms, by=0)
rownames(dds_tpm) <- dds_tpm$Row.names
dds_tpm <- dds_tpm[,-1]
cor(dds_tpm$log2FoldChange, dds_tpm$log2FC, method = c("pearson"))
head(dds_tpm[order(dds_tpm$padj),], 20)
withsymbols <- dds_tpm[- grep("[.]", dds_tpm$symbol),]
siggenes <- head(withsymbols |> arrange(desc(abs(log2FC))), 60)$symbol

EnhancedVolcano(dds_tpm, x = "log2FC", y = "padj", FCcutoff = 1, pCutoff = 0.05,
                lab = dds_tpm$symbol, selectLab = siggenes, labSize = 4.5, max.overlaps = 1000, drawConnectors = TRUE,
                legendPosition = 0, gridlines.major = FALSE, gridlines.minor = FALSE)
  ggsave("plots/volcanoplots.pdf")

 # ---------------------------------------------------------------------------------
 #  --------------------------------------------------------------------------------
 # DATA workshop 4

counts <- read.table("raw_data/Hi_PRJNA293882_counts.tsv", row.names = 1, header = TRUE)
featname <- read.table("raw_data/Hi_feature_names.tsv", row.names = 1, header = TRUE)
featlocs <- read.table("raw_data/Hi_feature_locations.bed", header = FALSE, col.names = c("chr","start","end","feat_ID","biotype","strand"))
rownames(featlocs) <- featlocs$feat_ID
gc <- read.table("raw_data/Hi_GC_1kb.bed", header=FALSE, col.names = c("chr", "start", "end", "GC"))
tpms <- read.table("proc_data/full_dataset_TPMs.tsv", row.names = 1, header = TRUE)
dds_tpm <- read.table("proc_data/MIV0vsMIV2_DEA_full_results.tsv", row.names = 1, header = TRUE)




dds_tpm$DEA <- "NO"
dds_tpm$DEA[dds_tpm$log2FC > 1 & dds_tpm$padj < 0.05] <- "UP"
dds_tpm$DEA[dds_tpm$log2FC < -1 & dds_tpm$padj < 0.05] <- "DOWN"
ggplot(dds_tpm, aes(x=log2FC, y=-log10(padj))) + 
  geom_point(aes(colour = DEA), show.legend = FALSE) + 
  scale_colour_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(-1,1), linetype = "dotted") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_text_repel(data=subset(withsymbols,abs(log2FC) > 3), aes(x=log2FC, y=-log10(padj), label=symbol), max.overlaps = 1000)
  ggsave("plots/volcanoplots2.pdf")


BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
organism  <- "org.EcK12.eg.db"
BiocManager::install(organism)
library(organism, character.only = TRUE)
symbols_only <- dds_tpm[- grep("[.]", dds_tpm$symbol.x),]
symbols_only <- symbols_only[order(symbols_only[,"symbol.x"],-symbols_only[,"log2FC"]),]
symbols_only <- symbols_only[!duplicated(symbols_only$symbol.x),]
symbols_only <- symbols_only[- grep("ribosomal_RNA", symbols_only$description.x),]
symbols_only <- symbols_only[- grep("tRNA", symbols_only$symbol.x),]
gene_list <- symbols_only$log2FC
names(gene_list) <- symbols_only$symbol.x
gene_list <- sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList = gene_list, ont = "ALL", keyType = "SYMBOL",
             minGSSize = 3, maxGSSize = 800, verbose = FALSE,
             pvalueCutoff = 0.05, OrgDb = organism, pAdjustMethod = "none")
gsearesults <- gse[,c("Description","setSize","qvalue","NES","core_enrichment")]
gsearesults_top <- head(gsearesults[order(gsearesults$qvalue),], 40)
gsearesults_top <- gsearesults_top[order(gsearesults_top$NES, decreasing = TRUE),]
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("plots/MIV0MIV2_GSEA_dotplot.pdf")
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

gsearesults[grep("gam", gsearesults$core_enrichment),]
gsearesults[grep("Competence|competence", gsearesults$Description),]

condcomp_tpms <- tpms[, grep("kw20.MIV0|kw20.MIV2|kw20.BHI3", colnames(tpms))]
condcomp_pca <- princomp(log10(condcomp_tpms+1))
condcomp_pcavar <- round((data.frame(condcomp_pca$sdev^2/sum(condcomp_pca$sdev^2)))*100,2)
colnames(condcomp_pcavar) <- "var_perc"
condcomp_pcavar
condcomp_pcacomps <- data.frame(condcomp_pca$loadings[, c("Comp.1", "Comp.2")])
pca_x <- paste(c("PC1 ("), condcomp_pcavar["Comp.1",], c("%)"))
pca_y <- paste(c("PC2 ("), condcomp_pcavar["Comp.2",], c("%)"))
ggplot(condcomp_pcacomps, aes(x=Comp.1, y=Comp.2)) + 
  geom_point() + 
  geom_text_repel(label=rownames(data.frame(condcomp_pca$loadings[, 1:2]))) +
  labs(x = pca_x, y = pca_y)
  ggsave("plots/PCA3.pdf")



MIV0BHI3_counts <- counts[, grep("kw20.MIV0|kw20.BHI3", colnames(counts))]
sample_info <- data.frame(colnames(MIV0BHI3_counts))
colnames(sample_info) <- c("sample")
rownames(sample_info) <- sample_info$sample
sample_info <- sample_info |> separate(sample, c("genotype", "condition", "replicate"))
MIV0BHI3_dds <- DESeqDataSetFromMatrix(countData = MIV0BHI3_counts, colData = sample_info, design = ~ condition)
MIV0BHI3_dds$condition <- relevel(MIV0BHI3_dds$condition, ref = "MIV0")
MIV0BHI3_dds <- DESeq(MIV0BHI3_dds)
MIV0BHI3_dds_results <- results(MIV0BHI3_dds)
MIV0BHI3_dds_results <- merge(as.data.frame(MIV0BHI3_dds_results), featname, by=0)
rownames(MIV0BHI3_dds_results) <- MIV0BHI3_dds_results$Row.names
MIV0BHI3_dds_results <- MIV0BHI3_dds_results[,-1]
MIV0BHI3_tpms <- condcomp_tpms[, -grep("MIV2", colnames(condcomp_tpms))]
MIV0BHI3_tpms <- mutate(MIV0BHI3_tpms, MIV0_avg = rowMeans(dplyr::select(MIV0BHI3_tpms, contains("MIV0"))), 
                        BHI3_avg = rowMeans(dplyr::select(MIV0BHI3_tpms, contains("BHI3"))))
MIV0BHI3_tpms <- mutate(MIV0BHI3_tpms, log2FC = log2((BHI3_avg+1)/(MIV0_avg+1)))
MIV0BHI3_dds_results <- merge(MIV0BHI3_dds_results, MIV0BHI3_tpms, by=0)
rownames(MIV0BHI3_dds_results) <- MIV0BHI3_dds_results$Row.names
MIV0BHI3_dds_results <- MIV0BHI3_dds_results[,-1]
MIV0BHI3_withsymbols <- MIV0BHI3_dds_results[- grep("[.]|HI_", MIV0BHI3_dds_results$symbol),]
MIV0BHI3_siggenes <- head(MIV0BHI3_withsymbols |> arrange(desc(abs(log2FC))), 50)$symbol
EnhancedVolcano(MIV0BHI3_dds_results, x = "log2FoldChange", y = "padj", FCcutoff = 1, pCutoff = 0.05,
                lab = MIV0BHI3_dds_results$symbol, selectLab = MIV0BHI3_siggenes, labSize = 4.5, max.overlaps = 1000, drawConnectors = TRUE,
                legendPosition = 0, gridlines.major = FALSE, gridlines.minor = FALSE)
  ggsave("plots/volcanoeplot3.pdf")

MIV0BHI3_withoutsymbols <- MIV0BHI3_dds_results[grep("[.]", MIV0BHI3_dds_results$symbol),]
MIV0BHI3_withoutsymbols_sig <- rownames(head(MIV0BHI3_withoutsymbols |> arrange(desc(abs(log2FC))), 10))


MIV0BHI3_DEA <- MIV0BHI3_dds_results[, c("log2FoldChange", "padj", "symbol")]
colnames(MIV0BHI3_DEA) <- c("MIV0BHI3_log2FC", "MIV0BHI3_padj", "symbol")
MIV0MIV2_DEA <- dds_tpm[, c("log2FoldChange", "padj")]
colnames(MIV0MIV2_DEA) <- c("MIV0MIV2_log2FC", "MIV0MIV2_padj")
dea_comp <- merge(MIV0BHI3_DEA, MIV0MIV2_DEA, by=0)
rownames(dea_comp) <- dea_comp$Row.names
dea_comp <- dea_comp[,-1]
dea_comp_symbols <- dea_comp[- grep("[.]", dea_comp$symbol),]
cor(dea_comp$MIV0BHI3_log2FC, dea_comp$MIV0MIV2_log2FC, method = c("pearson"))
ggplot(dea_comp, aes(x=MIV0BHI3_log2FC, y=MIV0MIV2_log2FC)) + 
  geom_point() + 
  geom_text_repel(data=subset(dea_comp_symbols, abs(MIV0MIV2_log2FC) > 3 & abs(MIV0BHI3_log2FC) < 1), aes(x=MIV0BHI3_log2FC, y=MIV0MIV2_log2FC, label=symbol), max.overlaps = 1000, colour="red") +
  geom_abline() + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave("plots/DEA.pdf")

tracklist <- BioCircosTextTrack("titletrack", "Hi kw20", opacity = 0.5, x = -0.2, y = 0)
tracklist <- tracklist + BioCircosArcTrack("prophage region", "L42023.1", 1558774, 1597183, 
                                           opacities = c(1), minRadius = 1.25, maxRadius = 1.4,
                                           labels = c("mu prophage"))
tracklist <- tracklist + BioCircosLineTrack("GC", "L42023.1", gc$start, gc$GC, 
                                            minRadius = 0.4, maxRadius = 0.68, color = "black",
                                                                                     labels = c("GC content"))

rownames(featlocs) <- featlocs$feat_ID
dea <- merge(dds_tpm, featlocs, by=0)
rownames(dea) <- dea$Row.names
dea <- dea[,c("symbol", "log2FoldChange", "chr", "start", "end", "strand")]
tracklist <- tracklist + BioCircosHeatmapTrack("DEA", "L42023.1", 
                                               dea$start, dea$end, dea$log2FC,
                                               minRadius = 0.8, maxRadius = 0.95,
                                               color = c("#FF0000", "#0000FF"),
                                               labels = c("MIV0vsMIV2 DEA"))
BioCircos(tracklist, genome = list("L42023.1" = 1830138), genomeLabelTextSize = 0,
          genomeTicksScale = 1e5, genomeTicksTextSize = 12)
ggsave("plots/circos.pdf")




t_tpms <- tpms |> rotate_df()
gam_tpm_cor <- cor(as.matrix(t_tpms[,c("gene-HI_1483")]), as.matrix(t_tpms))
gam_cors <- as.data.frame(gam_tpm_cor[,order(-gam_tpm_cor[1,])])
colnames(gam_cors) <- c("gam_cor")
gam_cors <- merge(gam_cors, featname, by=0)
rownames(gam_cors) <- gam_cors$Row.names
gam_cors <- gam_cors[,-1]
gam_cors <- gam_cors[order(-gam_cors$gam_cor),]
head(gam_cors, 25)
gam_cors <- merge(gam_cors, featlocs, by=0)
rownames(gam_cors) <- gam_cors$Row.names
gam_cors <- gam_cors[,-1]
gam_cors <- mutate(gam_cors, gam_circ_dist = ((start+1830138)-1565806), gam_noncirc_dist = abs(1565297-start))
gam_cors <- transform(gam_cors, min_gam_dist = pmin(gam_circ_dist, gam_noncirc_dist))
gam_cors_top <- (gam_cors[order(-gam_cors$gam_cor),])[,c("gam_cor", "symbol", "min_gam_dist")]
ggplot(gam_cors_top, aes(x=log10(min_gam_dist+1), y=gam_cor)) + 
  geom_point() + 
  geom_text_repel(data = gam_cors_top |> mutate(label = ifelse(gam_cor > 0.87, rownames(gam_cors_top), "")), aes(label = label), max.overlaps = 100) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave("plots/gamcors.pdf")


