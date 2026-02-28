############################################################
# Transcriptomic Analysis of Salt Stress Response
# Dataset: GSE52305
# Organism: Populus euphratica
############################################################

# ==========================
# 1. Load Packages
# ==========================

library(GEOquery)
library(limma)
library(stringr)
library(pheatmap)
library(ggplot2)

# ==========================
# 2. Download & Prepare Data
# ==========================

gset <- getGEO("GSE52305", GSEMatrix = TRUE)[[1]]

ex <- exprs(gset)
pheno <- pData(gset)

# ==========================
# 3. Define Groups
# ==========================

group <- pheno$`treatment:ch1`
group <- make.names(group)

group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# ==========================
# 4. Linear Model & Contrast
# ==========================

fit <- lmFit(ex, design)

contrast.matrix <- makeContrasts(
  Salt600_vs_Control = Seedlings.irrigated.with.600.mM.NaCl.solution -
    Control.irrigated.with.natural.water,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2,
                    adjust.method = "BH",
                    number = Inf)

# ==========================
# 5. DEG Summary
# ==========================

deg_up <- results[results$adj.P.Val < 0.05 & results$logFC > 0, ]
deg_down <- results[results$adj.P.Val < 0.05 & results$logFC < 0, ]

cat("Total DEG:", nrow(deg_up) + nrow(deg_down), "\n")
cat("Upregulated:", nrow(deg_up), "\n")
cat("Downregulated:", nrow(deg_down), "\n")

# ==========================
# 6. Volcano Plot
# ==========================

results$threshold <- "Not Significant"
results$threshold[results$adj.P.Val < 0.05 & results$logFC > 0] <- "Up"
results$threshold[results$adj.P.Val < 0.05 & results$logFC < 0] <- "Down"

png("Volcano_Salt600_vs_Control.png", width=2000, height=1500, res=300)

ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=threshold)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  ggtitle("Volcano Plot: Salt600 vs Control")

dev.off()

# ==========================
# 7. Heatmap Top 50 Genes
# ==========================

top50 <- head(results, 50)
genes_top50 <- rownames(top50)
heatmap_data <- ex[genes_top50, ]

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(heatmap_data)

png("Heatmap_Top50_Salt600_vs_Control.png", width=2000, height=2000, res=300)

pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes")

dev.off()

# ==========================
# 8. GO Biological Process (Descriptive)
# ==========================

annotation_data <- fData(gset)
sig_genes <- rownames(results[results$adj.P.Val < 0.05, ])

go_bp <- annotation_data[sig_genes, "Gene Ontology Biological Process"]
go_bp <- go_bp[!is.na(go_bp)]

go_list <- str_split(go_bp, "///")
go_list <- unlist(go_list)

go_clean <- str_extract(go_list, "(?<=// ).*?(?= //)")
go_clean <- go_clean[!is.na(go_clean)]

table_go <- sort(table(go_clean), decreasing = TRUE)

top10_go <- head(table_go, 10)

png("GO_Biological_Process.png", width=2000, height=1500, res=300)

par(mar=c(10,4,4,2))
barplot(top10_go,
        las=2,
        col="steelblue",
        main="Top Enriched GO Biological Processes")

dev.off()

############################################################
# END OF ANALYSIS
############################################################