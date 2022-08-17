install.packages("BiocManager", , repos = "http://cran.us.r-project.org")
BiocManager::install("RCAS", force = TRUE)
install.packages("kableExtra", , repos = "http://cran.us.r-project.org")
webshot::install_phantomjs()

library(RCAS)
library(kableExtra)

queryRegions <- importBed(filePath = "step11_downstream_analysis/peak_ucsc.bed", sampleN = 10000)
gff <- importGtf(filePath = "step11_downstream_analysis/ucsc_genome.gff3")
overlaps <- as.data.table(queryGff(queryRegions = queryRegions, gffData = gff))
biotype_col <- grep('gene_biotype', colnames(overlaps), value = T)
df <- overlaps[,length(unique(queryIndex)), by = biotype_col]
colnames(df) <- c("feature", "count")
df$percent <- round(df$count / length(queryRegions) * 100, 1)
df <- df[order(count, decreasing = TRUE)]
annotation <- ggplot2::ggplot(df, aes(x = reorder(feature, -percent), y = percent)) + 
  geom_bar(stat = 'identity', aes(fill = feature)) + 
  geom_label(aes(y = percent + 0.5), label = df$count) + 
  labs(x = 'transcript feature', y = paste0('percent overlap (n = ', length(queryRegions), ')')) + 
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 90))

png("step11_downstream_analysis/annotation.png")
print(annotation)
dev.off()

targetedGenes <- unique(overlaps$Dbxref)
targetedGenes <- as.numeric(gsub(".*?([0-9]+).*", "\\1", targetedGenes))  

res <- RCAS::findEnrichedFunctions(targetGenes = targetedGenes, species = 'hsapiens', significant = FALSE)
res <- res[order(res$p_value),]
resGO <- res[grep('GO:BP', res$source),]
knitr::kable(subset(resGO[1:10,], select = c('p_value', 'term_name', 'source'))) %>% kable_styling() %>% save_kable("step11_downstream_analysis/GO.png")
