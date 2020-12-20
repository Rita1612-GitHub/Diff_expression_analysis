# Diff_expression_analysis

Authors:
- Komarova Margarita
- Gainova Kristina
- Kvach Anna

This repositorium consist description of several pipelines of differential gene expression analysis.

**Aim**

In this project we tests several pipelines of differential gene expression analysis. 
The aim of this project was to compare different methods of differential gene expression analysis and identification of the most suitable for certain types of data. 
For the implementation of this project, several articles were taken and reanalyzed on open sequencing data.

# *Mus Musculus* (Komarova Margarita)
## Methods:
- STAR (v2.7.3a) alignment
- featureCounts (v2.0.1) (getting counts table)
- DESeq2 (v1.26.0) (gene expression analysis in Rstudio)
- also we use fgsea (v1.12.0) and ClusterProfiler (v3.14.3) packages to identify signal pathways. 

For the first pipline was taken article doi:10.3390/genes11091057. 
The object of the study in the article is the cell line of murine myoblasts C2C12:wild type and with a mutation in the LMNA gene (G232E). The aim of the study in this article was to assess the effect of the mutation on myogenic differentiation C2C12 (HS influence).

Samples (three replicas for each):
- wt_contr
- 232_contr
- wt_HS
- 232_HS

## The first step of all piplines was quality assessment of reads using FastQC program (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

Reads of this article were posted in open access and are available via the link https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP261257&o=acc_s%3Aa. The reads have of high quality and do not need trimming procedure.

## At the next step, the genome was indexed in the STAR program. Example of command for this step:
```
/home/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir /mnt4/transciptome_compare/komarova/STAR_genome_index --genomeFastaFiles /mnt4/transciptome_compare/komarova/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa --sjdbGTFfile /mnt4/transciptome_compare/komarova/Mus_musculus.GRCm38.101.gtf --sjdbOverhang 50 
```
Depending on the organism, this process can take 4-8 hours.

## Then you need to align the reads to the genome
```
/home/tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --genomeDir /mnt4/transciptome_compare/komarova/STAR_genome_index/ --runThreadN 10 --readFilesIn /mnt4/transciptome_compare/komarova/SRR11776837.fastq --outFileNamePrefix /mnt4/transciptome_compare/komarova/STAR_alignment/11776837 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
```
## Counts table
After STAR, counts are obtained in files with .tab resolution and aligned reads in .out.bam. But it is better to run featureCount to calculate counts.
```
/home/tools/subread-2.0.1-source/bin/featureCounts -T 4 -s 0 -a /mnt4/transciptome_compare/komarova/Mus_musculus.GRCm38.101.gtf -o /mnt4/transciptome_compare/komarova/STAR_alignment/featureCount_results/Counts.txt /mnt4/transciptome_compare/komarova/STAR_alignment/*.out.bam
```
Where is:
- -T number of cores
- -s these data are "reverse"ly stranded
- -a path to annotation
- -o output file. The output of this tool is 2 files, a count matrix and a summary file that tabulates how many the reads were “assigned” or counted and the reason they remained “unassigned”.

## "Cleaning" of output file
Unnecessary columns need to be removed.
```
cut -f1,7-18 /mnt4/transciptome_compare/komarova/STAR_alignment/featureCount_results/Counts.txt > /mnt4/transciptome_compare/komarova/STAR_alignment/featureCount_results/Counts.Rmatrix.txt
```
After that we need to add column names and remove unnecessary information.

## DESeq2 analysis
**1. Data pre-processing.**
Read the file with featureCounts results.
```
table_count <- read.table("Counts.Rmatrix2.txt", header=TRUE, row.names =1)
```
Let's remove the X before 232 in the title.
```
tags <- sapply(colnames(table_count), function (x) sub("X", "", basename(x)))
colnames(table_count) <- tags
```
**2. Condition assignment.**
Condition of the experiment:
```
cond <- data.frame(read.csv("Conditions_IB.tsv", header=T, sep = '\t', row.names = 1))
```
**3. Creation DESEq2 dataset.**
```
samples = names(table_count)
dds = DESeqDataSetFromMatrix(countData=table_count, colData=cond, design = ~ Condition)
#dds <- dds[rowSums(counts(dds)) > 20, ] # Filtering (not nessisary on this step).
dds <- DESeq(dds) #Normalisation
res <- results(dds)
```
**4. PCA plot.**
```
png("PCA_plot.png", 1000, 1000, pointsize=20)
vst_dds <- vst(dds)
counts.norm <- assay(vst_dds)
PCA_plot <- plotPCA(vst_dds,intgroup=c("Condition"))
PCA_plot
dev.off()
PCA_plot
```
**5. Volcano plot.**
```
gdata <- data.frame(
  x=res$log2FoldChange,
  y=-log10(res$padj)
)
png("Volcanoplot_before.png", 500, 500, pointsize=30)
ggplot(data=gdata, aes(x=x, y=y, color = (res$padj<0.01 & abs(res$log2FoldChange) >= 1) )) +
  geom_point(size=1) + theme_bw()  +
  xlab("Log fold change") +
  ylab("Adjusted p.value")
dev.off()
```
**6. Heatmap plot.**
```
counts.norm <- counts(dds, normalized=TRUE)
png("heatmap_large.png", width=6, height=20, units="in", res=300)
to_visualise <- counts.norm[rownames(res), order(cond[, 2])]
to_visualise <- t(apply(to_visualise, 1, function(r) {
  (r - min(r)) / (max(r) - min(r))
}))
to_visualise <- na.omit(to_visualise)
pheatmap(to_visualise, 
         show_rownames = F, cluster_rows = F,
         cluster_cols=F,
         annotation_col = cond)
dev.off()
```
## GO analysis
GO analysis allows to group differentially expressed genes by their molecular function, biological process or cellular component. We can identify which signal pathwats are upregulated or downregulated.

**1. Fast Gene Set Enrichment Analysis (GSEA) Using fgsea package.**
Download annotation
```
load("C:/R/IB/Projects/Data/mouse_c5_v5p2.rdata")
pathways <- Mm.c5 # http://bioinf.wehi.edu.au/software/MSigDB/
```
Let us carry out a comparative analysis of differentially expressed genes between the wt_contr and 232_contr groups. Let's see which genes were influenced by the 232 mutation in LMNA.
```
de <- results(dds, contrast = c("Condition", "wt_contr", "232_contr"))
de_Genes <- row.names(de) 
de$row_Genes <- de_Genes 
de <- as.data.table(de)
de <- na.omit(de) # NA omit
de$row_Genes <- as.character(de$row_Genes) 
de[, row_Genes2 := mapIds(org.Mm.eg.db, keys=row_Genes, 
                                        keytype="ENSEMBL", column="SYMBOL")] # Create new column with gene SYMBOL
de <- unique(de, by = "row_Genes") 
de <-de[order(de$padj),] #sorting padj
write.csv(de, "table_de.csv") # write in file
gene_list_de <- de$log2FoldChange
genes <- bitr(de$row_Genes,
              fromType = "ENSEMBL",
              toType = "ENTREZID", 
              OrgDb = org.Mm.eg.db) 
names(gene_list_de) <- genes$ENTREZID
#For fgsea need to create named vector with ENTREZID of genes and their log2FoldChange values.
fr_de <- fgsea(pathways, gene_list_de, nperm = 100000, nproc=4, minSize=15, maxSize=500)
#head(fr_de[order(pval), ])
fr_de[order(padj)]
```
Some graphs:
```
plotEnrichment(pathways[["GO_MUSCLE_TISSUE_DEVELOPMENT"]], gene_list_de) + 
  ggtitle("Myogenesis")
```
Upregulated and downregulated signal pathways:
```
topPathwaysUp <- fr_de[ES > 0][head(order(pval), n=5), pathway]
topPathwaysDown <- fr_de[ES < 0][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
```

```
plotGseaTable(pathways[topPathways], gene_list_de, fr_de, 
              gseaParam=0.5)
```

You can do the same for wt_HS vs 232_HS.

**2. Gene Set Enrichment Analysis with ClusterProfiler.**
Filtering:
```
sorted <- res[with(res, order(padj, -log2FoldChange)), ] # 16289 genes.
sorted_df <- data.frame("id"=rownames(sorted),sorted)
sorted_df <- na.omit(sorted_df) #удалили NA
p_adj_below_0.01 <- sorted_df %>% filter(padj<0.01 & abs(log2FoldChange)>1) #filtering to padj и log2FoldChange
```
Add some new column with ENTREZID and SYMBOL of genes.
```
p_adj_Genes <- row.names(p_adj_below_0.01)
p_adj_below_0.01$row_Genes <- p_adj_Genes
p_adj_below_0.01 <- as.data.table(p_adj_below_0.01)
p_adj_below_0.01 <- na.omit(p_adj_below_0.01)
p_adj_below_0.01$row_Genes <- as.character(p_adj_below_0.01$row_Genes)
p_adj_below_0.01[, row_Genes2 := mapIds(org.Mm.eg.db, keys=row_Genes, 
                                        keytype="ENSEMBL", column="SYMBOL")]
p_adj_below_0.01[, row_Genes3 := mapIds(org.Mm.eg.db, keys=row_Genes, 
                                        keytype="ENSEMBL", column="ENTREZID")]
p_adj_below_0.01 <-p_adj_below_0.01[order(p_adj_below_0.01$padj),]
p_adj_below_0.01$row_Genes3 <- as.character(p_adj_below_0.01$row_Genes3)
```
Drow some pictures:
```
de[, row_Genes3 := mapIds(org.Mm.eg.db, keys=row_Genes, 
                                        keytype="ENSEMBL", column="ENTREZID")]
ego_de <- enrichGO(gene          = de$row_Genes3,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
dotplot(ego_de, showCategory=30)
#BiocManager::install("Rgraphviz")
plot_GO <- plotGOgraph(ego_de)
```
2. Пайплайн Кристины Гайновой

3. Пайплайн Ани Квач

The first step of all piplines was quality assessment of reads using FastQC programm.

