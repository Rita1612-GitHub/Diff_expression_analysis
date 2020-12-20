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
# *Arabidopsis thaliana* (Gainova Kristina)
## Methods:
- TopHat2 (v2.1.1) alignment
- Cuffdiff (v2.2.1) to to find significant changes in transcript expression, splicing, and promoter use.
- fgsea (v1.12.0) and ClusterProfiler (v3.14.3) packages to identify signal pathways. 

To capture the dynamic transcriptional response of Arabidopsis plants to HL stress, scientists performed a time course RNA sequencing (RNA-seq) study (doi:10.1016/j.celrep.2019.11.051). Arabidopsis seedlings were treated with high light (HL: 1,200 mmol m-2 s-1, leaf temperature: 22C) for 0.5, 6, 12, 24, 48, and 72 h with corresponding growth light (GL: 60 mmol m-2 s-1, leaf temperature: 22C) treatments as control. 
In this study samples treated with HL for 24 h with corresponding growth light (highlighted in yellow) were chosen for following analysis with Tophat and Cuffdiff.

Samples (two replicas for each):
SRR6767639,SRR6767640	GL24h
SRR6767652,SRR6767653 HL24h

Genome (and bowtie2 index files) and annotation were downloaded from: (http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Arabidopsis_thaliana/Ensembl/TAIR10/Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz)

## FastQC
At the first, it was checked quality of raw reads using FastQC program (v0.11.5). The diagrams below demonstrate high quality of data. Thus read trimming is not required for mapping and quantification of RNA-seq reads.

## TopHat
Usage: **tophat [options]* <genome_index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]**
Arguments:

**genome_index_base** 	The basename of the genome index to be searched. The basename is the name of any of the index files up to but not including the first period. Bowtie first looks in the current directory for the index files, then looks in the indexes subdirectory under the directory where the currently-running bowtie executable is located, then looks in the directory specified in the BOWTIE_INDEXES  (or BOWTIE2_INDEXES) environment variable. **Please note that it is highly recommended that a FASTA file with the sequence(s) the genome being indexed be present in the same directory with the Bowtie index files** and having the name <genome_index_base>.fa. If not present, TopHat will automatically rebuild this FASTA file from the Bowtie index files.

**reads1_1[,...,readsN_1]**	A comma-separated list of files containing reads in FASTQ or FASTA format. When running TopHat with paired-end reads, this should be the *_1 ("left") set of files.

**[reads1_2,...readsN_2]**	A comma-separated list of files containing reads in FASTQ or FASTA format. Only used when running TopHat with paired end reads, and contains the *_2 ("right") set of files. The *_2 files MUST appear in the same order as the *_1 files.


Options:

**-i/--min-intron-length <int>**	The minimum intron length. TopHat will ignore donor/acceptor pairs closer than this many bases apart. The default is 70.
**-I/--max-intron-length <int>**	The maximum intron length. When searching for junctions ab initio, TopHat will ignore donor/acceptor pairs farther than this many bases apart, except when such a pair is supported by a split segment alignment of a long read. The default is 500000.

### Running TopHat

Working directory must contains: reference genome and bowtie index files, annotation (GTF) and reads in format fastq

I used the following commands to map the RNA-seq reads to the genome:

~$ tophat -p 2 -i 20 -I 5000 -G genes.gtf -o SRR6767639_tophat genome SRR6767639.fastq

~$ tophat -p 2 -i 20 -I 5000 -G genes.gtf -o SRR6767640_tophat genome SRR6767640.fastq

~$ tophat -p 2 -i 20 -I 5000 -G genes.gtf -o SRR6767652_tophat genome SRR6767652.fastq

~$ tophat -p 2 -i 20 -I 5000 -G genes.gtf -o SRR6767653_tophat genome SRR6767653.fastq


Here, SRR***_tophat represents the output directory for each run. If no output directory is specified by the –o option, TopHat automatically creates a directory called tophat_out in the working director, and stores the output files in it. The values of the intron lengths are adjusted as 20bp and 5000 bp respectively for Arabidopsis thaliana instead of the default values which are for mammals.

**Error occurs during the execution of the program**
samtools view: samtools view: writing to standard output failedwriting to standard output failed: Broken pipe
: Broken pipe
samtools view: error closing standard output: -1

The successful run creates a directory with the name specified by the user, containing the following files: **accepted_hits.bam**, **align_summary.txt**, **deletions.bed**, **insertions.bed**, **junctions.bed**, **prep_reads.info**, **unmapped.bam** and a directory with the name **logs**.

## Cuffdiff

Since it exists reference genome and annotation for Arabidopsis thaliana I skipped steps with Cufflinks and Cuffmerge (_de novo assembling_ step) and used Cuffdiff program.

Cuffdiff program is used to find significant changes in transcript expression, splicing, and promoter use. Manual for Cuffdiff you can find here <http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/>.


Cuffdiff options:

**-o/–output-dir <string>* Sets the name of the directory in which Cuffdiff will write all of its output. The default is “./”. From the command line, run cuffdiff as follows:

**-L/–labels <label1,label2,…,labelN>** Specify a label for each sample, which will be included in various output files produced by Cuffdiff.

**-p/–num-threads <int>** Use this many threads to align reads. The default is 1.

~$ cuffdiff --geometric-norm -p 16 -o cuffdiff_result_only genes.gtf -L GL_12h_control,HL_12h_treatment ./SRR6767639_tophat/accepted_hits.bam,./SRR6767640_tophat/accepted_hits.bam ./SRR6767652_tophat/accepted_hits.bam,./SRR6767653_tophat/accepted_hits.bam

Here parameter **--geometric-norm** was used for data normalization.

The successful run creates a directory with the name specified by the user, containing the following files:
bias_params.info, gene_exp.diff, isoforms.fpkm_tracking, tss_group_exp.diff, cds.count_tracking, genes.count_tracking,    isoforms.read_group_tracking, tss_groups.count_tracking, cds.diff, genes.fpkm_tracking, promoters.diff, tss_groups.fpkm_tracking,
cds_exp.diff, genes.read_group_tracking, read_groups.info, tss_groups.read_group_tracking, cds.fpkm_tracking, isoform_exp.diff, run.info, var_model.info, cds.read_group_tracking, isoforms.count_tracking, splicing.diff

Needed file called gene_exp.diff. Initially I viewed content of the file and remove rows, where encountered value "inf" using command line:

~$ cat gene_exp.diff | grep -v "inf"  > new_file.csv

## Data analysis in R 
Data analysis in R and visualization is located in branch “GainovaKristina”, Analysis of DEGs Arabidopsis thaliana.Rmd. This file contains main GSEA graphics (Dotplot, GseaPlot, Enricment map and volcano plot).

## Problems
1. Unfortunately normilization in cuffddif programs works awful! I find that for almost of genes q_value has same values. I decided to find amount of duplicated values in column "q_value" and it equals more 80 % of data!
2. Also using of pipeline cufflinks -> cuffmerge -> cuffdiff convert gene id to XLOC ** id. It makes difficult work with data.
3. Tophat2 v2.1.1 have some problems with error samtools view (samtools view: writing to standard output failedwriting to standard output failed: Broken pipe : Broken pipe samtools view: error closing standard output: -1) and unknown ways to solve this error.

## Conclusions:
Here it was conducted reanalysis of RNA-seq data to identify the core HL-responsive genes. It was found that 2620 DEGs (1260 upregulated and 1365 downregulated genes, padj < 0.05 and abs(log2FC) > 1). It was shown that upregulateg genes are related with different response to stress stimulus, while downregulated genes are involved in biosythesis and differentiation.

## References:
https://rstudio.com/wp-content/uploads/2015/03/rmarkdown-reference.pdf
https://www.bioconductor.org/packages/release/data/annotation/manuals/org.At.tair.db/man/org.At.tair.db.pdf 
https://cran.r-project.org/web/packages/ggridges/vignettes/introduction.html
https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
https://rdrr.io/github/GuangchuangYu/enrichplot/man/pairwise_termsim.html
Martin Morgan (2020). BiocVersion: Set the appropriate version of Bioconductor packages. R package version 3.12.0.
Ghosh S, Chan CK. Analysis of RNA-Seq Data Using TopHat and Cufflinks. Methods Mol Biol. 2016;1374:339-61. doi: 10.1007/978-1-4939-3167-5_18. PMID: 26519415.


3. Kvach Anna



