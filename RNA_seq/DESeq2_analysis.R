#!/usr/bin/env Rscript
#
#Author: Boulanger Mathias
#Last update 02-05-2019
#USAGE: Rscript path/DESeq2_analysis.R -m path/sampleMatrix -s path/Spike-inTable -d parametric -a 0.1
#

###packages
# check packages
packages <- c("optparse", "dplyr", "DESeq2", "ggplot2", "gplots", "ggrepel", "knitr", "magrittr", "affy", "lazyeval", "genefilter", "RColorBrewer", "VennDiagram")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("knitr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("affy"))
suppressPackageStartupMessages(library("lazyeval"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("VennDiagram"))

###Creation of script option and the working environment
#Make option list
option_list = list(
  make_option(c("-m", "--matrixfile"), type="character", default=NULL, help="path/matrixFile_to_analyse", metavar="character"),
  make_option(c("-s", "--spikeinmatrix"), type="character", help="path/spike-in_matrixFile_to_normalize", metavar="character"),
  make_option(c("-d", "--dispersion"), type="character", default="parametric", help="parameter to estimate the dispersion: 'parametric' or 'local' [default= %default]", metavar="character"),
  make_option(c("-a", "--alpha"), type="double", default=0.1, help="alpha factor to calculate the FDR, between 0 and 1 [default= %default]", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Check the ability to work
if (is.null(opt$matrixfile)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (-m input_file).\n", call.=FALSE)
}
if (!file.exists(opt$matrixfile)) {
	print_help(opt_parser)
	stop("The input file does not existe!\n", call.=FALSE)
}
if (!is.null(opt$spikeinmatrix) && !file.exists(opt$matrixfile)) {
	print_help(opt_parser)
	stop("The spike-in input file does not existe!\n", call.=FALSE)
}
if (opt$dispersion != "parametric" && opt$dispersion != "local") {
	print_help(opt_parser)
	stop("The argument in --dispersion is not good (-d parametric or -d local).\n", call.=FALSE)
}
if (class(opt$alpha) != "numeric" ) {
	print_help(opt_parser)
	stop("The argument in --alpha is not numeric (-a 0.1).\n", call.=FALSE)
}
if (opt$alpha > 1 | opt$alpha < 0 ) {
	print_help(opt_parser)
	stop("The number in --alpha is not between 0 and 1 (-a 0.1).\n", call.=FALSE)
}

### Preparing the data and environment
# load the data from HTSeq-count fusion matrix
cat("\n", "Loading data...", "\n")
countsTable_genome <- read.table(file=opt$matrixfile, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.name=1)
if (!is.null(opt$spikeinmatrix)){
	countsTable_ERCC92 <- read.table(file=opt$spikeinmatrix, sep="\t", header=TRUE, stringsAsFactors = FALSE, row.name=1)
}

#Working environment
cat("\n", "Preparing the work environment...", "\n")
matrix_name <- unlist(strsplit(basename(path = opt$matrixfile), "[.]"))
if (!is.null(opt$spikeinmatrix)){
	working_dir <-  as.character(paste0("DESeq2_analysis_", as.character(Sys.Date()), "_", matrix_name[1], "_spike-in_normalized_dispersion_", opt$dispersion, "_alpha_", opt$alpha))
} else {
	working_dir <-  as.character(paste0("DESeq2_analysis_", as.character(Sys.Date()), "_", matrix_name[1], "_dispersion_", opt$dispersion, "_alpha_", opt$alpha))
}
dir.create(working_dir)
setwd(dir = working_dir)
dir.create("Plots")
dir.create("Data")
dir.create("Data/Results")

#Recovery of the number total of read (count, uncount and ambiguous)
cat("\n", "Recovery of the total number of read...", "\n")
count.total=colSums(countsTable_genome)

#Suppression of HTseq-count informative lines from the matrix
cat("\n", "Preparing the matrix...", "\n")
countsTable_genome <- countsTable_genome[- grep("__", row.names.data.frame(countsTable_genome), value=FALSE),]
if (!is.null(opt$spikeinmatrix)){
	countsTable_ERCC92 <- countsTable_ERCC92[- grep("__", row.names.data.frame(countsTable_ERCC92), value=FALSE),]
}

# Experimental design (metadata)
cat("\n", "Creation of the metadata...", "\n", "\n")
label = as.factor(c("DMSO_R1", "DMSO_R2", "DMSO_R3", "DNR_R1", "DNR_R2", "DNR_R3", "ML_R1", "ML_R2", "ML_R3", "ML_DNR_R1", "ML_DNR_R2", "ML_DNR_R3"))
condition = as.factor(c(rep.int("DMSO", 3), rep.int("DNR", 3), rep.int("ML", 3), rep.int("ML_DNR", 3)))
replicats = as.factor(rep(1:3, 4))
col.strain <- c("DMSO"="green","DNR"="red","ML"="blue","ML_DNR"="orange") # Choose color for each replicat
metadata <- data.frame(label = label, condition = condition, replicats = replicats)
metadata$color <- col.strain[as.vector(metadata$condition)]
metadata
cat("\n")
write.table(metadata, "Data/metadata.txt", row.names=F,col.names=T,quote=F,sep="\t")

###Basic statistic
cat("\n", "Calculation of the basic statistics...", "\n")
#graph of counted reads per samples
cat("\n", "Plotting number of read per sample...", "\n")
pdf("Plots/Plot1_num_reads_per_sample.pdf")
barplot(colSums(countsTable_genome)/1000000, 
	main="Total number of genomic reads per sample (million)", 
	col=metadata$color, 
	names.arg =metadata$label, 
	las=1,  horiz=TRUE, 
	ylab="Samples", cex.names=0.5,
	xlab="Million counts")
invisible(dev.off())

# Dimensions
cat("\n", "Calculation of the data distribution...", "\n", "\n")
nb_samples <- ncol(countsTable_genome)
nb_genes <- nrow(countsTable_genome)
p <- paste0("The matrix present ", nb_samples, " samples and ", nb_genes, " genes.")
cat("\n", p, "\n")
stats.per.sample <- data.frame(t(do.call(cbind, lapply(countsTable_genome, summary))))
stats.per.sample$libsum <- apply(countsTable_genome, 2, sum)
stats.per.sample$perc05 <- apply(countsTable_genome, 2, quantile, 0.05)
stats.per.sample$perc10 <- apply(countsTable_genome, 2, quantile, 0.10)
stats.per.sample$perc90 <- apply(countsTable_genome, 2, quantile, 0.90)
stats.per.sample$perc95 <- apply(countsTable_genome, 2, quantile, 0.95)
stats.per.sample$zeros <- apply(countsTable_genome==0, 2, sum)
stats.per.sample$percent.zeros <- 100*stats.per.sample$zeros/nrow(countsTable_genome)
kable(stats.per.sample, caption = "** Table: statistics per sample **", "html") %>% cat(., file = "Data/basic_stats_raw_data.html") 
kable(stats.per.sample, caption = "** Table: statistics per sample **")
cat("\n")

#distribution of counts per gene
cat("\n", "Plotting the data distribution...", "\n")
epsilon <- 1 # pseudo-count to avoid problems with log(0)
pdf("Plots/Plot2_distribution_count_per_gene_raw_data.pdf")
par(mfrow=c(3,1))
hist(as.matrix(countsTable_genome), col="blue", border="white", breaks=100)
hist(as.matrix(countsTable_genome), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
hist(as.matrix(log2(countsTable_genome + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
invisible(dev.off())
# Boxplots
pdf("Plots/Plot3_distribution_count_per_sample_raw_data.pdf")
boxplot(log2(countsTable_genome  + epsilon), col=metadata$color, pch=".", 
	horizontal=TRUE, cex.axis=0.5, las=1, ylab="Samples", xlab="log2(Counts +1)")
invisible(dev.off())

## Density
pdf("Plots/Plot4_density_of_count_raw_data.pdf")
plotDensity(log2(countsTable_genome + epsilon), lty=1, col=metadata$color, lwd=2, xlab="log2(Counts +1)")
grid()
legend("topright", legend=names(col.strain), col=col.strain, lwd=2)
invisible(dev.off())

# take off negative gene in all conditions
cat("\n", "Pre-filtering: elimination of gene with 0 count in all samples...", "\n")
prop.null <- apply(countsTable_genome, 2, function(x) 100*mean(x==0))
pdf("Plots/Plot5_prop_of_null.pdf")
barplot(prop.null, main="Percentage of null counts per sample", horiz=TRUE, 
	cex.names=0.5, las=1, col=metadata$color, ylab='Samples', xlab='% of null counts', names.arg =metadata$label,)
invisible(dev.off())

countsTable_genome_withoutNeg <- countsTable_genome
countsTable_genome_withoutNeg <- countsTable_genome_withoutNeg[rowSums(countsTable_genome_withoutNeg) > 0 ,]

nb_genes_positive <- nrow(countsTable_genome_withoutNeg)
p2 <- paste0("The matrix present ", nb_genes_positive, " genes whith at least 1 count.")
cat("\n", p2, "\n")

#Distributions without neg
cat("\n", "Plotting the data distribution without negative genes...", "\n")
pdf("Plots/Plot6_distribution_count_per_gene_withoutNeg.pdf")
par(mfrow=c(3,1))
hist(as.matrix(countsTable_genome_withoutNeg), col="blue", border="white", breaks=100)
hist(as.matrix(countsTable_genome_withoutNeg), col="blue", border="white",
     breaks=20000, xlim=c(0,2000), main="Counts per gene",
     xlab="Counts (truncated axis)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
hist(as.matrix(log2(countsTable_genome_withoutNeg + epsilon)), breaks=100, col="blue", border="white",
     main="Log2-transformed counts per gene", xlab="log2(counts+1)", ylab="Number of genes", 
     las=1, cex.axis=0.7)
invisible(dev.off())
# Boxplots
pdf("Plots/Plot7_distribution_count_per_sample_withoutNeg.pdf")
boxplot(log2(countsTable_genome_withoutNeg + epsilon), col=metadata$color, pch=".", 
	horizontal=TRUE, cex.axis=0.5, las=1, ylab="Samples", xlab="log2(Counts +1)")
invisible(dev.off())

# Density
pdf("Plots/Plot8_density_of_count_withoutNeg.pdf")
plotDensity(log2(countsTable_genome_withoutNeg + epsilon), lty=1, col=metadata$color, lwd=2, xlab="log2(Counts +1)")
grid()
legend("topright", legend=names(col.strain), col=col.strain, lwd=2)
invisible(dev.off())

#Basic stat if spike-in
if (!is.null(opt$spikeinmatrix)){
	cat("\n", "Making statistics on Spike-in matrix...", "\n")
	pdf("Plots/Plot1.2_num_reads_ERCC92_per_sample.pdf")
	barplot(colSums(countsTable_ERCC92), 
		main="Total number of reads ERCC92 per sample", 
		col=metadata$color, 
		names.arg =metadata$label,  
		las=1,  horiz=TRUE, 
		ylab="Samples", cex.names=0.5,
		xlab="Counts")
	invisible(dev.off())

	#Pourcentage of spike dep of the count
	cat("\n", "Calculation of the proportion of counted reads...", "\n", "\n")
	count.ERCC92=colSums(countsTable_ERCC92)
	count.genome=colSums(countsTable_genome)

	Pour.ERCC92.allReads=((100*count.ERCC92)/count.total)
	norm1=Pour.ERCC92.allReads/Pour.ERCC92.allReads[1]
	Pour.genome.allReads=((100*count.genome)/count.total)
	norm2=Pour.genome.allReads/Pour.genome.allReads[1]

	Pour.ERCC92.countedReads=((100*count.ERCC92)/(count.genome+count.ERCC92))
	norm3=Pour.ERCC92.countedReads/Pour.ERCC92.countedReads[1]
	Pour.genome.countedReads=((100*count.genome)/(count.genome+count.ERCC92))
	norm4=Pour.genome.countedReads/Pour.genome.countedReads[1]

	Spike_stat <- data.frame(count.ERCC92=count.ERCC92, count.genome=count.genome, count.total=count.total, 
		Pour.ERCC92.allReads=Pour.ERCC92.allReads, Pour.ERCC92.allReads.norm.DMSO_R1=norm1, 
		Pour.ERCC92.countedReads=Pour.ERCC92.countedReads, Pour.ERCC92.countedReads.norm.DMSO_R1=norm3, 
		Pour.genome.allReads=Pour.genome.allReads, Pour.genome.allReads.norm.DMSO_R1=norm2, 
		Pour.genome.countedReads=Pour.genome.countedReads, Pour.genome.countedReads.norm.DMSO_R1=norm4,
		row.names=colnames(countsTable_genome))
	kable(Spike_stat, "html") %>% cat(., file = "Data/basic_stats_spike-in.html") 
	kable(Spike_stat)
	cat("\n")
	# take off negative gene in all conditions
	countsTable_genome_ERCC92_withoutNeg <- countsTable_ERCC92
	countsTable_genome_ERCC92_withoutNeg <- countsTable_genome_ERCC92_withoutNeg[rowSums(countsTable_genome_ERCC92_withoutNeg) > 0,]
}

### Analysis
cat("\n", "Starting to analyze data...", "\n")
#creation of dataset
cat("\n", "Creating DESeqDataSet...", "\n")
dds_genome <- DESeqDataSetFromMatrix(countData = countsTable_genome_withoutNeg, colData = metadata, design = ~ condition)
if (!is.null(opt$spikeinmatrix)){
	dds_ERCC92 <- DESeqDataSetFromMatrix(countData = countsTable_genome_ERCC92_withoutNeg, colData = metadata, design = ~ condition)
}

# Principal component analysis before normalization
cat("\n", "Plotting the PCA analysis...", "\n")
rld_dds_genome <- rlog(dds_genome, blind=FALSE)
pdf("Plots/plot9_PCA_analysis_genome.pdf")
plotPCA(rld_dds_genome, intgroup=c("condition")) + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
invisible(dev.off())
if (!is.null(opt$spikeinmatrix)){
	rld_dds_ERCC92 <- rlog(dds_ERCC92, blind=FALSE)
	pdf("Plots/plot9.2_PCA_analysis_spike-in.pdf")	
	plotPCA(rld_dds_ERCC92, intgroup=c("condition")) + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
 	invisible(dev.off())
}

#Find normalization factor with or without spike-in matrix
if (!is.null(opt$spikeinmatrix)){
	cat("\n", "Data normalization with the spike-in...", "\n")
	dds_ERCC92 <- estimateSizeFactors(dds_ERCC92)
	norm.factors <- sizeFactors(dds_ERCC92)	
} else {
	cat("\n", "Data normalization with the library size...", "\n")
	dds_genome <- estimateSizeFactors(dds_genome)
	norm.factors <- sizeFactors(dds_genome)
}

#Application of the normalization factor to genome matrix
cat("\n", "Application of the normalization factors on the matrix...", "\n")
sizeFactors(dds_genome) <- norm.factors
cat("\n", "The factors used are :", "\n", "\n", norm.factors, "\n")
cat("\n", "writing the factors in the metadata...", "\n")
metadata$norm.factors <- norm.factors
write.table(metadata, "Data/metadata.txt", row.names=F,col.names=T,quote=F,sep="\t")

#Checking the normalization
cat("\n", "Plotting the normalized data...", "\n")
pdf("Plots/Plot10_normalization.pdf")
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds_genome)+epsilon),  col=metadata$color, cex.axis=0.7, 
	las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts", names=metadata$label)
boxplot(log2(counts(dds_genome, normalized=TRUE)+epsilon),  col=metadata$color, cex.axis=0.7, 
	las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts", names=metadata$label) 
plotDensity(log2(counts(dds_genome)+epsilon),  col=metadata$color,
	xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
plotDensity(log2(counts(dds_genome, normalized=TRUE)+epsilon), col=metadata$color,
	xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
invisible(dev.off())

if (!is.null(opt$spikeinmatrix)){
	cat("\n", "Plotting the normalized spike-in data...", "\n")
	pdf("Plots/Plot10.2_normalization_ERCC92.pdf")
	par(mfrow=c(2,2),cex.lab=0.7)
	boxplot(log2(counts(dds_ERCC92)+epsilon),  col=metadata$color, cex.axis=0.7,
		las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts", names=metadata$label)
	boxplot(log2(counts(dds_ERCC92, normalized=TRUE)+epsilon),  col=metadata$color, cex.axis=0.7,
		las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts", names=metadata$label) 
	plotDensity(log2(counts(dds_ERCC92)+epsilon),  col=metadata$color,
		xlab="log2(counts)", cex.lab=0.7, panel.first=grid()) 
	plotDensity(log2(counts(dds_ERCC92, normalized=TRUE)+epsilon), col=metadata$color,
		xlab="log2(normalized counts)", cex.lab=0.7, panel.first=grid())
	invisible(dev.off())
}

#Write normalized matrix
cat("\n", "Writing normalized matrix...", "\n")
countsTable_genome.norm <- counts(dds_genome, normalized=TRUE)
write.table(countsTable_genome.norm, "Data/matrix_GRCh37p13_primaryAssembly_reverse_norm.txt", sep = "\t", quote= F, row.names=T, col.names = NA)
pdf("Plots/Plot11_pairs_normalized_matrix.pdf")
pairs(countsTable_genome.norm, upper.panel = NULL, panel=panel.smooth, main="Counts comparison among samples", pch=21, bg="dodgerblue")
pairs(log2(countsTable_genome.norm+1), upper.panel = NULL, panel=panel.smooth,main="Counts comparison among samples \n (log2)", pch=21, bg="dodgerblue")
par(xpd=TRUE)
invisible(dev.off())

if (!is.null(opt$spikeinmatrix)){
	countsTable_ERCC92.norm <- counts(dds_ERCC92, normalized=TRUE)
	write.table(countsTable_ERCC92.norm, "Data/matrix_ERCC92_reverse_norm.txt", sep = "\t", quote= F, row.names=T, col.names = NA)
	pdf("Plots/Plot11.2_pairs_normalized_spike-in_matrix.pdf")
	pairs(countsTable_ERCC92.norm, upper.panel = NULL, panel=panel.smooth, main="Counts comparison among samples", pch=21, bg="dodgerblue")
	pairs(log2(countsTable_ERCC92.norm+1), upper.panel = NULL, panel=panel.smooth,main="Counts comparison among samples \n (log2)", pch=21, bg="dodgerblue")
	par(xpd=TRUE)
	invisible(dev.off())
}

cat("\n", "Plotting the PCA analysis after normalization...", "\n")
rld_dds_genome2 <- rlog(dds_genome, blind=FALSE)
pdf("Plots/Plot12_PCA_analysis_genome_normalized.pdf")
plotPCA(rld_dds_genome2, intgroup=c("condition")) + geom_text_repel(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
invisible(dev.off())


# Computing mean and variance
cat("\n", "Checking the relation between the mean and the variance...", "\n")
norm.counts <- counts(dds_genome, normalized=TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)
norm.counts.stats <- data.frame(
  min=apply(norm.counts, 2, min),
  mean=apply(norm.counts, 2, mean),
  median=apply(norm.counts, 2, median),
  max=apply(norm.counts, 2, max),
  zeros=apply(norm.counts==0, 2, sum),
  percent.zeros=100*apply(norm.counts==0, 2, sum)/nrow(norm.counts),
  perc05=apply(norm.counts, 2, quantile, 0.05),
  perc10=apply(norm.counts, 2, quantile, 0.10),
  perc90=apply(norm.counts, 2, quantile, 0.90),
  perc95=apply(norm.counts, 2, quantile, 0.95)
)
if (!is.null(opt$spikeinmatrix)){
	kable(norm.counts.stats, caption = "** Table: statistics per sample normalized **", "html") %>% cat(., file = "Data/basic_stats_normalized_with_spike.html")
} else {
	kable(norm.counts.stats, caption = "** Table: statistics per sample normalized **", "html") %>% cat(., file = "Data/basic_stats_ERCC92_normalized_with_library_size.html")
}

# Mean and variance relationship
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
pdf("Plots/Plot13_mean_and_variance_relation_normalization.pdf")
plot(x=log2(mean.counts), y=log2(variance.counts), pch=16, cex=0.5, 
     col=mean.var.col, main="Mean-variance relationship",
     xlab="Mean log2(normalized counts) per gene",
     ylab="Variance of log2(normalized counts)",
     panel.first = grid())
abline(a=0, b=1, col="brown")
invisible(dev.off())

# Performing estimation of dispersion parameter
cat("\n", "Estimation of the data dispersion with the parameter: ", opt$dispersion, "...", "\n", "\n")
dds_genome.disp <- estimateDispersions(dds_genome, fitType=opt$dispersion)
pdf(paste0("Plots/Plot14_dispersion_estimation_", opt$dispersion, ".pdf"))
plotDispEsts(dds_genome.disp)
invisible(dev.off())

cat("\n", "Results calling through the negative binomial law...", "\n")
wald.test <- nbinomWaldTest(dds_genome.disp)
baseMeanDMSO <- rowMeans(counts(wald.test,normalized=TRUE)[,wald.test$condition == "DMSO"])
baseMeanDNR <- rowMeans(counts(wald.test,normalized=TRUE)[,wald.test$condition == "DNR"])
baseMeanML <- rowMeans(counts(wald.test,normalized=TRUE)[,wald.test$condition == "ML"])
baseMeanML_DNR <- rowMeans(counts(wald.test,normalized=TRUE)[,wald.test$condition == "ML_DNR"])

#Function to make volcano plot
plotVolcano <- function(res, plotGeneName=FALSE, geneName="normal", col="black"){
	res1 <- as.data.frame(res)
	res1 <- mutate(res1, gene=rownames(res))
	#cbind(rownames(res), as.data.frame(res))
	if(col == "density"){
		cols <- densCols(res1$log2FoldChange, -log10(res1$pvalue))
		volcano <- ggplot(as.data.frame(res1), aes(x = log2FoldChange, y = -log10(padj))) + 
		geom_point(color = cols, cex = 0.6)
	} else if (col == "DNR") {
		res.dnr <- as.data.frame(res1.DESeq2)
		res.dnr <- mutate(as.data.frame(res.dnr), sig=ifelse(log2FoldChange >= 1 & padj <= opt$alpha , "Up-regulate DNR vs DMSO", ifelse(log2FoldChange <= -1 & padj <= opt$alpha , "Down-regulate DNR vs DMSO", "Not Sig")))
		res1$sig <- as.factor(res.dnr$sig)
		volcano <- ggplot(as.data.frame(res1), aes(x = log2FoldChange, y = -log10(padj))) + 
		geom_point(aes(col = sig), cex = 0.6) +
		scale_colour_manual(values = c("Up-regulate DNR vs DMSO"= "red", "Down-regulate DNR vs DMSO"="green", "Not Sig"="blue"))

	} else if (col == "significative"){
		res1 <- mutate(as.data.frame(res1), sig=ifelse(abs(log2FoldChange) >= 1 & padj <= opt$alpha , "FDR < 0.05 and log2(FC) > |1|", "Not Sig"))
		volcano <- ggplot(as.data.frame(res1), aes(x = log2FoldChange, y = -log10(padj))) + 
		geom_point(aes(col = sig), cex = 0.6) +
		scale_colour_manual(values = c("FDR < 0.05 and log2(FC) > |1|"= "red", "Not Sig"="blue"))
	} else if (col == "black"){
		volcano <- ggplot(as.data.frame(res1), aes(x = log2FoldChange, y = -log10(padj))) + 
		geom_point(color = "black", cex = 0.6)
	}

	volcano <- volcano +
	scale_x_continuous(breaks = c(-4, -2, 0, 2, 4, 6, 8)) +
	scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5), legend.position="none", panel.grid.minor = element_line(linetype = "dotted")) +
	labs(title=paste0("Volcano plot ", main), x ="Effect size: log2(fold-change)", y = "-log10(adjusted p-value)") +
	geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
	geom_vline(xintercept = 1, col = "brown", linetype = "dashed") +
	geom_vline(xintercept = -1, col = "brown", linetype = "dashed") +
	geom_hline(yintercept = -log10(opt$alpha), col = "brown", linetype = "dashed")

	if(plotGeneName == TRUE){
		if(geneName=="normal"){
			res2 <- subset(res1, abs(log2FoldChange) >= 2 & padj <= opt$alpha)
		}
		if(geneName=="full"){
			res2 <- subset(res1, abs(log2FoldChange) >= 1 & padj <= opt$alpha)
		}
		if(geneName=="spe"){
			res2 <- subset(res1, abs(log2FoldChange) >= 1 & padj <= 1e-40 | 
				abs(log2FoldChange) >= 2 & padj <= 1e-10 |
				log2FoldChange >= 4 & padj <= opt$alpha | 
				log2FoldChange <= -3 & padj <= opt$alpha)
		}		 
		volcano + geom_text_repel(data = res2, aes(label = gene), color = "black", size=2, segment.size=0.2, segment.alpha=0.6)
	} else {
		volcano
	}
}



#DNR vs DMSO
cat("\n", "Results calling DNR vs DMSO", "\n")
res1.DESeq2 <- results(wald.test, alpha=opt$alpha, pAdjustMethod="BH", contrast=c("condition","DNR","DMSO"))
cat("\n", "Results with FDR < ", opt$alpha, " - DNR vs DMSO", "\n")
table(res1.DESeq2$padj < opt$alpha)
cat("\n", "Results with p-value < ", opt$alpha, " - DNR vs DMSO", "\n")
table(res1.DESeq2$pvalue < opt$alpha)

# Draw an histogram of the p-values adjusted
cat("\n", "Checking p-value and FDR distribution...", "\n")
pdf("Plots/1-DNRvsDMSO_plot_FDR_distribution.pdf")
hist(res1.DESeq2$padj, breaks=20, col="grey", main="DESeq2 FDR distribution DNR vs DMSO", xlab="DESeq2 FDR", ylab="Number of genes")
invisible(dev.off())
pdf("Plots/1-DNRvsDMSO_plot_pvalue_distribution.pdf")
hist(res1.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution DNR vs DMSO", xlab="DESeq2 P-value", ylab="Number of genes")
invisible(dev.off())

#Volcano plot
cat("\n", "Generating the volcano plot...", "\n")
main = as.character("DNR vs DMSO")

pdf("Plots/1-DNRvsDMSO_volcano_wo_name.pdf")
plotVolcano(res1.DESeq2, plotGeneName=FALSE, geneName="spe", col="black")
invisible(dev.off())
pdf("Plots/1-DNRvsDMSO_volcano.pdf")
plotVolcano(res1.DESeq2, plotGeneName=TRUE, geneName="spe", col="black")
invisible(dev.off())

pdf("Plots/1-DNRvsDMSO_volcano_sig_wo_name.pdf")
plotVolcano(res1.DESeq2, plotGeneName=FALSE, geneName="spe", col="significative")
invisible(dev.off())
pdf("Plots/1-DNRvsDMSO_volcano_sig.pdf")
plotVolcano(res1.DESeq2, plotGeneName=TRUE, geneName="spe", col="significative")
invisible(dev.off())

pdf("Plots/1-DNRvsDMSO_volcano_density_wo_name.pdf")
plotVolcano(res1.DESeq2, plotGeneName=FALSE, geneName="spe", col="density")
invisible(dev.off())
pdf("Plots/1-DNRvsDMSO_volcano_density.pdf")
plotVolcano(res1.DESeq2, plotGeneName=TRUE, geneName="spe", col="density")
invisible(dev.off())

pdf("Plots/1-DNRvsDMSO_volcano_dnrColor_wo_name.pdf")
plotVolcano(res1.DESeq2, plotGeneName=FALSE, geneName="spe", col="DNR")
invisible(dev.off())
pdf("Plots/1-DNRvsDMSO_volcano_dnrColor.pdf")
plotVolcano(res1.DESeq2, plotGeneName=TRUE, geneName="spe", col="DNR")
invisible(dev.off())

#MA plot
pdf("Plots/1-DNRvsDMSO_plotMA_dispersion.pdf")
plotMA(res1.DESeq2, colNonSig = "grey")
abline(h=c(-1:1), col="red")
invisible(dev.off())



#ML vs DMSO
cat("\n", "Results calling ML vs DMSO", "\n")
res2.DESeq2 <- results(wald.test, alpha=opt$alpha, pAdjustMethod="BH", contrast=c("condition","ML","DMSO"))
cat("\n", "Results with FDR < ", opt$alpha, " - ML vs DMSO", "\n")
table(res2.DESeq2$padj < opt$alpha)
cat("\n", "Results with p-value < ", opt$alpha, " - ML vs DMSO", "\n")
table(res2.DESeq2$pvalue < opt$alpha)

# Draw an histogram of the p-values adjusted
cat("\n", "Checking p-value and FDR distribution...", "\n")
pdf("Plots/2-MLvsDMSO_plot_FDR_distribution.pdf")
hist(res2.DESeq2$padj, breaks=20, col="grey", main="DESeq2 FDR distribution ML vs DMSO", xlab="DESeq2 FDR", ylab="Number of genes")
invisible(dev.off())
pdf("Plots/2-MLvsDMSO_plot_pvalue_distribution.pdf")
hist(res2.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution ML vs DMSO", xlab="DESeq2 P-value", ylab="Number of genes")
invisible(dev.off())

#Volcano plot
cat("\n", "Generating the volcano plot...", "\n")
main = as.character("ML vs DMSO")

pdf("Plots/2-MLvsDMSO_volcano_wo_name.pdf")
plotVolcano(res2.DESeq2, plotGeneName=FALSE, geneName="full", col="black")
invisible(dev.off())
pdf("Plots/2-MLvsDMSO_volcano.pdf")
plotVolcano(res2.DESeq2, plotGeneName=TRUE, geneName="full", col="black")
invisible(dev.off())

pdf("Plots/2-MLvsDMSO_volcano_sig_wo_name.pdf")
plotVolcano(res2.DESeq2, plotGeneName=FALSE, geneName="full", col="significative")
invisible(dev.off())
pdf("Plots/2-MLvsDMSO_volcano_sig.pdf")
plotVolcano(res2.DESeq2, plotGeneName=TRUE, geneName="full", col="significative")
invisible(dev.off())

pdf("Plots/2-MLvsDMSO_volcano_density_wo_name.pdf")
plotVolcano(res2.DESeq2, plotGeneName=FALSE, geneName="full", col="density")
invisible(dev.off())
pdf("Plots/2-MLvsDMSO_volcano_density.pdf")
plotVolcano(res2.DESeq2, plotGeneName=TRUE, geneName="full", col="density")
invisible(dev.off())

pdf("Plots/2-MLvsDMSO_volcano_dnrColor_wo_name.pdf")
plotVolcano(res2.DESeq2, plotGeneName=FALSE, geneName="full", col="DNR")
invisible(dev.off())
pdf("Plots/2-MLvsDMSO_volcano_dnrColor.pdf")
plotVolcano(res2.DESeq2, plotGeneName=TRUE, geneName="full", col="DNR")
invisible(dev.off())

#MA plot
pdf("Plots/2-MLvsDMSO_plotMA_dispersion.pdf")
plotMA(res2.DESeq2, colNonSig = "grey")
abline(h=c(-1:1), col="red")
invisible(dev.off())



#ML_DNR vs DMSO
cat("\n", "Results calling ML_DNR vs DMSO", "\n")
res3.DESeq2 <- results(wald.test, alpha=opt$alpha, pAdjustMethod="BH", contrast=c("condition","ML_DNR","DMSO"))
cat("\n", "Results with FDR < ", opt$alpha, " - ML_DNR vs DMSO", "\n")
table(res3.DESeq2$padj < opt$alpha)
cat("\n", "Results with p-value < ", opt$alpha, " - ML_DNR vs DMSO", "\n")
table(res3.DESeq2$pvalue < opt$alpha)

# Draw an histogram of the p-values adjusted
cat("\n", "Checking p-value and FDR distribution...", "\n")
pdf("Plots/3-ML_DNRvsDMSO_plot_FDR_distribution.pdf")
hist(res3.DESeq2$padj, breaks=20, col="grey", main="DESeq2 FDR distribution ML_DNR vs DMSO", xlab="DESeq2 FDR", ylab="Number of genes")
invisible(dev.off())
pdf("Plots/3-ML_DNRvsDMSO_plot_pvalue_distribution.pdf")
hist(res3.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution ML_DNR vs DMSO", xlab="DESeq2 P-value", ylab="Number of genes")
invisible(dev.off())

#Volcano plot
cat("\n", "Generating the volcano plot...", "\n")
main = as.character("ML DNR vs DMSO")

pdf("Plots/3-ML_DNRvsDMSO_volcano_wo_name.pdf")
plotVolcano(res3.DESeq2, plotGeneName=FALSE, geneName="spe", col="black")
invisible(dev.off())
pdf("Plots/3-ML_DNRvsDMSO_volcano.pdf")
plotVolcano(res3.DESeq2, plotGeneName=TRUE, geneName="spe", col="black")
invisible(dev.off())

pdf("Plots/3-ML_DNRvsDMSO_volcano_sig_wo_name.pdf")
plotVolcano(res3.DESeq2, plotGeneName=FALSE, geneName="spe", col="significative")
invisible(dev.off())
pdf("Plots/3-ML_DNRvsDMSO_volcano_sig.pdf")
plotVolcano(res3.DESeq2, plotGeneName=TRUE, geneName="spe", col="significative")
invisible(dev.off())

pdf("Plots/3-ML_DNRvsDMSO_volcano_density_wo_name.pdf")
plotVolcano(res3.DESeq2, plotGeneName=FALSE, geneName="spe", col="density")
invisible(dev.off())
pdf("Plots/3-ML_DNRvsDMSO_volcano_density.pdf")
plotVolcano(res3.DESeq2, plotGeneName=TRUE, geneName="spe", col="density")
invisible(dev.off())

pdf("Plots/3-ML_DNRvsDMSO_volcano_dnrColor_wo_name.pdf")
plotVolcano(res3.DESeq2, plotGeneName=FALSE, geneName="spe", col="DNR")
invisible(dev.off())
pdf("Plots/3-ML_DNRvsDMSO_volcano_dnrColor.pdf")
plotVolcano(res3.DESeq2, plotGeneName=TRUE, geneName="spe", col="DNR")
invisible(dev.off())

#MA plot
pdf("Plots/3-ML_DNRvsDMSO_plotMA_dispersion.pdf")
plotMA(res3.DESeq2, colNonSig = "grey")
abline(h=c(-1:1), col="red")
invisible(dev.off())



#ML_DNR vs ML
cat("\n", "Results calling ML_DNR vs ML", "\n")
res4.DESeq2 <- results(wald.test, alpha=opt$alpha, pAdjustMethod="BH", contrast=c("condition","ML_DNR","ML"))
cat("\n", "Results with FDR < ", opt$alpha, " - ML_DNR vs ML", "\n")
table(res4.DESeq2$padj < opt$alpha)
cat("\n", "Results with p-value < ", opt$alpha, " - ML_DNR vs ML", "\n")
table(res4.DESeq2$pvalue < opt$alpha)

# Draw an histogram of the p-values adjusted
cat("\n", "Checking p-value and FDR distribution...", "\n")
pdf("Plots/4-ML_DNRvsML_plot_FDR_distribution.pdf")
hist(res4.DESeq2$padj, breaks=20, col="grey", main="DESeq2 FDR distribution ML_DNR vs ML", xlab="DESeq2 FDR", ylab="Number of genes")
invisible(dev.off())
pdf("Plots/4-ML_DNRvsML_plot_pvalue_distribution.pdf")
hist(res4.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution ML_DNR vs ML", xlab="DESeq2 P-value", ylab="Number of genes")
invisible(dev.off())

#Volcano plot
cat("\n", "Generating the volcano plot...", "\n")
main = as.character("ML DNR vs ML")

pdf("Plots/4-ML_DNRvsML_volcano_wo_name.pdf")
plotVolcano(res4.DESeq2, plotGeneName=FALSE, geneName="spe", col="black")
invisible(dev.off())
pdf("Plots/4-ML_DNRvsML_volcano.pdf")
plotVolcano(res4.DESeq2, plotGeneName=TRUE, geneName="spe", col="black")
invisible(dev.off())

pdf("Plots/4-ML_DNRvsML_volcano_sig_wo_name.pdf")
plotVolcano(res4.DESeq2, plotGeneName=FALSE, geneName="spe", col="significative")
invisible(dev.off())
pdf("Plots/4-ML_DNRvsML_volcano_sig.pdf")
plotVolcano(res4.DESeq2, plotGeneName=TRUE, geneName="spe", col="significative")
invisible(dev.off())

pdf("Plots/4-ML_DNRvsML_volcano_density_wo_name.pdf")
plotVolcano(res4.DESeq2, plotGeneName=FALSE, geneName="spe", col="density")
invisible(dev.off())
pdf("Plots/4-ML_DNRvsML_volcano_density.pdf")
plotVolcano(res4.DESeq2, plotGeneName=TRUE, geneName="spe", col="density")
invisible(dev.off())

pdf("Plots/4-ML_DNRvsML_volcano_dnrColor_wo_name.pdf")
plotVolcano(res4.DESeq2, plotGeneName=FALSE, geneName="spe", col="DNR")
invisible(dev.off())
pdf("Plots/4-ML_DNRvsML_volcano_dnrColor.pdf")
plotVolcano(res4.DESeq2, plotGeneName=TRUE, geneName="spe", col="DNR")
invisible(dev.off())

#MA plot
pdf("Plots/4-ML_DNRvsML_plotMA_dispersion.pdf")
plotMA(res4.DESeq2, colNonSig = "grey")
abline(h=c(-1:1), col="red")
invisible(dev.off())


#ML_DNR vs DNR
cat("\n", "Results calling ML_DNR vs DNR", "\n")
res5.DESeq2 <- results(wald.test, alpha=opt$alpha, pAdjustMethod="BH", contrast=c("condition","ML_DNR","DNR"))
cat("\n", "Results with FDR < ", opt$alpha, " - ML_DNR vs DNR", "\n")
table(res5.DESeq2$padj < opt$alpha)
cat("\n", "Results with p-value < ", opt$alpha, " - ML_DNR vs DNR", "\n")
table(res5.DESeq2$pvalue < opt$alpha)

# Draw an histogram of the p-values adjusted
cat("\n", "Checking p-value and FDR distribution...", "\n")
pdf("Plots/5-ML_DNRvsDNR_plot_FDR_distribution.pdf")
hist(res5.DESeq2$padj, breaks=20, col="grey", main="DESeq2 FDR distribution ML_DNR vs DNR", xlab="DESeq2 FDR", ylab="Number of genes")
invisible(dev.off())
pdf("Plots/5-ML_DNRvsDNR_plot_pvalue_distribution.pdf")
hist(res5.DESeq2$padj, breaks=20, col="grey", main="DESeq2 p-value distribution ML_DNR vs DNR", xlab="DESeq2 P-value", ylab="Number of genes")
invisible(dev.off())

#Volcano plot
cat("\n", "Generating the volcano plot...", "\n")
main = as.character("ML DNR vs DNR")

pdf("Plots/5-ML_DNRvsDNR_volcano_wo_name.pdf")
plotVolcano(res5.DESeq2, plotGeneName=FALSE, geneName="full", col="black")
invisible(dev.off())
pdf("Plots/5-ML_DNRvsDNR_volcano.pdf")
plotVolcano(res5.DESeq2, plotGeneName=TRUE, geneName="full", col="black")
invisible(dev.off())

pdf("Plots/5-ML_DNRvsDNR_volcano_sig_wo_name.pdf")
plotVolcano(res5.DESeq2, plotGeneName=FALSE, geneName="full", col="significative")
invisible(dev.off())
pdf("Plots/5-ML_DNRvsDNR_volcano_sig.pdf")
plotVolcano(res5.DESeq2, plotGeneName=TRUE, geneName="full", col="significative")
invisible(dev.off())

pdf("Plots/5-ML_DNRvsDNR_volcano_density_wo_name.pdf")
plotVolcano(res5.DESeq2, plotGeneName=FALSE, geneName="full", col="density")
invisible(dev.off())
pdf("Plots/5-ML_DNRvsDNR_volcano_density.pdf")
plotVolcano(res5.DESeq2, plotGeneName=TRUE, geneName="full", col="density")
invisible(dev.off())

pdf("Plots/5-ML_DNRvsDNR_volcano_dnrColor_wo_name.pdf")
plotVolcano(res5.DESeq2, plotGeneName=FALSE, geneName="full", col="DNR")
invisible(dev.off())
pdf("Plots/5-ML_DNRvsDNR_volcano_dnrColor.pdf")
plotVolcano(res5.DESeq2, plotGeneName=TRUE, geneName="full", col="DNR")
invisible(dev.off())

#MA plot
pdf("Plots/5-ML_DNRvsDNR_plotMA_dispersion.pdf")
plotMA(res5.DESeq2, colNonSig = "grey")
abline(h=c(-1:1), col="red")
invisible(dev.off())


#Generation of the result tables
cat("\n", "Generating the final table...", "\n")
colnames <- as.character(c("genes", "baseMeanDMSO", "baseMeanDNR", "baseMeanML", "baseMeanML_DNR", "baseMean", 
	"DNRvsDMSO.log2FoldChange", "DNRvsDMSO.lfcSE", "DNRvsDMSO.stat", "DNRvsDMSO.pvalue", "DNRvsDMSO.padj", 
	"MLvsDMSO.log2FoldChange", "MLvsDMSO.lfcSE", "MLvsDMSO.stat", "MLvsDMSO.pvalue", "MLvsDMSO.padj", 
	"ML_DNRvsDMSO.log2FoldChange", "ML_DNRvsDMSO.lfcSE", "ML_DNRvsDMSO.stat", "ML_DNRvsDMSO.pvalue", "ML_DNRvsDMSO.padj", 
	"ML_DNRvsML.log2FoldChange", "ML_DNRvsML.lfcSE", "ML_DNRvsML.stat", "ML_DNRvsML.pvalue", "ML_DNRvsML.padj", 
	"ML_DNRvsDNR.log2FoldChange", "ML_DNRvsDNR.lfcSE", "ML_DNRvsDNR.stat", "ML_DNRvsDNR.pvalue", "ML_DNRvsDNR.padj"))
res.final = cbind(genes=rownames(res1.DESeq2), baseMeanDMSO, baseMeanDNR, baseMeanML, baseMeanML_DNR, as.data.frame(res1.DESeq2), as.data.frame(res2.DESeq2[,-1]), as.data.frame(res3.DESeq2[,-1]), as.data.frame(res4.DESeq2[,-1]), as.data.frame(res5.DESeq2[,-1]))
names(res.final) <- colnames
write.table(as.data.frame(res.final),"Data/Results/results.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res.final[,1]),"Data/Results/results.list", row.names=F, col.names=F, quote=F, sep="\t")
write.csv2(as.data.frame(res.final), file = "Data/Results/results.csv", row.names=F, quote=F)


res1.sig <- subset(res.final, abs(DNRvsDMSO.log2FoldChange) >= 1 & DNRvsDMSO.padj <= opt$alpha)
res1.sig.up <- subset(res.final, DNRvsDMSO.log2FoldChange >= 1 & DNRvsDMSO.padj <= opt$alpha)
res1.sig.down <- subset(res.final, DNRvsDMSO.log2FoldChange <= -1 & DNRvsDMSO.padj <= opt$alpha)
write.table(as.data.frame(res1.sig),"Data/Results/results_DNRvsDMSO_sig.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res1.sig.up),"Data/Results/results_DNRvsDMSO_sig_upReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res1.sig.down),"Data/Results/results_DNRvsDMSO_sig_downReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res1.sig[,1]),"Data/Results/results_DNRvsDMSO_sig.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res1.sig.up[,1]),"Data/Results/results_DNRvsDMSO_sig_upReg.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res1.sig.down[,1]),"Data/Results/results_DNRvsDMSO_sig_downReg.list", row.names=F, col.names=F, quote=F, sep="\t")

res2.sig <- subset(res.final, abs(MLvsDMSO.log2FoldChange) >= 1 & MLvsDMSO.padj <= opt$alpha)
res2.sig.up <- subset(res.final, MLvsDMSO.log2FoldChange >= 1 & MLvsDMSO.padj <= opt$alpha)
res2.sig.down <- subset(res.final, MLvsDMSO.log2FoldChange <= -1 & MLvsDMSO.padj <= opt$alpha)
write.table(as.data.frame(res2.sig),"Data/Results/results_MLvsDMSO_sig.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res2.sig.up),"Data/Results/results_MLvsDMSO_sig_upReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res2.sig.down),"Data/Results/results_MLvsDMSO_sig_downReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res2.sig[,1]),"Data/Results/results_MLvsDMSO_sig.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res2.sig.up[,1]),"Data/Results/results_MLvsDMSO_sig_upReg.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res2.sig.down[,1]),"Data/Results/results_MLvsDMSO_sig_downReg.list", row.names=F, col.names=F, quote=F, sep="\t")

res3.sig <- subset(res.final, abs(ML_DNRvsDMSO.log2FoldChange) >= 1 & ML_DNRvsDMSO.padj <= opt$alpha)
res3.sig.up <- subset(res.final, ML_DNRvsDMSO.log2FoldChange >= 1 & ML_DNRvsDMSO.padj <= opt$alpha)
res3.sig.down <- subset(res.final, ML_DNRvsDMSO.log2FoldChange <= -1 & ML_DNRvsDMSO.padj <= opt$alpha)
write.table(as.data.frame(res3.sig),"Data/Results/results_ML_DNRvsDMSO_sig.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res3.sig.up),"Data/Results/results_ML_DNRvsDMSO_sig_upReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res3.sig.down),"Data/Results/results_ML_DNRvsDMSO_sig_downReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res3.sig[,1]),"Data/Results/results_ML_DNRvsDMSO_sig.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res3.sig.up[,1]),"Data/Results/results_ML_DNRvsDMSO_sig_upReg.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res3.sig.down[,1]),"Data/Results/results_ML_DNRvsDMSO_sig_downReg.list", row.names=F, col.names=F, quote=F, sep="\t")

res4.sig <- subset(res.final, abs(ML_DNRvsML.log2FoldChange) >= 1 & ML_DNRvsML.padj <= opt$alpha)
res4.sig.up <- subset(res.final, ML_DNRvsML.log2FoldChange >= 1 & ML_DNRvsML.padj <= opt$alpha)
res4.sig.down <- subset(res.final, ML_DNRvsML.log2FoldChange <= -1 & ML_DNRvsML.padj <= opt$alpha)
write.table(as.data.frame(res4.sig),"Data/Results/results_ML_DNRvsML_sig.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res4.sig.up),"Data/Results/results_ML_DNRvsML_sig_upReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res4.sig.down),"Data/Results/results_ML_DNRvsML_sig_downReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res4.sig[,1]),"Data/Results/results_ML_DNRvsML_sig.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res4.sig.up[,1]),"Data/Results/results_ML_DNRvsML_sig_upReg.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res4.sig.down[,1]),"Data/Results/results_ML_DNRvsML_sig_downReg.list", row.names=F, col.names=F, quote=F, sep="\t")

res5.sig <- subset(res.final, abs(ML_DNRvsDNR.log2FoldChange) >= 1 & ML_DNRvsDNR.padj <= opt$alpha)
res5.sig.up <- subset(res.final, ML_DNRvsDNR.log2FoldChange >= 1 & ML_DNRvsDNR.padj <= opt$alpha)
res5.sig.down <- subset(res.final, ML_DNRvsDNR.log2FoldChange <= -1 & ML_DNRvsDNR.padj <= opt$alpha)
write.table(as.data.frame(res5.sig),"Data/Results/results_ML_DNRvsDNR_sig.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res5.sig.up),"Data/Results/results_ML_DNRvsDNR_sig_upReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res5.sig.down),"Data/Results/results_ML_DNRvsDNR_sig_downReg.table", row.names=F, col.names=T, quote=F, sep="\t")
write.table(as.data.frame(res5.sig[,1]),"Data/Results/results_ML_DNRvsDNR_sig.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res5.sig.up[,1]),"Data/Results/results_ML_DNRvsDNR_sig_upReg.list", row.names=F, col.names=F, quote=F, sep="\t")
write.table(as.data.frame(res5.sig.down[,1]),"Data/Results/results_ML_DNRvsDNR_sig_downReg.list", row.names=F, col.names=F, quote=F, sep="\t")

#plot data for a specific gene
cat("\n", "Generating plot for NFkB2...", "\n")
pdf("Plots/Plot15_NFKB2_all_replicats.pdf")
gn.most.diff.val <- counts(wald.test, normalized=T)["NFKB2",]
barplot(gn.most.diff.val, col=metadata$color, main="NFKB2", las=2, cex.names=0.5)
invisible(dev.off())

cat("\n", "Generating plot for IER3...", "\n")
pdf("Plots/Plot16_IER3_all_replicats.pdf")
gn.most.diff.val <- counts(wald.test, normalized=T)["IER3",]
barplot(gn.most.diff.val, col=metadata$color, main="IER3", las=2, cex.names=0.5)
invisible(dev.off())

cat("\n", "Generating plot for TNFRSF1A...", "\n")
pdf("Plots/Plot17_TNFRSF1A_all_replicats.pdf")
gn.most.diff.val <- counts(wald.test, normalized=T)["TNFRSF1A",]
barplot(gn.most.diff.val, col=metadata$color, main="TNFRSF1A", las=2, cex.names=0.5)
invisible(dev.off())

cat("\n", "Generating plot for CEBPE...", "\n")
pdf("Plots/Plot18_CEBPE_all_replicats.pdf")
gn.most.diff.val <- counts(wald.test, normalized=T)["CEBPE",]
barplot(gn.most.diff.val, col=metadata$color, main="CEBPE", las=2, cex.names=0.5)
invisible(dev.off())
#heatmap
#1
cat("\n", "Generating heatmap 1 ...", "\n")
pdf("Plots/heatmap1_topVar_50.pdf")
topVarGenes <- head(order(rowVars(assay(rld_dds_genome2)), decreasing=TRUE), 50)
heatmap.2(assay(rld_dds_genome2)[topVarGenes,], scale="row", trace="none", dendrogram="none", Colv=FALSE ,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())
pdf("Plots/heatmap1_topVar_500.pdf")
topVarGenes <- head(order(rowVars(assay(rld_dds_genome2)), decreasing=TRUE), 500)
heatmap.2(assay(rld_dds_genome2)[topVarGenes,], scale="row", trace="none", dendrogram="none", Colv=FALSE ,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())


rld_dds_genome3=cbind(rld_dds_genome2[,c(1:3,7:9,4:6,10:12)])
pdf("Plots/heatmap1_topVar_50_2.pdf")
topVarGenes <- head(order(rowVars(assay(rld_dds_genome3)), decreasing=TRUE), 50)
heatmap.2(assay(rld_dds_genome3)[topVarGenes,], scale="row", trace="none", dendrogram="none", Colv=FALSE ,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())


#recovery of gene up ans down in DNR treatment
res1.keptup <- res1.DESeq2[res1.DESeq2$padj <= opt$alpha & !is.na(res1.DESeq2$padj) & res1.DESeq2$log2FoldChange > 1,]
res1.keptdown <- res1.DESeq2[res1.DESeq2$padj <= opt$alpha & !is.na(res1.DESeq2$padj) & res1.DESeq2$log2FoldChange < -1,]
res1.keptup <- res1.keptup[order(res1.keptup[,2], decreasing=T),]
res1.keptdown <- res1.keptdown[order(res1.keptdown[,2], decreasing=F),]

gene.kept <- c(rownames(head(res1.keptup, 25)), rownames(head(res1.keptdown, 25)))
gene.kept1 <- c(rownames(head(res1.keptup, 50)), rownames(head(res1.keptdown, 50)))
gene.keptup <- rownames(head(res1.keptup, 50))
gene.keptdown <- rownames(head(res1.keptdown, 50))

pdf("Plots/heatmap1_DNR_25up_25down.pdf")
heatmap.2(assay(rld_dds_genome3)[gene.kept,], scale="row", trace="none", dendrogram="none", Colv=FALSE , Rowv=FALSE,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())
pdf("Plots/heatmap1_DNR_50up_50down.pdf")
heatmap.2(assay(rld_dds_genome3)[gene.kept1,], scale="row", trace="none", dendrogram="none", Colv=FALSE , Rowv=FALSE,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())
pdf("Plots/heatmap1_DNR_50up.pdf")
heatmap.2(assay(rld_dds_genome3)[gene.keptup,], scale="row", trace="none", dendrogram="none", Colv=FALSE , Rowv=FALSE,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())
pdf("Plots/heatmap1_DNR_50down.pdf")
heatmap.2(assay(rld_dds_genome3)[gene.keptdown,], scale="row", trace="none", dendrogram="none", Colv=FALSE , Rowv=FALSE,
	col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), cexRow = 0.5 , cexCol = 0.5)
invisible(dev.off())


#2
cat("\n", "Generating heatmap 2 ...", "\n")
## We select gene names based on FDR (1%) DNR vs DMSO
gene.kept <- rownames(res1.DESeq2)[res1.DESeq2$padj <= opt$alpha & !is.na(res1.DESeq2$padj)]

## We retrieve the normalized counts for gene of interest
count.table.kept <- log2(countsTable_genome.norm + epsilon)[gene.kept, ]
dim(count.table.kept)

## Perform the hierarchical clustering with a distance based on Pearson-correlation coefficient and average linkage clustering as agglomeration criteria
pdf("Plots/heatmap2_hierarchical_clustering.pdf")
heatmap.2(as.matrix(count.table.kept), scale="row", hclust=function(x) hclust(x,method="average"), 
	distfun=function(x) as.dist((1-cor(t(x)))/2), trace="none", density="none", labRow="", cexCol=0.5)
invisible(dev.off())


save.image("imageDESeq2_analysis.RData")
cat("\n", "Analysis finished with success! ", "\n")
