#!/usr/bin/env Rscript
#
#Author: Boulanger Mathias
#Derived from JCAL team (IGMM Montpellier) original script 
#Last update 23-10-2019
#
#This script is used to merge several WIGs together
#REQUIREMENT: optparse and Pasha   
#USAGE: Rscript MergeWigs.R -m path/sampleMatrix -s path/Spike-inTable -d parametric -a 0.1
#

###packages
# check packages
# packages <- c("optparse", "Pasha")
# if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
# 	install.packages(setdiff(packages, rownames(installed.packages())))
# }

#Library loading
suppressPackageStartupMessages(library("optparse"))
#suppressPackageStartupMessages(library("Pasha"))

###Creation of script option and the working environment
#Make option list
option_list = list(
  make_option(c("-f", "--filesToMerge"), type="character", default=NULL, help="Comma separated list of WIG file paths to merge [Example: path_to_my/wigFile1.wig,path_to_my/wigFile2.wig].", metavar="character"),
  make_option(c("-b", "--binSize"), type="integer", default=50, help="Size of the bins of the WIG files", metavar="numeric"),
  make_option(c("-n", "--outputFileName"), type="character", default="out_MergeWig.wig", help="Name of the output file [default= %default]", metavar="character"),
  make_option(c("-o", "--outputFolder"), type="character", default="./", help="Path to the folder were the merge wig will be put [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Check the ability to work
if (is.null(opt$filesToMerge)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (-f file1,file2).\n", call.=FALSE)
}
if (class(opt$binSize) != "numeric" ) {
	print_help(opt_parser)
	stop("The argument in --binSize sould be numeric (-b 20).\n", call.=FALSE)
}


if (!file.exists(opt$matrixfile)) {
	print_help(opt_parser)
	stop("The input file does not exist!\n", call.=FALSE)
}
if (!is.null(opt$spikeinmatrix) && !file.exists(opt$matrixfile)) {
	print_help(opt_parser)
	stop("The spike-in input file does not exist!\n", call.=FALSE)
}
if (opt$dispersion != "parametric" && opt$dispersion != "local") {
	print_help(opt_parser)
	stop("The argument in --dispersion is not good (-d parametric or -d local).\n", call.=FALSE)
}

if ( > 1 | opt$alpha < 0 ) {
	print_help(opt_parser)
	stop("The number in --alpha is not between 0 and 1 (-a 0.1).\n", call.=FALSE)
}

#Start to work
cat( "\nCreating WIG file list\n")
wigFileList = vector(unlist(strsplit(opt$filesToMerge, ",")))
if (length(wigFileList) < 2){
	print_help(opt_parser)
	stop("At least two files must be supplied (-f file1,file2).\n", call.=FALSE)
}
for (i in 0:length(wigFileList)) {
	if (!file.exists(path=wigFileList[i])) {
		print_help(opt_parser)
		stop(paste0("The input file (", wigFileList[i], ") does not exist!\n"), call.=FALSE)
	}
}

cat("----- Setup parameters -----\nList of WIG files to merge:\n")
print(wigFileList)
cat("\nBin size of ", opt$binSize, ".\n")

# Test outputFolder finish by a "/"
if( !grepl(opt$outputFolder, ".*/$")) {
	outputFolder = paste( outputFolder, "/",sep="")
} else {
	outputFolder = opt$outputFolder
}
cat("Output folder = ", outputFolder, "\n", "Output file name = ", opt$outputFileName)

#Test if the file and the folder already exist.
if (dir.exists(path=outputFolder) == TRUE) {
	if (file.exists(path=paste(outputFolder, opt$outputFileName))) {
		cat("WARNING: The output file already present in he outup folder has been overwritten!\n")
	}
}

#Merge calling
files=list()
files[[opt$outputFileName]]= wigFileList
mergeWigs(files, binSize=opt$binSize, outputFolder=outputFolder)
