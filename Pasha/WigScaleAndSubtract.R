##############################################################
# This script is used to Scale Wig files and subtract their
# corresponding input WIG
##############################################################

library(Pasha)

library(Rargs)

##############################################################
### Parameters from command line using RIO
##############################################################
# Parameters format definition
paramsDefinition=list()

# mandatory params
paramsDefinition[["--wigFileList"]]=list(variableName="wigFileRegexList", numeric=F, mandatory=T, description="Space separated list of WIG file names (could be regular expressions).");
paramsDefinition[["--binSize"]]=list(variableName="binSize", numeric=T, mandatory=T, description="Size of the bins of the WIG files", default=50);

# non mandatory params
inputFileRegex = NULL
paramsDefinition[["--bgWigFile"]]=list(variableName="inputFileRegex", numeric=F, mandatory=F, description="Wig file of the input experiment", default=NULL);
paramsDefinition[["--rescaleInput"]]=list(variableName="rescaleInput", numeric=F, mandatory=F, description="Indicates if wether the input must be rescaled or not.", default=FALSE);
paramsDefinition[["--subtractInput"]]=list(variableName="subtractInput", numeric=F, mandatory=F, description="Indicates if whether the input must be subtracted from the WIG files or not.", default=TRUE);

# Retreives the parameters
getParams(paramsDefinition);

##############################################################
# Retrieve WIG files from regular expressions
##############################################################
cat( "\nRetrieving WIG files\n")
wigFileList = vector()
for( i in 1:length(wigFileRegexList))
{
	cat( "|-- Searching for file corresponding to ", wigFileRegexList[[i]],"\n")
	file_path = dirname(wigFileRegexList[[i]])
	file_list=list.files( path=file_path, pattern=basename(wigFileRegexList[[i]]))
	cat( "|---|-- File(s) found : \n")
	print( file_list)
	if( length(file_list)!=1){
		warning( paste( "|---|-- Could not locate unique file with regex :", wigFileRegexList[[i]]))
		stop()
	}
	file_full_path = paste( file_path, file_list[[1]],sep="/")
	wigFileList = append( wigFileList, file_full_path) 
}

cat("List of WIG files to manage :\n")
print( wigFileList)

inputFile = NULL
if( !is.null( inputFileRegex))
{
	cat( "\nRetrieving input WIG files\n")
	cat( "|-- Searching for input file corresponding to ", inputFileRegex, "\n")
	file_path = dirname(inputFileRegex)
	file_list=list.files( path=file_path, pattern=basename(inputFileRegex))
	cat( "|---|-- File(s) found : \n")
	print( file_list)
	if( length(file_list)!=1){
		warning( paste( "|---|-- Could not locate unique file with regex :", inputFileRegex))
		stop()
	}
	inputFile = paste( file_path, file_list[[1]],sep="/")
}


##############################################################
### Call the function to execute Normalization and subtraction
##############################################################
normAndSubtractWIG( wigFileList, inputFile=inputFile, rescaleInput=rescaleInput, meanBGLevelRescale=c(0,10), subtractInput=subtractInput, binSize=binSize)
