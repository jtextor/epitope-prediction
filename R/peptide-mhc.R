#' @importFrom methods getPackageName
#' @importFrom grDevices rainbow
#' @importFrom graphics axis barplot strheight text
#' @importFrom stats quantile var
#' @importFrom utils read.table tail
NULL

.smm.cache <- new.env()

#' Supported MHC Molecules and Peptide Lengths
#'
#' @param l integer giving the desired peptide length. If not given (default), then all
#' possible combinations of peptide lengths and MHC molecules are returned.
#' @return A data frame containing the combinations of supported MHCs and peptide lengths.
#' @examples
#' ## Which MHC molecules are supported?
#' levels( supportedMHCs()$mhc )
#' ## Which peptide lengths are supported for HLA-A02:01?
#' with( supportedMHCs(), l[mhc=="HLA-A-02:01"] )
#'
#' @export
supportedMHCs <- function( l=NULL ){
	print( getPackageName() )
	r <- read.table( system.file( "extdata", "model_list.txt", package=getPackageName() ),
		as.is=TRUE )$V1
	mhcname <- gsub( '-[0-9][0-9]?$', '', r )
	mhcname <- gsub( '-([0-9][0-9])([0-9][0-9])', '-\\1:\\2', mhcname )
	peptidelength <- as.integer(gsub( ".*-([0-9][0-9]?)$", "\\1", r ))
	if( !is.null(l) && is.finite(l) ){
		mhcname <- mhcname[peptidelength==l]
		return(data.frame(mhc=mhcname,l=rep(l,length(mhcname))))
	} else {
		data.frame( mhc=mhcname, l=peptidelength )
	}
}

#' Get Prediction Matrix
#'
#' Loads the prediction matrix for a given MHC molecule at the given peptide length.
#' Matrices are loaded from the package folder upon first load and are then cached in
#' memory for further uses.
#' The SMM method predicts binding by summing the matrix entries for each amino acid
#' per position and then adding and MHC-dependant number. The resulting number is the
#' log IC50 value (IC50=half maximal inhibitory concentration). Hence, low numbers mean
#' strong binding and high numbers mean weak binding.
#'
#' @param mhc name of the MHC molecule.
#' @param l the peptide length.
#' @return A list with components \code{M} (the matrix) and \code{c} (a number to be
#' added to each prediction made with this matrix).
#' @examples
#' ## load prediction matrix for HLA-A-02:01 9-mers
#' M <- smmMatrix( "HLA-A-02:01", 9 )$M
#' ## How well does Leucine bind at each position?
#' M["L",]
#' ## Which amino acid is most preferred at the C-terminal position?
#' names( which.min( M[,9] ) )
#'
#' ## Naming of the MHC molecule is permissive
#' for( mhc in c("HLA-A0201","HLA-A-0201","HLA-A*0201","HLA-A*02:01") ){
#'    names( which.min( smmMatrix( "HLA-A-02:01", 9 )$M[,9] ) )
#' }
#'
#' @export
smmMatrix <- function( mhc="HLA-A-02:01", l=9 ){
	mhc <- gsub( "HLA-([AB])-([0-9])", "HLA-\\1\\2", as.character(mhc) )
	matrixname <- paste(gsub("[:*]","",mhc),"-",l,sep="")
	if( exists( matrixname, .smm.cache ) ){
		return(
			list( M=get( matrixname, envir=.smm.cache ),
				c=get( paste(matrixname,"-c",sep=""), envir=.smm.cache ) )
		)
	}
	matrixfile <- system.file( "extdata", paste(matrixname,".txt",sep=""),
		package=getPackageName() )
	if( matrixfile == "" ){
		stop( paste( "SMM matrix for MHC",mhc,"and peptide length",l,"not found!") )
	} else {
		M <- as.matrix( read.table(matrixfile, skip=1, nrows=20, row.names=1) )
		colnames(M) <- NULL
		Mc <- as.numeric( read.table(matrixfile, skip=21) )
		assign( paste(matrixname,"-c",sep=""), Mc,
			.smm.cache )
		assign( matrixname, M, .smm.cache )
		return( list(M=M, c=Mc) )
	}
}

#' Find MHC Binders in Protein Sequence
#'
#' Performs binding predictions for all l-mers in the given protein sequence, extracts
#' the peptide whose IC50 values are below the given thresholds, and returns the results
#' as a data frame.
#'
#' @param x string, a protein sequence given in single-letter coding. Only the 20 common
#' amino acids are supported.
#' @param mhc string identifying the MHC molecule.
#' @param l the peptide length.
#' @param ic50.threshold peptides with a predicted IC50 value lower than this will be
#' considered binders. A threshold of 500 nM is common. Use \code{Inf} to show predictions
#' for all peptides of length \code{l}.
#' @param quantile.threshold a number between 0 and 1. If this is not \code{NULL}, then
#' the parameter \code{ic50.threshold} is ignored and the peptides for which the predicted
#' IC50 falls within the given quantile will be returned. For instance, a value of .02
#' return the peptides whose binding strength is in the top 2\% for the given protein.
#' @param include.peptide logical, whether to include the actual peptide in the output
#' data frame. This may not be desired in some circumstances, e.g. for very long
#' proteins or for converting the result to a matrix.
#' @param method string defining which prediction method to use. Currently the only
#' implemented method is the stabilized matrix method (SMM), such that this setting is
#' ignored. But further methods may be implemented in the future.
#' @return a data frame containing the peptide (if \code{include.peptide=TRUE}),
#' start position, end position, and predicted IC50 for
#' every peptide below the threshold.
#'
#' @examples
#' ## This is the CORE protein from the Hepatitis C virus reference sequence available
#' ## at https://hcv.lanl.gov/content/sequence/LOCATE/locate.html
#' hcv.core <- paste("MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRR",
#'	"QPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGA",
#'  "PLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA",sep="")
#' binders( hcv.core )
#' @export
binders <- function( x, mhc="HLA-A-02:01", l=9, ic50.threshold=500,
	quantile.threshold=NULL,
	include.peptide=TRUE, method="smm" ){
	x <- as.character(x)
	start <- seq_len( nchar(x)-l+1 )
	end <- start+l-1
	ic50 <- smm( substring(x,start,end), mhc, l )
	if( is.null( quantile.threshold ) || !is.finite( quantile.threshold ) ||
		quantile.threshold < 0 || quantile.threshold > 1 ){
		i <- which( ic50 < ic50.threshold )
	} else {
		i <- which( ic50 <= quantile( ic50, quantile.threshold ) )
	}
	start <- start[i]
	end <- end[i]
	if( include.peptide ){
		data.frame( peptide=substring(x,start,end),
			start=start, end=end, ic50=ic50[i] )
	} else {
		data.frame( start=start, end=end, ic50=ic50[i] )
	}
}

#' Peptide-MHC Binding Prediction
#'
#' Predicts peptide-MHC binding using the stabilized matrix method (SMM)
#' algorithm with a specifically constructed amino acid substitution matrix.
#' See \code{citation(package='EpitopePrediction')} for the reference.
#'
#' @param x vector of strings containing the peptides for which to predict MHC
#' binding.
#' @param mhc string or vector of strings identifying the MHC molecules. See
#' \code{supportedMHCs} for allowed values. If a vector is given, it must be of the same
#' length as \code{x}.
#' @param output.IC50 whether to output the IC50 value itself (default) or its
#' base10-logarithm.
#' @examples
#' ## Predict IC50 binding values for two famous peptides
#' smm( c("SLYNTVATL","SYFPEITHI"), "HLA-A-02:01" )
#' @export
smm <- function( x=c("SLYNTVATL","SYFPEITHI"), mhc="HLA-A-02:01",
	output.IC50=TRUE ){
	if( length(mhc)>1 && ( length(mhc) != length(x) ) ){
		stop( "If 'mhc' is a vector, it must have the same length as 'peptides'!" )
	}
	if( length(mhc) == 1 && length(x) > 1 ){
		mhc <- rep.int( mhc, length(x) )
	}
	pred <- function( i ){
		l <-  nchar(x[i])
		M <- smmMatrix( mhc[i], l )
		M$c + sum( sapply( seq_len(l), function(j) M$M[substring(x[i],j,j),j] ) )
	}
	v <- sapply( seq_along(x), pred )
	if( output.IC50 ){
		return( 10^v )
	} else {
		return( v )
	}
}

#' Determine Anchor Positions
#'
#' Uses the information in the SMM matrices to predict anchor positions for a given MHC
#' molecule. Specifically, SMM matrix columns are sorted by their variance, and the
#' \code{k} columns with lowest variance are defined as anchor positions.
#'
#' @param mhc the mhc molecules.
#' @param l the peptide length.
#' @param k how many anchor positions to return.
#' @examples
#' ## Investigate anchor positions for the HLA molecules A02:01 and B27:05
#' anchorPositions( "HLA-A-02:01" )
#' anchorPositions( "HLA-A-02:01", k=3 )
#' anchorPositions( "HLA-B-27:05" )
#' @export
anchorPositions <- function( mhc="HLA-A-02:01", l=9, k=2 ){
	M <- smmMatrix( mhc, l )$M
	return(sort(tail(order(apply( M, 2, var )),k)))
}

#' Plot an MHC Binding Motif
#'
#' Visualizes the set of peptides which are predicted to bind to a certain MHC molecule.
#'
#' @param mhc name of the MHC molecule.
#' @param l length of the peptide.
#' @param motif.matrix a matrix. If this is not \code{NULL}, then the parameters
#' \code{mhc} and \code{l} are ignored and the supplied matrix is used instead. This
#' makes it possible to use this function for plotting sequence logos unrelated to
#' MHC-peptide binding. The supplied matrix must have row names to indicate the
#' letters out of which the sequence logo is made.
#' @param width bar width for profile per position.
#' @param space amount of space left before each bar, as a fraction of the
#' bar width.
#' @param col vector of colors to use for the amino acids (given in alphabetical order
#' in single-letter coding).
#' @param main title for the plot.
#' @param ... further options to be passed on to \code{\link[graphics]{barplot}}
#' @examples
#' ## Compare binding motifs of HLA-A02 at all supported peptide lengths
#' par( mfrow=c(2,2) )
#' for( l in 8:11 ){
#'    plotBindingMotif( "HLA-A-02:01", l )
#' }
#' @export
plotBindingMotif <- function( mhc="HLA-A-02:01",
		l=9, motif.matrix=NULL, width=.5,
		space=1, col=rainbow(20), main=paste(mhc,", ",l,"-mers",sep=""),
		... ){
	if( is.null( motif.matrix ) ){
		M <- 10^(-smmMatrix( mhc, l )$M)
	} else {
		M <- motif.matrix
		l <- ncol(M)
	}
	M <- scale(M,center=FALSE,scale=colSums(M))
	colEntropies <- apply(M, 2, function(x){
		x <- x*log(x)
		x[is.nan(x)] <- 0
		log(20) + sum(x)
  })
	M <- scale(M,center=FALSE,scale=1/colEntropies)
	barplot(M,col=col,width=width,space=space,main=main,...)
	for( i in seq_len( l ) ){
		ypos <- cumsum(M[,i])-M[,i]/2
		showl <- M[,i]>1.5*strheight("M")
		if( sum(showl) > 0 ){
			text( i*width*(1+space)-width/2, ypos[showl], rownames(M)[showl] )
		}
	}
	axis(1,at=(1:l)*(width*(1+space))-width/2,labels=1:l)
}
