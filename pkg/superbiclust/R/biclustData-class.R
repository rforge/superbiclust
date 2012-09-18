setClass("BiclustSet",
		representation(GenesMembership="matrix",ColumnMembership="matrix",
				Number="numeric")
)
#constructor for the class BiclustSet
BiclustSet <-function(x){myBiclustSet = new("BiclustSet",GenesMembership=x@RowxNumber, ColumnMembership=x@NumberxCol, 
			Number=x@Number)}
		
setGeneric("BiclustSet")
setMethod("BiclustSet",signature(x="Biclust"),
		function(x){myBiclustSet = new("BiclustSet",GenesMembership=x@RowxNumber, ColumnMembership=x@NumberxCol, 
				Number=x@Number)})
setMethod("BiclustSet",signature("Factorization"),function(x){
	require(fabia)
	N <-  x@p1
	tmp <- extractBic(x)
	X <- x@X
	FabiaBicl <- tmp$numn
	DataFabiaColB <- c()
	DataFabiaRowB <- c()
	for (j in 1:N){
		colMem <- vector(length= ncol(X), mode="logical")
		rowMem <- vector (length = nrow(X), mode = "logical")
		
		if(length(unlist(FabiaBicl[j,2])) > 0){
			for (i in unlist(FabiaBicl[j,2])) {colMem[i] <- TRUE} 
			for (i in unlist(FabiaBicl[j,1])) {rowMem[i] <- TRUE}
			DataFabiaColB <- rbind(DataFabiaColB, colMem)
			DataFabiaRowB <- cbind(DataFabiaRowB,  rowMem)
		}else next
	}
	myBiclustSet = new("BiclustSet",GenesMembership=DataFabiaRowB , ColumnMembership=DataFabiaColB, 
			Number=N)
} )

setMethod("BiclustSet",signature("list"),
		function(x){
			biclust1Number<- ncol(x$rows)
			bicArows <- matrix(as.logical(x$rows),ncol=biclust1Number,byrow=T)
			bicAcols <- matrix(as.logical(x$columns),nrow=biclust1Number)
			myBiclustSet = new("BiclustSet",GenesMembership=bicArows, ColumnMembership=bicAcols, 
					Number=biclust1Number)
		})
