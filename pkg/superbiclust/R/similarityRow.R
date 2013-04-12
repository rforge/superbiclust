similarity <-
		function (x, index="jaccard")
{
	options <- c("jaccard", "sorensen", "ochiai","kulczynski","sensitivity","specificity")
	if(any(!index %in% options)) {
		stop("index must be one of jaccard, sorensen, ochiai,kulczynski,sensitivity, or specificity")
	}		
	output <- switch(index,
			jaccard = jaccardMat(x, x), 
			sorensen = sorensenMat(x,x), 
			ochiai = ochiaiMat(x,x),
			kulczynski = kulczynskiMat(x,x),
			sensitivity = sensitivityMat(x,x),
			specificity = specificityMat(x,x))
	
	class(output)<- "similarity"
	output
}

jaccardMat <- function(x,y){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	jamat <- matrix(0,pA,pB)	
	
	for (i in 1:pA) {
		biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
		apos <- biclA > 0
		sizeBiclA <- sum(apos) #size of bicluster A
		for (j in 1:pB) {
			biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
			bpos <- biclB > 0
			sizeBiclB <- sum(bpos)
			
			biclAB <- biclA + biclB
			abpos <- biclAB > 0
			sab <- sum(abpos)
			
			if (sab>0) {
				jamat[i,j] <- (sizeBiclA + sizeBiclB)/sab - 1
			} else {
				jamat[i,j] <- 0
			}
			
		}
	}
	return(jamat)
}

sorensenMat <- function(x,y){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	somat <- matrix(0,pA,pB)	
	
	for (i in 1:pA) {
		biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
		apos <- biclA > 0
		sizeBiclA <- sum(apos) #size of bicluster A
		for (j in 1:pB) {
			biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
			bpos <- biclB > 0
			sizeBiclB <- sum(bpos)
			
			biclAB <- biclA + biclB
			abpos <- biclAB > 0
			sab <- sum(abpos)
			
			if (sab>0) {
				somat[i,j] <- 2.0-2.0*sab/(sizeBiclA+sizeBiclB)
			} else {
				somat[i,j] <- 0
			}
			
		}
	}
	return(somat)	
}

kulczynskiMat <- function(x,y){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	kumat <- matrix(0,pA,pB)	
	
	for (i in 1:pA) {
		biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
		apos <- biclA > 0
		sizeBiclA <- sum(apos) #size of bicluster A
		for (j in 1:pB) {
			biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
			bpos <- biclB > 0
			sizeBiclB <- sum(bpos)
			
			biclAB <- biclA + biclB
			abpos <- biclAB > 0
			sab <- sum(abpos)
			
			if ((sizeBiclA>0)&&(sizeBiclB>0))
			{
				kumat[i,j] <- 1.0+0.5*( (sizeBiclA-sab)/sizeBiclB + (sizeBiclB -sab)/sizeBiclA )				
			}else {
				kumat[i,j] <- 0
			}			
		}
	}
	return(kumat)	
}

ochiaiMat <- function(x,y){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	ocmat <- matrix(0,pA,pB)
	
	for (i in 1:pA) {
		biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
		apos <- biclA > 0
		sizeBiclA <- sum(apos) #size of bicluster A
		for (j in 1:pB) {
			biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
			bpos <- biclB > 0
			sizeBiclB <- sum(bpos)
			
			biclAB <- biclA + biclB
			abpos <- biclAB > 0
			sizeBiclAB <- sum(abpos)
			
			if ((sizeBiclA>0)&&(sizeBiclB>0))
			{
				ocmat[i,j] <- (sizeBiclA+sizeBiclB-sizeBiclAB)/sqrt(sizeBiclB*sizeBiclA)				
			}else {
				ocmat[i,j] <- 0
			}			
		}
	}
	return(ocmat)	
}

sensitivityMat <- function(x,y){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	senmat <- matrix(0,pA,pB)
	
	for (i in 1:pA) {
		biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
		apos <- biclA > 0
		sizeBiclA <- sum(apos) #size of bicluster A
		for (j in 1:pB) {
			biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
			bpos <- biclB > 0
			sizeBiclB <- sum(bpos)
			
			biclAB <- biclA + biclB
			abpos <- biclAB > 0
			sizeBiclAB <- sum(abpos)
			
			if ((sizeBiclA>0)&&(sizeBiclB>0))
			{
				senmat[i,j] <- 1+ (sizeBiclB-sizeBiclAB)/sizeBiclA
			}else {
				senmat[i,j] <- 0
			}			
		}
	}
	return(senmat)	
}

specificityMat <- function(x,y){
	pA = x@Number
	pB = y@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	bicBrows = y@GenesMembership
	bicBcols = y@ColumnMembership
	
	n <- length(bicArows[1,]) #number of rows in data 
	l <- length(bicAcols[1,]) #number of columns in data
	u <- n*l
	spemat <- matrix(0,pA,pB)
	
	for (i in 1:pA) {
		biclA <- tcrossprod(bicArows[i,],bicAcols[i,])
		apos <- biclA > 0
		sizeBiclA <- sum(apos) #size of bicluster A
		for (j in 1:pB) {
			biclB <- tcrossprod(bicBrows[,j],bicBcols[j,])
			bpos <- biclB > 0
			sizeBiclB <- sum(bpos)
			
			biclAB <- biclA + biclB
			abpos <- biclAB > 0
			sizeBiclAB <- sum(abpos)
			
			if ((sizeBiclA>0)&&(sizeBiclB>0))
			{
				spemat[i,j] <- 1+ (sizeBiclA-sizeBiclAB)/sizeBiclB
			}else {
				spemat[i,j] <- 0
			}			
		}
	}
	return(spemat)	
}