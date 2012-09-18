similarity <-
function (x, index="jaccard")
{
	options <- c("jaccard", "sorensen", "ochiai","kulczynski","sensitivity","specificity")
	if(any(!index %in% options)) {
		stop("index must be one of jaccard, sorensen, ochiai,kulczynski,sensitivity, or specificity")
	}
#    require(truecluster)
	require(Matrix)
	numA = x@Number
	bicArows = t(x@GenesMembership)
	bicAcols =x@ColumnMembership
	
	numB = x@Number
	bicBrows =x@GenesMembership
	bicBcols =x@ColumnMembership
	
    pA <- numA
    pB <- numB
	options <- c("jaccard", "sorensen", "ochiai","kulczynski","sensitivity","specificity")
        n <- length(bicArows[1,])
        l <- length(bicAcols[1,])
        u <- n*l
        jamat <- matrix(0,pA,pB)
        kumat <- matrix(0,pA,pB)
        ocmat <- matrix(0,pA,pB)
        somat <- matrix(0,pA,pB)
	senmat <- matrix(0,pA,pB)
	spemat <- matrix(0,pA,pB)
	
        for (i in 1:pA) {
            bcA <- tcrossprod(bicArows[i,],bicAcols[i,])
            apos <- bcA > 0
            sa <- sum(apos)
            if (sa > 0.5*u) {
                bcA[apos] <- 0
                 sa <- 0
            }
            for (j in 1:pB) {
                bcB <- tcrossprod(bicBrows[,j],bicBcols[j,])
                bpos <- bcB > 0
                sb <- sum(bpos)
                if (sb > 0.5*u) {
                    bcB[bpos] <- 0
                   sb <- 0
                }
                bcAB <- bcA + bcB
                abpos <- bcAB > 0
                sab <- sum(abpos)

                if (sab>0) {
                    jamat[i,j] <- (sa + sb)/sab - 1
		    somat[i,j] <- 2.0-2.0*sab/(sa+sb)
                } else {
                    jamat[i,j] <- 0
		    somat[i,j] <- 0
                }
		 if ((sa>0)&&(sb>0))
                {
                    kumat[i,j] <- 1.0+0.5*( (sa-sab)/sb + (sb -sab)/sa )
                    ocmat[i,j] <- (sa+sb-sab)/sqrt(sb*sa)
		    senmat[i,j] <- 1+ (sb-sab)/sa
		    spemat[i,j] <- 1+ (sa-sab)/sb
		
                }else {
                    kumat[i,j] <- 0
                    ocmat[i,j] <- 0
		    senmat[i,j] <- 0
		    spemat[i,j] <- 0
                }
		
            }
    }

	
if(index==options[1]){
	output <- jamat
} else if (index==options[2]){
	output <- somat	
} else if (index==options[3]){
	output <-ocmat
}else if (index==options[4]){
	output <- kumat
} else if (index==options[5]){
	output <- senmat
} else if (index==options[6]){
	output <- spemat
}
class(output)<- "similarity"
output
}

