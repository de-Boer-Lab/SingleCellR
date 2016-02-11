

#' Calculate TMM normalization factors
#' @param  x: a matrix of TPM values, for which normalization factors are calculated
#' @param  normFactors: a vector of normalization factors to divide by if provided, otherwise runs calcTMMFactors
#' @param  ...: other parameters to be passed to calcTMMFactors
#' @return the TMM normalized data
#' @examples
#' allTMM = normalizeTMM(allTPM, verbose=0);
normalizeTMM = function(x, normFactors=NULL, ...){
  if (is.null(normFactors)){
    normFactors = calcTMMFactors(x, ...);
  }
  for (c in 1: ncol(x)){
    x[,c] = x[,c]/normFactors[c];
  }
  return(x);
}

#' Calculate TMM normalization factors
#' @param  x: a matrix of TPM values, for which normalization factors are calculated
#' @param  trimByExpr: the fraction of data points to remove from the lowly expressed and highly expressed (respectively) genes
#' @param  trimByRatio: the fraction of data points to remove from either side of the log ratios of the cell/reference
#' @param  referenceSample: the vector of gene expression values that are used as a reference sample. Defaults to the median of non-0 expression values
#' @param  zeroExpression: the expression value is considered 0
#' @param  verbose: verbosity of messaging; will report every nth cell for verbose>0
#' @return the normalization factors that each sample should be divided by
#' @examples
#' normFactors = calcTMMFactors(tpmMatrix);
calcTMMFactors = function(x, trimByExpr = c(0.3,0.1), trimByRatio = c(0.25,0.25), referenceSample=NULL, zeroExpression=0, verbose=0){
  #x is a matrix of TPMs (genes x samples)
  #get the median expression of each gene
  if (is.null(referenceSample)){ # take the reference sample as the median of non-zero expression values in the data (TPM usually)
    referenceSample = apply(x,1,function(y){median(y[y>zeroExpression], na.rm=T)});
    referenceSample[is.na(referenceSample)]=0;
  }
  #for each cell, take only those genes that 1) are expressed 2) are somewhere in the "middle" of expression values and 3) have intermediate ratios
  keepReference = referenceSample>zeroExpression; # non-zero
  referenceBounds = quantile(referenceSample[keepReference], c(trimByExpr[1], 1-trimByExpr[2]), na.rm=T); #trim extremes of expression
  keepReference = keepReference & referenceSample>=referenceBounds[1] & referenceSample<=referenceBounds[2];
  
  #p = ggplot(data.frame(refExpr = c(referenceSample[referenceSample>zeroExpression],referenceSample[keepReference]),class=c(rep("all",sum(referenceSample>zeroExpression)), rep("used",sum(keepReference)))),aes(x=log2(refExpr+1),fill=class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = log2(1+referenceBounds)); print(p)
  
  if (verbose>0){
    message(sprintf("Keeping %f%% of reference expression values.",100*mean(keepReference)))
  }
  normFactors = rep(1, ncol(x));
  for (c in 1: ncol(x)){
    keepCurCell = x[,c]>zeroExpression; #non zero for this cell
    curExprBounds = quantile(x[keepCurCell,c], c(trimByExpr[1], 1-trimByExpr[2]), na.rm=T); #define extremes of expression
    keepCurCell = keepCurCell & x[,c]>=curExprBounds[1] & x[,c]<=curExprBounds[2]; 
    #p = ggplot(data.frame(cellExpr = c(x[x[,c]>zeroExpression,c],x[keepCurCell,c]),class=c(rep("all",sum(x[,c]>zeroExpression)), rep("used",sum(keepCurCell)))),aes(x=log2(cellExpr+1),fill=class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = log2(1+curExprBounds)); print(p)
    curLogRatios = log2(x[,c]/referenceSample); # get current log ratios for all datapoints
    curRatioBounds = quantile(curLogRatios[keepCurCell & keepReference], c(trimByRatio[1], 1-trimByRatio[2]), na.rm=T); #define extremes of fold change
    keepCurCell = keepCurCell & keepReference & curLogRatios>=curRatioBounds[1] & curLogRatios<=curRatioBounds[2];
    keepCurCell[is.na(keepCurCell)]=F;
    curOffset = mean(curLogRatios[keepCurCell], na.rm=T) # calculate the offset of this sample with the outlier-free data
    #p = ggplot(data.frame(cellOverRef = c(curLogRatios,curLogRatios[keepCurCell]),class=c(rep("all",length(curLogRatios)), rep("used",sum(keepCurCell)))),aes(x=cellOverRef,fill=class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = curRatioBounds); print(p)
    if (verbose>0 && c %% verbose ==0){
      message(sprintf("Cell %i is being centered by %i genes and is offset from the reference by %g (log2)",c,sum(keepCurCell), curOffset))
    }
    if (sum(keepCurCell)<50){
      message(sprintf("WARNING: fewer than 50 genes used to center sample %i (#genes = %i)",c,sum(keepCurCell)));
    }
    if (sum(keepCurCell)>0){
      normFactors[c] = (2^curOffset);
    }
  }
  return(normFactors);
}
