
#' Calculate TMM normalization factors
#' @param  x: a matrix of TPM values, for which normalization factors are calculated
#' @param  normFactors: a vector of normalization factors to divide by if provided, otherwise runs calcTMMFactors
#' @param  ...: other parameters to be passed to calcTMMFactors
#' @return the TMM normalized data
#' @examples
#' allTMM = normalizeTMM(allTPM, verbose=0);
normalizeTMM = function(x, normFactors=NULL, isLog=F, ...){
  if (is.null(normFactors)){
    normFactors = calcTMMFactors(x, isLog, ...);
  }
  for (c in 1: ncol(x)){
    if (isLog){
      x[,c] = x[,c] - normFactors[c];
    }else{
      x[,c] = x[,c]/normFactors[c];
    }
  }
  return(x);
}

#' Calculate TMM normalization factors
#' @param  x: a matrix of TPM values, for which normalization factors are calculated
#' @param  trimByExpr: the fraction of data points to remove from the lowly expressed and highly expressed (respectively) genes
#' @param  trimByRatio: the fraction of data points to remove from either side of the log ratios of the cell/reference
#' @param  referenceSample: the vector of gene expression values that are used as a reference sample. Defaults to the median of non-0 expression values
#' @param  verbose: verbosity of messaging; will report every nth cell for verbose>0
#' @param  printGraphs: will print graphs at each stage of normalization
#' @param  isLog: data is already log-scale
#' @return the normalization factors that each sample should be divided by
#' @examples
#' normFactors = calcTMMFactors(tpmMatrix);
calcTMMFactors = function(x, trimByExpr = c(0.1,0.1), trimByRatio = c(0.25,0.25), referenceSample=NULL, verbose=0, printGraphs=F, isLog=F){
  #x is a matrix of TPMs (genes x samples)
  #get the median expression of each gene
  if (is.null(referenceSample)){ # take the reference sample as the median expression values in the data (TPM usually)
    referenceSample = apply(x,1,function(y){median(y, na.rm=T)});
    referenceSample[is.na(referenceSample)]=0;
  }
  #for each sample, take only those genes that 1) are expressed 2) are somewhere in the "middle" of expression values and 3) have intermediate ratios
  referenceBounds = quantile(referenceSample, c(trimByExpr[1], 1-trimByExpr[2]), na.rm=T); #trim extremes of expression
  keepReference = referenceSample>=referenceBounds[1] & referenceSample<=referenceBounds[2];
  keepReference[is.na(keepReference)]=F;
  if (printGraphs){
    if (isLog){
      p = ggplot(data.frame(refExpr = c(referenceSample,referenceSample[keepReference]), class = c(rep("all",length(keepReference)), rep("used",sum(keepReference)))),aes(x=refExpr, fill=class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = referenceBounds)+ggtitle("reference"); print(p)
    }else{
      p = ggplot(data.frame(refExpr = c(referenceSample,referenceSample[keepReference]), class = c(rep("all",length(keepReference)), rep("used",sum(keepReference)))),aes(x=log2(refExpr+1), fill=class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = log2(1+referenceBounds))+ggtitle("reference"); print(p)
    }
    cat ("Press [enter] to continue")
    line <- readline()
  }
  if (verbose>0){
    message(sprintf("Keeping %f%% of reference expression values.",100*mean(keepReference)))
  }
  normFactors = rep(1, ncol(x));
  for (c in 1: ncol(x)){
    curExprBounds = quantile(x[,c], c(trimByExpr[1], 1-trimByExpr[2]), na.rm=T); #define extremes of expression
    keepCurSample = x[,c]>=curExprBounds[1] & x[,c]<=curExprBounds[2]; 
    keepCurSample[is.na(keepCurSample)]=F;
    if (printGraphs){
      if (isLog){
        p = ggplot(data.frame(cellExpr = c(x[,c],x[keepCurSample,c]), class = c(rep("all",length(keepCurSample)), rep("used",sum(keepCurSample)))),aes(x=cellExpr, fill = class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = curExprBounds)+ggtitle(colnames(x)[c]); print(p)
      }else{
        p = ggplot(data.frame(cellExpr = c(x[,c],x[keepCurSample,c]), class = c(rep("all",length(keepCurSample)), rep("used",sum(keepCurSample)))),aes(x=log2(cellExpr+1), fill = class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = log2(1+curExprBounds))+ggtitle(colnames(x)[c]); print(p)
      }
      cat ("Press [enter] to continue")
      line <- readline()
    }
    if (isLog){
      curLogRatios = x[,c] - referenceSample; # get current log ratios for all datapoints
    }else{
      curLogRatios = log2(x[,c]/referenceSample); # get current log ratios for all datapoints
    }
    curRatioBounds = quantile(curLogRatios[keepCurSample & keepReference], c(trimByRatio[1], 1-trimByRatio[2]), na.rm=T); #define extremes of fold change
    keepCurSample = keepCurSample & keepReference & curLogRatios>=curRatioBounds[1] & curLogRatios<=curRatioBounds[2];
    keepCurSample[is.na(keepCurSample)]=F;
    curOffset = mean(curLogRatios[keepCurSample], na.rm=T) # calculate the offset of this sample with the outlier-free data
    if (printGraphs){
      p = ggplot(data.frame(cellOverRef = c(curLogRatios,curLogRatios[keepCurSample]),class=c(rep("all",length(curLogRatios)), rep("used",sum(keepCurSample)))),aes(x=cellOverRef,fill=class))  +geom_histogram(alpha=0.25, position="identity", bins=100) + theme_bw()+geom_vline(xintercept = curRatioBounds)+ggtitle(colnames(x)[c]); print(p)
      cat ("Press [enter] to continue")
      line <- readline()
    }
    if (verbose>0 && c %% verbose ==0){
      message(sprintf("Cell %i is being centered by %i genes and is offset from the reference by %g (log2)",c,sum(keepCurSample), curOffset))
    }
    if (sum(keepCurSample)<50){
      message(sprintf("WARNING: fewer than 50 genes used to center sample %i (#genes = %i)",c,sum(keepCurSample)));
    }
    if (sum(keepCurSample)>0){
      if (isLog){
        normFactors[c] = (curOffset);
      }else{
        normFactors[c] = (2^curOffset);
      }
    }
  }
  return(normFactors);
}
