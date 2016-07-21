randomRows = function(x, n=10){ x[sample(nrow(x),n),]}

na.rm = function(x){ x[!is.na(x)]}

naIsFalse = function(x){x[is.na(x)]=F; return(x)}

corner = function(x, n=10){x[1:(min(c(nrow(x),n))),1:(min(c(ncol(x),n)))]}

combinePFisher = function(x){pchisq( -2*sum(log(x)), 2*length(x), lower.tail=FALSE)}

combineZStouffer = function(x){sum(x, na.rm=T)/sqrt(sum(!is.na(x)))}

clusterDataFrame = function(df, formula, value, method="euclidean", NA.set=NA, ...){ #formula must have only one y variable
  if (is.null(hcluster)){
    error("This function requires amap:hcluster.   library(amap) first.")
  }
  df = cast(df, formula = formula, value=value);
  dfIDs = as.character(df[[1]]);
  df[[1]]=NULL
  df[is.na(df)]=NA.set;
  myClust=hcluster(as.matrix(df), method=method, ...) # get clustering
  myClust$order[myClust$order<0] = (1:nrow(df))[!1:nrow(df) %in% myClust$order]
  return(dfIDs[myClust$order]);
}

ranksumROC = function(x,y,na.rm=T,...){
  if (na.rm){
    x=na.rm(x);
    y=na.rm(y);
  }
  curTest = wilcox.test(x,y,...);
  curTest$AUROC = curTest$statistic/(length(x)*length(y))
  return(curTest)
}
