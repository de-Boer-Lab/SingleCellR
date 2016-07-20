randomRows = function(x, n=10){ x[sample(nrow(x),n),]}

na.rm = function(x){ x[!is.na(x)]}

naIsFalse = function(x){x[is.na(x)]=F; return(x)}

corner = function(x, n=10){x[1:(min(c(nrow(x),n))),1:(min(c(ncol(x),n)))]}

combinePFisher = function(x){pchisq( -2*sum(log(x)), 2*length(x), lower.tail=FALSE)}


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
