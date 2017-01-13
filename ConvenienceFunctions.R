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
  curTest$AUROC = (curTest$statistic/length(x))/length(y)
  return(curTest)
}

inputDataFromFilesAsMeltedTable = function(d, fileCols, formatStr="%s", ...){
  i=1;
  curData = read.table(do.call(sprintf, as.list(c(formatStr,d[i,fileCols]))), ...);
  
  inputData = data.frame(i = 1:(nrow(curData) * nrow(d)));
  for (c in names(curData)){ # add col names from curData
    inputData[c] = NA
  }
  if (all(inputData$i == 1:(nrow(curData) * nrow(d)))){
    inputData$i=NULL;
  }
  nrow(inputData)
  for (c in names(d)){ # add col names from d (list of inputs)
    inputData[c]=d[i,c];
  }
  
  z=1;
  for (i in 1:nrow(d)){
    message(i/nrow(d));
    curData = read.table(do.call(sprintf, as.list(c(formatStr,d[i,fileCols]))), ...);
    for (c in names(d)){
      curData[c]=d[i,c];
    }
    
    if ((nrow(curData)+z-1) > nrow(inputData)){
      inputData = inputData[1:(z-1),];
      inputData = rbind(inputData,curData)
    }else{
      inputData[z:(nrow(curData)+z-1),] = curData;
    }
    z=z+nrow(curData);
  }
  inputData = inputData[1:(z-1),];
  return(inputData);
}
