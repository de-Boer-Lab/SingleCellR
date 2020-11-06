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


#Compares the original eigen values to permuted matrix values (done ntimes) for each PC, each time removing the last nPCs significant PCs, stopping when one is not significant. 
findSigPCs = function(myPCA, nfolds=10, alpha=0.01){
  permPCAStats = data.frame();
  stillSig=T;
  nPCs = 0
  rePCA = prcomp(t(myPCA$rotation) * matrix(rep(myPCA$sdev,times=nrow(myPCA$rotation)),nrow=ncol(myPCA$rotation), ncol=nrow(myPCA$rotation)),center = F, scale. = F);
  while(stillSig){
    curPCAStats = data.frame();
    for (i in 1:nfolds){
      remainder = ncol(myPCA$rotation) - nPCs;
      x1 = as.vector(t(myPCA$rotation[,(nPCs+1):ncol(myPCA$rotation)])* matrix(rep(myPCA$sdev[(nPCs+1):ncol(myPCA$rotation)],times=nrow(myPCA$rotation)),nrow=remainder,ncol=nrow(myPCA$rotation)))
      x1 = x1[sample(length(x1),length(x1))]
      x1 = matrix(x1,nrow=nrow(myPCA$rotation),ncol=remainder)
      curPCA = prcomp(t(x1),center = F, scale. = F);
      curPCAStats = rbind(curPCAStats, data.frame(sdev = curPCA$sdev, perm=i, eig = curPCA$sdev^2, ranks = 1:length(curPCA$sdev),nPCs=nPCs));
    }
    curPCAStats=curPCAStats[order(curPCAStats$sdev,decreasing = T)[1:nfolds],]
    message(sprintf("PC=%i; P=%g; sdev_act=%g, sdev_max_perm=%g",nPCs+1, mean(rePCA$sdev[nPCs+1] < curPCAStats$sdev),rePCA$sdev[nPCs+1], max(curPCAStats$sdev)))
    
    if (mean(rePCA$sdev[nPCs+1] > curPCAStats$sdev) < (1-alpha)){
      stillSig=F
    }else{
      nPCs = nPCs+1;
    }
    if(nPCs >= ncol(myPCA$rotation)){
      stillSig=F;
    }
    permPCAStats = rbind(permPCAStats, curPCAStats)
  }
  myPCA$permPCA = permPCAStats;
  myPCA$rePCA_sdev = rePCA$sdev;
  myPCA$nPCs = nPCs
  return(myPCA)
}

proportionalDownsample = function(x, nbins=100, desiredN=100000){
  minX = min(x);
  maxX = max(x);
  range = maxX-minX;
  window = range/nbins;
  expectedPerBin = desiredN/nbins;
  keepVector = rep(F, length(x))
  randVec = runif(length(x));
  for (i in 1:nbins){
    thesePoints = x > (minX + window*(i-1)) & x <= (minX + window*(i))
    numHere = sum(thesePoints)
    keepP = expectedPerBin/numHere; #prob of keeping
    keepVector = keepVector | (thesePoints & (randVec < keepP));
  }
  return (keepVector)
}

#allMotifStimCorrs = m2apply(log(eps+allProbeActivity[grepl("_A(_RC)?$", allProbeActivity$ID),names(allProbeActivity) %in% goodMotifIDs]), allProbeActivity[grepl("_A(_RC)?$", allProbeActivity$ID),grepl("_ratio",names(allProbeActivity))], FUN = function(x,y){xmin = as.numeric(quantile(x,probs=1-0.035)); cor(x[x>xmin],y[x>xmin], use="pairwise.complete.obs")})
m2apply = function(x,y,FUN, ...){
  results = matrix(nrow=ncol(x), ncol=ncol(y));
  row.names(results)=colnames(x);
  colnames(results)=colnames(y)
  for(xi in 1:ncol(x)){
    for(yi in 1:ncol(y)){
      results[xi,yi]=FUN(x[,xi],y[,yi], ...);
    }
  }
  return(results)
}

fixExcelGenesHuman = function(x){
  x = gsub("([0-9]+)-Mar","MARCH\\1", x)
  x = gsub("([0-9]+)-Sep","SEPT\\1", x)
  x = gsub("1-Dec","DELEC1", x)
}
