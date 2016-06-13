randomRows = function(x, n=10){ x[sample(nrow(x),n),]}

na.rm = function(x){ x[!is.na(x)]}
naIsFalse = function(x){x[is.na(x)]=F; return(x)}

corner = function(x, n=10){x[1:(min(c(nrow(x),n))),1:(min(c(ncol(x),n)))]}
