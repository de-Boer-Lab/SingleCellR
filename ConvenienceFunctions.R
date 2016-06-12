randomRows = function(x, n=10){ x[sample(n,nrow(x)),]}

na.rm = function(x){ x[!is.na(x)]}
naIsFalse = function(x){x[is.na(x)]=F; return(x)}

corner = function(x){t(head(t(head(x))))}
