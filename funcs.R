euc_norm<-function(row){row/sqrt(sum(row*row))}

euc_norm_rows<-function(mat){t(apply(mat,1,euc_norm))}

cos_dist<-function(dframe){as.dist(crossprod(t(euc_norm_rows(dframe))))}