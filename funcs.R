norm<-function(data){ t(t(data)/rowSums(t(data))) }

euc_norm<-function(row){ row/sqrt(sum(row*row)) }

euc_norm_rows<-function(mat){ t(apply(mat,1,euc_norm)) }

cos_dist<-function(dframe){ as.dist(acos(crossprod(t(euc_norm_rows(norm(dframe)))))) }

get_monkey<-function(ensembl_ls, outfile){ library('biomaRt'); mart<-useMart('ensembl'); mart<-useMart('ensembl',dataset='mmulatta_gene_ensembl'); seq<-getSequence(id=ensembl_ls, type='ensembl_peptide_id', seqType='peptide', mart=mart); exportFASTA(seq, outfile)}

go_bp_5<-function(ls,outfile){ library('biomaRt'); mart<-useMart('ensembl'); mart<-useMart('ensembl',dataset='hsapiens_gene_ensembl'); bm<-getBM(attributes=c('entrezgene','go_biological_process_id'), filter='entrezgene',values=ls,mart=mart); write.table(bm,file=outfile,sep='\t',append=FALSE,col.names=FALSE,row.names=FALSE,quote=FALSE)}