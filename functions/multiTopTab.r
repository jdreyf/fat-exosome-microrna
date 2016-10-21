##jmd
##12.26.13

multiTopTab <- function(fit, cols=c('P.Value', 'adj.P.Val', 'logFC'), adjust.method='BH'){
  contrasts <- gsub(' ', '', colnames(fit$contrasts))
  #get gene order
  #limma 3.16 has row.names=row number & 'ID' column; limma 3.18 has row.names=ID
  ttf <- topTableF(fit, number=Inf)
  if ('ID' %in% colnames(ttf)){ genes <- ttf$ID } else { genes <- rownames(ttf) }
	#go thru contrasts
	for (i in 1:length(contrasts)){
	  mtt.tmp <- my.toptab(fit, coef=i, cols=cols, adjust.method=adjust.method)
	  mtt.tmp <- mtt.tmp[genes,]
	  colnames(mtt.tmp) <- paste(contrasts[i], colnames(mtt.tmp), sep='.')
    if (i==1){
      mtt <- mtt.tmp
    } else {
      mtt=cbind(mtt, mtt.tmp)
    }
	}
	return(mtt)
}