##jmd
##5.27.14

library('limma')

#sort by p
#assume that if 'logFC' in cols, then want 'FC'
my.toptab <- function(fit, cols=c('P.Value', 'adj.P.Val', 'logFC'), adjust.method='BH', prefix='', coef=NULL){
  tt <- topTable(fit, number=Inf, sort.by='P', adjust.method=adjust.method, coef=coef)
  #limma 3.16 has row.names=row number & 'ID' column; limma 3.18 has row.names=ID
  if ('ID' %in% colnames(tt)){ rownames(tt) <- tt$ID }
  #FC
  if ('logFC' %in% cols){
    tt$FC <- logfc2fc(tt$logFC)
    cols <- c(cols, 'FC')
  }    
  tt <- tt[,cols]
  colnames(tt) <- sub('P.Value', 'p', sub('adj.P.Val', 'FDR', colnames(tt)))
  if (prefix!=''){ colnames(tt) <- paste(prefix, colnames(tt), sep='.') }
  return(tt)
}