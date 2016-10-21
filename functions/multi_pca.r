##jmd
##4.10.14

library('ggplot2')

#pca for mult dimensions of pheno
#attributes: color, shape, size (default: 2), fill
#all.size gives size for all points w/o appearing in legend
#inc. facet as proof of principle
#qplot really hard to use -- can't send "col=time" or "col='time'" to ID column of dat
#ex use: multi.pca(mat, pheno, name=NA, color='time', shape='infusion', facet='.~tissue')
multi.pca <- function(mat, pheno, name='pca', scale.=FALSE, alpha=1, all.size=NULL, facet=NULL, ret.pcs=FALSE, 
                      rm.leg.title=FALSE, add.labels=FALSE, ...){
  stopifnot(ncol(mat)==nrow(pheno), colnames(mat)==rownames(pheno))
  
  pca <- prcomp(t(mat[rowSums(is.na(mat))==0,]), scale.=scale.)
  pve <- signif(summary(pca)$importance['Proportion of Variance', 1:2]*100, 2)
  
  dat <- data.frame(pca$x[rownames(pheno), 2:1], pheno)
  
  #qp <- qplot(x=PC1, y=PC2, data=dat, xlab=paste('PC1 (', pve[1], '%)', sep=''), ylab=paste('PC2 (', pve[2], '%)', sep=''), ...)
  #need to set alpha/all.size in geom_point, else it appears in legend
  qp <- ggplot(dat, aes_string(y='PC1', x='PC2', ...)) + theme_bw() + theme(panel.grid=element_line(color='black'))
  if (!is.null(all.size)){ qp <- qp + geom_point(size=all.size, alpha=alpha) } else { qp <- qp + geom_point(alpha=alpha) }
  if (!is.null(facet)){ qp <- qp + facet_grid(facet) }
  qp <- qp + ylab(paste0('PC1 (', pve[1], '%)')) + xlab(paste0('PC2 (', pve[2], '%)', sep=''))
  if (rm.leg.title){ qp <- qp + theme(legend.title=element_blank()) }
  if (add.labels){
    dat2 <- dat
    dat2$row_names <- rownames(pheno)
    qp <- qp + geom_text(data=dat2, mapping=aes_string(y='PC1', x='PC2', label='row_names'), size=2, vjust=-.7)
  }
  
  if (!is.na(name)){ ggsave(filename=paste0(name, '.png'), plot=qp) } else { print(qp) }
  
  if (ret.pcs){
    return(dat)
  } else {
    return(head(dat))
  }
}