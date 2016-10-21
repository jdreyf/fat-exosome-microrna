##jmd+hp
##3.17.16

library('limma')

##run limma with contrasts and multiTopTab
#dat is a matrix or EList from limma voom
#grp isn't needed if provide design & add.means is FALSE
#contrasts.v is a vector
#can make contrasts.v based on design matrix not defined by grp
# dat=mat; grp=pheno$grp; contrasts.v=contr.v; add.means=TRUE; adjust.method='BH'; weights=NULL; design=NULL; 
# trend=FALSE; cols=c('P.Value', 'adj.P.Val', 'logFC')
limma.contrasts <- function(dat, grp=NULL, contrasts.v, add.means=TRUE, adjust.method='BH', weights=NULL, design=NULL, 
                            trend=FALSE, cols=c('P.Value', 'adj.P.Val', 'logFC'), block = NULL, correlation = NULL){
  if (is.null(design)|add.means) stopifnot(ncol(dat)==length(grp), colnames(dat)==names(grp))

  if (is.null(design)){
    design <- model.matrix(~0+grp)
    colnames(design) <- sub('grp', '', colnames(design), fixed=TRUE)
  }
  
  #can't set weights=NULL in lmFit when using voom, since lmFit only assigns weights "if (missing(weights) && !is.null(y$weights))"
  if (!missing(weights)){
    if (!is.matrix(dat) && !is.null(dat$weights)){ cat('dat$weights are being ignored\n') }
    fit <- lmFit(dat, design, block = block, correlation = correlation, weights=weights)
  } else {
    fit <- lmFit(dat, design, block = block, correlation = correlation)
  }
  
  contr.mat <- makeContrasts(contrasts=contrasts.v, levels=design)
  fit2 <- contrasts.fit(fit, contr.mat)
  fit2 <- eBayes(fit2, trend=trend)
  #limma ignores names of contrasts.v when it's given as vector
  if (!is.null(names(contrasts.v))){
    stopifnot(colnames(fit2$contrasts)==contrasts.v)
    colnames(fit2$contrasts) <- names(contrasts.v)
  }
  mtt <- multiTopTab(fit2, cols=cols, adjust.method=adjust.method)
  
  #cbind grp means
  if (add.means){
    grp.means <- t(apply(dat, 1, FUN=function(v) tapply(v, grp, mean, na.rm=TRUE)))
    colnames(grp.means) <- paste(colnames(grp.means), 'avg', sep='.')
    mtt <- cbind(grp.means[rownames(mtt),], mtt)
  }
  return(mtt)
}