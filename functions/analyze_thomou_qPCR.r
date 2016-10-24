##jmd
##10.4.16

##also return ndct so that Rmd can make heatmaps for supp
analyze_thomou_qPCR <- function(ct, filt.thresh=34, norm.nm="Mouse U6 snRNA", contrasts.v=c(WTvsKO="WT-KO"), quiet=FALSE){
  stopifnot(norm.nm %in% c("globalmean", rownames(ct)))
  
  ##group/phenotype info, taken from regular sample names
  grp <- setNames(gsub('(\\.|)[0-9]$', '', colnames(ct)), colnames(ct))
  pheno <- data.frame(sample=names(grp), grp=grp)
  rownames(pheno) <- names(grp)
  
  ##filter
  #want 50% of smallest group size present
  #can't use group info, as per bourgon, gentleman & huber (pmid:20460310)
  if (!is.na(filt.thresh)){
    #get min n per group for filtering
    min.n.per.grp <- min(table(as.factor(grp)))
    mirs.keep <- rownames(ct)[rowSums(ct<=filt.thresh) > 0.5*min.n.per.grp]
    #but don't cut out normalization probe
    ct.filt <- ct[mirs.keep,]
  }
  
  ##tell number of microRNAs not counting u6
  if (!quiet){
    cat("There were", nrow(ct)-1, "miRNAs profiled")
    if (!is.na(filt.thresh)){
      cat(" of which", nrow(ct.filt)-1, "were detected.\n")
      #min number of times a mir must be detectable
      cat("Min n detectable:", floor(0.5*min.n.per.grp)+1, "\n")
    } else {
      cat("\n") 
    }
  }
  
  ##norm
  #row needs to be coerced to vector
  #Mean normalization from Mestdagh et al. (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2009-10-6-r64).
  if (norm.nm=="globalmean"){
    #filt.thresh is expression threshold
    if (is.na(filt.thresh)) filt.thresh <- max(ct.filt)
    #normalizer is mean value of expressed miRNAs
    norm.v <- apply(ct, MARGIN=2, FUN=function(v){ mean(v[v<filt.thresh]) })
  } else {
    norm.v <- as.numeric(ct[norm.nm,])
  }
  #negative delta Ct
  ndct0 <- -1*scale(ct, center=norm.v, scale=FALSE)
  ndct <- ndct0[setdiff(rownames(ndct0), norm.nm),]
  if (!is.na(filt.thresh)){ ndct.filt <- ndct[intersect(mirs.keep, rownames(ndct)),] } else { ndct.filt <- ndct  }

  ##stats
  if (!is.na(contrasts.v)[1]){
    toptab <- limma.contrasts(ndct, grp, contrasts.v = contrasts.v)
    if (!is.na(filt.thresh)){ toptab.filt <- limma.contrasts(ndct.filt, grp, contrasts.v = contrasts.v) } else { toptab.filt <- toptab  }
  } else {
    toptab <- NULL
  }
  
  return(list(pheno=pheno, ndct=ndct, ndct.filt=ndct.filt, toptab=toptab, toptab.filt=toptab.filt))
}