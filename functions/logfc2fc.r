##jmd
##10.21.16

logfc2fc <- function(logFC){
  #sign is -1 if logFC<0; 1 if logFC>=0
  sgn <- (-1)^(1+as.numeric(logFC>=0))
  fc <- sgn*2^abs(logFC)
  return(fc)
}