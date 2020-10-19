function(knodes, x1,y1,x2,y2) {
  #k-nodes - a dataframe with node data
  #x1,y1
  #x2,y2 - two points that form a line in lcc coords
  
  y <- c(y1,y2)
  x <- c(x1,x2)
  ln <- lm(y~x)
  side = knodes$y - predict.lm(ln, newdata = data.frame(x=knodes$x)) 
  if(length(which(side<0)) > length(which(side>0)) ) {
    ramside = which(side>0)
    noramside = which(side<0)
    noramside = c(noramside,which(side==0))
  } else {
    ramside = which(side<0)
    noramside = which(side>0)
    noramside = c(noramside,which(side==0))
  }
  return(knodes$`$#   nid`[noramside])
  #return(knodes$keep <- knodes$`$#   nid` %in% ramside)
}
