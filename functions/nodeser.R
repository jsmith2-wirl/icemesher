  nodeser = function(kfilemesh, elev) {
  #takes the desired amount of node layers and outputs dfs of these for processing concave hulls
  #kfilemesh - the 3-D mesh K-file
  #elev - the elevation that you want to select (there and above)
  #returns a dataframe
  
  myNodes = readk(kfilemesh, type = "nodes")
  top = myNodes[which(myNodes$z >= elev),] #I've just been changing the values in-function rather than making them an argument
  xyz=cbind(top[,1],top$x,top$y, top$z, top$tc, top$rc)  
  dfz <- as.data.frame(xyz)
  names(dfz) <- c("nid", "x", "y", "z", "r", "t")
  return(dfz)
} 
