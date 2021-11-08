vertex_kfile = function(ci2d3_spdf, row, buff=0) {
  #function to write a K-file of the vertices of an ice island polygon. It also buffers the points outward to account for the ram.
  #input: ci2d3_spdf -- a spatial data frame with ice island polygons and attributes as per ci2d3
  #       row -- the row number within that data frame that contains the ice island of interest
  #       buff -- a buffer value in metres - defaults to 0 (no buffer)
  #output: writes a kfile to the current directory - using the inst of the ice island
  
  #Buffer Function
  #Widens the area of the polygon to account for underwater ram
  
  #Buffer variable (m)
  #Assumes a ram of x length
  
  newpoly = buffer(ci2d3_spdf[row,], width=buff)
  coords = newpoly@polygons[[1]]@Polygons[[1]]@coords
  coords = cbind(1:dim(coords)[1], coords,rep(0,dim(coords)[1]))
  fname = paste0(ci2d3_spdf$inst[row],"_", buff, '.k')
  write.table(coords, 'temp.csv', row.names = FALSE, col.names = FALSE, sep=",")
  
  # read that file back in as one block of text and add the header + footer and overwrite it ##removed line breaks where necessary.
  vertices = readLines('temp.csv')
  header = '*Keyword\n*Node\n$ Node,X,Y,Z'
  footer = '*End'
  outfile = file(fname)
  writeLines(c(header, vertices, footer), outfile)
  close(outfile)
}
