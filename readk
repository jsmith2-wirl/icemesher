readk = function(kfile_mesh, type='nodes') {
  # function to read a kfile_mesh file and return a dataframes for one of: nodes, elements and segments
  # kfile_mesh = a 3-D mesh K-file
  # type = either "nodes", "elements" or "segments" or "segment_set" 
  
  kfile = readLines(kfile_mesh)
  #hdrs = grep(kfile, pattern="*",fixed = T)
  #kfile[hdrs]
  if (type=='nodes') {
    node1 = grep(kfile, pattern="*NODE",fixed = T)
    node2 = grep(kfile, pattern="*END",fixed = T)
    knodes = read.table(kfile_mesh, header = FALSE, sep="", skip = node1+1, nrows = node2-2-node1, fill = TRUE)
    names(knodes) =c('$#   nid','x','y', 'z', 'tc', 'rc')
    return(knodes)
    
  }  else if (type=='elems') {    
    elem1 = grep(kfile, pattern="*ELEMENT_SOLID",fixed = T)
    elem2 = grep(kfile, pattern="*NODE",fixed = T)
    kelems = read.table(kfile_mesh, header = FALSE, sep="", fill= TRUE, skip = elem1+1, nrows = elem2-2-elem1)
    names(kelems) =c('eid','pid', 'n1','n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8')
    return(kelems)
    
  } else if (type=='segments') {
    seg1 = grep(kfile, pattern="*SET_SEGMENT",fixed = T)
    seg2 = grep(kfile, pattern="*DAMPING_GLOBAL",fixed = T)
    if(length(seg1) == 0 | length(seg2) == 0) {return("SET_SEGMENT Not Found")}
    ksegs = read.table(kfile_mesh, header = FALSE, sep="", skip = seg1+3, nrows = seg2-4-seg1)
    names(ksegs) =c('n1','n2','n3', 'n4', 'a1', 'a2', 'a3','a4')
    return(ksegs)
    
  } else {
    return("Error")
  }
}
