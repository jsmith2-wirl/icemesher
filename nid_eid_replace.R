nid_eid_replace = function(kfilemesh, eidtable, nidtable, inst_names, buff, output_dir) {
  #takes in the fixed node and element files and pastes them to their respective positions in the K-files
  #kfilemesh - kfile in question
  #eidtable/nidtable - the fixed-width format output (.txt) to read in
  #inst_names - name of the ice island inst
  #buff - equal to the vertex buffer/ram length of the ice island
  #output_dir - location to deposit output
  
  kfilemesh <- readLines(kfilemesh)
  element = which(kfilemesh == "*ELEMENT_SOLID") #new variable for elements
  node = which(kfilemesh == "*NODE") #identifies the range of lines to delete (nodes in this case)
  end = which(kfilemesh == "*END")
  
  filecon = file(paste0(inst_names,"_", buff, ".k"))
  eidtable <- readLines(con = eidtable)
  nidtable <- readLines(con = nidtable)
  
  result <- c(kfilemesh[1:element], eidtable, kfilemesh[node], nidtable, kfilemesh[end])
  
  setwd(output_dir)
  writeLines(result, filecon)
  close(filecon)

  } 
