function(kfile, x1, y1, x2, y2, sailheight, concave_iterations) {
  
  ##kfile should be the flist vector
  
  knodes <- readk(kfile_mesh = kfile, type = "nodes") #all the nodes in the mesh
  noramside <- data.frame(isolate_a_ram(knodes = knodes, x1 = x1, y1 = y1, x2 = x2, y2 =y2)) #find nodes NOT on the edge of this inst
  names(noramside) <- "nid"
  names(knodes)[1] <- "nid"
  knodes$dontkeep <-knodes$nid %in% noramside$nid #flag nonramside for deletion
  knodes$dontkeep[knodes$z >= sailheight & knodes$dontkeep == 'FALSE'] <- "TRUE" #flag ramside points >== 90 m for deletion
  nodes_to_delete <- knodes[which(knodes$dontkeep == 'TRUE'),] #now subset the ones that need to be deleted.
  
  
  return(concaved_nodes <- nick_concave(nodes_to_delete, concave_iterations, kfile_mesh = kfile)) # draw a 20 m concave hull, and delete those nodes.
  
}
