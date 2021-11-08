nick_concave <- function(nodestable, iterations, kfile_mesh) {
  # takes a data frame of node (upper slices) and draws a concave hull, deletes intersecting points 
  # interations - number of times to shave off the ram which should be based on the resolution of your elements in x/y directions
  # nodestable - table of nodes you want "shaved"
  # kfile_mesh - the mesh you're interested in
  points <- data.table(nodestable[,2:3])
  
  for (i in 1:iterations){  
    points <- ramstep(points)
  }
  
  
  all_nodes <- readk(kfile_mesh, type = "nodes")
  inner_points <- dplyr::semi_join(x = nodestable, y = points, by = c("x" = "x", "y" = "y"))
  concave_hull <- dplyr::anti_join(x = nodestable, y = inner_points, by = "nid")
  new_node_points <- dplyr::anti_join(x = all_nodes, y = concave_hull, by = c("$#   nid" = "nid")) 
  
  return(new_node_points)
} 
