delete_elements_improved <- function(kfilemesh, new_nodes, ht_cutoff) { 
  #deletes elements that reference non-existent nodes (so ones that have been subjected to concave hull funct.)
  #kfilemesh - 3-D mesh K-files to extract element tables from
  #new_nodes - table of nodes with the concave hull/buffer zone removed
  #ht_cutoff - only work above this z value to speed up the function
  
  elementdf <- readk(kfilemesh, type = 'elems') #all the elems
  names(new_nodes) <- c("nid", "x", "y", "z", "tc", "rc")
  upper_elements <- elementdf %>% 
    gather(source, nid, -eid) %>% 
    inner_join(new_nodes %>% filter(z>=ht_cutoff)) %>% 
    semi_join(elementdf, .)
  elements_deleted <- upper_elements[apply(upper_elements[,3:10], 1, function(x) all(x %in% new_nodes[,1])),] #result of slicing upper
  deleted_elements <- anti_join(x = upper_elements, y = elements_deleted, by = "eid") 
  final_elements <- dplyr::anti_join(x= elementdf, y = deleted_elements, by = "eid") 
}
###SUPER IMPORTANT: THE Z HERE MUST BE LESS THAN THE ELEVATION USED IN NODESER
