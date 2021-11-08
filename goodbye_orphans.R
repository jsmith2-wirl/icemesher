goodbye_orphans <- function(new_nodes, deleted_elements_df) {
  #takes the csv with deleted elements and eliminates the unref'd nodes
  #new_nodes - the ramcut nodes df
  #deleted_elements_df - df of elements that have been deleted
  

nid <- new_nodes[,1]
included <- array(data=NA,dim = c(length(nid), 8))
included[,1] <- nid %in% deleted_elements_df$n1  # if true then the node is in the first column of this element
included[,2] <- nid %in% deleted_elements_df$n2
included[,3] <- nid %in% deleted_elements_df$n3
included[,4] <- nid %in% deleted_elements_df$n4
included[,5] <- nid %in% deleted_elements_df$n5
included[,6] <- nid %in% deleted_elements_df$n6
included[,7] <- nid %in% deleted_elements_df$n7
included[,8] <- nid %in% deleted_elements_df$n8
ind <-which(rowSums(included)==0) 

included_df <- data.frame(included)
included_df <- add_column(included_df, nid = 1:nrow(included_df), .before = 1)
included_df <- new_nodes[-ind,]

}
