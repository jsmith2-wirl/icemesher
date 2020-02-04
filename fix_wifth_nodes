fix_width_nodes <- function(row_nid, col_nid, node_dataframe) { 
  #takes the node df and then saves it
  #row_nid - a vector for the fixed-width spacing of each element in a row
  #col_nid - a vector for the fixed-width spacing of each element in a column
  #node_dataframe - a .csv of the node to output to the text file
  
  names(node_dataframe) <- c("$#   nid", "x", "y", "z", "t", "c")
  colnames(node_dataframe) <- formatC(colnames(node_dataframe), width = col_nid, flag = " ")
  write.fwf(node_dataframe, file = "temp_fwf_nodes.txt", width = row_nid, sep="")
  
}
