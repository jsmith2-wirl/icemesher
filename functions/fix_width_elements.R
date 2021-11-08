fix_width_elements <- function(row_elements, col_elements, element_dataframe) { 
  #takes the element .csvs and then writes a .txt file that respects the fixed-width format that LS DYNA can read
  #row_elements - a vector for the fixed-width spacing of each element in a row
  #col_elements - a vector for the fixed-width spacing of each element in a column
  #element_dataframe - dataframe to fix width of
names(element_dataframe) <- c("$#   eid", "pid", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8")
colnames(element_dataframe) <- formatC(colnames(element_dataframe), width = col_elements, flag = " ")
write.fwf(element_dataframe, file = "temp_fwf_elements.txt", width = row_elements, sep="")

}
