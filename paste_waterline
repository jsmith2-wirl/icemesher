paste_waterline <- function(vertex_directory, ci2d3_df, t_sail, t_keel, vertex_kfile, kfile_input_directory, kfile_output_directory, kfile_mesh) {
  #vertex_directory - location of original vertex-only K-files
  #ci2d3_df - 
  #t_sail - thickness of freeboard (m)
  #t_keel thickness of submerged part of ice island (m)
  #2-D mesh K-file (vertices only)
  #kfile_input_directory - location of most up-to-date 3-D K-files
  #kfile_output_directory - where you want to put the output 3-D K-Files
  
  setwd(vertex_directory)  
  refz <- find_waterline(ci2d3_df, t_sail, t_keel, vertex_kfile)
  refz <- signif(refz, 2)
  
  setwd(kfile_input_directory)
  kmesh <- readLines(kfile_mesh)
  line1 = grep(kmesh, pattern = "trise = 0.1; refz = 40; rho = 1024; grav = 9.81;",fixed = T)
  water_line <- gsub("40", refz, kmesh[line1])
  kmesh_final = c(kmesh[1:line1-1], water_line, kmesh[line1+1:length(kmesh)])
  end = grep(kmesh, pattern = "*END",fixed = T)
  
  setwd(kfile_output_directory)
  writeLines(kmesh_final[1:end], kfile_mesh)
  
}
