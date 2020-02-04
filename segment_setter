segment_setter = function(kfile_in, inputdir, outputdir, lysdynadir) {
  # identifies the bottom surface of a mesh (for determining hydrostatic forces) 
  # Input: kfile_in - a 3-D mesh K-file
  # inputdir - the directory to take the meshes from
  # outputdir - the directory to put the updated meshes in (should be the same as input)
  # lsdynadir - directory w/ lsdyna .exec files
  # Assumes that working version of lsprepost4.3_x64.exe is in lsdynadir
  # assumes that meshoutline.scl is in the inputdir
  
  input <- paste0('open keyword \"', inputdir,'\\',kfile_in,'\"')
  output <- paste0('save keyword  \"', outputdir,'\\',kfile_in,'\"')
  fileConn<-file("segment_setter.cfile")
  two = paste0("Message 1")  
  three = paste0("save keywordoutversion 9") 
  four = paste0("Message 0") 
  five = paste0("Message 1")  
  six = paste0("save keywordoutversion 9") 
  seven = paste0("Message 0") 
  eight = paste0("genselect target part")  
  nine = paste0("setsegment") 
  ten = paste0("genselect target segment") 
  eleven = paste0("genselect clear") 
  twelve = paste0("genselect element add plane in 100.0 100.0 0.0 0 0 1 1 0") 
  thirteen = paste0("setsegment createset 1 1 0 0 0 0 ")  
  fourteen = paste0("genselect clear") 
  fifteen = paste0("save keywordabsolute 0")
  sixteen = paste0("save keywordbylongfmt 0") 
  seventeen = paste0("save outversion 7")
  
  
  writeLines(c(input, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, thirteen, fourteen, fifteen, sixteen, seventeen, output), fileConn)
  close(fileConn)
  system(paste0(lsdynadir,"\\lsprepost4.3_x64.exe c=segment_setter.cfile -nographics"))
  #segment_mesh <- readLines(kfile_in)
}
