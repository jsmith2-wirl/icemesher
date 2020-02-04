
mesh_kfile = function(kfile_in, inputdir, outputdir, lysdynadir, mesh_size=60.0, z_length=50.0, z_elem=10) {
  #Takes the vertex kfile (kfile_in) from an input directory and runs LS DYNA to make a mesh with the inst. name in the outputdir
  #Assumes that working version of lsprepost4.3_x64.exe is in lsdynadir
  #assumes that meshoutline.scl is in the inputdir
  #mesh_size and z_length must be floating points
  #inputdir - directory for vertex kfiles
  #outputdir - output location for this function
  #kfile_in - the vertex kfile
  #lsdynadir - location of lsdyna exec. file
  #mesh_size - length/width dimensions of elem (sq. root) 
  #z_length - Height of z in m
  #z_elem - amount of elems in Z dimension
  
  input <- paste0('parameter input \"', inputdir,'\\',kfile_in,'\"')
  output <- paste0('parameter output \"', outputdir,'\\',kfile_in,'\"')
  fileConn<-file("meshmaker.cfile")
  mesh_size = paste0("parameter meshsize ", format(mesh_size, nsmall=1))
  z_length = paste0("parameter zlength ", format(z_length, nsmall=1))
  z_elem = paste0("parameter zelem ", zelem)
  writeLines(c(input, output, mesh_size,z_length, z_elem, 
               'runscript "meshoutline.scl"  &input &output &meshsize &zlength &zelem'), fileConn)
  close(fileConn)
  system(paste0(lsdynadir,"\\lsprepost4.3_x64.exe c=meshmaker.cfile -nographics"))
  setwd(outputdir)
  flist[i] <- readLines(kfile_in)
}
