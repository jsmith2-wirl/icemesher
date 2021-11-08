#create solutions for meshes

for (i in 1:length(flist)) {
  
  setwd(paste0(solutions_dir, "\\", flist[i])) 
  
  shell(paste0(lsdyna_solver,"\\ls-dyna_smp_d_R10.0_winx64_ifort160.exe I=", kfiles_complete_isolated, "\\", flist[i], " O=",solutions_dir, "\\", flist[i], "\\solved   NCPU =  8 MEMORY = 200000000"))
  shell("exit")
  

}
