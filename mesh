#Packages
install.packages(c("igraph","RPostgreSQL","rgeos","sp","plyr","raster", "concaveman", "gdata", "raster", "tibble", "concaveman", "sf", "spData", "spDataLarge", "rcpp", "ggplot2", "splitstackshape"))
install.packages("tibble", type = 'binary')
install.packages("units", type='binary')
install.packages("pillar", type='binary')
install.packages("spDataLarge")
install.packages("sf", dependencies = TRUE)

#Load libraries
library(igraph)
library(RPostgreSQL)
library(rgeos)
library(sp)
library(plyr)
library(raster) 
library(concaveman) 
library(gdata) 
library(raster) 
library(tibble) 
library(concaveman) 
library(sf) 
library(data.table) 
library(spData) 
library(spDataLarge)
library(pillar)
library(tibble)
library(Rcpp)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(parallel)
library(foreach)
library(splitstackshape)
library(pracma)
library(rlang)

#variables
kfile_vertex <- "dir-for-vertex-buffers"
kfile_vertex_iso <- "dir-for-isolated-ram-buffers"
kfile_mesh <- "mesh-dir"    
lsdynadir <- "ls-prepost-dir"
lsdyna_solver= "solver-dir"
concavehulls <- "concave-hulls-dir"
pre_concave_hulls <- "node-coords-dir"                  
solutions_dir_iso <- "solutions-dir" 
kfiles_complete_isolated <- "completed-mesh-dir"
goodelements <- ""                                               # dir for 'corrected' elements (ram that results from deletions)
topviewdir_af <- "after-fracture-plots-dir"             
topviewdir_bf <- "before-fracture-plots-dir"            
buff = 20                                                        # ram/bench buffer value (m)
mesh_size <- 20.0                                                # mesh resolution ~5x5 m
z_length <- 80.0                                                 # thickness of ice island (m)
zelem <- 8                                                       # no.of elements in the z direction = 5 m
ram_ht <- 80
mat <- readLines("mat_paste.txt")                                   
nidwidth <- c(8, 16, 16, 16, 8, 8)                               # vector of fixed-width row values for nodes section of k files
nidcols <- c(8, 16, 16, 16 ,8, 8)                                # vector of fixed-width column values for nodes section of k file
eidwidth <- rep(8, 10)                                           # vector of fixed-width row values for elements section of k file
eidcols <- rep(8,10)                                          
rho_ice <- 900                                                   # density of ice
rho_water <- 1024                                                # density of water
g <- 9.8                                                         # accelration of gravity on Earth
t_sail <- 10                                                     # sail thickness
t_keel <- 70                                                     # keel thickness 

#make text files to define parameters for various parts of the kfile meshes:

setwd(kfiles_complete_isolated)
cat("*DAMPING_GLOBAL",
    "$#    lcid    valdmp       stx       sty       stz       srx       sry       srz",
    "        0       uuu       0.0       0.0       0.0       0.0       0.0       0.0", sep="\n", file = "damping_coef.txt", append=TRUE)
damp <- readLines("damping_coef.txt")

#mat properties, erosion
cat("*MAT_ADD_EROSION
$#     mid      excl    mxpres     mneps    effeps    voleps    numfip       n
         1       0.0       0.0       0.0       0.0       0.0       1.0       1.0
$#  mnpres     sigp1     sigvm     mxeps     epssh     sigth   impulse    failtm
       0.0  500000.0       0.0       0.0       0.0       0.0       0.0       0.0
$#    idam    dmgtyp     lcsdg     ecrit    dmgexp     dcrit    fadexp    lcregd
         0       0.0         0       0.0       1.0       0.0       1.0         0
$#   lcfld             epsthin    engcrt    radcrt      
         0         0       0.0       0.0       0.0", sep="\n", file="erosion_paste.txt", append=TRUE)
erosion_lines <- readLines("erosion_paste.txt")

### Make a text string of material properties that you can paste into the k files. Not very elegant.######
cat("*CONTROL_ENERGY",
    "$#    hgen      rwen    slnten     rylen   ",
    "         2         2         1         1",
    "*CONTROL_TERMINATION",
    "$#  endtim    endcyc     dtmin    endeng    endmas   ",
    "      200.0         0       0.0       0.01.000000E8", 
    "*CONTROL_TIMESTEP", 
    "$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st",
    "       0.0       0.6         0       0.0       0.0         0         0   ", 
    "$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl    ",
    "       0.0         0         0                           0.0",
    "*DATABASE_BINARY_D3PLOT",
    "$#      dt      lcdt      beam     npltc    psetid  ",
    "      0.10         0         0         0         0",
    "$#   ioopt",
    "         0", 
    "*LOAD_BODY_Z",
    "$#    lcid        sf    lciddr        xc        yc        zc       cid   ",
    "         1       1.0         0       0.0       0.0       0.0         0", 
    "*LOAD_SEGMENT_SET",
    "$#    ssid      lcid        sf        at  ",
    "         1        10       1.0       0.0",
    "*PART",
    "$#",
    "ice island",
    "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid",
    "       2         1         1         0         0         0         0         0",
    "*SECTION_SOLID_TITLE",
    "Solid section",
    "$#   secid    elform       aet",
    "         1         2         0",
    "*MAT_ELASTIC_TITLE",
    "Elastic material", 
    "$#     mid        ro         e        pr        da        db  not used", 
    "         1     900.08.999999E9      0.33       0.0       0.0         0",
    "*DEFINE_FUNCTION",
    "$#     fid                                                               heading",
    "        10        ",
    "$#                                                                      function",
    "float hpres(float t, float x, float y, float z, float x0, float y0, float z0)",
    "{",
    "float fac, trise, refz, rho, grav;",
    "trise = 0.1; refz = 40; rho = 1024; grav = 9.81;",
    "fac = 1.0;",
    "if(t<=trise) fac = t/trise;",
    "return fac*rho*grav*(refz-z);",
    "}",
    "*DEFINE_CURVE_TITLE",
    "Gravity curve",
    "$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint",
    "         1         0       1.0       1.0       0.0       0.0         0         0",
    "$#                a1                  o1 ",
    "                 0.0                 0.0",
    "                 0.1                9.81",
    "               200.0                9.81", sep="\n", file="mat_paste.txt", append=TRUE)

setwd(kfiles_complete_isolated)
mat <- readLines("mat_paste.txt")

#Build the meshes via for loop:   

for (i in 1:length(flist)) {
  setwd(kfile_vertex)
  vertex_kfile(ft_suspects_sp, i ,buff)
  mesh_kfile(kfile_in = flist[i], kfile_vertex, kfilemesh,lsdynadir,mesh_size,z_length, zelem)
  segment_setter(flist[i], kfilemesh, kfilemesh,lsdynadir)
  setwd(kfilemesh)
  mat_paster(flist[i], mat, input_dir = kfilemesh, output_dir = kfilemesh)
  toplevel = nodeser(kfilemesh = flist[i], elev = ram_ht) ##intermediate step here where only the 'right' nodes are supplied for isolated ramget
  ramcut = nick_concave(toplevel, iterations = 1, flist[i])
  setwd(kfilemesh) ##start here
  elements_deleted <- delete_elements_improved(kfilemesh = flist[i], new_nodes = ramcut, ht_cutoff = ram_ht-10) #
  orphans_departed <- goodbye_orphans(new_nodes = ramcut, elements_deleted)
  nodes_width_fixed <- fix_width_nodes(nidwidth, nidcols, orphans_departed)
  elements_width_fixed <- fix_width_elements(eidwidth, eidcols,elements_deleted)
  final_kfile <- nid_eid_replace(kfilemesh = flist[i], nidtable = "temp_fwf_nodes.txt", eidtable = "temp_fwf_elements.txt", inst_names = inst_names_meshes[i], buff = buff, output_dir = kfile_complete_100)
  paste_waterline(vertex_directory = kfile_vertex, ci2d3_df = ft_suspects_sp, t_sail = 10, t_keel = 90, vertex_kfile = flist[i], kfile_input_directory = kfile_complete_100, kfile_output_directory = kfile_complete_100, kfile_mesh = flist[i])
  
}

#add damping coefficient and erosion limit:

for (i in 1:length(flist[i]))
  {
  
  setwd(kfile_vertex_iso) 
  damp_coef <- predict_damp(inst = inst_names_meshes[i], buff_length = 60, thickness = 100, ft_suspect_df = ft_suspects_sp)
  
  setwd(kfiles_complete_isolated_100)
  erode_damp(kfile_mesh = flist[i], damp_coef = damp_coef)
  
}
