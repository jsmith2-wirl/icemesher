### Setting things up ###
# This script generates simple meshes of icebergs and ice islands from vertices of polygons that Luc and co. digitized for CI2D3.

# It then runs these meshes through a finite element model developed by Dr. Mahmud Sazidy to analyze their 1st principal stress distributions and determine
# if/where a fracture will occur in the iceberg/island.


###get your workspace back####
### BE CAREFUL!!! depending on where you are in this script you will
### be picking up from different directories. 


getwd()
setwd("J:/SCRATCH/jsmith/Solutions")
setwd(scratchdrivemesh) #most mesh stuff happens from here/leave hashed out so you have to think before running it
setwd("C:/Users/jsmith/Desktop")
save.image(file='.RData') #always save environment to scratchdrivemesh
setwd()

# Packages ####
#Necessary packages
install.packages(c("igraph","RPostgreSQL","rgeos","sp","plyr","raster", "concaveman", "gdata", "raster", "tibble", "concaveman", "sf", "spData", "spDataLarge", "rcpp", "ggplot2", "splitstackshape"))

install.packages("tibble", type = 'binary')
                 
                 ###new packages!!!! (06/02/2019)

install.packages("units", type='binary')
install.packages("pillar", type='binary')
install.packages("spDataLarge")
install.packages("sf", dependencies = TRUE)

"igraph"# Load libraries
library(igraph)
library(RPostgreSQL)
library(rgeos)
library(sp)
library(plyr)
library(raster) 
library(concaveman) #new library to load (15/01/2019)
library(gdata) #new library to load (21/01/2019)
library(raster) #new library (24/01/2019)
library(tibble) #new library (05/02/2019)
library(concaveman) #new library (05/02/2019)
library(sf) #new library (05/02/2019)
library(data.table) #new library (06/02/2019)
library(spData) #new library (06/02/2019)
library(spDataLarge) #new library (06/02/2019)
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



setwd(dir = "C:/Users/jsmith/Desktop/all_thesis_stuff/CI2D3Stuff")
source('f.fract.R')
source('f.fract_af.R')
source('f.fract_bf.R')
source('f.Spatialdf.R')
source('f.subquery.R')

# Parameters ####

# Need to set up these directories properly.  Make sure you have the correct path for the OS use double \\ for windows
kfile_vertex="C:\\Users\\jsmith\\Desktop\\KFiles1"                  # dir for vertex buffers
kfile_vertex_iso="C:\\Users\\jsmith\\Desktop\\KFiles1_Iso"
kfile_mesh="C:\\Users\\jsmith\\Desktop\\KFiles2"  #                 # old dir for meshes. decommissioned.
lsdynadir="C:\\LSTC\\LS-PrePost\\4.3-x64"                           # dir for ls prepost
lsdyna_solver= "C:\\LSDYNA\\program"
scratchdrivemesh <- "J:\\SCRATCH\\jsmith\\KFiles2"                  # dir for meshes (scratch drive bec files are huge)
scratchdrivemesh2 <- "J:\\SCRATCH\\jsmith\\KFiles2_2"
concavehulls = "C:\\Users\\jsmith\\Desktop\\Concave_Hulls"          # dir for concave hull/node deletions
pre_concave_hulls = "C:\\Users\\jsmith\\Desktop\\Pre_Concave_Hulls" # dir for df with node coordinates for runninin concave hull
solutions_dir <- "J:\\SCRATCH\\jsmith\\Solutions"                   #dir for solutions
solutions_dir_600 <- "J:\\SCRATCH\\jsmith\\Solutions_600"  
solutions_dir_700 <- "J:\\SCRATCH\\jsmith\\Solutions_700"  
solutions_dir2 <- "J:\\SCRATCH\\jsmith\\Solutions_2" 
solutions_dir_iso <- "J:\\SCRATCH\\jsmith\\Solutions_Isolated" 
solutions_dir_iso_600 <- "J:\\SCRATCH\\jsmith\\Solutions_Isolated_600" 
solutions_dir_iso_700 <- "J:\\SCRATCH\\jsmith\\Solutions_Isolated_700"
solutions_dir_no_max <- "J:\\SCRATCH\\jsmith\\Solutions_No_Max"
kfile_complete <- "J:\\SCRATCH\\jsmith\\KFiles_Complete"            #dir for finished meshes
kfile_complete2 <- "J:\\SCRATCH\\jsmith\\KFiles_Complete_2"
kfiles_complete_isolated <- "J:\\SCRATCH\\jsmith\\KFiles_Complete_Isolated"
test_complete2 <- "J:\\SCRATCH\\jsmith\\Test_Complete_2" 
goodelements = "C:\\Users\\jsmith\\Desktop\\Good_Elements/"         # dir for 'corrected' elements (ram that results from deletions)
topviewdir_af = "C:\\Users\\jsmith\\Desktop\\AFKFiles"              # dir for after-fracture plots
topviewdir_bf = "C:\\Users\\jsmith\\Desktop\\BFKFiles"             # dir for before-fracture plots
buff = 0                                                           # ram/bench buffer value (m)
mesh_size=20.0 # must be a floating point number                     # mesh resolution ~5x5 m
z_length=80.0 # must be a floating point number                     # thickness of ice island (m)
zelem=8                                                            # no.of elements in the z direction = 5 m
ram_ht =80
mat <- readLines("mat_paste.txt")                                   
nidwidth = c(8, 16, 16, 16, 8, 8)                                   # vector of fixed-width row values for nodes section of k files
nidcols = c(8, 16, 16, 16 ,8, 8)                                    # vector of fixed-width column values for nodes section of k file
eidwidth = rep(8, 10)                                                # vector of fixed-width row values for elements section of k file
eidcols = rep(8,10)  
rho_ice <- 900
rho_water <- 1024
g <- 9.8
t_sail <- 10 
t_keel <- 70                                                      # vector of fixed-width column values for elements section of k file

# Functions ####

vertex_kfile = function(ci2d3_spdf, row, buff=0) {
  #function to write a k file of the vertices of a polygon. also buffers the points 50 m outward to account for ram.
  #input: ci2d3_spdf -- a spatial data frame with ice island polygons and attributes as per ci2d3
  #       row -- the row number within that data frame that contains the ice island of interest
  #       buff -- a buffer value in metres - defaults to 0 (no buffer)
  #output: writes a kfile to the current directory - using the inst of the ice island
  
  #Buffer Function
  #Widens the area of the polygon to account for underwater ram
  
  #Buffer variable (m)
  #Assumes a ram of x length
  
  newpoly = buffer(ci2d3_spdf[row,], width=buff)
  coords = newpoly@polygons[[1]]@Polygons[[1]]@coords
  coords = cbind(1:dim(coords)[1], coords,rep(0,dim(coords)[1]))
  fname = paste0(ci2d3_spdf$inst[row],"_", buff, '.k')
  write.table(coords, 'temp.csv', row.names = FALSE, col.names = FALSE, sep=",")
  
  # read that file back in as one block of text and add the header + footer and overwrite it ##removed line breaks where necessary.
  vertices = readLines('temp.csv')
  header = '*Keyword\n*Node\n$ Node,X,Y,Z'
  footer = '*End'
  outfile = file(fname)
  writeLines(c(header, vertices, footer), outfile)
  close(outfile)
}

mesh_kfile = function(kfile_in, inputdir, outputdir, lysdynadir, mesh_size=60.0, z_length=50.0, z_elem=10) {
  #Takes the vertex kfile (kfile_in) from inputdir and runs lsdyna to make a mesh with the inst. name in the outputdir
  #Assumes that working version of lsprepost4.3_x64.exe is in lsdynadir
  #assumes that meshoutline.scl is in the inputdir
  # mesh_size and z_length must be floating points
  
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

segment_setter = function(kfile_in, inputdir, outputdir, lysdynadir) {
  # identifies the segment set (faces of elements for hydrostatic force) of each mesh 
  # Input: kfile_in - a k file that has been meshed
  # inputdir - the directory to take the meshes from
  # outputdir - the directory to put the updated meshes in (should be the same as input)
  # lsdynadir - directory w/ lsdyna .exec files
  # Assumes that working version of lsprepost4.3_x64.exe is in lsdynadir
  # assumes that meshoutline.scl is in the inputdir
  # mesh_size and z_length must be floating points
  
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

readk = function(kfile_mesh, type='nodes') {
  # function to read a kfile_mesh file and return a dataframes for one of: nodes, elements and segments
  # kfile_mesh = a valid kfile (3D)
  # type = either "nodes", "elements" or "segments" or "segment_set" 
  
  kfile = readLines(kfile_mesh)
  #hdrs = grep(kfile, pattern="*",fixed = T)
  #kfile[hdrs]
  if (type=='nodes') {
    node1 = grep(kfile, pattern="*NODE",fixed = T)
    node2 = grep(kfile, pattern="*END",fixed = T)
    knodes = read.table(kfile_mesh, header = FALSE, sep="", skip = node1+1, nrows = node2-2-node1, fill = TRUE)
    names(knodes) =c('$#   nid','x','y', 'z', 'tc', 'rc')
    return(knodes)
    
  }  else if (type=='elems') {    
    elem1 = grep(kfile, pattern="*ELEMENT_SOLID",fixed = T)
    elem2 = grep(kfile, pattern="*NODE",fixed = T)
    kelems = read.table(kfile_mesh, header = FALSE, sep="", fill= TRUE, skip = elem1+1, nrows = elem2-2-elem1)
    names(kelems) =c('eid','pid', 'n1','n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8')
    return(kelems)
    
  } else if (type=='segments') {
    seg1 = grep(kfile, pattern="*SET_SEGMENT",fixed = T)
    seg2 = grep(kfile, pattern="*DAMPING_GLOBAL",fixed = T)
    if(length(seg1) == 0 | length(seg2) == 0) {return("SET_SEGMENT Not Found")}
    ksegs = read.table(kfile_mesh, header = FALSE, sep="", skip = seg1+3, nrows = seg2-4-seg1)
    names(ksegs) =c('n1','n2','n3', 'n4', 'a1', 'a2', 'a3','a4')
    return(ksegs)
    
  } else {
    return("Error")
  }
}

mat_paster = function(kfilemesh, mat, input_dir, output_dir) {
  # pastes the material property lines to the appropriate section of each k file
  # kfilemesh - k file (that has already been given a segment set)
  # mat - read-in file of material properties
  # outputs a k file with material properties
  
  setwd(input_dir)
  kmesh <- readLines(kfilemesh)
  kmesh = c(kmesh[1:6], mat, kmesh[12:length(kmesh)])
  setwd(output_dir)
  writeLines(kmesh, kfilemesh)
  
} 

nodeser = function(kfilemesh, elev) {  
  #takes the desired amount of node layers and outputs dfs of these for processing concave hulls
  #kfile mesh - the k file(s)
  #elev - the elevation that you want to select (there and above)
  #returns a dataframe
  
  myNodes = readk(kfilemesh, type = "nodes")
  top = myNodes[which(myNodes$z >= elev),] #I've just been changing the values in-function rather than making them an argument
  xyz=cbind(top[,1],top$x,top$y, top$z, top$tc, top$rc)  
  dfz <- as.data.frame(xyz)
  names(dfz) <- c("nid", "x", "y", "z", "r", "t")
  return(dfz)
} 

ramstep <- function(dt) {
  #supply a data table with nodes coords x and y and outer concave hull will be removed
  DT_sf = sf::st_as_sf(dt, coords = c("x", "y"),agr = "constant")
  point_coords <- as.data.frame(st_coordinates(st_as_sf(DT_sf)))
  names(point_coords) <- c("x", "y")
  #split_coords <- read.table(text=gsub('[c()]', '', DT_sf$geometry), 
  #sep=",", col.names=c('x', 'y'))                                      #Useful if you need the nids
  #all_points <- data.frame(DT_sf[,1], split_coords[,1:2], DT_sf[,-1])
  #all_points$geometry <- NULL
  #all_points$geometry.1 <- NULL
  
  c_hull <- concaveman(DT_sf)
  c_hull_coords <- st_coordinates(c_hull$polygons)
  c_hull_coords <- as.data.frame(c_hull_coords[, -c(3:4)])
  names(c_hull_coords) <- c("x", "y")
  
  hulled_points <- dplyr::anti_join(point_coords, c_hull_coords, by = c("x" = "x", "y" = "y")) 
  return(as.data.table(hulled_points))
} 
 
nick_concave <- function(nodestable, iterations, kfile_mesh) {
  # takes a data frame of node (upper slices) and draws a concave hull, deletes intersecting points 
  # interations - number of times to shave off the ram which should be based on the resolution of your elements in x/y directions
  # nodestable - table of nodes you want "shaved"
  #kfilemesh - the mesh you're interested in
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

delete_elements <- function(kfilemesh, new_nodes) { 
  #deletes elements that reference non-existent nodes (so ones that have been subjected to concave hull funct.)
  #kfilemesh - k files to extract element tables from
  #concavehulls - the directory to put the updated element in (.csv)
  #nodedf - node .csvs after concave hull function
  
  elementdf <- readk(kfilemesh, type = 'elems')
  elements_deleted <- elementdf[apply(elementdf[,3:10], 1, function(x) all(x %in% new_nodes[,1])),] ####kind of works :)))))

} 

delete_elements_improved <- function(kfilemesh, new_nodes, ht_cutoff) { 
  #deletes elements that reference non-existent nodes (so ones that have been subjected to concave hull funct.)
  #kfilemesh - k files to extract element tables from
  #concavehulls - the directory to put the updated element in (.csv)
  #nodedf - node .csvs after concave hull function
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
###SUPER IMPORTANT: THE Z HERE MUST BE LESS THAN THE ELEV IN NODESER

goodbye_orphans <- function(new_nodes, deleted_elements_df) {
  #takes the csv with deleted elements and eliminates the unref'd nodes
  #lesselements_csv - the csv that has been both concaved, and element-shrunk
  #output_dir - where you want the result. preferably 'concavehulls'
  

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

} ##good 02/15/2019

fix_width_nodes <- function(row_nid, col_nid, node_dataframe) { 
  #takes the node df and then saves it
  #row_nid - a vector for the fixed-width spacing of each element in a row
  #col_nid - a vector for the fixed-width spacing of each element in a column
  #node_dataframe - a .csv of the node to output to the text file
  
  names(node_dataframe) <- c("$#   nid", "x", "y", "z", "t", "c")
  colnames(node_dataframe) <- formatC(colnames(node_dataframe), width = col_nid, flag = " ")
  write.fwf(node_dataframe, file = "temp_fwf_nodes.txt", width = row_nid, sep="")
  
}  #good 02/15/2019

fix_width_elements <- function(row_elements, col_elements, element_dataframe) { 
  #takes the element .csvs and then writes a .txt file that respects the fixed-width format that LS DYNA can read???
  #row_elements - a vector for the fixed-width spacing of each element in a row
  #col_elements - a vector for the fixed-width spacing of each element in a column
  #element_dataframe - dataframe to fix width of
names(element_dataframe) <- c("$#   eid", "pid", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8")
colnames(element_dataframe) <- formatC(colnames(element_dataframe), width = col_elements, flag = " ")
write.fwf(element_dataframe, file = "temp_fwf_elements.txt", width = row_elements, sep="")

} ##good 02/15/2019

nid_eid_replace = function(kfilemesh, eidtable, nidtable, inst_names, buff, output_dir) {
  
  #takes in the fixed node and element files and pastes them to their respective positions in the kfiles
  #kfilemesh - kfile in question
  # eidtable/nidtable - the fixed-width format output (.txt) to read in
  kfilemesh <- readLines(kfilemesh)
  element = which(kfilemesh == "*ELEMENT_SOLID") #new variable for elements
  node = which(kfilemesh == "*NODE") #identifies the range of lines to delete (nodes in this case)
  end = which(kfilemesh == "*END")
  
  filecon = file(paste0(inst_names,"_", buff, ".k"))
  eidtable <- readLines(con = eidtable)
  nidtable <- readLines(con = nidtable)
  
  result <- c(kfilemesh[1:element], eidtable, kfilemesh[node], nidtable, kfilemesh[end])
  
  setwd(output_dir)
  writeLines(result, filecon)
  close(filecon)

  } 

get_stress_table <- function(solutions_dir, flist, lsdynadir, states) {
  ###gets a table of max principal stresses throughout all elements/states from solution 1
  ###solutins_dir - directory to grab solution file from
  ###flist - list of kfile names
  ###lsdynadir - directory fort lspp 10.exe file
  ###state - number of states to look at
  
  #numbs = seq(1:100)
  #states <- capture.output(cat(numbs, sep=":"))
  biggest_stresses <- data.frame(eid= numeric(),
                                 stress = numeric(),
                                 stringsAsFactors=FALSE) 
  #second_stresses <- data.frame(eid= numeric(),
                                #stress = numeric(),
                                #stringsAsFactors=FALSE)
 
  for (j in 1:states) {
    fileconn <- file("get_stresses.cfile")
    line_one <- paste0("open d3plot ", solutions_dir, "\\", flist, "\\", "d3plot")
    line_two <- "ac"
    line_three <- "fringe 14"
    line_four <- "pfringe"
    line_five <- "anim forward"
    line_six <- "anim stop; state 100;"
    line_seven <- paste0("output ", solutions_dir, "\\", flist, "\\", flist, " ", j, " 1 1 1 0 0 0 0 1 0 0 0 0 0 0 1.000000")
    writeLines(c(line_one, line_two, line_three, line_four, line_five, line_six, line_seven), fileconn)
    close(fileconn)
    system(paste0(lsdynadir,"\\lsprepost4.3_x64.exe c=get_stresses.cfile -nographics"))
    
  
    stresses <- readLines(flist)
    start <- grep(stresses, pattern="*KEYWORD",fixed = T)
    stop <- grep(stresses, pattern="$Interpreted from averaged nodal data",fixed = T)
    stresses <- stresses[-seq(start, stop, by = 1)]
    writeLines(stresses, flist)
    stresses <- read.table(flist, header = FALSE)
    names(stresses) <- c("eid", "stress")
    max_stress <- which(stresses$eid == which.max(stresses$stress))
    #second_stress <- subset(stresses,sort(z<-rank(stress),T)[2]==z)[,1]
    
    
    
    biggest_stresses <- rbind(biggest_stresses, stresses[max_stress,])
    #biggest_stresses[nrow(biggest_stresses)+1,] = list(stresses[max_stress,1], stresses[max_stress,2])
    #second_stresses <- rbind(second_stresses, stresses[second_stress,])
    
  }
  
  #return(biggest_stresses)
  return(biggest_stresses[which.max(biggest_stresses$stress),1:2])
  #return(second_stresses[which.max(second_stresses$stress),1])

  
}

get_max_stress_eid <- function(solutions_dir, sub_directory, kfile, lsdynadir, states) {
  ###gets a table of max principal stresses throughout all elements/states from a solution
  ###solutins_dir - directory to grab solution file fro.m
  ###kfile - kfile mesh in question
  ###lsdynadir - directory fort lspp 10.exe file
  ###state - number of states to look at
  ###output is eid and max princip stress
  
  #numbs = seq(1:100)
  #states <- capture.output(cat(numbs, sep=":"))
  biggest_stresses <- data.frame(eid= numeric(),
                                 stress = numeric(),
                                 stringsAsFactors=FALSE) 
  
  
  for (j in 900:states) {
    fileconn <- file("get_stresses.cfile")
    line_one <- paste0("open d3plot ", solutions_dir, "\\", sub_directory, "\\", "d3plot")
    line_two <- "ac"
    line_three <- "fringe 14"
    line_four <- "pfringe"
    line_seven <- paste0("output ", solutions_dir, "\\", sub_directory, "\\", kfile, " ", j, " 1 1 1 0 0 0 0 1 0 0 0 0 0 0 1.000000")
    writeLines(c(line_one, line_two, line_three, line_four, line_seven), fileconn)
    close(fileconn)
    system(paste0(lsdynadir,"\\lsprepost4.3_x64.exe c=get_stresses.cfile -nographics"))
    
    setwd(paste0(solutions_dir, "\\", sub_directory))
    
    stresses <- readLines(kfile)
    start <- grep(stresses, pattern="*KEYWORD",fixed = T)
    stop <- grep(stresses, pattern="$Interpreted from averaged nodal data",fixed = T)
    stresses <- stresses[-seq(start, stop, by = 1)]
    writeLines(stresses, kfile)
    stresses <- read.table(kfile, header = FALSE)
    names(stresses) <- c("eid", "stress")
    #max_stress <- which(stresses$eid == which.max(stresses$stress))
    max_stress <- stresses[which(stresses$stress == max(stresses$stress)), ]
    #second_stress <- subset(stresses,sort(z<-rank(stress),T)[2]==z)[,1]
    print(j)
    
    
    biggest_stresses <- rbind(biggest_stresses, max_stress)
    #biggest_stresses[nrow(biggest_stresses)+1,] = list(stresses[max_stress,1], stresses[max_stress,2])
    #second_stresses <- rbind(second_stresses, stresses[second_stress,])
    
  }
  
  #return(biggest_stresses)
  max_stresses<- (biggest_stresses[which.max(biggest_stresses$stress),1:2])
  write.csv(x = max_stresses, file = paste0(stats_list[i],"_max_stress"), row.names = FALSE)
  return(max_stresses)
  #return(biggest_stresses)
  #return(second_stresses[which.max(second_stresses$stress),1])
  
  
}

element_time_series <- function(stressed_eid, solutions_dir, kfile, lsdynadir) {
  ###once you know the EID that has max stress, this will grab the time-series table for that element
  ###stressed_eid <- EID number you determined
  ###input_dir - whereever the solution is stored
  ###output_dir - where you want the time-series table saved
  ###lsdynadir - lspp .exe file location
  
  fileConn <- file("get_stress_plots.cfile")
  line_zero <- paste0("open d3plot ", solutions_dir, "\\", kfile, "\\", "d3plot")
  line_one <- paste0('genselect element add solid ', stressed_eid,'/', 'F1/0')
  line_three <- paste0("xyplot 1 savefile xypair ", solutions_dir, "\\", kfile, "\\", "file_xy ", "1 all")
  
  writeLines(c(line_zero, "genselect target element", line_one, 'etime 14', line_three, "xyplot 1 donemenu"), fileConn)
  close(fileConn)
  
  system(paste0(lsdynadir,"\\lsprepost4.3_x64.exe c=get_stress_plots.cfile -nographics"))
  
} 

make_small_ram <- function(kfile, x1, y1, x2, y2, sailheight, concave_iterations) {
  
  ##kfile should be the flist vector
  
  knodes <- readk(kfile[i], type = "nodes") #all the nodes in the mesh
  noramside <- data.frame(isolate_a_ram(knodes = knodes, x1 = x1, y1 = y1, x2 = x2, y2 =  y2)) #find nodes NOT on the edge of this inst
  names(noramside) <- "nid"
  names(knodes)[1] <- "nid"
  knodes$dontkeep <-knodes$nid %in% noramside$nid #flag nonramside for deletion
  knodes$dontkeep[knodes$z >= sailheight & knodes$dontkeep == 'FALSE'] <- "TRUE" #flag ramside points >== 90 m for deletion
  nodes_to_delete <- knodes[which(knodes$dontkeep == 'TRUE'),] #now subset the ones that need to be deleted.
  
  
  return(concaved_nodes <- nick_concave(nodes_to_delete, concave_iterations, flist[i])) # draw a 20 m concave hull, and delete those nodes.
  
}

erode_damp <- function(kfile_mesh, damp_coef) {
  ### pastes the damp coeff. back into the kfile in the appropriate spot
  ### erosion_lines - character string of erosion card for kfile
  ### damp_coef - damping coefficient
  ### kfile_mesh - the mesh you're working on
  
  
  kmesh <- readLines(kfile_mesh)
  #line1 = grep(kmesh, pattern="*DEFINE_FUNCTION",fixed = T)
  #kmesh = c(kmesh[1:line1-1], erosion_lines, kmesh[line1:length(kmesh)])
  
  damped <- gsub("uuu", damp_coef, damp)
  line1 = grep(kmesh, pattern="*ELEMENT_SOLID",fixed = T)
  kmesh = c(kmesh[1:line1-1], damped, kmesh[line1:length(kmesh)])
  
  writeLines(kmesh, paste0(kfile_mesh))
  
}

mass_calculator <- function(inst, suffix, thickness) {
  #calculates the mass of an ice island taking into account the ram and and all that
  #this is to try and find set damping coefficients for different mass ranges
  
  x_y <- read.table(paste0(inst, suffix), sep = ",", dec = ".", comment.char = "*", skip = 3)
  names(x_y) <- c("1", "x", "y", "4")
  keel_area <- abs(polyarea(x_y$x, x_y$y))
  
  volume_keel <- keel_area * thickness
  volume_sail <- (ft_suspects_df[,"poly_area"][which(ft_suspects_df$inst == inst)] * 1e6) * thickness 
  volume <- volume_keel + volume_sail
  
  mass <- (volume * 900)/ 1e9
  print(mass)
  
}
ci2d3
predict_damp <- function(inst, buff_length, thickness, ft_suspect_df) {
  
  v1 <- c(89.75714, 62.83, 62.83, 25.137, 15.7075, 8.975714, 5.98, 4.33)
  v2 <- c(717737.5, 5332638, 29553931, 1.61E+08, 4.47E+08, 9.41E+08, 2.4E+09, 3.6E+09)
  df <- data.frame(v1, v2)
  names(df) <- c("damp", "volume")
  #plot the data
  #plot(df$damp~df$volume)
  #plot(fit)
  
  #fit log model
  fit <- lm(damp~log(volume), data = df)
  #Results of the model
  #summary(fit)
  
  x_y <- read.table(paste0(inst, "_", buff_length, ".k"), sep = ",", dec = ".", comment.char = "*", skip = 3)
  #x_y <- readk(paste0(inst, "_", buff_length, ".k"), 'nodes')
  names(x_y) <- c("1", "x", "y", "4")
  keel_area <- abs(polyarea(x_y$x, x_y$y))
  
  volume_keel <- keel_area * thickness
  volume_sail <- (ft_suspects_df[,"poly_area"][which(ft_suspects_df$inst == inst)] * 1e6) * thickness 
  volume <- volume_keel + volume_sail
  
  volume <-data.frame(volume=volume)
  damping_coefficient <- predict(fit, newdata = volume)
  damping_coefficient <- signif(damping_coefficient,digits=3)
  
  
  return(damping_coefficient)
  
  
}

volume_finder <- function(inst, buff_length, thickness) {
  
  v1 <- c(89.75714, 62.83, 62.83, 25.137, 15.7075, 8.975714)
  v2 <- c(717737.5, 5332638, 29553931, 1.61E+08, 4.47E+08, 9.41E+08)
  df <- data.frame(v1, v2)
  names(df) <- c("damp", "volume")
  #plot the data
  #plot(df$damp~df$volume)
  #plot(fit)
  
  #fit log model
  fit <- lm(damp~log(volume), data = df)
  #Results of the model
  summary(fit)
  
  x_y <- read.table(paste0(inst, "_", buff_length, ".k"), sep = ",", dec = ".", comment.char = "*", skip = 3)
  names(x_y) <- c("1", "x", "y", "4")
  keel_area <- abs(polyarea(x_y$x, x_y$y))
  
  volume_keel <- keel_area * thickness
  volume_sail <- (ft_suspects_df[,"poly_area"][which(ft_suspects_df$inst == inst)] * 1e6) * thickness 
  volume <- volume_keel + volume_sail
  
  
  return(volume)
  
  
}

find_waterline <- function(ft_suspects_df, t_sail, t_keel, vertex_kfile) {
  
  x_y <- read.table(vertex_kfile, sep = ",", dec = ".", comment.char = "*", skip = 3)
  names(x_y) <- c("1", "x", "y", "4")
  keel_area <- abs(polyarea(x_y$x, x_y$y))
  
  sail_area <- ft_suspects_df$poly_area[i] * 1e6
  v_sail <- sail_area * t_sail 
  
  v_keel <- (keel_area) * t_keel 
  v <- v_sail + v_keel
  m_ice <- v * rho_ice
  w_ice <- g * m_ice
  w_water <- w_ice
  m_water <- w_water/g 
  v_water <- m_water/rho_water
  
  if (v_water > v_keel) {
    
    print(((v_water - v_keel)/v_sail)*t_sail + t_keel)  # from the bottom of the keel
    
  } else {
    
    print ((v_water/v_keel)*t_keel)  # from the bottom of the keel
    
  }
  
}

paste_waterline <- function(vertex_directory, ci2d3_df, t_sail, t_keel, vertex_kfile, kfile_input_directory, kfile_output_directory, kfile_mesh) {
  
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

make_small_ram <- function(kfile, x1, y1, x2, y2, sailheight, concave_iterations) {
  
  ##kfile should be the flist vector
  
  knodes <- readk(kfile_mesh = kfile, type = "nodes") #all the nodes in the mesh
  noramside <- data.frame(isolate_a_ram(knodes = knodes, x1 = x1, y1 = y1, x2 = x2, y2 =y2)) #find nodes NOT on the edge of this inst
  names(noramside) <- "nid"
  names(knodes)[1] <- "nid"
  knodes$dontkeep <-knodes$nid %in% noramside$nid #flag nonramside for deletion
  knodes$dontkeep[knodes$z >= sailheight & knodes$dontkeep == 'FALSE'] <- "TRUE" #flag ramside points >== 90 m for deletion
  nodes_to_delete <- knodes[which(knodes$dontkeep == 'TRUE'),] #now subset the ones that need to be deleted.
  
  
  return(concaved_nodes <- nick_concave(nodes_to_delete, concave_iterations, kfile_mesh = kfile)) # draw a 20 m concave hull, and delete those nodes.
  
}

isolate_a_ram <- function(knodes, x1,y1,x2,y2) {
  #k-nodes - a dataframe with node data
  #x1,y1
  #x2,y2 - two points that form a line in lcc coords
  
  y <- c(y1,y2)
  x <- c(x1,x2)
  ln <- lm(y~x)
  side = knodes$y - predict.lm(ln, newdata = data.frame(x=knodes$x)) 
  if(length(which(side<0)) > length(which(side>0)) ) {
    ramside = which(side>0)
    noramside = which(side<0)
    noramside = c(noramside,which(side==0))
  } else {
    ramside = which(side<0)
    noramside = which(side>0)
    noramside = c(noramside,which(side==0))
  }
  return(knodes$`$#   nid`[noramside])
  #return(knodes$keep <- knodes$`$#   nid` %in% ramside)
}

delete_lines <- function(kfile) { 
  
  file <- readLines(kfile)
  #start = grep(file, pattern="*SET_SEGMENT",fixed = T)
  #stop = grep(file, pattern="*NODE",fixed = T)
  #my_new_file = file[-(start:stop-1)]
  
  start2 = grep(file, pattern="*DAMPING_GLOBAL",fixed = T)
  stop2 = grep(file, pattern="*ELEMENT_SOLID",fixed = T)
  my_newer_file = file[-(start2:(stop2-1))]
  
  filecon = file(kfile)
  writeLines(my_newer_file, filecon)
  close(filecon)
  
}
#### for loop/function to do sensitivity analysis/collect stats on each mesh

step_clean_model <- function(mesh_dir, kfile_mesh) {
  
  fileConn<-file("step_cleaner.cfile")
  
  one <- paste0('open keyword \"', mesh_dir, '\\',kfile_mesh,'\"')      
  two <- paste0('save keywordoutversion 8')
  three <- paste0('modelcheck checkgeneral')
  four <- paste0('selectentity getalllabel')
  five <- paste0('selectentity getalllabel')
  six <- paste0('modelcheck autocleaninfo')
  seven <- paste0('modelcheck applyclean')
  eight <- paste0('selectentity getalllabel')
  nine <- paste0('selectentity getalllabel')
  ten <- paste0('genselect clear')
  eleven <- paste0('save keywordabsolute 0')
  twelve <- paste0('keywordbylongfmt 0')
  fourteen <- paste0('save outversion 7')
  fifteen <- paste0('save keyword \"', mesh_dir, '\\',kfile_mesh,'\"')
  fourteen <- paste0('save outversion 7')
  sixteen <- paste0('selectentity getalllabel')
  seventeen <- paste0('genselect clear')
  
  
  writeLines(c(one, two, three, four, five, six, seven, eight, nine, ten, eleven, twelve, fourteen, fifteen, sixteen, seventeen), fileConn)
  close(fileConn)
  system(paste0(lsdynadir,"\\lsprepost4.3_x64.exe c=step_cleaner.cfile -nographics"))
}

make_ram_picks <- function(kfile_mesh, iterations) {
  
  setwd(scratchdrivemesh)
  nodes <- readk(kfile_mesh, type = 'nodes')
  colnames(nodes)[1] <- "nid"
  nodes <- nick_concave(nodes, iterations = iterations, kfilemesh = kfile_mesh)
  

  setwd("C:/Users/jsmith/Desktop/make_ram_picks")
  write.csv(nodes[1:3], file = paste0(kfile_mesh, "_nodes.csv"), row.names = 'FALSE')
  
}

get_extent_coords_of_iso_rams <- function(kfiles_complete_isolated_dir, kfile, top_elev, kfile_vertex_iso_dir) {
  
## gets the extent coords of a mesh with iso ram to get volume to find damping coefficient.
  
  setwd(kfiles_complete_isolated_dir)
  dt <- nodeser(kfile, top_elev)
  
  DT_sf = sf::st_as_sf(dt, coords = c("x", "y"),agr = "constant")
  point_coords <- as.data.frame(st_coordinates(st_as_sf(DT_sf)))
  names(point_coords) <- c("x", "y")
  
  c_hull <- concaveman(DT_sf)
  c_hull_coords <- st_coordinates(c_hull$polygons)
  c_hull_coords <- as.data.frame(c_hull_coords[, -c(3:4)])
  names(c_hull_coords) <- c("x", "y")
  c_hull_coords$z <- 0
  setwd(kfile_vertex_iso_dir)
  write.table(c_hull_coords, 'temp.csv', row.names = TRUE, col.names = FALSE, sep=",", quote=FALSE)
  
  # read that file back in as one block of text and add the header + footer and overwrite it ##removed line breaks where necessary.
  vertices = readLines('temp.csv')
  header = '*Keyword\n*Node\n$ Node,X,Y,Z'
  footer = '*End'
  outfile = file(kfile)
  writeLines(c(header, vertices, footer), outfile)
  close(outfile)
  
}

paste_flex_strength <- function(kfile_mesh, flex_strength, kfile_input_directory, kfile_output_directory) {
  
  
  setwd(kfile_input_directory)
  kmesh <- readLines(kfile_mesh)
  line1 = grep(kmesh, pattern = "       0.0  1900000.0       0.0       0.0       0.0       0.0       0.0       0.0",fixed = T)
  flex_line <- gsub("1900000.0", format(flex_strength, scientific=F), kmesh[line1])
  kmesh_final = c(kmesh[1:line1-1], flex_line, kmesh[line1+1:length(kmesh)])
  end = grep(kmesh, pattern = "*END",fixed = T)
  
  setwd(kfile_output_directory)
  print(getwd())
  writeLines(kmesh_final[1:end], kfile_mesh)
  
}

paste_term_time <- function(kfile_mesh, term_time, kfile_input_directory, kfile_output_directory) {
  
  
  setwd(kfile_input_directory)
  kmesh <- readLines(kfile_mesh)
  line1 = grep(kmesh, pattern = "      60.0         0       0.0       0.01.000000E8",fixed = T)
  time_line <- gsub("60.0", format(term_time, scientific=F), kmesh[line1])
  kmesh_final = c(kmesh[1:line1-1], time_line, kmesh[line1+1:length(kmesh)])
  end = grep(kmesh, pattern = "*END",fixed = T)
  
  setwd(kfile_output_directory)
  print(getwd())
  writeLines(kmesh_final[1:end], kfile_mesh)
  
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))


#given an EID and a radius, find the stress within that distance from an element
# Arguments:
## k - a valid k file in the current folder
## eid - a valid EID number from the mesh                                         
## radius - the distance in mesh units (m) to seek for other elements

# Output: 
## a dataframe with 2 columns: EID and stress

#Pseudocode
k = "10_10_10.k"
eid = 10
radius = 75 

get_element_buffer_stress <- function(k, eid, radius) {

# Find all the nodes in the element [x]
kfile <- readLines(k)
start <- grep(kfile, pattern="*ELEMENT_SOLID",fixed = T)
stop <- grep(kfile, pattern = "*NODE")
elementdf <- read.table(k, header = FALSE, sep="", skip = start+1, nrows = stop-2-start, fill = TRUE)
elementnodes <- elementdf[which(elementdf$V1 == eid),-c(1,2)]


# Find all the coords of all these nodes
start_coord <- grep(kfile, pattern="*NODE",fixed = T)
stop_coord <- grep(kfile, pattern = "*END")
coords <- read.table(k, header = FALSE, sep="", skip = start_coord+1, nrows = stop_coord-2-start_coord, fill = TRUE)
names(coords) <- c("nid", "x", "y","z")
eid_coords <-coords[which(coords$nid %in% elementnodes),-c(1,5,6)]

# Find the centroid - average all the x, y and z to get one point
centroid <- c(mean(eid_coords$x), mean(eid_coords$y), mean(eid_coords$z))

# Find all nodes within the search volume
ind = which(coords$x >= centroid[1] - radius & coords$x <= centroid[1] + radius &
            coords$y >= centroid[2] - radius & coords$y <= centroid[2] + radius &
            coords$z >= centroid[3] - radius & coords$z <= centroid[3] + radius)

subsetnodes <- coords$nid[ind]

# Determine which EIDs they are 
subsetelements <- elementdf[apply(elementdf[,3:10], 1, function(x) all(x %in% subsetnodes)),1]
# Find stress in each of those
# Output
}


stresses <- get_stress_table("J:/SCRATCH/jsmith/Solutions_2/sens_test_edit/", "40_smaller", lsdynadir = lsdynadir, states = 5)



library()
# Begin Code here ####

# Connect to database
drv <- dbDriver('PostgreSQL')

# Create a connection to the postgres database
con <- dbConnect(drv, dbname = "estewart",
                 host = "pollux-new.wirl.carleton.ca", port = 5432,
                 user = "estewart", password = 'dbpass') 

# -----------------------------------------------------------------------------------------------------------------------------

### Instances before and after fracturing ###
table_f <- f.fract(con, calvingyr = NULL, calvingloc = NULL, wk_num = NULL)

table_bf <- f.fract_bf(con, calvingyr = NULL, calvingloc = NULL) 
table_af <- f.fract_af(con, calvingyr = NULL, calvingloc = NULL) 

### Create spatial polygons datframe
#df_sp_f <- f.Spatialdf(table_f)
df_sp_bf <- f.Spatialdf(table_bf)
df_sp_af <- f.Spatialdf(table_af)



### Make a text string for damp coeff.
#####

setwd(kfile_complete)
cat("*DAMPING_GLOBAL",
    "$#    lcid    valdmp       stx       sty       stz       srx       sry       srz",
    "        0       uuu       0.0       0.0       0.0       0.0       0.0       0.0", sep="\n", file = "damping_coef.txt", append=TRUE)
damp <- readLines("damping_coef.txt")
###Make a text string for element erosion data
####
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

setwd(scratchdrivemesh)
mat <- readLines("mat_paste.txt")


##Paste mat properties at the appropriate lines in the kfile meshes:

setwd(scratchdrivemesh)
for (i in 16:length(flist)) {
  mat_paster(kfilemesh = flist[i], mat)
}

## Export x/y tables of top nodes and their plots as PNGs/JPEGs to new directories for fracture references

setwd(topviewdir_bf)
for (i in 1:length(df_sp_bf$inst)) {
  # i is the index of bf polygons
  # j is the index of associated af (daughter) polygons
  bf <-  df_sp_bf$inst[i]
  j <- which(df_sp_af$motherinst == bf)
  j <- j[order(df_sp_af$poly_area[j],decreasing = T)]  # order by size
  
  cols <- length(j)+1
  if (cols > 4) { 
    print(bf)
    next
    }
  png(paste0(bf,'.png'))
  par(mfrow=c(1,cols))
  x_size <- bbox(df_sp_bf[i,])[1,2]-bbox(df_sp_bf[i,])[1,1]
  y_size <- bbox(df_sp_bf[i,])[2,2]-bbox(df_sp_bf[i,])[2,1]
  size <- max(x_size,y_size)
  x_lim <- c(df_sp_af$centroid_x[af]-size/2, df_sp_af$centroid_x[af]+size/2) 
  y_lim <- c(df_sp_af$centroid_y[af]-size/2, df_sp_af$centroid_y[af]+size/2)
  plot(df_sp_bf[i,], col='darkgrey')
  mtext(bf, 1, cex=.5)
  for (k in 1:(cols-1)) {
    af <- j[k]
    x_lim <- c(df_sp_af$centroid_x[af]-size/2, df_sp_af$centroid_x[af]+size/2) 
    y_lim <- c(df_sp_af$centroid_y[af]-size/2, df_sp_af$centroid_y[af]+size/2)
    plot(df_sp_af[af,], xlim=x_lim, ylim=y_lim, col='lightgrey')
    mtext(df_sp_af$inst[af], 1, cex=.5)
    } # end plotting daughters
dev.off()
}  #end for loop  

par(mfrow=c(1,1))

install.packages("pillar", type = 'binary')
library(sf)

####################
## make a ram###
#####################


###Buld the meshes

for (i in 1:length(flist)) {
  setwd(kfile_vertex)
  vertex_kfile(ft_suspects_sp, i ,buff)
  mesh_kfile(kfile_in = flist[i], kfile_vertex, scratchdrivemesh,lsdynadir,mesh_size,z_length, zelem)
  segment_setter(flist[i], scratchdrivemesh, scratchdrivemesh,lsdynadir)
  setwd(scratchdrivemesh)
  mat_paster(flist[i], mat, input_dir = scratchdrivemesh, output_dir = scratchdrivemesh) #####
  toplevel = nodeser(kfilemesh = flist[i], elev = ram_ht) ##intermediate step here where only the 'right' nodes are supplied for isolated ramget
  ramcut = nick_concave(toplevel, iterations = 1, flist[i])
  setwd(scratchdrivemesh) ##start here
  elements_deleted <- delete_elements_improved(kfilemesh = flist[i], new_nodes = ramcut, ht_cutoff = ram_ht-10) #
  orphans_departed <- goodbye_orphans(new_nodes = ramcut, elements_deleted)
  nodes_width_fixed <- fix_width_nodes(nidwidth, nidcols, orphans_departed)
  elements_width_fixed <- fix_width_elements(eidwidth, eidcols,elements_deleted)
  final_kfile <- nid_eid_replace(kfilemesh = flist[i], nidtable = "temp_fwf_nodes.txt", eidtable = "temp_fwf_elements.txt", inst_names = inst_names_meshes[i], buff = buff, output_dir = kfile_complete_100)
  paste_waterline(vertex_directory = kfile_vertex, ci2d3_df = ft_suspects_sp, t_sail = 10, t_keel = 90, vertex_kfile = flist[i], kfile_input_directory = kfile_complete_100, kfile_output_directory = kfile_complete_100, kfile_mesh = flist[i])
  
}



####if you mess up the above, run the code below and try again

for (i in 14:length(flist)) {
  
  delete_lines(kfile = flist[i])
  
}

###if you end up with a bunch of NAs at the bottom of the file.....
for (i in 1:length(flist_meshes)) {
  
  bad_file <- readLines(flist_meshes[i])
  good_file <- bad_file[-grep("^NA$", bad_file)]
  writeLines(good_file, con = flist_meshes[i])
  
}


### make a place to place solutions###
solutions_list<- paste0(inst_names_meshes, "_20.k") 
solutions_list <- paste0(inst_names_meshes, "_40.k") 
solutions_list <- paste0(inst_names_meshes, "_60.k") 
solutions_list <- paste0(inst_names_meshes, "_80.k") 

flist <- solutions_list
setwd(solutions_dir_no_max_100)
for (i in 1:length(flist)){
  dir.create(path = flist[i])
}

setwd(solutions_dir2)
for (i in 1:length(flist_meshes))
{
 
  unlink(flist_meshes[i], recursive = TRUE)
}

setwd(kfiles_complete_isolated)

for (i in 1:length(flist)) {
  
  delete_lines(kfile = flist[i])
  
}


delete_lines <- function(kfile) { 
  
  file <- readLines(kfile)
  start = grep(file, pattern="*SET_SEGMENT",fixed = T)
  stop = grep(file, pattern="*ELEMENT_SOLID",fixed = T)
  my_new_file = file[-(start:stop-1)]
  
  #start2 = grep(my_new_file, pattern="*DAMPING_GLOBAL",fixed = T)
  #stop2 = grep(my_new_file, pattern="*ELEMENT_SOLID",fixed = T)
  #my_newer_file = my_new_file[-(start2:(stop2-1))]
  
  filecon = file(kfile)
  writeLines(my_new_file, filecon)
  close(filecon)
  
}

flist <- paste(inst_names_meshes, "_60.k")

for (i in 1:length(flist)) {
  
  paste_flex_strength(flist[i], 600000.0, kfiles_complete_isolated, kfiles_complete_isolated)
  
}

## perhaps this for loop will give you max stresses for each element in the limitless yield strength models....
##this funct/forloop works. Try to make it so that it begins checking at states near end of sim

stats_list <- flist[c(2, 3, 9, 13, 15, 16, 18, 19, 20)]
max_stresses <- data.frame(eid= numeric(),
                               stress = numeric(),
                               stringsAsFactors=FALSE) 
setwd(kfile_complete)
for (i in 9:length(stats_list)) {
 
index <- get_max_stress_eid(solutions_dir, sub_directory =  stats_list[i], kfile = stats_list[i], lsdynadir = lsdynadir, states = 1501)
max_stresses <- rbind(max_stresses, index[1,])
}

write.csv(max_stresses, file = "max_sym_stresses.csv", row.names = FALSE)

### Below is a for loop to get solutions. 
setwd(paste0(solutions_dir, "\\", flist[1])) 


for (i in 1:length(flist)) {
  
  setwd(paste0(solutions_dir_no_max_100, "\\", flist[i])) 
  
  shell(paste0(lsdyna_solver,"\\ls-dyna_smp_d_R10.0_winx64_ifort160.exe I=", kfiles_complete_isolated_100, "\\", flist[i], " O=",solutions_dir_no_max_100, "\\", flist[i], "\\solved   NCPU =  8 MEMORY = 200000000"))
  shell("exit")
  

}

### for non-isolated rams (recall that your workspace is now saved in solutions_dir):
for (i in 3:length(stats_list)) {
  
  setwd(paste0(solutions_dir, "\\", stats_list[4])) 
  
  shell(paste0(lsdyna_solver,"\\ls-dyna_smp_d_R10.0_winx64_ifort160.exe I=", kfile_complete, "\\", stats_list[4], " O=",solutions_dir, "\\", stats_list[4], "\\solved   NCPU =  8 MEMORY = 200000000"))
  shell("exit")
  
}

#for loops to add delete old data and add new data to kfiles and run solution #2
#######################################################################
for (i in 1:length(flist)) {
  
  step_clean_model(mesh_dir = kfiles_complete_isolated, kfile_mesh = flist[i])
  
  
}

for (i in 13:length(flist[13]))
  {
  
  setwd(kfile_vertex_iso) 
  damp_coef <- predict_damp(inst = inst_names_meshes[i], buff_length = 60, thickness = 100, ft_suspect_df = ft_suspects_sp)
  
  setwd(kfiles_complete_isolated_100)
  erode_damp(kfile_mesh = flist[i], damp_coef = damp_coef)
  
}

for (i in 1:length(flist))
{
  
  setwd(kfile_vertex_iso) 
  damp_coef <- predict_damp(inst = inst_names_meshes[i], buff_length = 60, thickness = 100, ft_suspect_df = ft_suspects_sp)
  
  print(damp_coef)
  
}

####################################################
######################################
## import table of footloose suspects:

setwd("C:/Users/jsmith/Desktop")
ft_suspects <- read.csv("footloose_suspects.csv")
library(rgeos)
colnames(ft_suspects)[1] <- "Inst"
ft_suspects[,2:4] <- NULL
ft_suspects$X.2 <- NULL
ft_suspects_df <- dplyr::semi_join(x = table_bf, y = ft_suspects, by = c("inst" = "Inst"))
ft_suspects_sp <- f.Spatialdf(ft_suspects_df)
save.image(file='myEnvironment.RData')
flist <- paste0(ft_suspects_sp$inst, "_50.k")
setwd(scratchdrivemesh)
save.image(file='myEnvironment.RData')

ft_suspects_df <- ft_suspects_df[ft_suspects_df$inst %in% inst_names_meshes,]

flist <- paste0(ft_suspects_df$inst, "_50.k") 



###run code in parallel?
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

foreach(i=1:length(ft_suspects_sp)) %dopar% {
  
  for (i in 1:length(ft_suspects_sp)) {
    setwd(kfile_vertex)
    vertex_kfile(ft_suspects_sp, i ,buff)
    mesh_kfile(kfile_in = flist[i], kfile_vertex, scratchdrivemesh,lsdynadir,mesh_size,z_length, z_elem)
    segment_setter(flist[i], scratchdrivemesh, scratchdrivemesh,lsdynadir)
    setwd(scratchdrivemesh)
    mat_paster(flist[i], mat)
    toplevel = nodeser(kfilemesh = flist[i], elev = ram_ht)
    ramcut = nick_concave(toplevel, 0, flist[i])
    setwd(scratchdrivemesh)
    elements_deleted <- delete_elements_improved(kfilemesh = flist[i], new_nodes = ramcut, ht_cutoff = ram_ht - mesh_size)
    orphans_departed <- goodbye_orphans(new_nodes = ramcut, elements_deleted)
    nodes_width_fixed <- fix_width_nodes(nidwidth, nidcols, orphans_departed)
    elements_width_fixed <- fix_width_elements(eidwidth, eidcols,elements_deleted)
    final_kfile <- nid_eid_replace(kfilemesh = flist[i], nidtable = "temp_fwf_nodes.txt", eidtable = "temp_fwf_elements.txt")
    
  }
}
  
#stop cluster
stopCluster(cl)

for (i in 1:length(flist_meshes)) {

bad_file <- readLines(flist_meshes[i])
good_file <- bad_file[-grep("^NA$", bad_file)]
writeLines(good_file, con = flist_meshes[i])

}



##### get my mahfuckin segment sets back beeyotch

for (i in 1:length(flist)){
  

segment_setter(flist[i], kfile_complete2, kfile_complete2,lsdynadir)

}


#igraph stuff <-
required_df <- table_f[table_f$motherinst %in% ft_suspects_sp$inst,]
required_df <- as.data.frame(required_df)
required_df_test <- required_df[,c(12,11,8,13,20)]
required_df_test <- required_df[,c(12,11, 8)]
g <- graph_from_data_frame(required_df)

required_df_test <- rbind(required_df_test, grandmothers)
required_df_test <- required_df_test[!duplicated(required_df_test$inst), ]

pre_edgelist <- data.frame(table_f$motherinst, table_f$inst)
colnames(pre_edgelist) <- c("motherinst", "inst")
edgelist_df <- as.data.frame(pre_edgelist[pre_edgelist$motherinst %in% ft_suspects_sp$inst,]) ### This code matches the inst and motherinst together [google "value matching in r"]
#orphans = df[-which(df$motherinst %in% df$inst),]  #motherless instances
#orphans_1generation = which(edgelist_df$inst %in% orphans$inst)

edgelist_df <- as.data.frame(na.omit(edgelist_df))  
edgelist_df <- cbind(edgelist_df, required_df$poly_area)
names(edgelist_df)[3] <- "poly_area"
names(grandmothers) <- c("motherinst", "inst", "poly_area")
edgelist_df_test <- rbind(edgelist_df, grandmothers)

g <- graph_from_data_frame(edgelist_df_test, vertices = required_df_test, directed = T)
g <-simplify(g, remove.multiple = F, remove.loops = T)

V(g)$color <- randomColor(72, hue = "random")

colors <- colors[as.numeric(as.factor(vertex_attr(g, "motherinst")))]
V(g)$color <- colors[V(g)$motherinst]
V(g)$color <

# create color pallet based on unique values for vertex attribute
colors <- randomColor(23, hue="random")[f]
vcols <- 
plot(g, vertex.color = color[as.numeric(as.factor(vertex_attr(g, "motherinst")))])
colors <- colors[as.numeric(as.factor(vertex_attr(g, "motherinst")))]
V(g)$color <- colors[as.numeric(as.factor(vertex_attr(g, "motherinst")))] ##this is the key to gettin labels to match plot
V(g)$inst <- stri_sub(V(g)$, -3,-1)
              
V(g)$color
             
pt.cex = 1.22  
setwd("C:/Users/jsmith/Desktop")
tiff("igraph_ice_islands_alt_labs.tif", units="in", width=14, height=10, res=300)
plot(g, vertex.size = (V(g)$poly_area + 6.65), arrow.size = 0.5, arrow.width = 0.1, vertex.label = stri_sub(V(g)$name, -3,-1), vertex.label.font= 2, vertex.label.color= "black", vertex.label.dist = 0, vertex.label.cex = 1.22, vertex.frame.color = "black", edge.color = "black", alpha = 2, curved = TRUE, edge.arrow.size=0.08, edge.arrow.width=0.009, directed = T)
legend(1.4, 1.2, legend = stri_sub(levels(f), -3,-1), col = colors, pch = 16, bty = "n", title = "MOTHERS", pt.cex = 1.75, cex = 1.75) 
dev.off()

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

substrRight(x, 6)

