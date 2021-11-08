mat_paster = function(kfilemesh, mat, input_dir, output_dir) {
  # pastes the material property lines to the appropriate section of each k file
  # kfilemesh - k file (that has already been given a segment set)
  # mat - read-in file of material properties/see below this function code
  # outputs a k file with material properties
  
  setwd(input_dir)
  kmesh <- readLines(kfilemesh)
  kmesh = c(kmesh[1:6], mat, kmesh[12:length(kmesh)])
  setwd(output_dir)
  writeLines(kmesh, kfilemesh)
  
} 


##MAT CARD BELOW ##
cat("*CONTROL_ENERGY",
    "$#    hgen      rwen    slnten     rylen   ",
    "         2         2         1         1",
    "*CONTROL_TERMINATION",
    "$#  endtim    endcyc     dtmin    endeng    endmas   ",
    "      60.0         0       0.0       0.01.000000E8", 
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
    "         1         1         0",
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
    "               100.0                9.81", sep="\n", file="mat_paste.txt", append=TRUE)
