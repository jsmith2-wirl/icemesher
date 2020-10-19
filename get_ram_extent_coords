get_ram_extent_coords <- function(kfiles_complete_isolated_dir, kfile, top_elev, kfile_vertex_iso_dir) {
  
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
