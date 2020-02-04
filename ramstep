ramstep <- function(dt) {
  # supply a data table with nodes coords x and y and outer concave hull will be removed. This function is 
  # embedded within "nick_cocave"
  # dt - datatable of node coords
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
