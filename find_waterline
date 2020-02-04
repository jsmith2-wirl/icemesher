
find_waterline <- function(ft_suspects_df, t_sail, t_keel, vertex_kfile) {
#ft_suspects_df - spatial data df of ice islands
#t_sail - thickness of the freeboard (m)
#t_keel - thickness of the ram (m)
#vertex_kfile - 2-D vertices-only mesh K-file
  
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
    
  }polyarea(DF$Strain, DF$Stress)
  
}
