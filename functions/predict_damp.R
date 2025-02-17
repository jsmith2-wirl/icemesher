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
