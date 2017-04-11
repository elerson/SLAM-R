#!/usr/bin/r

#.libPaths(c("/home/lg/anaconda3/envs/my-r-env/lib/R/library", .libPaths()))
#install.packages("jpeg", repos="http://cran.rstudio.com/")

#print(.libPaths())
library("plyr")
library("matrixStats")

##Parameters
#load("database_real.RData")
#intialize particles and other variables


##
##  Simulation or Real data selection
simulation = TRUE 
if(simulation){
  n_particles = 10
  n_perception_particles = 50
   
  load("slam_database.RData")
  alphas = c(2,  16, 0.5,  8)#3,  10, 2,  8 #2,10,2,10 ## motion model variables
  seq_step = 1

  gamma=3 # 3 -- **60
  map_particles_weights=20
  lookup_resolution_param = 1
} else{
  n_particles = 40
  n_perception_particles = 200

  load("database_real.RData")
  alphas = c(3/10,  1/10, 2/10,  8/10)#3,  10, 2,  8 #2,10,2,10 ## motion model variables
  seq_step = 5

  gamma=10# 3 -- **60
  map_particles_weights=30
  lookup_resolution_param = 10
}
##Map variables

l0 = 0.5
loc = 0.0
lfree = 0.6

resolution = database_info$map_metadata$resolution
beta = 0.7

#############################################################################################
##				Motion Model
#############################################################################################


sample_motion_model <- function(particles, xt_curr, xt_prev, alpha, resolution){ # xt, xt' -> [x, y, theta]
  num_particles = nrow(particles)
  #calcula theta1 -> angulo de rotacao no movimento do robo
  #distance -> distancia do movimento
  #theta2 -> restante do angulo do movimento
  theta1 = atan2(xt_curr[2]-xt_prev[2], xt_curr[1]-xt_prev[1]) - xt_prev[3]
  distance = sqrt((xt_curr[2]-xt_prev[2])**2 + (xt_curr[1]-xt_prev[1])**2)
  theta2 = xt_curr[3]  - theta1  - xt_prev[3]

  #generate samples for theta1, theta2 and distance
  theta1_sample = theta1 - rnorm(num_particles, 0, alpha[1]*(abs(theta1)) + alpha[2]*(distance))
  distance_sample = distance - rnorm(num_particles, 0, alpha[4]*(distance) + alpha[3]*(abs(theta1) + abs(theta2)))
  theta2_sample = theta2 - rnorm(num_particles, 0, alpha[1]*abs(theta2) + alpha[2]*(distance))
 
  #update the particles pose
  particles[,"x"] = particles[,"x"] + (distance_sample*cos(particles[,"theta"]+theta1_sample))/resolution
  particles[,"y"] = particles[,"y"] + (distance_sample*sin(particles[,"theta"]+theta1_sample))/resolution
  particles[,"theta"] = particles[,"theta"] + theta1_sample + theta2_sample

  return (particles)

}

#############################################################################################
##				Perception Model
#############################################################################################


###
### BINARY SCAN MATCHING SENSOR MODEL
###


scan_match <- function(particle, particle_map, local_map_update, map_aux_data){
  

  global_xy = c(particle_map[[3]][1] + (particle["x"] - particle_map[[2]][1]), particle_map[[3]][2] + (-particle["y"] + particle_map[[2]][2]))
  ##
  ##  ROTATE THE LOCAL MAP
  ##
  theta <- particle["theta"]
  rot <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)), 2, 2)
  ##matrix rotation, faster than using %*%
  xy_loc <- local_map_update[[2]]-dim(map_aux_data[[1]])[1]/2

  xy_loc_ <-  (cbind(c(global_xy[1]+xy_loc[,1]*rot[1,1]+xy_loc[,2]*rot[2,1]),c(global_xy[2]+xy_loc[,1]*rot[1,2]+xy_loc[,2]*rot[2,2]))) 
  
  return ((exp(gamma*sum(particle_map[[4]][xy_loc_/particle_map[[5]]])/1080)))

}



#############################################################################################
##				Map Functions
#############################################################################################

plot_map = function(map){
  map = t(map)
  res = dim(map)[1:2] # get the resolution
  plot(1,1,xlim=c(1,res[2]),ylim=c(1,res[1]),asp=1,type='n',xaxs='i',yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(map,1,1,res[2],res[1])
}

create_map <- function(map_dimension_x, map_dimension_y, l0, lookup_resolution){
  map <- matrix(l0, map_dimension_x, map_dimension_y)
  scan_lookup <- matrix(0, map_dimension_x/lookup_resolution, map_dimension_y/lookup_resolution)
  initial_point <- c(0, 0)
  map_central_point <- ceiling(c(map_dimension_x/2,map_dimension_y/2))
  return (list(map, initial_point, map_central_point, scan_lookup, lookup_resolution,c()))
}

update_map_size <- function(map, increase_sizes, l0) {
  new_dimension <- dim(map[[1]]) + c(increase_sizes[2]+increase_sizes[4], increase_sizes[1]+increase_sizes[3])
  
  #update map
  new_map <- matrix(l0, new_dimension[1], new_dimension[2])
  new_map[(increase_sizes[2]+1):(increase_sizes[2]+dim(map[[1]])[1]), (increase_sizes[1]+1):(increase_sizes[1]+dim(map[[1]])[2])] = map[[1]]

  #Update lookup map
  lookup_resolution = map[[5]]
  new_scan_lookup <- matrix(0, new_dimension[1]/lookup_resolution, new_dimension[2]/lookup_resolution)  
  new_scan_lookup[((increase_sizes[2]/lookup_resolution)+1):((increase_sizes[2]/lookup_resolution)+dim(map[[4]])[1]), ((increase_sizes[1]/lookup_resolution)+1):((increase_sizes[1]/lookup_resolution)+dim(map[[4]])[2])] = map[[4]]

  #
  map_central_point <- ceiling(new_dimension/2)
  initial_point <- map[[2]] + round(c(increase_sizes[4]-increase_sizes[2], increase_sizes[3]-increase_sizes[1]))/2

  return (list(new_map, initial_point, map_central_point, new_scan_lookup, lookup_resolution,0))
}

transform_angles2 <- function(angle_vector){ 
  angle_vector[angle_vector<(-pi)] <- 2*pi + angle_vector[angle_vector<(-pi)]
  angle_vector[angle_vector> pi] <- angle_vector[angle_vector>2*pi]-2*pi
  return(angle_vector)
}

get_map_aux_data <- function(sensor_max_range, resolution){
  index_ranges <- 2*(sensor_max_range/resolution) + 1
  y_index = (t(matrix(1:index_ranges,index_ranges,index_ranges)) - (sensor_max_range/resolution)-1)*resolution
  x_index = ((matrix(1:index_ranges,index_ranges,index_ranges)) - (sensor_max_range/resolution)-1)*resolution
  angle_matrix = atan2(y_index, x_index)
  dist_matrix = sqrt(y_index**2 + x_index**2)
  matrix_index_to_consider = dist_matrix<(sensor_max_range-1.5)
 
  y_index_ = (t(matrix(1:index_ranges,index_ranges,index_ranges)) - (index_ranges/2)-1)
  x_index_ = ((matrix(1:index_ranges,index_ranges,index_ranges)) - (index_ranges/2)-1)
  xy_ = cbind(as.vector(y_index_[matrix_index_to_consider]), as.vector(x_index_[matrix_index_to_consider]))


  return(list(x_index, y_index, angle_matrix, dist_matrix, matrix_index_to_consider, xy_))
}

##
##  Create Local Map Update
##
create_local_map_update <- function(sensor_data, map_aux_data, l0, loc, lfree, alpha, beta){
  ranges = sensor_data$ranges[1:length(sensor_data$ranges)]
  ranges[1] = -1000 #discard the first measurement -- just for convenience
  

  sensor_measurement_matrix <- matrix(ranges[mat_index], dim(mat_index)[1], dim(mat_index)[2])

  #local_map <- matrix(l0, dim(mat_index)[1], dim(mat_index)[2])
  #
  sensor_measurement_matrix[!map_aux_data[[5]]] = -1
  sensor_measurement_matrix[!valid_angles] = -1
  #
  measurement_diff = sensor_measurement_matrix-map_aux_data[[4]]
  loc_position = (measurement_diff>0)*(measurement_diff<(alpha/2))
  #update the local map
  free <- which(map_aux_data[[4]] < sensor_measurement_matrix, arr.ind=T)
  occ <- which(loc_position==1, arr.ind=T)

  local_map_update = list(free, occ)

  return (local_map_update)
}


###
### Occupancy Grid update
###

update_map <- function(particle, particle_map, local_map_update, map_aux_data){
  
  ##
  ##  Increase map size
  ##
  x = particle_map[[3]][1] + (particle["x"] - particle_map[[2]][1])
  y = particle_map[[3]][2] + (-particle["y"] + particle_map[[2]][2])
 
  inc_map_size = dim(map_aux_data[[1]])/2
  increase_size = round(inc_map_size)
  inc_vec <- c(increase_size*((x-inc_map_size)<1), increase_size*((y-inc_map_size)<1), increase_size*((x+inc_map_size)>dim(particle_map[[1]])[1]), increase_size*((y+inc_map_size)>dim(particle_map[[1]])[1]))

  colnames(inc_vec) <- NULL
  if(sum(inc_vec)>0){
    particle_map <- update_map_size(particle_map, inc_vec, l0)
  }


  global_xy = c(particle_map[[3]][1] + (particle["x"] - particle_map[[2]][1]), particle_map[[3]][2] + (-particle["y"] + particle_map[[2]][2]))
  ##
  ##  ROTATE THE LOCAL MAP
  ##
  theta <- particle["theta"]
  rot <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)), 2, 2)
  ##matrix rotation, faster than using %*%
  xy_free <- local_map_update[[1]]-dim(map_aux_data[[1]])[1]/2
  xy_loc <- local_map_update[[2]]-dim(map_aux_data[[1]])[1]/2



  xy_free_ <- (cbind(c(global_xy[1]+xy_free[,1]*rot[1,1]+xy_free[,2]*rot[2,1]),c(global_xy[2]+xy_free[,1]*rot[1,2]+xy_free[,2]*rot[2,2])))
  xy_loc_ <-  (cbind(c(global_xy[1]+xy_loc[,1]*rot[1,1]+xy_loc[,2]*rot[2,1]),c(global_xy[2]+xy_loc[,1]*rot[1,2]+xy_loc[,2]*rot[2,2]))) 

  
  #UPDATE MAP
  particle_map_ = particle_map
  
  particle_map_[[1]][xy_free_] = particle_map[[1]][xy_free_] +(lfree-l0)
  particle_map_[[1]][xy_loc_] = particle_map[[1]][xy_loc_] +(loc-l0)
  
  
  lookup_resolution = particle_map[[5]]
  particle_map_[[4]][particle_map[[6]]] =  particle_map_[[4]][particle_map[[6]]]*0.99
  particle_map_[[4]][xy_free_/lookup_resolution] = 0
  particle_map_[[4]][xy_loc_/lookup_resolution] = 1
  particle_map_[[6]] = xy_loc_/lookup_resolution
  

  return (particle_map_)
}



##CREATE PARTICLES
sensor_data = database[[1]]$base_scan_message
particles <- matrix(c(rep(0, n_particles),rep(-500, n_particles),rep(0, n_particles)), ncol=3)
colnames(particles) <- c("x","y","theta")

particle_maps = lapply(1:n_particles, function(idx) create_map(1500, 2000, l0, lookup_resolution_param))
particles_weights = rep(1, n_particles)
map_aux_data = get_map_aux_data(sensor_data$range_max, resolution)


##
#  GLOBAL
#
##


angle_min <- sensor_data$angle_min
angle_max <- sensor_data$angle_max
angle_increment <- sensor_data$angle_increment
num_angles <- (angle_max-angle_min)/angle_increment

angle2index = (num_angles/(angle_max - angle_min))
  
mat_angles <- transform_angles2(map_aux_data[[3]])
valid_angles <- (mat_angles < angle_max)+(mat_angles < angle_min)
mat_index <- -((mat_angles-angle_max))*angle2index + 1
  #considering the beam thickness
valid_angles = valid_angles*(abs(mat_index - round(mat_index))<beta/2)
mat_index[mat_index > num_angles] = 1#discard the first measurement -- just for convenience
mat_index[mat_index < 1] = 1 #discard the first measurement -- just for convenience

sensor_angles = seq(angle_min, angle_max, angle_increment)
###
###
###
##



###
### FAST SLAM
###
start.time <- Sys.time()

for(i in seq(1+seq_step,length(database),seq_step)){
#sample new particles for motion model
  x_curr = c(database[[i]]$odom_message$x, database[[i]]$odom_message$y, database[[i]]$odom_message$z)
  x_prev = c(database[[i-seq_step]]$odom_message$x, database[[i-seq_step]]$odom_message$y, database[[i-seq_step]]$odom_message$z)
  base_scan_message =  database[[i]]$base_scan_message
  local_map_update <- create_local_map_update(base_scan_message, map_aux_data, l0, loc, lfree, 0.5, 1)

  ##Motion Model and Perception Model
  perception_particles = lapply(1:n_particles, function(idx){
     print(idx)
     particle = particles[idx, ]
     perception_particles <- matrix(c(rep(particle[1], n_perception_particles),rep(particle[2], n_perception_particles),rep(particle[3], n_perception_particles)), ncol=3)
     colnames(perception_particles) <- c("x","y","theta")

     perception_particles = sample_motion_model(perception_particles, x_curr, x_prev, alphas, database_info$map_metadata$resolution)
     ws=unlist(lapply(1:n_perception_particles, function(perc_idx) scan_match(perception_particles[perc_idx, ], particle_maps[[idx]], local_map_update, map_aux_data)))
   
     ws = (ws/(sum(ws)+0.01))
     pidx = which.max(ws)
     return(list(perception_particles[pidx, ],ws[pidx]))
  })

  particles = matrix(unlist(lapply(perception_particles, function(p) p[[1]])), ncol=3, byrow = TRUE)
  colnames(particles) <- c("x","y","theta")
  
  #MAP Update - update all
  particle_maps = lapply(1:n_particles, function(idx) update_map(particles[idx,], particle_maps[[idx]], local_map_update, map_aux_data))
  
  ws = unlist(lapply(perception_particles, function(p) p[[2]]))
  ##Update parcicles using the motion model
  if(sum(ws) > 0){
    new_particles_idx = sample(1:n_particles, n_particles, replace = TRUE, ws**map_particles_weights)
  
    particle_maps = particle_maps[new_particles_idx]
    particles <- particles[new_particles_idx,]
 
  ##
    particles_weights =  particles_weights + log(ws)
    particles_weights = particles_weights[new_particles_idx]
  
  }
  max_idx = which.max(particles_weights)

  ##Plot best particle
  m = particle_maps[[max_idx]][[1]]
  m[which(m>l0)] = 1
  m[which(m<l0)] = 0
  plot_map(m)
  global_xy = cbind(particle_maps[[max_idx]][[3]][1] + (particles[,"x"] - particle_maps[[max_idx]][[2]][1]), particle_maps[[max_idx]][[3]][2] + (particles[,"y"] + particle_maps[[max_idx]][[2]][2]))
   
  plot_map(m)
  points(global_xy[,1],global_xy[,2])
  
  
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("Tempo de Execucao:", time.taken))
#pdf("slam6_10_200.pdf")
#plot_map(m)
#dev.off()






 
