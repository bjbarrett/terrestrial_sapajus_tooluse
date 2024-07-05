library(lubridate)
library(rethinking)
library(units)
library(dplyr)

# collate dataframe
# number of tool use events in a cell (22m or 110m) per individual per day
# hours spent following an indivdual ona a given day
# denisty of stones, density of palms, 
# currently with aborted follows, the last GPS point uses the first sampling point 
# to do heights integration
mapview(grid_22m ) + mapview(focals_data_cropped , cex=0.5) +  mapview(nutcracking_data_cropped, col.regions ="red", cex=0.5) + mapview(palm_trees, col.regions = "green", cex=0.5)
mapview(grid_110m )
str(grid_22m)
str(focals_data_cropped)
focals_data_cropped$date <- date(focals_data_cropped$timestamp) #unique date

# number of nutcracks per individual per day across all of data
daily_crack <- aggregate(nutcrackin ~ individual + date + corresponding_cell_id , focals_data_cropped, sum) 
hist(focals_data_cropped$nutcrackin)
hist(daily_crack$nutcrackin)
str(daily_crack)
unique(daily_crack$individual)
unique(focals_data_cropped$individual)
str(focals_data_cropped)
str(heights_data) ## this is where i will get sampling time by summing up minutes per day

#need unique number of minutes per follow 
heights_data$date <- date(heights_data$timestamp) #create dates
temp <- aggregate(minute ~ individual + date + focal, heights_data, length) # number of nutcracks per individual across all of data
temp$minutes_followed <- temp$minute - 1 # ges=ts rid of zero issues
hist(temp$minutes_followed) #there is a duplicated follow, check up on it
### STOPPING HERE need to get dataframes to talk and run first model

# temp[which(temp$minutes_followed > 5),]


#WERE MONKEYS IN A CELL, AND IF THEY WERE DOES WAS NUT CRACKING OBSERVED THERE, IRRESPECTIVE OF TIME
# HOW MUCH TIME TO MONKEYS SPEND IN CELLS
# GIVEN A MONKEY WAS IN A CELL FOR X TIME, WHAT WAS THEIR DAILY RATE OF TOOL USE THEIR

# intercepts only poisson
# data_list_daily <- list(
#   nutcrackin = daily_crack$nutcrackin
# )
# ##ZIP intercepts only of 1/0
# m0ZIP <- ulam(
#   alist(
#     nutcrackin ~ dzipois( p , lambda ),
#     logit(p) <- ap,
#     log(lambda) <- al,
#     ap ~ dnorm( 0 , 1 ),
#     al ~ dnorm( 0 , 2 )
#   ) , data=data_list_daily , chains=4)
# 
# precis( m0ZIP )
# mean(daily_crack$nutcrackin)
# exp(-0.24)
# 
# post <- extract.samples(m0ZIP)

## bernouli count of mean
# craete a list
data_list <- list(
  nutcrackin = focals_data_cropped$nutcrackin
)

rpois(n=1000 , lambda=2 )
#simplehist(rpois(n=10000 , lambda=mean(exp(post$al)) ) )

###what cuses terrestriality
str(grid_22m)


str(focals_data_cropped)
focals_data_cropped$corresponding_cell_id
grid_22m$cell_id

str(focals_data_cropped)
str(grid_22m)

# joined_df <- merge(focals_data_cropped,grid_22m, by.x = "corresponding_cell_id", 
#                   by.y = "cell_id", all.x = TRUE, all.y = FALSE)
# ab<-union(focals_data_cropped, grid_22m)
g22df <- as.data.frame(grid_22m)
fdcdf <- as.data.frame(focals_data_cropped)
joined_df <- merge(fdcdf,g22df, by.x = "corresponding_cell_id", 
                   by.y = "cell_id", all.x = TRUE, all.y = FALSE)
str(joined_df)
joined_df$id <- as.integer(as.factor(joined_df$individual))
joined_df$log_average_stone <- log(joined_df$average_stone)

sort(unique(joined_df$timestamp)) # this gives timestamps which we need to link to the stones
##### add in raw stones

raw_stones <- read.csv("data/raw_stones_subplots.csv")
raw_stones$rawstones_count_subplot
str(raw_stones)
dens(raw_stones$rawstones_count_subplot)
raw_stones$rawstones_count_subplot
mean(raw_stones$rawstones_count_subplot)
for(i in 1:max(raw_stones$grid_id)){
  dens(raw_stones$rawstones_count_subplot[raw_stones$grid_id==i])
  #points(rep(0,3) , raw_stones$rawstones_count_subplot[raw_stones$grid_id==i])
  points( raw_stones$rawstones_count_subplot[raw_stones$grid_id==i] , rep(0,3) )
}

# set.seed(942)
# 
# mstones<- ulam(
#   alist(
#     rawstones_count_subplot ~ dpois( lambda ), #define poisson likelihood
#     log(lambda) <-  a_bar + a_grid[grid_id],
#     #priors
#     a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
#     a_grid[grid_id] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
#     sigma_stone ~ dexp(1) #p
#   ), data=raw_stones , chains=4, cores=4, control=list(adapt_delta=0.999) , iter=10000)
# 
# # mstones2<- ulam(
# #   alist(
# #     rawstones_count_subplot ~ dpois( lambda ), #define poisson likelihood
# #     log(lambda) <-  a_bar + a_grid[grid_id],
# #     #priors
# #     a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
# #     a_grid[grid_id] ~ dnorm( 0 , 1 ) #priovarying effects centered on mean
# #     #sigma_stone ~ dexp(1) #p
# #   ), data=raw_stones , chains=4, cores=4, control=list(adapt_delta=0.999) , iter=10000)
# 
# #look at prior predictive simulation of mean
# prior_samps <- extract_prior_ulam(mstones)
# dens(exp(prior_samps$a_bar))
# dens(exp(prior_samps$a_bar) , xlim=c(0,10000))
# 
# #summary of ouputs
# precis(mstones , depth=2 , pars="a_grid")
# precis(mstones2 , depth=2 , pars="a_grid")
# 
# plot(precis(mstones , depth=2 , pars='a_grid'))
# #plots preds per plot
# post <- extract.samples(mstones)
# 
# str(post)
# #real scale
# for(i in 1:40){
#   dens(exp(post$a_grid[,i]) , xlim=range(raw_stones$rawstones_count_subplot) , main=i)
#   points( raw_stones$rawstones_count_subplot[raw_stones$grid_id==i] , rep(0,3) )
# }
# 
# #link (log) scale
# for(i in 1:40){
#   dens(post$a_grid[,i] , main=i)
#   points( log(raw_stones$rawstones_count_subplot[raw_stones$grid_id==i] ), rep(0,3) )
# }

#in stan i would model this jointly, but i will extract posterior mean and SD on log scale instead and feed that into the measurement error model
# 
# stone_precis <- as.data.frame(precis(mstones2 , depth=2 , pars="a_grid"))
# stone_precis <- stone_precis[,1:2]
# stone_precis$stone_post_index <- 1:40
# names(stone_precis)[1:2] <- c("stone_post_mean" , "stone_post_sd")
# merge
# 
# joined_df2 <- merge(joined_df,stone_precis, by.x = "grid_id", 
#                    by.y = "stone_post_index", all.x = TRUE, all.y = FALSE)
# plot(joined_df2$stone_post_mean ~ joined_df2$log_average_stone )


####causes of terrestriality

str(heights_data)
heights_data$height
heights_data$id <- as.integer(as.factor(heights_data$individual))
heights_data$follow_index <- as.integer(as.factor(heights_data$focal))
heights_data$day_index <- as.integer(as.factor(date(heights_data$timestamp)))
heights_data$follow_length<-NA
for(i in 1:max(heights_data$follow_index)){
  heights_data$follow_length[i] <- length(which(heights_data$follow_index==heights_data$follow_index[i]))
}
nrow(heights_data)
# heights_data <- heights_data[heights_data$follow_length>1,]
heights_data$follow_length <- heights_data$follow_length - 1 #change to minutes

which(heights_data$follow_length==0)
which(heights_data$follow_length=='NA')

nrow(heights_data)

temp2 <- heights_data[,c("timestamp" ,"follow_length")]

joined_df$timestamp <- substr(joined_df$timestamp,1,nchar(joined_df$timestamp)-4)#remove weird trailing zeros
#make dates compatible  
temp2$timestamp <- ymd_hms(temp2$timestamp)
joined_df$timestamp <- ymd_hms(joined_df$timestamp)
temp2 <- temp2[!duplicated(temp2), ]

joined_df2<- merge(temp2,joined_df, by="timestamp")
str(joined_df2)
#lets drop the rows where follows end
joined_df3 <- joined_df2[joined_df2$focal_begi==0,]
joined_df4 <- joined_df3[which(is.na(joined_df3$follow_length)==FALSE),]
joined_df4 <- joined_df4[joined_df4 $follow_length>0,]
str(joined_df4)

data_list_daily <- list(
  nutcrackin = joined_df4$nutcrackin,
  id = joined_df4$id,
  grid_id_follow = joined_df4$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
  grid_id_follow_dummy = as.integer(as.factor(joined_df4$grid_id)),# had to reassign ids so indexing works finwe with ulam
  avg_stone = joined_df4$average_stone, #average stones in big grid
  count_palm = joined_df4$count_palm_trees, #palms in 22m grid
  avg_stone_std = standardize(joined_df4$average_stone), #average stones in big grid
  log_avg_stone_std = standardize(joined_df4$log_average_stone), #log average stones in big grid
  count_palm_std = standardize(joined_df4$count_palm_trees), #palms in 22m grid
  follow_length=joined_df4$follow_length,
  log_follow_length=log(joined_df4$follow_length),
  sex_index = as.integer(as.factor(joined_df4$sex)),
  # stone_post_mean = stone_precis$stone_post_mean,
  # stone_post_sd = stone_precis$stone_post_sd
  # stone_post_mean = joined_df2$stone_post_mean,
  # stone_post_sd = joined_df2$stone_post_sd,
  # N=nrow(joined_df2)
  stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
  stone_grid_index = raw_stones$grid_id # list of grid ids associated with observations
 
)


dens(joined_df4$count_palm_trees)
dens(joined_df4$average_stone)

m1ZIP <- ulam(
  alist(
    nutcrackin ~ dzipois( p , lambda ),
    logit(p) <- ap + a_id[id,1],
    log(lambda) <- al + a_id[id,2],
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    al ~ dnorm( 0 , 2 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)

plot(precis( m1ZIP , depth=2 ))
plot(precis( m1ZIP , depth=3 , pars='a_id' ))
plot(precis( m1ZIP , depth=3 , pars='L_Rho_id' ))

str(heights_data)
str(joined_df)

joined_df4$ts2  <- ymd_hms(joined_df4$timestamp)
heights_data$ts2  <- ymd_hms(heights_data$timestamp)
colnames(joined_df)
heights_data <- heights_data
joined_df_blah <- joined_df4[,c("corresponding_cell_id","location.l","location_1"           
                               ,"focal_begi","geometry.x","grid_id","count_palm_trees","average_stone","rawstones_average",
                            "geometry.y","log_average_stone","ts2")]
blah <- merge(select(heights_data, -c(follow_length)),joined_df_blah, by="ts2")
str(blah)
sort(unique(blah$ts2))
ddd <- sort(unique(blah$ts2))
as.POSIXct(heights_data$timestamp)
##terrest data
data_list_height <- list(
  
    height_dm=as.integer(heights_data$height*10),
  height_m=heights_data$height,
  nutcrackin=heights_data$nutcracking,
  id = heights_data$id,
  follow_index = heights_data$follow_index,
  day_index = heights_data$day_index
  #follow_length = heights_data$follow_length,
  #log_follow_length = log(heights_data$follow_length)
)

blah_tu <- blah[blah$nutcracking==1,]
blah_ntu <- blah[blah$nutcracking==0,]

###. The euclidian distance matrix:

#reorder grid so cell 1 is first
grid_110m_subset2 <- grid_110m_subset %>% 
    arrange(grid_id)
#drop the cell no monkeys visit (15) and calculate grid centroids
mtx_distance_110 <- st_centroid(grid_110m_subset2[grid_110m_subset2$grid_id!=15,]) 
mtx_distance_110_full <- st_centroid(grid_110m_subset2) 

#calculate matrix distance
Dmat_110_dm <- drop_units(st_distance(mtx_distance_110,mtx_distance_110))/100
Dmat_110_m <- drop_units(st_distance(mtx_distance_110,mtx_distance_110))
Dmat_110_dm_full <- drop_units(st_distance(mtx_distance_110_full,mtx_distance_110_full))/100

grid_22m <- grid_22m %>% 
    arrange(cell_id)
obs_22_monkey <- sort(unique(blah$corresponding_cell_id))
grid_22m_subset <- grid_22m[which(grid_22m$cell_id %in% obs_22_monkey),]
mtx_distance_22_full <- st_centroid(grid_22m$geometry) 
Dmat_22_dm_full <- drop_units(st_distance(mtx_distance_22_full,mtx_distance_22_full))/100
Dmat_22_m_full <- drop_units(st_distance(mtx_distance_22_full,mtx_distance_22_full))

mtx_distance_22 <- st_centroid(grid_22m$geometry) 
data_list_height2 <- list(
    height_dm=as.integer(blah$height*10),
    height_m=blah$height,
    nutcrackin=blah$nutcracking,
    id = blah$id,
    follow_index = blah$follow_index,
    day_index = blah$day_index,
    grid_id_follow = blah$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    grid_id_follow_dummy = as.integer(as.factor(blah$grid_id)),# had to reassign ids so indexing works finwe with ulam
    avg_stone = blah$average_stone, #average stones in big grid
    count_palm = blah$count_palm_trees, #palms in 22m grid
    avg_stone_std = standardize(blah$average_stone), #average stones in big grid
    log_avg_stone_std = standardize(blah$log_average_stone), #log average stones in big grid
    count_palm_std = standardize(blah$count_palm_trees), #palms in 22m grid
    sex_index = as.integer(as.factor(blah$sex)),
    stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
    stone_grid_index = raw_stones$grid_id, # list of grid ids associated with observations
    #follow_length = heights_data$follow_length,
    #log_follow_length = log(heights_data$follow_length)
    cell_id=blah$corresponding_cell_id,
    cell_id_index=as.integer(as.factor(blah$corresponding_cell_id)),
    nc=blah$nutcracking + 1,
    Dmat110=Dmat_110_dm,
    Dmat110_full=Dmat_110_dm_full
)

data_list_height3 <- list(
    height_dm=as.integer(blah$height*10),
    height_m=blah$height,
    nutcrackin=blah$nutcracking,
    id = blah$id,
    follow_index = blah$follow_index,
    day_index = blah$day_index,
    grid_id_follow = blah$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    grid_id_follow_dummy = as.integer(as.factor(blah$grid_id)),# had to reassign ids so indexing works finwe with ulam
    avg_stone = blah$average_stone, #average stones in big grid
    count_palm = blah$count_palm_trees, #palms in 22m grid
    avg_stone_std = standardize(blah$average_stone), #average stones in big grid
    log_avg_stone_std = standardize(blah$log_average_stone), #log average stones in big grid
    count_palm_std = standardize(blah$count_palm_trees), #palms in 22m grid
    sex_index = as.integer(as.factor(blah$sex)),
    stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
    stone_grid_index = raw_stones$grid_id, # list of grid ids associated with observations
    #follow_length = heights_data$follow_length,
    #log_follow_length = log(heights_data$follow_length)
    cell_id=blah$corresponding_cell_id,
    cell_id_index=as.integer(as.factor(blah$corresponding_cell_id)),
    nc=blah$nutcracking + 1,
    Dmat110=Dmat_22_dm
    
)
data_list_height_tu <- list(
    height_dm=as.integer(blah_tu$height*10),
    height_m=blah_tu$height,
    nutcrackin=blah_tu$nutcracking,
    id = as.integer(as.factor(blah_tu$id)),
    follow_index = blah_tu$follow_index,
    day_index = blah_tu$day_index,
    grid_id_follow = blah_tu$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    grid_id_follow_dummy = as.integer(as.factor(blah_tu$grid_id)),# had to reassign ids so indexing works finwe with ulam
    avg_stone = blah_tu$average_stone, #average stones in big grid
    count_palm = blah_tu$count_palm_trees, #palms in 22m grid
    avg_stone_std = standardize(blah_tu$average_stone), #average stones in big grid
    log_avg_stone_std = standardize(blah_tu$log_average_stone), #log average stones in big grid
    count_palm_std = standardize(blah_tu$count_palm_trees), #palms in 22m grid
    sex_index = as.integer(as.factor(blah_tu$sex)),
    stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
    stone_grid_index = raw_stones$grid_id, # list of grid ids associated with observations
    #follow_length = heights_data$follow_length,
    #log_follow_length = log(heights_data$follow_length)
    cell_id=blah_tu$corresponding_cell_id,
    cell_id_index=as.integer(as.factor(blah_tu$corresponding_cell_id)),
    nc=blah_tu$nutcracking + 1,
    Dmat110=mtx_distance_110
)

data_list_height_ntu <- list(
    height_dm=as.integer(blah_ntu$height*10),
    height_m=blah_ntu$height,
    nutcrackin=blah_ntu$nutcracking,
    id = blah_ntu$id,
    follow_index = blah_ntu$follow_index,
    day_index = blah_ntu$day_index,
    grid_id_follow = blah_ntu$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    grid_id_follow_dummy = as.integer(as.factor(blah_ntu$grid_id)),# had to reassign ids so indexing works finwe with ulam
    avg_stone = blah_ntu$average_stone, #average stones in big grid
    count_palm = blah_ntu$count_palm_trees, #palms in 22m grid
    avg_stone_std = standardize(blah_ntu$average_stone), #average stones in big grid
    log_avg_stone_std = standardize(blah_ntu$log_average_stone), #log average stones in big grid
    count_palm_std = standardize(blah_ntu$count_palm_trees), #palms in 22m grid
    sex_index = as.integer(as.factor(blah_ntu$sex)),
    stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
    stone_grid_index = raw_stones$grid_id, # list of grid ids associated with observations
    #follow_length = heights_data$follow_length,
    #log_follow_length = log(heights_data$follow_length)
    cell_id=blah_ntu$corresponding_cell_id,
    cell_id_index=as.integer(as.factor(blah_ntu$corresponding_cell_id)),
    nc=blah_ntu$nutcracking + 1,
    Dmat110=mtx_distance_110
    
)

data_list_daily_2 <- list(
    nutcrackin = joined_df4$nutcrackin,
    id = joined_df4$id,
    grid_id_follow = joined_df4$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    grid_id_follow_dummy = as.integer(as.factor(joined_df4$grid_id)),# had to reassign ids so indexing works finwe with ulam
    avg_stone = joined_df4$average_stone, #average stones in big grid
    count_palm = joined_df4$count_palm_trees, #palms in 22m grid
    avg_stone_std = standardize(joined_df4$average_stone), #average stones in big grid
    log_avg_stone_std = standardize(joined_df4$log_average_stone), #log average stones in big grid
    count_palm_std = standardize(joined_df4$count_palm_trees), #palms in 22m grid
    follow_length=joined_df4$follow_length,
    log_follow_length=log(joined_df4$follow_length),
    sex_index = as.integer(as.factor(joined_df4$sex)),
    # stone_post_mean = stone_precis$stone_post_mean,
    # stone_post_sd = stone_precis$stone_post_sd
    # stone_post_mean = joined_df2$stone_post_mean,
    # stone_post_sd = joined_df2$stone_post_sd,
    # N=nrow(joined_df2)
    stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
    stone_grid_index = raw_stones$grid_id, # list of grid ids associated with observations
    cell_id=joined_df4$corresponding_cell_id,
    Dmat110=Dmat_110_dm_full,
    Dmat22=Dmat_22_dm_full,
    cell_id_all=grid_22m$cell_id,
    count_palm_all=grid_22m$count_palm_trees
)


###arrive and leave
str(blah)
blah_zero <- blah[blah$minute==0,]

data_list_height_zero <- list(
    height_dm=as.integer(blah_zero$height*10),
    height_m=blah_zero$height,
    nutcrackin=blah_zero$nutcracking,
    id = blah_zero$id,
    follow_index = blah_zero$follow_index,
    day_index = blah_zero$day_index,
    grid_id_follow = blah_zero$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    grid_id_follow_dummy = as.integer(as.factor(blah_zero$grid_id)),# had to reassign ids so indexing works finwe with ulam
    avg_stone = blah_zero$average_stone, #average stones in big grid
    count_palm = blah_zero$count_palm_trees, #palms in 22m grid
    avg_stone_std = standardize(blah_zero$average_stone), #average stones in big grid
    log_avg_stone_std = standardize(blah_zero$log_average_stone), #log average stones in big grid
    count_palm_std = standardize(blah_zero$count_palm_trees), #palms in 22m grid
    sex_index = as.integer(as.factor(blah_zero$sex)),
    stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
    stone_grid_index = raw_stones$grid_id, # list of grid ids associated with observations
    #follow_length = heights_data$follow_length,
    #log_follow_length = log(heights_data$follow_length)
    cell_id=blah_zero$corresponding_cell_id,
    cell_id_index=as.integer(as.factor(blah_zero$corresponding_cell_id)),
    nc=blah_zero$nutcracking + 1,
    Dmat110=Dmat_110_dm,
    Dmat110_full=Dmat_110_dm_full
)


str(blah_zero)
