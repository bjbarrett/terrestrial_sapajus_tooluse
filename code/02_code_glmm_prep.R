library(lubridate)
library(rethinking)
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
data_list_daily <- list(
  nutcrackin = daily_crack$nutcrackin
)
##ZIP intercepts only of 1/0
m0ZIP <- ulam(
  alist(
    nutcrackin ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- al,
    ap ~ dnorm( 0 , 1 ),
    al ~ dnorm( 0 , 2 )
  ) , data=data_list_daily , chains=4)

precis( m0ZIP )
mean(daily_crack$nutcrackin)
exp(-0.24)

post <- extract.samples(m0ZIP)

## bernouli count of mean
# craete a list
data_list <- list(
  nutcrackin = focals_data_cropped$nutcrackin
)

rpois(n=1000 , lambda=2 )
simplehist(rpois(n=10000 , lambda=mean(exp(post$al)) ) )

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

set.seed(942)

mstones<- ulam(
  alist(
    rawstones_count_subplot ~ dpois( lambda ), #define poisson likelihood
    log(lambda) <-  a_bar + a_grid[grid_id],
    #priors
    a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
    a_grid[grid_id] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
    sigma_stone ~ dexp(1) #p
  ), data=raw_stones , chains=4, cores=4, control=list(adapt_delta=0.999) , iter=10000)

# mstones2<- ulam(
#   alist(
#     rawstones_count_subplot ~ dpois( lambda ), #define poisson likelihood
#     log(lambda) <-  a_bar + a_grid[grid_id],
#     #priors
#     a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
#     a_grid[grid_id] ~ dnorm( 0 , 1 ) #priovarying effects centered on mean
#     #sigma_stone ~ dexp(1) #p
#   ), data=raw_stones , chains=4, cores=4, control=list(adapt_delta=0.999) , iter=10000)

#look at prior predictive simulation of mean
prior_samps <- extract_prior_ulam(mstones)
dens(exp(prior_samps$a_bar))
dens(exp(prior_samps$a_bar) , xlim=c(0,10000))

#summary of ouputs
precis(mstones , depth=2 , pars="a_grid")
precis(mstones2 , depth=2 , pars="a_grid")

plot(precis(mstones , depth=2 , pars='a_grid'))
#plots preds per plot
post <- extract.samples(mstones)

str(post)
#real scale
for(i in 1:40){
  dens(exp(post$a_grid[,i]) , xlim=range(raw_stones$rawstones_count_subplot) , main=i)
  points( raw_stones$rawstones_count_subplot[raw_stones$grid_id==i] , rep(0,3) )
}

#link (log) scale
for(i in 1:40){
  dens(post$a_grid[,i] , main=i)
  points( log(raw_stones$rawstones_count_subplot[raw_stones$grid_id==i] ), rep(0,3) )
}

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

data_list_daily <- list(
  nutcrackin = joined_df$nutcrackin,
  id = joined_df$id,
  grid_id_follow = joined_df$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
  avg_stone = joined_df$average_stone, #average stones in big grid
  count_palm = joined_df$count_palm_trees, #palms in 22m grid
  avg_stone_std = standardize(joined_df$average_stone), #average stones in big grid
  log_avg_stone_std = standardize(joined_df$log_average_stone), #log average stones in big grid
  count_palm_std = standardize(joined_df$count_palm_trees), #palms in 22m grid
  # stone_post_mean = stone_precis$stone_post_mean,
  # stone_post_sd = stone_precis$stone_post_sd
  # stone_post_mean = joined_df2$stone_post_mean,
  # stone_post_sd = joined_df2$stone_post_sd,
  # N=nrow(joined_df2)
  stone_raw_grid = raw_stones$rawstones_count_subplot, # list of all stones, 3 observations per grid
  stone_grid_index = raw_stones$grid_id # list of grid ids associated with observations
)


dens(joined_df$count_palm_trees)
dens(joined_df$average_stone)

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

####causes of terrestriality

str(heights_data)
heights_data$height
heights_data$id <- as.integer(as.factor(heights_data$individual))
heights_data$follow_index <- as.integer(as.factor(heights_data$focal))
heights_data$day_index <- as.integer(as.factor(date(heights_data$timestamp)))

data_list_height <- list(
  height_dm=heights_data$height*10,
  height_m=heights_data$height,
  nutcrackin=heights_data$nutcracking,
  id = heights_data$id,
  follow_index = heights_data$follow_index,
  day_index = heights_data$day_index
)