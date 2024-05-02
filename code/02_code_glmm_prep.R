library(lubridate)
library(rethinking)
set_ulam_cmdstan(TRUE) #set cmdstan to true
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
m0 <- ulam(
  alist(
    nutcrackin ~ binomial(1,p),
    logit(p) <- ap ,
    ap ~ dnorm( 0 , 1 )
    ) , data=data_list , chains=4)

# to go from probability to log-odds use logit link 
logit(.99)
# to go from log-odds to proability/real scale use
logistic(.99)
inv_logit(.99)

post <- extract.samples(m0)
str(post)
dens(post$ap) #log odds scale
dens(logistic(post$ap)) #log odds scale
abline(v=mean(logistic(post$ap)))
precis( m0 ) # real probability scale
inv_logit(-1.27)
mean(focals_data_cropped$nutcrackin)


###what cuses terrestriality
str(grid_22m)
str(focals_data_cropped)
focals_data_cropped$corresponding_cell_id
grid_22m$cell_id

joined_df <- merge(focals_data_cropped,grid_22m, by.x = "corresponding_cell_id", 
                  by.y = "cell_id", all.x = TRUE, all.y = FALSE)

