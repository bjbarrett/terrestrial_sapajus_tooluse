library(lubridate)
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

temp <- aggregate(nutcrackin ~ individual + date + corresponding_cell_id , focals_data_cropped, sum) # number of nutcracks per individual per day across all of data
hist(focals_data_cropped$nutcrackin)
hist(temp$nutcrackin)
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
m0 <- ulam(
  alist(
    y ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- al,
    ap ~ dnorm( -1.5 , 1 ),
    al ~ dnorm( 1 , 0.5 )
  ) , data=list(y=y) , chains=4 )
precis( m12.3 )