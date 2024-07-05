#Load packages

library(sf)
library(terra)
library(mapview)
library(tmaptools)
#library(sp)
library(raster)
library(tidyverse)
library(leafem)
library(dplyr)
library(units)

#git config --global user.name "tati-valenca"


###Prepare the grids

#Obs: grid_id correspond to original grid numbers


#Read the shapefile
grid_110m <- st_read("data/grid_110m/grid_110m.shp")

#View the first few rows of the shapefile data
head(grid_110m)

#Erase unnecessary variables (keep just the polygon)
grid_110m_subset <- grid_110m[, !(names(grid_110m) == "top")]
grid_110m_subset <- grid_110m_subset[, !(names(grid_110m_subset) == "left")]
grid_110m_subset <- grid_110m_subset[, !(names(grid_110m_subset) == "right")]
grid_110m_subset <- grid_110m_subset[, !(names(grid_110m_subset) == "bottom")]
grid_110m_subset <- grid_110m_subset[, !(names(grid_110m_subset) == "id")]

head(grid_110m_subset)
plot(grid_110m_subset)


#Adjusting a rectangle of 22m cells to the grid 110m,
rectangle_22m <- st_make_grid(grid_110m_subset, square = T, cellsize = c(22,22), n=c(60,25) )%>%
  st_sf()


#Crop the 22m rectangle to the size of the grid 110m
grid_22m_size_110 <- crop_shape(rectangle_22m, grid_110m_subset, polygon = TRUE)


#Filter cells without area
grid_22m_size_110$area <- st_area(grid_22m_size_110)
minimum_area <- set_units(0.001, "m^2")
grid_22m_filtered <- grid_22m_size_110[grid_22m_size_110$area > minimum_area, ]

mapview(grid_22m_filtered)


#Erase an unnecessary variable (area)
grid_22m <- grid_22m_filtered[, !(names(grid_22m_filtered) == "area")]

# give cells unique identifier
grid_22m$cell_id <- c(1:nrow(grid_22m))

mapview(grid_22m)



#Extracting grid_id from grid_110m to grid_22m

  # Step 1: give id columns in both grids
  grid_110m$id <- grid_110m$grid_id
  grid_22m$id <- grid_22m$cell_id

  # Step 2: Compute centroids of the smaller cells
  small_centroids <- st_centroid(grid_22m)

  # Step 3: Match centroids to larger grid cells, using st_within to find centroids within larger cells
  matches <- st_within(small_centroids, grid_110m, sparse = FALSE)

  # Step 4: Extract the big grid ID for each small cell based on centroid location
  big_cell_ids <- vector("list", length = nrow(small_centroids))
  for (i in 1:nrow(small_centroids)) {
    if (any(matches[i,])) {
      big_cell_ids[[i]] <- grid_110m$id[which(matches[i, ])]
    } else {
     big_cell_ids[[i]] <- NA  # For centroids that do not match, though this should not happen
    }
  }
  
  
  # Step 5: Flatten the list to a vector (since each centroid matches exactly one big cell)
  grid_id <- unlist(big_cell_ids)

  # Step 6: Add the big grid IDs to the small grid data
  grid_22m$grid_id <- grid_id

    # Step 7: Erase unnecessary variables (id)
  grid_22m <- grid_22m[, !(names(grid_22m) == "id")]
  
  mapview(grid_22m)

  # Checking whether the previous steps worked correctly
    # count number of small cells in each big cell, should be 25
    
    temp_data <- st_drop_geometry(grid_22m)
    grouped_data <- group_by(temp_data, grid_id)
    summary_data <- summarise(grouped_data, n = n())
    print(summary_data, n = 40)

    
### Extract number of palm trees for each cell

  # Set the working directory
  # setwd("C:/Users/tativ/OneDrive/Doutorado_Spatial data/Terrestriality and tool use/palmeiras_UB/")

  # give ids unique identifier
  palm_trees_all <- st_read("data/palmeiras_UB_spatial/Palm_nuts_UB_2023_reprojected.shp")
  plot(palm_trees_all)
  mapviewOptions(fgb = FALSE)
  mapview(palm_trees_all)
  
  palm_trees <- crop_shape(palm_trees_all, grid_22m, polygon = TRUE)
  palm_trees$id <- c(1:nrow(palm_trees))
  mapview(palm_trees)


  #Checking whether projections match
  crs(grid_22m)
  crs(palm_trees)

  #Extract number of points for each cell_id

  intersection <- st_intersection(x = grid_22m, y = palm_trees)
  plot(intersection)
  table(intersection$cell_id) 
  
  int_result <- intersection %>% 
   group_by(cell_id) %>% 
    count()

  as.data.frame(int_result)[,-3] # table with the values for each cell_id (remove the third column - the geometry)
  
  #Put palm tree data into grid cell table
  
  all_cells <- data.frame(cell_id = 1:1000) #create a complete dataframe of all cell ids
  int_result <- merge(all_cells, int_result, by = "cell_id", all.x = TRUE) #fill the missing cell ids
  int_result$n[is.na(int_result$n)] <- 0 #fill all missing values with zero
  int_result <- int_result[order(int_result$cell_id),] #order by cell id results
  grid_22m$count_palm_trees <- int_result$n #both have the same cell id order
  
  mapview(grid_22m) + mapview(palm_trees)
  

### Reading average availability of stones for each plot
  
  # Set the working directory
  # setwd("C:/Users/tativ/OneDrive/Doutorado_Spatial data/Terrestriality and tool use/data/")
  # Read data
  average_stones_plot <- read.csv("data/raw_stones_average_plots.csv")
  
  all_cell_ids <- st_drop_geometry(grid_22m)
  average_stones_plot <- merge(all_cell_ids, average_stones_plot, by = "grid_id", all.x = TRUE) #fill the missing cell ids
  average_stones_plot <- average_stones_plot[order(average_stones_plot$cell_id),] #order by cell id results
  grid_22m$average_stone <- average_stones_plot$rawstones_average #both have the same cell id order
  

### Extract information from grids to the data
    
  # Set the working directory
  # setwd("C:/Users/tativ/OneDrive/Doutorado_Spatial data/Terrestriality and tool use/data/")
    
  # Read data
    heights_data <- read.csv("data/sertao_heights_2022-2023.csv")
    focals_data <- st_read("data/sertao_focals_spatial/sertao_focals_2022-2023_shape.shp")
    nutcracking_data <- st_read("data/sertao_nutcracking_spatial/sertao_nutcracking_2022-2023_partial_shape.shp")

    mapview(focals_data)
    mapview(nutcracking_data)
    mapview(grid_22m) + mapview(focals_data) +  mapview(nutcracking_data, col.regions ="red")
    
     
    # Preparing focals_data
    
      # Step 1: Crop the focals data to the size of the grid 22m
      focals_data_cropped <- crop_shape(focals_data, grid_22m, polygon = TRUE)
      mapview(focals_data_cropped) + mapview(grid_22m)
    
      # Step 2: Match focal points to grid cells
      matches <- st_within(focals_data_cropped, grid_22m, sparse = FALSE)
    
      # Step 3: Extract the grid_cell for each focal point
      cell <- vector("list", length = nrow(focals_data_cropped))
      for (i in 1:nrow(focals_data_cropped)) {
        if (any(matches[i,])) {
          cell[[i]] <- grid_22m$cell_id[which(matches[i, ])]
        } else {
          grid_22m[[i]] <- NA
        }
      }
    
      # Step 4: Flatten the list to a vector
      cell <- unlist(cell)
      corresponding_cell_id <- data.frame(Value= cell)
      corresponding_cell_id$New_Column <- c()
    
      # Step 5: Add the grid_id to the focal points
      focals_data_cropped$corresponding_cell_id <- cell

      #Step 6: Check whether the focal points are showing the right corresponding grid cell
      mapview(grid_22m) + mapview(focals_data_cropped)
    
    
     # Preparing nutcracking_data
      
      # Step 1: Crop the nutcracking data to the size of the grid 22m
      nutcracking_data_cropped <- crop_shape(nutcracking_data, grid_22m, polygon = TRUE)
      mapview(nutcracking_data_cropped) + mapview(grid_22m)
      
      # Step 2: Match nutcracking points to grid cells
      matches <- st_within(nutcracking_data_cropped, grid_22m, sparse = FALSE)
      
      # Step 3: Extract the grid_cell for each nutcracking point
      cell_2 <- vector("list", length = nrow(nutcracking_data_cropped))
      for (i in 1:nrow(nutcracking_data_cropped)) {
        if (any(matches[i,])) {
          cell_2[[i]] <- grid_22m$cell_id[which(matches[i, ])]
        } else {
          grid_22m[[i]] <- NA
        }
      }
      
      # Step 4: Flatten the list to a vector
      cell_2 <- unlist(cell_2)
      corresponding_cell_id_2 <- data.frame(Value= cell_2)
      corresponding_cell_id_2$New_Column <- c()
      
      # Step 5: Add the grid_id to the nutcracking points
      nutcracking_data_cropped$corresponding_cell_id_2 <- cell_2
      
      #Step 6: Check whether the nutcracking points are showing the right corresponding grid cell
      mapview(grid_22m) + mapview(nutcracking_data_cropped)
      
      # Step 7: add mean tree and stone counts to 110m grid and 22m grid
      grid_110m <- merge(grid_110m,average_stones_plot, duplicateGeoms = T)
      grid_22m <- merge(grid_22m,average_stones_plot, duplicateGeoms = T)
      
      
      
mapview(grid_22m ) + mapview(focals_data_cropped , cex=0.5) +  mapview(nutcracking_data_cropped, col.regions ="red", cex=0.5) + mapview(palm_trees, col.regions = "green", cex=0.5)


mapview(grid_110m ) +  mapview(nutcracking_data_cropped, col.regions ="red", cex=0.5)
mapview(grid_110m ) +  mapview(focals_data_cropped , cex=0.5)
