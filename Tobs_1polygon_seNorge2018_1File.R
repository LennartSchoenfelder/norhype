library(dplyr)
library(ncdf4)
library(tictoc)
library(sf)
library(xts)
library(stringr)

# #this script clips gridded Precipitation from netCDF files to subbasins in HYPE
# NC files can be downloaded from kss.no
# output is a "Tobs.txt" for direct use in a HYPE model setup
##  hist_CNRM_CCLM_RR_daily_mm_1990.nc is an example filename
dir.gis = "C:/Users/lennarts/OneDrive - SINTEF/SINTEF/PHD/HYPE_T/GIS/"
dir.nc = "E:/ClimateData/seNorge/"


#sets work directory to path where script is located
setwd(dir.nc)
shell.exec(getwd())



# import shapefile (polygons) that delineate the subcatchments (exportable from WHIST)
#see WHIST from SMHI
#Polygons can be in any projection, projection will be changed to  (according to netcdf file)
#shapefile
#subbasins <- readOGR(dsn = "C:/Users/lennarts/OneDrive - SINTEF/SINTEF/PHD/HYPE_T/GIS", layer= "Gaula_delineation.shp", encoding = "UTF8")
subbasins = st_read(dsn = paste0(dir.gis,"Gaula_catchment.gpkg"), layer = "Gaula_catchment")



st_crs(subbasins)
subbasins <- st_transform(subbasins, crs = 32633)

subbasins$ROWNR = seq(1:nrow(subbasins))





# # #load first NC PREC data to create spatial structure of dataset
setwd(dir.nc)

#filelist <- list.files(path = dir.nc, pattern = "*.nc")
nc <- nc_open(filename = "./seNorge2018_TG_1957_2017.nc", verbose = F, write = FALSE)
nc$var$tg
#check if projection is correct and find variable name
show(nc)
nc$dim$time
#select variable
varname = "tg"

#lonlat by values for wgs84 utm 33 projection from VALS
lat <- ncvar_get(nc,"Y")
lon <- ncvar_get(nc,"X")
lonlat <- as.data.frame(merge(lon,lat)) 
rm(lat, lon)

# #check the cornerpoints of your domain (optional)
# #left top
# lonlat[1,]
# #right top
# lonlat[nc$dim$X$len,]
# #left bottom
# lonlat[(nc$dim$Y$len-1)*nc$dim$X$len+1,]
# #right bottom
# lonlat[nc$dim$X$len*nc$dim$Y$len,]

#create matrix of indices, this will be needed to access the NC file rapidly
lonind <- as.double(rep(1:nc$dim$X$len, times = nc$dim$Y$len))
latind <- as.double(rep(1:nc$dim$Y$len, each = nc$dim$X$len))
lonlatind <- cbind (lonind, latind)
lonlatind <- as.data.frame(lonlatind)
coords_ind_mat <- cbind(lonlat,lonlatind)
rm(lonind, latind)

# ## check coordinates
# # #left top
# coords_ind_mat[1,]
# #right top
# coords_ind_mat[nc$dim$X$len,]
# #left bottom
# coords_ind_mat[(nc$dim$Y$len-1)*nc$dim$X$len+1,]
# #right bottom
# coords_ind_mat[nc$dim$X$len*nc$dim$Y$len,]

lonlat_mat = as.matrix(lonlat)


# # # create & write dataframe
sf_Prec_points = st_as_sf(coords_ind_mat, coords = c("x","y"), crs = 32633)
st_crs(sf_Prec_points)
plot(head(sf_Prec_points))
class(sf_Prec_points)



# #DEPRECATED
# Prec_points_old <- SpatialPointsDataFrame(coords = lonlat, data = lonlatind,
#                                       proj4string  = CRS("+init=epsg:32633"))
# Prec_points <- spTransform(Prec_points, CRS("+init=epsg:32633")) #transform to same system as sub-catchments
# coordinates(Prec_points)
# crs(Prec_points)


 # # # writes shapefile, for control in GIS tool
# plot(subbasins)
# plot(sf_Prec_points[1], add = TRUE, col = 'red')
# plot(sf_Prec_points[1])
# #write_sf(sf_Prec_points, dsn="Pre_nodes_seNorge2018.gpkg")
shell.exec(getwd())


#create function to retrieve dataframe of indices of temp/prec points within a subbasin from the subbasin ID
fn_getIdx <- function(subbasinID){
  spdf_idx <- sf_Prec_points[subbasins[subbasinID,], ]
  idx_list <- cbind(spdf_idx$lonind, spdf_idx$latind)
  colnames(idx_list) <- c("lon", "lat")
  return(idx_list)
}



# #testing
# single_polygon = subbasins[1,]
# plot(single_polygon)
# sf_Prec_points[single_polygon, ]
# 
# 
# subbasinID = 1
# spdf_idx <- Prec_points[subbasins[subbasinID,], ]
# rm(subbasinID)


# # #create list of internal prec/temp points by indexes for all subbasins ( took 400 sec on i7 CPU with 250 subbasins and 1km grid)
tic()
indexlist <- lapply(subbasins$ROWNR, fn_getIdx)
toc()




save.image(file = "Gaula_1polygon_seNorge2018_line140.R")

names(indexlist) <- subbasins$SUBID
rm(subbasins, coords_ind_mat)

#check if all subbasins have at least one prec/temp point within them, also assign closest point
#if ever used when a raster has a very coarse resolution compared to the catchment size, an interpolation scheme should be implemented here
#get list entries with no lat/lon entries or no gridnodes within subbasin polygon respectively
cond_small <- sapply(indexlist, function(x) length(x) < 1)
small_subID <- indexlist[cond_small]
small_subID <- names(small_subID)
rm(cond_small)

# #manual fixing of subcatchments that do not contain any precipitaton grid points
# #check in GIS system first if small subcatchments are meaningful to have (maybe only sliver polygons due to cutting)
# #add entries such as "indexlist$`119`" according to subbasinIDs in vector "small_subID"
# this is poorly hardcoded for the example of two points replacing the empty indexlist of subID 119
# #1 (example)
# indexlist$`119` <- indexlist$`1` #give entry of list same structure as first functionable entry
# indexlist$`119`[1,] <- c(132,261)
# indexlist$`119`[2,] <- c(132,260)
# indexlist$`119` <- indexlist$`119`[1:2,] #give manual length of how many entries you've given in lines before
# rm(small_subID)


#function to get mean of each subcatchment for given timestep by position in indexlist
fn_get_T_avg <- function(subbasinindex){
  avg_T_idx <- round(mean(T_ts[indexlist[[subbasinindex]]]), digits = 2)
     return(avg_T_idx )
}


#index vector from 1 to number of subbasins
index <- c(1:length(indexlist))
#save.image()
# load(".RData")


# # TEST average first subbasin of indexlist into one value
nc$var$rr
ts=3
T_ts = ncvar_get(nc, varid= varname, start = c(1,1,ts), count=c(-1,-1, 1), verbose=FALSE,
                    signedbyte=TRUE, collapse_degen=TRUE)
tic()
avg_T_idx <- round(mean(T_ts[indexlist[[1]]]), digits = 2)
toc()



tic()
avg_all_T_idx <- round(mean(T_ts[indexlist[[1]]]), digits = 2)
toc()

rm(ts)


ts = 1
tic()
for (ts in 1:nc$dim$time$len){#loop through timestep, for each ts a new row with subbasin averaged T is added
    if (ts == 1){
      T_ts <- ncvar_get(nc, varid= varname, start = c(1,1,ts), count=c(-1,-1, 1), verbose=FALSE,
                      signedbyte=TRUE, collapse_degen=TRUE)
      df_T <- sapply(index, fn_get_T_avg)
      df_T <- as.data.frame(t(df_T))
      colnames(df_T) <- names(indexlist)}
    else {
      T_ts <- ncvar_get(nc, varid= varname, start = c(1,1,ts), count=c(-1,-1,1), verbose=FALSE,
                    signedbyte=TRUE, collapse_degen=TRUE)
      new_row <- sapply(index, fn_get_T_avg)
      df_T <- rbind(df_T, new_row)
    }
  if (ts%%100 == 0){
    print(ts)
    }
  }
toc()
rm(new_row, ts)
save.image()
#load(".RData")

#retrieve dates
#check time format
nc$dim$time


timeseries_temp <- as.Date("1900-01-01") + nc$dim$time$vals
tail(timeseries_temp)




timeseries <- as.character(timeseries_temp)

first(timeseries)
last(timeseries)

testts = as.data.frame(timeseries)

df_T <- cbind(timeseries, df_T)
colnames(df_T)[1] <- "SUBID"

#check nc file name
nc$filename
filename =  paste0("Tobs_",str_remove(nc$filename,"./"), "_Gaula_1catchment.txt")
write.table(df_T, file = filename, quote = F, row.names = F, sep = "\t")
file.show(filename)

shell.exec(getwd())

# #extract first two timesteps, and all x-y- values
# tic()
#Prec_ts_test <- ncvar_get(nc, varid= "precipitation", count=c(-1,-1,2), start = c(1,1,1), verbose=FALSE,
#                 signedbyte=TRUE, collapse_degen=TRUE)
# toc()


## open textfiles fast 
# #write meta data description file
# sink("nc_description.txt")
# cat(print(nc))
# sink()
# file.show("nc_description.txt")

# #Later needed for code improvement: automatic fixing of too small subbasins
#domain_extent <- readOGR(dsn = getwd(), layer= "ModelExtent_WGS84", encoding = "UTF8")
#domain_extent <- spTransform(domain_extent, CRSobj = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# centroids_small_subID <- gCentroid(small_subbasins, byid = T)
# centroids_small_subID_zone33 <- spTransform(centroids_small_subID, CRSobj = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# rm(centroids_small_subID)
# 
# T_points_zone33 <- spTransform(T_points, CRSobj = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
# T_points_zone33_clip <- gIntersection(domain_extent, T_points_zone33, byid = TRUE, drop_lower_td = TRUE)
