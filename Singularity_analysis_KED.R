########################################################################################################################
##########                    This scripts performs KED merging with singularity analysis                       ########
########################################################################################################################
rm(list = ls())  # clean memory
gc()
########################################           Input              ###################################################
# Rain gauge file
RGfile <- "/home/fc14509/UK/P3/Dataset/EA_RG/EA"
RGdates <- "/home/fc14509/UK/P3/Dataset/EA_RG/dates"
RGcoord <- "/home/fc14509/UK/P3/Dataset/EA_RG/coord"

# Radar folder
Radfolder <- "/home/fc14509/UK/P3/Dataset/Radar/"
Radcoord <- "/home/fc14509/UK/P3/Dataset/Radar/alex_catchment_xy.txt"

# Projection parameters in proj4 format
myproj <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"

# Transformed variable folder
Transformed <- "/home/fc14509/UK/P3/Transformed/Singularity Analysis/"

# Result folder
Results <- "/home/fc14509/UK/P3/Results/Singularity Analysis/"

###################################### load libraries  ####################################
library(sp)
library(maptools)
library(gstat)
library(spacetime)
library(raster)
library(minpack.lm)
library(lubridate)
library(plyr)
library(gdata)
library(ggplot2)
library(plotKML)
library(RSAGA)
library(abind)
source("Negentropy.R")
###########################################################################################
# load rain gauge data, coord, and dates
load(RGfile)
RG <- EA

RG_points <- dim(RG)[2]
RG_steps <- dim(RG)[1]
RG[RG<0] <- 0
RG_ns <- RG

load(RGdates)
RGdates <- dates

load(RGcoord)
RGcoord <- coord

toremove <- which(RGcoord[,1]<320000 | RGcoord[,1]>520000 | RGcoord[,2]<350000 | RGcoord[,2]>550000)
RG <- RG[,-toremove]
RGcoord <- RGcoord[-toremove,]

rm(EA,dates,coord)

# loop on the monthly radar files
for(i in 1:12){
  
  #load the radar data
  radfile <- paste0(Radfolder, "radar2009", sprintf("%02d",i), ".txt")
  rad <- read.table(radfile)
  
  raddates <- rad[,1:6]
  rad <- rad[,7:dim(rad)[2]]
  
  radcoord <- read.table(Radcoord)
  
  rad_steps <- dim(rad)[1]
  
  rad_ns <- rad
  ked_ns <- rad
  
  #loop on the available time steps
  for (j in 1:rad_steps) {
    
    # Prepare the spatial data
    R <- data.frame(t(rad[j,]), radcoord)
    names(R) <- c("rain","x","y")
    dat <- raddates[j,]
    RG_position <- which(RGdates[,1]==as.numeric(dat[1]) & RGdates[,2]==as.numeric(dat[2]) & RGdates[,3]==as.numeric(dat[3]) & RGdates[,4]==as.numeric(dat[4]))
    G <- data.frame(t(RG[RG_position,]),RGcoord)
    names(G) <- c("rain","x","y") 
    G <- G[complete.cases(G),]
    coordinates(G) <- ~x+y
    G@proj4string@projargs <- myproj
    
    coordinates(R) <- ~x+y
    bb <- R@bbox
    bb[,2] <- bb[,2]+1000
    R <- vect2rast(R, cell.size=1000, bbox=bb)
    R@proj4string@projargs <- myproj
    
    #loop on each pixel
    X <- seq(from=R@bbox[1,1]+500, to=R@bbox[1,2]-500, by=1000)
    Y <- seq(from=R@bbox[2,1]+500, to=R@bbox[2,2]-500, by=1000)
    res <- c(1,3,5,7,9)
    eps <- res/max(res)
    rho <- matrix(0,1,length(res))
    
    R_ns <- R
    cr <- 0
    xside <- R_ns@grid@cells.dim[1]
    yside <- R_ns@grid@cells.dim[2]
    l <- matrix(0,xside, yside)
    for (r in res){
      cr <- cr+1
      
      # to make the areal average more efficient, we define a larger matrix with more cells on the external frame
      # and we move the original matrix along x and y summing the cells at each step 
      # then we divide by the total mumber of steps taken
      bigm <- matrix(0, xside+2*(r-1),yside+2*(r-1))
      movex <- xside+2*(r-1) - xside +1
      movey <- yside+2*(r-1) - yside +1
      for (x in 1:movex){
        for (y in 1:movey){
          bigm[x:(x-1+xside),y:(y-1+yside)] <- bigm[x:(x-1+xside),y:(y-1+yside)] + as.matrix(R_ns)
        }
      }
      bigm <- bigm/(movex*movey)
      l<- abind(l,bigm[r:(xside+(r-1)),r:(yside+(r-1))], along=3)
    }
    l <- l[,,2:dim(l)[3]] # just to exclude the first layer, defined as zeros only to create the stack
    cp <- 0
    
    #loop on each pixel
    for (x in 1:xside){
      for (y in 1:yside){
        cp <- cp+1
        rho <- log(as.numeric(l[x,y,]))
        df <- data.frame(log(eps),as.numeric(rho))
        names(df) <- c("eps", "rho")
        regression <- lm(rho~eps, data=df)
        R_ns$rain[cp] <- coefficients(regression)[1]
      }
    }
  # Save teh results
  rad_ns[j,] <- R_ns$rain
  
  # Apply the correction to Rain Gauges as well
  ratio <- (over(G,R_ns)+0.0001)/(over(G,R)+0.0001)   #simple trick to avoid the division by 0 and excessive ratio values, it shouldn't really affect the other results
  G$rain <- G$rain*ratio
  G$radar <- over(G,R_ns)$rain
  RG_ns[RG_position,is.na(RG_ns[RG_position,])==F] <- RG_ns[RG_position,is.na(RG_ns[RG_position,])==F]*t(ratio)
  
  # Check for zeros in the radar at RG locations (indeterminate system)
  checkzeros <- sum(G$rain)+sum(G$radar, na.rm=T)  
  rgcheckzeros <- sum(G$rain)
  radcheckzeros <- sum(G$radar, na.rm=T)
  if (checkzeros==0) { #Case in which both the radar and the gauges say it didn't rain: pred = 0, var = NA
    ked <- R_ns
    names(ked) <- c("pred")
    ked$pred <- ked$pred*0
  } else if (!(rgcheckzeros==0) & radcheckzeros==0) { # case in which it is impossible to do UK and we do OK
    
    # variogram calculation
    v <- variogram(rain~1, G)
    v <- fit.variogram(v, vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp"))
    if(v[2,2]==0){  # If the fitted sill is zero, we change it to the mean of the observed variogram
      v <- variogram(rain~1, G)
      v=vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp")
    }
    if (v[1,2]>v[2,2]){    # if the nugget is larger than the sill, we set them equal
      v[1,2] <- v[2,2]
    }
    
    # Ordinary kriging
    pred_grid <- R_ns
    names(pred_grid) <- "radar"
    ked <- krige(rain~1, G, newdata=pred_grid, v, na.action = na.pass)
    names(ked) <- c("pred", "var")
  } else {
    
    # residual scaling for the KED variogram
    if (sum((G$radar - mean(G$radar))^2)==0){
      r2 <- 1
    } else {
      r2 <- 1-sum((G$radar-G$rain)^2)/sum((G$radar - mean(G$radar))^2)
    }
    
    # variogram calculation
    v <- variogram(rain~1, R_ns)
    v <- fit.variogram(v, vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp"))
    if(v[2,2]==0){
      v <- variogram(rain~1, R_ns)
      v=vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp")
    }
    if (v[1,2]>v[2,2]){
      v[1,2] <- v[2,2]
    }
    
    v[2,2] <- v[2,2]*(1-r2)
    
    # Universal kriging
    pred_grid <- R_ns
    names(pred_grid) <- "radar"
    ked <- krige(rain~radar, G, newdata=pred_grid, v, na.action = na.pass)
    names(ked) <- c("pred", "var")
  }

  ked_ns[j,] <- ked$pred
  }
  save(rad_ns, file=paste0(Transformed, "radar2009", sprintf("%02d",i), "_ns"))
  save(ked_ns, file=paste0(Results, "ked2009", sprintf("%02d",i), "_ns"))
}

save(RG_ns, file=paste0(Transformed, "EA_ns"))

