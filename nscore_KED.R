#########################################################################################################################
#####              This scripts transforms radar and rain gauge data with nscore transformations               #########
#########################################################################################################################
rm(list = ls())  # clean memory
gc()
########################################           Input              ###################################################

# Rain gauge files
RGfile <- "/home/fc14509/UK/P3/Dataset/EA_RG/EA"
RGdatefile <- "/home/fc14509/UK/P3/Dataset/EA_RG/dates"
RGcoordfile <- "/home/fc14509/UK/P3/Dataset/EA_RG/coord"

# Radar folder
Radfolder <- "/home/fc14509/UK/P3/Dataset/Radar/"
Radcoord <- "/home/fc14509/UK/P3/Dataset/Radar/alex_catchment_xy.txt"

# Transformed data folder
Transformed <- "/home/fc14509/UK/P3/Transformed/"

# KED folder
Results <- "/home/fc14509/UK/P3/Results/Nscore/"

# Projection parameters in proj4 format
myproj <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"

######################################            load libraries           #############################################
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
source("nscore.R")
#################################################################################################################################

# Load Rain Gauge data
load(RGfile)
RG <- EA

load(RGdatefile)
RGdates <- dates

load(RGcoordfile)
RGcoord <- coord

# Remove points outside the study area
toremove <- which(RGcoord[,1]<320000 | RGcoord[,1]>520000 | RGcoord[,2]<350000 | RGcoord[,2]>550000)
RG <- RG[,-toremove]
RGcoord <- RGcoord[-toremove,]

RG_points <- dim(RG)[2]
RG_steps <- dim(RG)[1]
RG[RG<0] <- 0

rm(EA,dates,coord)

# Load radar coordinates
radcoord <- read.table(Radcoord)

# Loop on the months
for (m in 1:12){
  radfile <- paste0(Radfolder, "radar2009", sprintf("%02d",m), ".txt")
  Rad <- read.table(radfile)
  
  R <- Rad[,7:dim(Rad)[2]]
  raddates <- Rad[,1:6]
  
  n_points <- dim(R)[2]
  n_steps <- dim(R)[1]

  R_ns <- Rad
  ked_ns <- Rad
  RG_ns <- matrix(NA, dim(Rad)[1], dim(RG)[2])
  
  for (t in 1:n_steps){
    
    # Prepare the spatial data
    r <- R[t,]
    dat <- raddates[t,]
    
    # calculates nscores
    vec <- as.vector(t(r))                  
    ns <- qqnorm(vec, plot.it = FALSE)
    ns <- data.frame(ns)
    names(ns) <- c("nscore", "rain")
    
    # eliminate all the duplicates
    x <- sort(unique(ns$rain))
    y <- matrix(0,length(x), 1)
    for (i in 1:length(x)){
      ns_val <- ns[which(ns[,"rain"]==x[i]),"nscore"]           # Select only the scores associated with the repeated value
      y[i] <- median(ns_val, na.rm=T)    # Use the median as value for all the repeated values, following Lien, 2013
    }
    
    # create a lookup table for this time step
    lookup <- data.frame(rain=x,nscore=y)
    
    # convert the radar
    r_ns <- lookup[match(r, lookup[,1]),2]
    R_ns[t,7: dim(R_ns)[2]] <- r_ns

    # Select corresponding rain gauges
    RG_position <- which(RGdates[,1]==as.numeric(dat[1]) & RGdates[,2]==as.numeric(dat[2]) & RGdates[,3]==as.numeric(dat[3]) & RGdates[,4]==as.numeric(dat[4]))
    g <- RG[RG_position,]
    c <- RGcoord[is.na(g)==F,]
    g <- g[is.na(g)==F]
    
    # create a lookup table for rain gauges only at this time step
    g_value <- sort(unique(g))
    g_score <- g_value
    for (i in 1:length(g_value)){
      if (is.na(match(g_value[i], lookup[,1]))==F){
        g_score[i] <- lookup[match(g_value[i], lookup[,1]),2]
      } else if (sum(g_value[i]>lookup[,2])==0){
        g_score[i] <- lookup[1,2] + (lookup[2,2]-lookup[1,2])/(lookup[2,1]-lookup[1,1])*(g_value[i]-lookup[1,1])
      } else if (sum(g_value[i]>lookup[,2])==length(g_value)){
        n <- dim(lookup)[1]
        g_score[i] <- lookup[n,2] + (lookup[(n-1),2]-lookup[n,2])/(lookup[(n-1),1]-lookup[n,1])*(g_value[i]-lookup[n,1])
      } else {
        n <- sum(lookup[,1]<g_value[i])
        g_score[i] <- lookup[(n+1),2] + (lookup[n,2]-lookup[(n+1),2])/(lookup[n,1]-lookup[(n+1),1])*(g_value[i]-lookup[(n+1),1])
      }
    }
    lookupRG <- data.frame(rain=g_value,nscore=g_score)
    
    # convert the rain gauge data
    g_ns <- lookupRG[match(g, lookupRG[,1]),2]
    RG_ns[t,which(is.na(RG[t,])==F)] <- g_ns
    
    # Prepare data for merging
    g_ns <- data.frame(g_ns, c)
    names(g_ns) <- c("rain","x","y")
    coordinates(g_ns) <- ~x+y
    g_ns@proj4string@projargs <- myproj
    
    r_ns <- data.frame(r_ns, radcoord)
    names(r_ns) <- c("rain","x","y")
    coordinates(r_ns) <- ~x+y
    bb <- r_ns@bbox
    bb[,2] <- bb[,2]+1000
    r_ns <- vect2rast(r_ns, cell.size=1000, bbox=bb)
    r_ns@proj4string@projargs <- myproj
    
    g_ns$radar <- over(g_ns,r_ns)$rain
    
    # Check for zeros in the radar at RG locations (indeterminate system)
    zeroval <- lookup[lookup[,1]==0, 2]              # nscore of 0
    checkzeros <- mean(g_ns$rain, na.rm=T)+mean(g_ns$radar, na.rm=T)  
    rgcheckzeros <- mean(g_ns$rain, na.rm=T)
    radcheckzeros <- mean(g_ns$radar, na.rm=T)
    if (checkzeros==zeroval) { # Case in which both the radar and the gauges say it didn't rain: pred = 0
      ked <- r_ns
      names(ked) <- c("pred")
      ked$pred <- ked$pred*0
    } else if (!(rgcheckzeros==zeroval) & radcheckzeros==zeroval) { # case in which it is impossible to do UK and we do OK
      
      # variogram calculation
      v <- variogram(rain~1, g_ns)
      v <- fit.variogram(v, vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp"))
      if(v[2,2]==0){
        v <- variogram(rain~1, g_ns)
        v=vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp")
      }
      if (v[1,2]>v[2,2]){
        v[1,2] <- v[2,2]
      }
      
      # Ordinary kriging
      pred_grid <- r_ns
      names(pred_grid) <- "radar"
      ked <- krige(rain~1, g_ns, newdata=pred_grid, v, na.action = na.pass)
      names(ked) <- c("pred", "var")
    } else {
      
      # residual scaling for the KED variogram
      if (sum((g_ns$radar - mean(g_ns$radar))^2)==0){
        r2 <- 1
      } else {
        r2 <- 1-sum((g_ns$radar-g_ns$rain)^2)/sum((g_ns$radar - mean(g_ns$radar))^2)
      }
      
      # variogram calculation
      v <- variogram(rain~1, r_ns)
      v <- fit.variogram(v, vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp"))
      if(v[2,2]==0){
        v <- variogram(rain~1, r_ns)
        v=vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp")
      }
      if (v[1,2]>v[2,2]){
        v[1,2] <- v[2,2]
      }
      
      v[2,2] <- v[2,2]*(1-r2)
      
      # Universal kriging
      pred_grid <- r_ns
      names(pred_grid) <- "radar"
      ked <- krige(rain~radar, g_ns, newdata=pred_grid, v, na.action = na.pass)
      names(ked) <- c("pred", "var")
    }
    ked_ns[t,7:dim(ked_ns)[2]] <- ked$pred
  }
  save(R_ns, file=paste0(Transformed,"radar2009", sprintf("%02d",m),"_nscore"))
  save(RG_ns, file=paste0(Transformed,"rg_2009", sprintf("%02d",m),"_nscore"))
  save(ked_ns, file=paste0(Results,"ked_2009", sprintf("%02d",m),"_nscore"))
}
