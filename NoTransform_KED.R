#########################################################################################################################
#####               This scripts Merges untransformed rain gauge with BoxCox transformedradar data              #########
#########################################################################################################################
rm(list = ls())  # clean memory
gc()
########################################           Input              ###################################################

# Radar folder
Folder <- "/home/fc14509/UK/P3/Dataset/Radar/"

# RG date and coord files
RGfile <- "/home/fc14509/UK/P3/Dataset/EA_RG/EA"
RGcoordfile <- "/home/fc14509/UK/P3/Dataset/EA_RG/coord"
RGdatesfile <- "/home/fc14509/UK/P3/Dataset/EA_RG/dates"

#Radar coord file
Radcoordfile <- "/home/fc14509/UK/P3/Dataset/Radar/alex_catchment_xy.txt"

# Projection parameters in proj4 format
myproj <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"

# Result folder
Results <- "/home/fc14509/UK/P3/Results/NoTransform/"

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
###########################################################################################

#Load radar coordinates
Radcoord <- read.table(Radcoordfile)

# Load the rain gauge file
load(RGfile)
RG <- EA
RG[RG<0] <- NA
load(RGcoordfile)
load(RGdatesfile)
RGdates <- dates
RGcoord <- coord

# remove RG out of the box
RG <- RG[,-which(RGcoord[,1]<320000 | RGcoord[,1]>520000 | RGcoord[,2]<350000 | RGcoord[,2]>550000)]
RGcoord <- RGcoord[-which(RGcoord[,1]<320000 | RGcoord[,1]>520000 | RGcoord[,2]<350000 | RGcoord[,2]>550000),]

# loop on the radar months
for (m in 1:12) {
  
  #Load radar file with dates
  radfile <- paste0(Folder, "Rad_2009",sprintf("%02d", m), ".txt")
  R <- read.tavle(radfile)
  ked_bc <- R
  Raddates <- R[,1:6]
  R <- R[,7:dim(R)[2]]
  t_steps <- dim(R)[1]
  
  # Loop on the available time steps
  for (t in 1:t_steps){
    
    # Select radar and rain gauges for the right time step
    r <- R[t,]
    dat <- Raddates[t,]
    RG_pos <- which(RGdates[,1]==as.numeric(dat[1]) & RGdates[,2]==as.numeric(dat[2]) & RGdates[,3]==as.numeric(dat[3]) & RGdates[,4]==as.numeric(dat[4]))
    g <- RG[RG_pos,]
    
    # prepare the data for merging
    r <- data.frame(t(r), Radcoord)
    names(r) <- c("rain", "x", "y")
    coordinates(r) <- ~x+y
    bb <- r@bbox
    bb[,2] <- bb[,2]+1000
    r <- vect2rast(r, cell.size=1000, bbox=bb)
    r@proj4string@projargs <- myproj
    
    g <- data.frame(g, RGcoord)
    g <- g[complete.cases(g),]
    names(g) <- c("rain", "x", "y")
    coordinates(g) <- ~x+y
    g@proj4string@projargs <- myproj
    
    # attach the radar to the rain gauge data frame
    g$radar <- over(g,r)$rain
    
    # prediction grid
    pred_grid <- r
    names(pred_grid) <- "radar"
    
    # Check for zeros in the radar at RG locations (indeterminate system or no solution)
    checkzeros <- (mean(g$rain)+mean(g$radar, na.rm=T))/2
    rgcheckzeros <- mean(g$rain)
    radcheckzeros <- mean(g$radar, na.rm=T)
    if (checkzeros==0) { #Case in which both the radar and the gauges say it didn't rain: pred = 0, var = NA
      ked <- pred_grid
      names(ked) <- c("pred")
      ked$pred <- ked$pred*0
    } else if (!(rgcheckzeros==0) & radcheckzeros==0) {  # case in which it is impossible to do UK and we do OK. Variogram is calculated on rg
      
      # variogram calculation
      v <- variogram(rain~1, g)
      v <- fit.variogram(v, vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp"))
      if(v[2,2]==0){
        v <- variogram(rain~1, g)
        v=vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp")
      }
      if (v[1,2]>v[2,2]){
        v[1,2] <- v[2,2]
      }
      
      # Finally Merge!!!
      ked <- krige(rain~1, g, newdata=pred_grid, v, na.action = na.pass)
      names(ked) <- c("pred", "var")
      
    } else {    # this is the proper ked
      
      # residual scaling for the KED variogram
      if (sum((g$radar - mean(g$radar))^2)==0){
        r2 <- 1
      } else {
        r2 <- 1-sum((g$radar-g$rain)^2)/sum((g$radar - mean(g$radar))^2)
      }
      
      # variogram calculation
      v <- variogram(rain~1, r)
      v <- fit.variogram(v, vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp"))
      if(v[2,2]==0){
        v <- variogram(rain~1, r)
        v=vgm(nugget=0, psill=mean(v$gamma), range=30000, model="Exp")
      }
      if (v[1,2]>v[2,2]){
        v[1,2] <- v[2,2]
      }
      
      v[2,2] <- v[2,2]*(1-r2)
      
      # Finally Merge!!!
      ked <- krige(rain~radar, g, newdata=pred_grid, v, na.action = na.pass)
      names(ked) <- c("pred", "var")
    }
    ked_bc[t,7:dim(ked_bc)[2]] <- ked$pred
  }
  save(ked, file=paste0(Results, "ked_2009",sprintf("%02d", m), "_NoTransf"))
}

