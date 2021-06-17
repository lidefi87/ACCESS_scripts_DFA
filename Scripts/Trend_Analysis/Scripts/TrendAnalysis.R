################################################################################
# Trend Analyses
# Author: Denisse Fierro Arcos
# Based on script shared by Amelie Meyer
################################################################################

# Loading libraries -------------------------------------------------------
# library(ncdf4)
# library(sf)
library(tidyverse)
library(nlme)
library(magrittr)

# Setting up data into R --------------------------------------------------
#Getting list of csv files
filelist <- list.files(path = "Data/", pattern = ".csv$", full.names = T)[-5]

vars_int <- list()

#Reading files
for(i in seq_along(filelist)){
  vars_int[[i]] <- read.csv(filelist[i])
  names(vars_int)[i] <- str_match(filelist[i], pattern = "\\/(.*)\\.csv")[2]
}

#Months into which data will be divided
months_int <- list(May_Sep = month.abb[5:9],
                   Oct_Dec = month.abb[10:12],
                   Jan_Apr = month.abb[1:4])

for(i in seq_along(vars_int)){
  for(j in seq_along(months_int)){
    data <- vars_int[[i]] %>% 
      rename("MLD" = "X0") %>% 
      mutate(month = month.abb[month]) %>% 
      filter(month %in% months_int[[j]]) %>% 
      rowid_to_column("obs") %>% 
      mutate(group = (obs-1)%/%length(months_int[[j]]))
  
    fit_LM <- data %$% 
      lm(MLD ~ obs)
    
    coef_LM <- summary(fit_LM)$coefficients
    
    data_season <- data %>% 
      group_by(group) %>%
      summarise(obs = mean(obs), MLD = mean(MLD))
    
    fit_LM_season <- data_season %$% 
      lm(MLD~obs)
    
    coef_LM_season <-summary(fit_LM_season)$coefficients
    
    data.confInt <- predict(fit_LM_season, interval = "confidence", data = data_season) %>% 
      as_tibble() %>% 
      bind_cols(data_season, .)
    
    fit_AR1 <- gls(MLD ~ obs, correlation = corAR1(), data = data)
    coeff_AR1 <- summary(fit_AR1)$tTable
    
    data_SLA.LM <- data.frame(obs = data$obs) %>% 
      mutate(SLA = predict(fit_LM, newdata = .))
    
    data_SLA.MB_LM <- data.frame(obs = data$obs) %>% 
      mutate(SLA = predict(fit_LM_season, newdata = .))
    
    data_SLA.AR1 <- data.frame(obs = data$obs) %>% 
      mutate(SLA = predict(fit_AR1, newdata = .))
    
    p <- data %>% 
      ggplot(aes(obs, MLD))+
      geom_point()+
      geom_line(col = "azure3")+
      # ylim(0, max(obs*1.5))+
      geom_ribbon(data = data.confInt, aes(x = obs, ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.4, 
                  show.legend = F)+
      geom_line(data = data_SLA.LM, aes(x = obs, y = SLA, 
                                        colour = sprintf("Linear regression (LM), p-value = %.3f", coef_LM[2,4])))+
      geom_line(data = data_SLA.MB_LM, aes(x = obs, y = SLA, 
                                           colour = sprintf("Moving Block LM, p-value = %.3f", coef_LM_season[2,4])))+
      geom_line(data = data_SLA.AR1, aes(x = obs, y = SLA, 
                                         colour = sprintf("Autoregressive (AR1), p-value = %.3f", coeff_AR1[2,4])))+
      geom_point(data = data_season, aes(x = obs, y = MLD), col = "chocolate1", show.legend = F)+
      scale_color_manual(values = c("dodgerblue4", "chocolate1", "chartreuse"))+
      labs(x = "Timesteps", y = "Mixed layer depth (MLD, m)",
           caption = sprintf("Trend = %.2f m/year ± %.2f", coef_LM_season[2,1]*12, coef_LM_season[2,2]*12))+
      guides(colour = guide_legend(title = element_blank(), ncol = 1))+
      theme_bw()+
      theme(panel.grid = element_blank(), legend.position = "top", legend.box.spacing = unit(0.1, "cm"),
            plot.caption = element_text(hjust = 0))
    ggsave(paste0("Outputs/", names(vars_int)[i], "-", names(months_int)[j], ".png"), p, device = "png")
}}



# Sea ice seasonality -----------------------------------------------------
library(lubridate)
library(cowplot)

#Location of sea ice seasonality data
Inputs <- "C:/Users/ldfierro/OneDrive - University of Tasmania/ACCESS_Outputs/Calculations/SeaIceSeasonality/Yearly/Sectors/"
#Location where outputs will be saved
Outputs <- "C:/Users/ldfierro/OneDrive - University of Tasmania/ACCESS_Outputs/Figures/TimeSeries/SeaIceSeasonality/PolarZones/"

#Get lists of files for sea ice seasons
adv <- list.files(Inputs, pattern = "advance", full.names = T)
ret <- list.files(Inputs, pattern = "retreat", full.names = T)
dur <- list.files(Inputs, pattern = "duration", full.names = T)

#Defining function that will calculate time series trends
Trend_Calculation <- function(filepath, group_yrs, var_name){
  #Load csv file containing sea ice data
  data <- read_csv(filepath) %>% 
    #Rename variable of interest to standarise code
    rename("var" = var_name) %>% 
    #Create a column of unique identifiers per row
    rowid_to_column("obs") %>% 
    #Create groups based on the number of years set as input
    mutate(group = (obs-1)%/%group_yrs)
  
  #Calculate a simple linear model
  fit_LM <- lm(var ~ obs, data)
  #Extract coefficients of linear model
  coef_LM <- summary(fit_LM)$coefficients
  
  ###
  #Group data by group_yrs
  data_season <- data %>% 
    group_by(group) %>%
    summarise(obs = mean(obs), var_grp = mean(var))
  
  #Calculate LM on grouped data
  fit_LM_season <- lm(var_grp ~ obs, data_season)
  #Extract coefficient of LM on grouped data
  coef_LM_season <-summary(fit_LM_season)$coefficients
  
  ###
  #Calculating predictions based on grouped LM results
  data.confInt <- predict(fit_LM_season, interval = "confidence", data = data_season) %>% 
    as_tibble() %>% 
    bind_cols(data_season, .)
  
  ###
  #Calculate Generalised Linear Squares with AR1
  fit_AR1 <- gls(var ~ obs, correlation = corAR1(), data = data)
  #Extract coefficients of GLS
  coeff_AR1 <- summary(fit_AR1)$tTable
  
  ###
  #Calculating predictions based on all model defined above
  data.LM <- data.frame(obs = data$obs) %>% 
    mutate(pred = predict(fit_LM, newdata = .))
  
  data.MB_LM <- data.frame(obs = data$obs) %>% 
    mutate(pred = predict(fit_LM_season, newdata = .))
  
  data.AR1 <- data.frame(obs = data$obs) %>% 
    mutate(pred = predict(fit_AR1, newdata = .))
  
  return(list(raw = data, LM_raw = fit_LM, coef_LM_raw = coef_LM, 
              group_data = data_season, LM_group = fit_LM_season, coef_LM_group = coef_LM_season, 
              confInt = data.confInt, AR1 = fit_AR1, coef_AR1 = coeff_AR1,
              pred_LM = data.LM, pred_group_LM = data.MB_LM, pred_AR1 = data.AR1))
}


#####
#Defining function to plot results
plot_Trends <- function(list, y_lab, units_change, group_yrs, title){
  #Creating figure
  ggplot()+
    geom_point(data = list$raw, aes(obs, var))+
    geom_line(data = list$raw, aes(obs, var), col = "azure3")+
    geom_ribbon(data = list$confInt, aes(x = obs, ymin = lwr, ymax = upr), fill = "grey70", alpha = 0.4, 
                show.legend = F)+
    geom_line(data = list$pred_LM, aes(x = obs, y = pred,
                                       colour = sprintf("Linear regression (LM), p-value = %.3f",
                                                        list$coef_LM_raw["obs", "Pr(>|t|)"])))+
    geom_line(data = list$pred_group_LM, aes(x = obs, y = pred, colour = sprintf("Moving Block LM, p-value = %.3f",
                                                                         list$coef_LM_group["obs", "Pr(>|t|)"])))+
    geom_line(data = list$pred_AR1, aes(x = obs, y = pred, colour = sprintf("Autoregressive (AR1), p-value = %.3f", 
                                                                       list$coef_AR1["obs", "p-value"])))+
    geom_point(data = list$group_data, aes(x = obs, y = var_grp), col = "chocolate1", show.legend = F)+
    scale_color_manual(values = c("dodgerblue4", "chocolate1", "chartreuse"))+
    labs(x = "Timesteps", y = y_lab, title = title,
         caption = sprintf("Trend = %.2f ± %.2f %s", list$coef_LM_group["obs", "Estimate"]*group_yrs, 
                           list$coef_LM_group["obs", "Std. Error"]*group_yrs, units_change))+
    guides(colour = guide_legend(title = element_blank(), nrow = 2))+
    theme_bw()+
    theme(panel.grid = element_blank(), legend.position = "top", legend.box.spacing = unit(0.1, "cm"),
          plot.caption = element_text(hjust = 0, size = 12), legend.text = element_text(size = 12), 
          axis.text = element_text(size = 11), axis.title = element_text(size = 12), title = element_text(size = 14))
}


####
#Advance
#Creating an empty list to save figures
figs <- list()
#Create an empty vector to standarise y axis limits of all plots
lims <- vector()

#Applying functions
for(i in seq_along(adv)){
  out <- Trend_Calculation(adv[i], group_yrs = 10, var_name = "aice")
  sector <- str_match(adv[i], ".*_(.*)_\\d{4}-\\d{4}")[,2]
  figs[[sector]] <- plot_Trends(out, "Sea ice advance (Days from February 15)", "days/decade", group_yrs = 10, 
                                title = paste0(sector, " sector"))
  lims <- append(lims, layer_scales(figs[[i]])$y$range$range)
}

#Get y limits to standarise them across figures before saving
lims <- c(min(lims), max(lims))
for(i in seq_along(figs)){
  figs[[i]] <- figs[[i]]+ylim(lims)
  fileout <- paste0(Outputs, "SeaIceAdvance_", str_match(adv[i], ".*_(.*_\\d{4}-\\d{4})")[,2], ".png")
  ggsave(fileout, figs[[i]], device = "png")
  }

####
#Retreat
#Creating an empty list to save figures
figs <- list()
#Create an empty vector to standarise y axis limits of all plots
lims <- vector()

#Applying functions
for(i in seq_along(ret)){
  out <- Trend_Calculation(ret[i], group_yrs = 10, var_name = "aice")
  sector <- str_match(ret[i], ".*_(.*)_\\d{4}-\\d{4}")[,2]
  figs[[sector]] <- plot_Trends(out, "Sea ice retreat (Days from February 15)", "days/decade", group_yrs = 10, 
                                title = paste0(sector, " sector"))
  lims <- append(lims, layer_scales(figs[[i]])$y$range$range)
}

#Get y limits to standarise them across figures before saving
lims <- c(min(lims), max(lims))
for(i in seq_along(figs)){
  figs[[i]] <- figs[[i]]+ylim(lims)
  fileout <- paste0(Outputs, "SeaIceRetreat_", str_match(ret[i], ".*_(.*_\\d{4}-\\d{4})")[,2], ".png")
  ggsave(fileout, figs[[i]], device = "png")
}


####
#Duration
#Creating an empty list to save figures
figs <- list()
#Create an empty vector to standarise y axis limits of all plots
lims <- vector()

#Applying functions
for(i in seq_along(dur)){
  out <- Trend_Calculation(dur[i], group_yrs = 10, var_name = "aice")
  sector <- str_match(dur[i], ".*_(.*)_\\d{4}-\\d{4}")[,2]
  figs[[sector]] <- plot_Trends(out, "Sea ice duration (Total number of days)", "days/decade", group_yrs = 10, 
                                title = paste0(sector, " sector"))
  lims <- append(lims, layer_scales(figs[[i]])$y$range$range)
}

#Get y limits to standarise them across figures before saving
lims <- c(min(lims), max(lims))
for(i in seq_along(figs)){
  figs[[i]] <- figs[[i]]+ylim(lims)
  fileout <- paste0(Outputs, "SeaIceDuration_", str_match(dur[i], ".*_(.*_\\d{4}-\\d{4})")[,2], ".png")
  ggsave(fileout, figs[[i]], device = "png")
}


# Loading data ------------------------------------------------------------
SeaIceAdv <- nc_open("Yearly/SeaIceAdv_1965-1966.nc")
print(SeaIceAdv)

#Get latitude, longitude and time
lat <- ncvar_get(SeaIceAdv, "yt_ocean")
lat

lon <- ncvar_get(SeaIceAdv, "xt_ocean")
lon
i
# time <- ncvar_get(SeaIceAdv, "time")
time <- str_extract(filelist[1], pattern = "[0-9]{4}-[0-9]{4}")
time

#Get variable values
ice_season <- ncvar_get(SeaIceAdv, "aice")
fill_value <- ncatt_get(SeaIceAdv, "aice", "_FillValue")


#Closing netcdf connection
nc_close(SeaIceAdv)


# Data manipulation -------------------------------------------------------
#Setting fill values to NAs
ice_season[ice_season == paste(fill_value$value)] <- NA

as.data_frame(ice_season) %>% 
  ggplot(aes())

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(ylim = c(-90, -50), xlim = c(-180, 180), expand = FALSE)#, crs = "EPSG:3976")




library(ncdf4)
library(raster)
library(tidyverse)

setwd("C:/Users/ldfierro/OneDrive - University of Tasmania/ACCESS_Outputs/Calculations/SeaIceSeasonality")
filelist <- list.files(pattern = ".nc$")

x <- nc_open(filelist[1])
x
lat <- ncvar_get(x, "yt_ocean", verbose = F)
lon <- ncvar_get(x, "xt_ocean")
SST <- ncvar_get(x, "aice")
time <- 2018

fillvalue <- ncatt_get(x, "aice", "_FillValue")
fillvalue

SST[SST == "NaN"] <- NA
SST



