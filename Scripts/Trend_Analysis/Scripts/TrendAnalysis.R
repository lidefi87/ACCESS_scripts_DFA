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
filelist <- list.files(path = "Data/", pattern = ".csv$", full.names = T)

#Reading files
for(i in filelist){
  var_int <- str_match(i, pattern = "\\/([A-Z]*)\\.csv")[2]
  assign(var_int, read.csv(filelist))
}

Atl_MLD <- read.csv(filelist[1])
CentInd_MLD <- read.csv(filelist[2])
EastInd_MLD <- read.csv(filelist[3])
EP_MLD <- read.csv(filelist[4])
WP_MLD <- read.csv(filelist[5])

#Months into which data will be divided
months1 <- month.abb[5:9]
months2 <- month.abb[10:12]
months3 <- month.abb[1:4]


# data <- Atl_MLD %>% 
# data <- CentInd_MLD %>% 
# data <- EastInd_MLD %>% 
# data <- EP_MLD %>% 
data <- WP_MLD %>% 
  rename("MLD" = "X0") %>% 
  mutate(month = month.abb[month]) %>% 
  filter(month %in% months1) %>% 
  rowid_to_column("obs") %>% 
  mutate(group = (obs-1)%/%length(months1))

fit_LM <- data %$% 
  lm(MLD ~ obs)

summary(fit_LM)
coef_LM <- summary(fit_LM)$coefficients

data_season <- data %>% 
  group_by(group) %>%
  summarise(obs = mean(obs), MLD = mean(MLD))

fit_LM_season <- data_season %$% 
  lm(MLD~obs)

summary(fit_LM_season)
coef_LM_season <-summary(fit_LM_season)$coefficients

data.confInt <- predict(fit_LM_season, interval="confidence", data = data_season) %>% 
  as_tibble() %>% 
  bind_cols(data_season, .)

fit_AR1 <- gls(MLD~obs, correlation = corAR1(), data = data)
coeff_AR1 <- summary(fit_AR1)$tTable

data_SLA.LM <- data.frame(obs = data$obs) %>% 
  mutate(SLA = predict(fit_LM, newdata = .))

data_SLA.MB_LM <- data.frame(obs = data$obs) %>% 
  mutate(SLA = predict(fit_LM_season, newdata = .))

data_SLA.AR1 <- data.frame(obs = data$obs) %>% 
  mutate(SLA = predict(fit_AR1, newdata = .))


plot(MLD~obs,data = data, xlab="Timesteps", ylab = "Mixed layer depth (MLD, m) ",
     pch=19,cex=0.5,col="black", ylim =c(0,350))
polygon(c(data.confInt$obs, rev(data.confInt$obs)), c(data.confInt$lwr, rev(data.confInt$upr)),
        col=adjustcolor("grey80",alpha.f=0.5),border=NA)
lines(MLD~obs,data=data, col="azure3", lwd=0.5)

lines(SLA~obs,data=data_SLA.LM, lwd=1,col="dodgerblue4")
lines(SLA~obs,data=data_SLA.MB_LM, lwd=1,col="chocolate1")
lines(SLA~obs,data=data_SLA.AR1, lwd=1,col="chartreuse")

points(MLD~obs,data=data_season,pch=19,cex=0.8,col="chocolate1")

legend(1, 70,c("Linear regression (LM)","Moving Block LM","Autoregressive (AR1)"), lwd=c(1,1,1),
       col=c("dodgerblue4","chocolate1","chartreuse"), y.intersp=0.8,box.lty=0,bg="transparent",cex=0.8)

text(40,42, paste0("p-value = ", round(coef_LM[2,4], 3)),cex=0.8)
# text(97,35, "p-value < 0.001", cex=0.8)
text(40,20, paste0("p-value = ", round(coef_LM_season[2,4], 3)),cex=0.8)
# text(97,22, "p-value < 0.001", cex=0.8)
text(40,0, paste0("p-value = ", round(coeff_AR1[2,4], 3)),cex=0.8)
# text(15,300, paste0("Trend = ", round(coef_LM_season[2,1]*12, 2),"m/year ±", 
#                      round(coef_LM_season[2,2]*12, 2)),cex=0.8,font=2)
# text(40,340, paste0("Trend = ", round(coef_LM_season[2,1]*12, 2),"m/year ±", 
#                     round(coef_LM_season[2,2]*12, 2)),cex=0.8,font=2)
text(40,70, paste0("Trend = ", round(coef_LM_season[2,1]*12, 2),"m/year ±", 
                    round(coef_LM_season[2,2]*12, 2)),cex=0.8,font=2)
# title("1965-2018 MLD timeseries Atlantic Sector (May-Sep) with fits and estimated trend",cex.main=0.9)
# title("1965-2018 MLD timeseries Central Indian Sector (May-Sep) \n with fits and estimated trend",cex.main=0.9)
# title("1965-2018 MLD timeseries East Indian Sector (May-Sep) \n with fits and estimated trend",cex.main=0.9)
# title("1965-2018 MLD timeseries East Pacific Sector (May-Sep) \n with fits and estimated trend",cex.main=0.9)
title("1965-2018 MLD timeseries West Pacific Sector (May-Sep) \n with fits and estimated trend",cex.main=0.9)


dev.off()

ggplot()+
  geom_point(data = data, aes(obs, MLD))+
  geom_line(data = data, aes(obs, MLD))
  labs(x = "Time (months)", y = "Mixed layer depth (MLD, m)", 
     title = "1965-2018 MLD timeseries with fits and estimated trend: location 1")+
  geom_ribbon(data = data.confInt, aes(ymin = lwr, ymax = upr, x = obs), fill = "grey80", alpha = 0.5)+
  theme_bw()
  # ylim(0, 275)+
  # scale_y_reverse()+
  # theme(panel.grid = element_blank())
  # 




# xRas <- raster(SST, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat))
# plot(xRas)
# 
# data <- data.frame(aice = SST[!is.na(SST)][1:24])
# data$obs <- 1:nrow(data)

# fit_LM <- lm(aice~obs,data=data)
# summary(fit_LM)
# coef_LM<-summary(fit_LM)$coefficients

# %/% divides two numbers and gives you the result as an integer
# data$group <- (data$obs-1)%/%3
# data_season <- data %>%
#   group_by(group) %>%
#   summarise(obs = mean(obs), aice = mean(aice))

# fit_LM_season <- lm(aice~obs, data_season)
# summary(fit_LM_season)

# coef_LM_season <-summary(fit_LM_season)$coefficients

# data.confInt <- data_season
# data.confInt <- cbind(data.confInt, predict(fit_LM_season, interval="confidence", data=data.confInt))


# fit_AR1 <- gls(aice~obs, correlation = corAR1(), data = data)
# coeff_AR1 <- summary(fit_AR1)$tTable

# plot(aice~obs,data = data, xlab="Time (months) ", ylab="Sea level anomaly (cm) ",pch=19,cex=0.5,col="black",ylim=c(-70,70))
# polygon(c(data.confInt$obs, rev(data.confInt$obs)), c(data.confInt$lwr, rev(data.confInt$upr)),
#         col=adjustcolor("grey80",alpha.f=0.5),border=NA)
# lines(aice~obs,data=data, col="azure3", lwd=0.5)
# 
# 
# data_SLA.LM <- data.frame(obs=1:24)
# data_SLA.LM$SLA <- predict(fit_LM,newdata=data_SLA.LM)
# data_SLA.MB_LM <- data.frame(obs=1:24)
# data_SLA.MB_LM$SLA <- predict(fit_LM_season,newdata=data_SLA.MB_LM)
# data_SLA.AR1 <- data.frame(obs=1:24)
# data_SLA.AR1$SLA <- predict(fit_AR1,newdata=data_SLA.AR1)

# lines(SLA~obs,data=data_SLA.LM, lwd=1,col="dodgerblue4")
# lines(SLA~obs,data=data_SLA.MB_LM, lwd=1,col="chocolate1")
# lines(SLA~obs,data=data_SLA.AR1, lwd=1,col="chartreuse")
# 
# points(aice~obs,data=data_season,pch=19,cex=0.8,col="chocolate1")
# 
# legend(10,-50,c("Linear regression (LM)","moving block LM","Autoregressive (AR1)"), lwd=c(1,1,1),
#        col=c("dodgerblue4","chocolate1","chartreuse"), y.intersp=0.8,box.lty=0,bg="transparent",cex=0.8)
# 
# text(17,-55, paste0("p-value=", round(coef_LM[2,4], 3)),cex=0.8)
# text(17,-59, paste0("p-value=", round(coef_LM_season[2,4], 3)),cex=0.8)
# text(17,-63, paste0("p-value=", round(coeff_AR1[2,4], 3)),cex=0.8)
# text(2.75,65, paste0("Trend =", round(coef_LM_season[2,1]*12, 2),"cm/year  ±", round(coef_LM_season[2,2]*12, 2)),cex=0.8,font=2)
# title("1993-2014 SLA timeseries with fits and estimated trend: location 1",cex.main=0.9)
# 
# 
# dev.off()




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



