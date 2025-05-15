
rm(list=ls())


library(ggplot2)
library(dplyr) 
library(viridis) 
library(Interpol.T) 
library(lubridate) 
library(ggExtra) 
library(tidyr) 
library(tidyverse)
library(timetk)


setwd("YourPath")


#-----------------------------------------------------------------------------------------


data<- read.csv("WaterLevel19.csv",header=T)

data$Day <- as.POSIXct((paste(data$Date, data$Time)), format="%Y-%m-%d %H:%M:%S")

data <- na.omit(data)

# Setup for the plotly charts (# FALSE returns ggplots)
interactive <- FALSE

p1 <-data %>% 
  plot_time_series(Day,WaterLevel,
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (30-min intervals)",
                   .y_lab = "Water level (m) ",)
b1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "WaterLevel.svg", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")

p1<-data %>% 
  plot_time_series(Day,Flow,
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .facet_scales = "free", 
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (30-min intervals)",
                   .y_lab = "Flow rate (l/s) ",)
b1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "Flow.png", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")

#------------------------ conductivity/DO------------------------

data<- read.csv("Conductivity19.csv",header=T)

data$Day <- as.POSIXct((paste(data$Date, data$Time)), format="%Y-%m-%d %H:%M:%S")

data <- na.omit(data)

interactive <- FALSE

b1<-data %>% 
  plot_time_series(Day,Conductivity,
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .facet_ncol = 2, 
                   .facet_scales = "free", 
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (2-min intervals)",
                   .y_lab = "Conductivity (mS/cm) ",)

b1 <- b1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "Conductivity.png", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")

p1 <- data %>% 
  plot_time_series(Day,(DO_mgl),
                   .color_var = year(Day),      # Color applied to Week transformation
                   # Facet formatting
                   .facet_scales = "free", 
                   .interactive = interactive,
                   .title = NULL,
                   .legend_show = FALSE,
                   .x_lab = "Date (2-min intervals)",
                   .y_lab = "Dissolved oxygen (mg/l) ",)
b1 <- p1 + theme_bw() + theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
  theme(legend.position="none")+
  theme(text=element_text(family="Calibri", face="bold", size=11))

b1

ggsave(filename = "DO.png", plot = b1, width = 14, height = 10, dpi = 300, units = "cm")