# Map figure
library(ggplot2)
library(dplyr)
library(mapdata)
library(RColorBrewer)
library(maptools)

setwd("~/Dropbox/Projects/population_genomics/maps")
sampleData <- read.csv("metadata.csv", header = TRUE, 
                       na.strings = c("", "NA"))
source("scalebars.R")
sampleData <- filter(sampleData, Continent=="N America",
                     Environment=="Copepod", 
                     !(is.na(Code)), 
                     !(Sample %in% c("VE", "V2E", "MIE", "TBE", "POE", "BBE1", "BBE2")))

redfresh <- filter(sampleData, Region=="Northeast", WaterType=="Fresh")
rf <- brewer.pal(6, "Paired")[5]
redsalt <- filter(sampleData, Region=="Northeast", WaterType=="Salt")
rs <- brewer.pal(6, "Paired")[6]
#change SJPE point for clearer mapping. Will add line pointing to exact location
redsalt$Latitude[2] <- 48
redsalt$Longitude[2] <- -70.8

theme_nothing <- function(base_size = 12, base_family = "Helvetica")
{
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            line             = element_blank(),
            text             = element_blank()
        )
}

usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")
Mexico <- map_data("worldHires", "Mexico")

shape <- readShapeLines("ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
shape@data$id <- rownames(shape@data)
watershedPoints <- fortify(shape, region = "id")
watershedDF <- merge(watershedPoints, shape@data, by = "id")
watershedDF <- filter(watershedDF, name %in% c("Mississippi", "Missouri", "Ohio", "St. Lawrence", "N. Fk. Red"))

NAmap <- ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "white", 
                                 color="black") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "white", color="black") + 
    geom_polygon(data = Mexico, aes(x=long, y = lat, group = group), fill = "white", color="black") +
    geom_path(data=watershedDF, aes(x = long, y = lat, group = group), color = 'black', size=0.5) +
    coord_fixed(1.3) + coord_fixed(xlim = c(-100, -65),  ylim = c(25, 50), ratio = 1.3)
NAmap <- NAmap + geom_point(data=redfresh, aes(x=Longitude, y=Latitude), fill=rf, color = "black", shape=21, size=8.9) +
    geom_point(data=redsalt, aes(x=Longitude, y=Latitude), fill=rs, color = "black", shape=21, size=8.9) +
    geom_point(data=greenfresh, aes(x=Longitude, y=Latitude), fill=gf, color = "black", shape=21, size=8.9) +
    geom_point(data=greensalt, aes(x=Longitude, y=Latitude), fill=gs, color = "black", shape=21, size=8.9)
pdf("NAmerica.pdf",width=9,height=9, useDingbats = FALSE)
NAmap + theme_nothing() + theme(panel.background = element_rect(fill = "grey88")) +
    scaleBar(lon = -95.5, lat = 25, distanceLon = 500, distanceLat = 100, distanceLegend = 200, dist.unit = "km", legend.size = 7.5, orientation = FALSE)
dev.off()


purpfresh <- filter(europedata, WaterType=="Fresh")
pf <- brewer.pal(10, "Paired")[9]
purpsalt <- filter(europedata, WaterType=="Salt")
ps <- brewer.pal(10, "Paired")[10]

europe <- map_data("worldHires")
EUmap <- ggplot() + geom_polygon(data = europe, aes(x=long, y = lat, group = group), 
                                 fill = "white", color="black") +
    coord_fixed(1.3) + coord_fixed(xlim = c(3, 8.6),  ylim = c(50, 54), ratio = 1.3)
EUmap <- EUmap + geom_point(data=purpfresh, aes(x=Longitude, y=Latitude), fill=pf, color = "black", shape=21, size=8.9) +
    geom_point(data=purpsalt, aes(x=Longitude, y=Latitude), fill=ps, color = "black", shape=21, size=8.9)
EUmap <- EUmap + theme_nothing() + theme(panel.background = element_rect(fill = "grey88")) +
    scaleBar(lon = 3, lat = 50.25, distanceLon = 50, distanceLat = 20, distanceLegend = 40, dist.unit = "km", legend.size = 7.5, orientation = FALSE)
pdf("Europe.pdf",width=9,height=9, useDingbats = FALSE)
EUmap
dev.off()

EUfull <- ggplot() + geom_polygon(data = europe, aes(x=long, y = lat, group = group), 
                                  fill = "white", color="black") +
    coord_fixed(1.3) + coord_fixed(xlim = c(-10, 23),  ylim = c(36, 58), ratio = 1.3) + 
    theme_nothing() + theme(panel.background = element_rect(fill = "grey88")) +
    annotate("rect", xmin = 3, xmax = 8.6, ymin = 49.3, ymax = 53.5, size=2, color = "red", fill="transparent")
pdf("EUFullMap.pdf", width=9, height=9, useDingbats = FALSE)
EUfull
dev.off()
