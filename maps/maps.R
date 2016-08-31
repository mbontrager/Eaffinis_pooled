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
                     !(Sample %in% c("VE", "V2E", "MIE", "TBE", "POE", "BBE1", "BBE2", "SJE")))

redfresh <- filter(sampleData, Region=="Northeast", WaterType=="Fresh")
rf <- brewer.pal(6, "Paired")[5]
redsalt <- filter(sampleData, Region=="Northeast", WaterType=="Salt")
rs <- brewer.pal(6, "Paired")[6]

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

shape <- readShapeLines("ne_50m_rivers_lake_centerlines/ne_50m_rivers_lake_centerlines.shp")
shape@data$id <- rownames(shape@data)
watershedPoints <- fortify(shape, region = "id")
watershedDF <- merge(watershedPoints, shape@data, by = "id")
watershedDF <- filter(watershedDF, name %in% c("St. Lawrence"))

NAmap <- ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "white", 
                                 color="black") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "white", color="black") + 
    geom_path(data=watershedDF, aes(x = long, y = lat, group = group), color = 'black', size=0.5) +
    coord_fixed(1.3) + coord_fixed(xlim = c(-92.5, -65),  ylim = c(40, 50), ratio = 1.3)
NAmap <- NAmap + geom_point(data=redfresh, aes(x=Longitude, y=Latitude), fill=rf, color = "black", shape=21, size=8.9) +
    geom_point(data=redsalt, aes(x=Longitude, y=Latitude), fill=rs, color = "black", shape=21, size=8.9)
pdf("NAmerica.pdf",width=9,height=4.5, useDingbats = FALSE)
NAmap + theme_nothing() + theme(panel.background = element_rect(fill = "grey88")) +
    scaleBar(lon = -77.5, lat = 41, distanceLon = 300, distanceLat = 50, distanceLegend = 100, dist.unit = "km", legend.size = 7.5, orientation = FALSE)
dev.off()

