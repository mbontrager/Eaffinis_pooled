## TAKEN FROM
#http://egallic.fr/scale-bar-and-north-arrow-on-a-ggplot2-map/

#
# Result #
#--------#
# Return a list whose elements are :
# 	- rectangle : a data.frame containing the coordinates to draw the first rectangle ;
# 	- rectangle2 : a data.frame containing the coordinates to draw the second rectangle ;
# 	- legend : a data.frame containing the coordinates of the legend texts, and the texts as well.
#
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distanceLon : length of each rectangle ;
# distanceLat : width of each rectangle ;
# distanceLegend : distance between rectangles and legend texts ;
# dist.units : units of distance "km" (kilometers) (default), "nm" (nautical miles), "mi" (statute miles).
createScaleBar <- function(lon,lat,distanceLon,distanceLat,distanceLegend, dist.units = "km"){
    # First rectangle
    bottomRight <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon, dist.units = dist.units, model = "WGS84")
    
    topLeft <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLat, dist.units = dist.units, model = "WGS84")
    rectangle <- cbind(lon=c(lon, lon, bottomRight[1,"long"], bottomRight[1,"long"], lon),
                       lat = c(lat, topLeft[1,"lat"], topLeft[1,"lat"],lat, lat))
    rectangle <- data.frame(rectangle, stringsAsFactors = FALSE)
    
    # Second rectangle t right of the first rectangle
    bottomRight2 <- gcDestination(lon = lon, lat = lat, bearing = 90, dist = distanceLon*2, dist.units = dist.units, model = "WGS84")
    rectangle2 <- cbind(lon = c(bottomRight[1,"long"], bottomRight[1,"long"], bottomRight2[1,"long"], bottomRight2[1,"long"], bottomRight[1,"long"]),
                        lat=c(lat, topLeft[1,"lat"], topLeft[1,"lat"], lat, lat))
    rectangle2 <- data.frame(rectangle2, stringsAsFactors = FALSE)
    
    # Now let's deal with the text
    onTop <- gcDestination(lon = lon, lat = lat, bearing = 0, dist = distanceLegend, dist.units = dist.units, model = "WGS84")
    onTop2 <- onTop3 <- onTop
    onTop2[1,"long"] <- bottomRight[1,"long"]
    onTop3[1,"long"] <- bottomRight2[1,"long"]
    
    legend <- rbind(onTop, onTop2, onTop3)
    legend <- data.frame(cbind(legend, text = c(0, distanceLon, distanceLon*2)), stringsAsFactors = FALSE, row.names = NULL)
    return(list(rectangle = rectangle, rectangle2 = rectangle2, legend = legend))
}


#
# Result #
#--------#
# This function enables to draw a scale bar on a ggplot object, and optionally an orientation arrow #
# Arguments : #
#-------------#
# lon, lat : longitude and latitude of the bottom left point of the first rectangle to draw ;
# distanceLon : length of each rectangle ;
# distanceLat : width of each rectangle ;
# distanceLegend : distance between rectangles and legend texts ;
# dist.units : units of distance "km" (kilometers) (by default), "nm" (nautical miles), "mi" (statute miles) ;
# rec.fill, rec2.fill : filling colour of the rectangles (default to white, and black, resp.);
# rec.colour, rec2.colour : colour of the rectangles (default to black for both);
# legend.colour : legend colour (default to black);
# legend.size : legend size (default to 3);
# orientation : (boolean) if TRUE (default), adds an orientation arrow to the plot ;
# arrow.length : length of the arrow (default to 500 km) ;
# arrow.distance : distance between the scale bar and the bottom of the arrow (default to 300 km) ;
# arrow.North.size : size of the "N" letter (default to 6).
scaleBar <- function(lon, lat, distanceLon, distanceLat, distanceLegend, dist.unit = "km", rec.fill = "white", rec.colour = "black", rec2.fill = "black", rec2.colour = "black", legend.colour = "black", legend.size = 3, orientation = TRUE, arrow.length = 500, arrow.distance = 300, arrow.North.size = 6){
    laScaleBar <- createScaleBar(lon = lon, lat = lat, distanceLon = distanceLon, distanceLat = distanceLat, distanceLegend = distanceLegend, dist.unit = dist.unit)
    # First rectangle
    rectangle1 <- geom_polygon(data = laScaleBar$rectangle, aes(x = lon, y = lat), fill = rec.fill, colour = rec.colour)
    
    # Second rectangle
    rectangle2 <- geom_polygon(data = laScaleBar$rectangle2, aes(x = lon, y = lat), fill = rec2.fill, colour = rec2.colour)
    
    # Legend
    scaleBarLegend <- annotate("text", label = paste(laScaleBar$legend[,"text"], dist.unit, sep=""), x = laScaleBar$legend[,"long"], y = laScaleBar$legend[,"lat"], size = legend.size, colour = legend.colour, fontface="bold")
    
    res <- list(rectangle1, rectangle2, scaleBarLegend)
    
    if(orientation){# Add an arrow pointing North
        coordsArrow <- createOrientationArrow(scaleBar = laScaleBar, length = arrow.length, distance = arrow.distance, dist.unit = dist.unit)
        arrow <- list(geom_segment(data = coordsArrow$res, aes(x = x, y = y, xend = xend, yend = yend)), annotate("text", label = "N", x = coordsArrow$coordsN[1,"x"], y = coordsArrow$coordsN[1,"y"], size = arrow.North.size, colour = "black"))
        res <- c(res, arrow)
    }
    return(res)
}