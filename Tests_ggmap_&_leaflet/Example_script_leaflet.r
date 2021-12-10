library(leaflet)

mydata = read.csv("mydata.turkeyLLv2.csv")

map = leaflet() %>%
	addProviderTiles(providers$Esri.OceanBasemap) %>% # add ocean basemap
	addTiles() %>% # add another layer with place names (default = OpenStreetMap tiles)
	setView(lng = 4.1, lat = 50.8, zoom = 8) # focus map in a certain area and zoom level
map

map = map %>% addCircleMarkers(data=mydata, ~longitude,~latitude, popup=mydata$Beslag.Sanitelnr., color="blue", stroke=F, fillOpacity=0.7, radius="2.5")
map

map = map %>%
	addScaleBar(position = "topright") %>% # add a map scalebar
	addMeasure( # add measurement tool
    primaryLengthUnit="kilometers",
    secondaryLengthUnit="miles", 
    primaryAreaUnit="hectares",
    secondaryAreaUnit="acres", 
    position="topleft") %>%
    addLegend(position="bottomright", colors= c("blue"), label=c("turkey"))
map
