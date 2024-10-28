devtools::install_github("dkahle/ggmap")

library(ggmap)
library(raster)

register_google("AIzaSyBk2NSduWGCf1jZAqxpzFPy39Dpyb-2uEc")
	# https://console.cloud.google.com/google/maps-apis/credentials?project=gis-visualisation

ggmap_rast = function(map, google=T)
	{
  		outputs = list(); map_bbox = attr(map,"bb")
		extent_map = extent(as.numeric(map_bbox[c(2,4,1,3)]))
		my_map = raster(extent_map, nrow=nrow(map), ncol=ncol(map), crs=CRS("+init=epsg:4326"))
		rgb_cols = setNames(as.data.frame(t(col2rgb(map))), c("red","green","blue"))
		red = my_map; values(red) = rgb_cols[["red"]]
		green = my_map; values(green) = rgb_cols[["green"]]
		blue = my_map; values(blue) = rgb_cols[["blue"]]
		vals = my_map; val1 = 1; val2 = length(values(vals))
		values(vals) = seq(val1,val2,1)-1
		cols = rgb(red[], green[], blue[], maxColorValue=255)
		if (google == TRUE) cols = cols[(dim(map)[1]+1):length(cols)] # bug R ??
		cols[(length(cols)-50*dim(map)[1]):length(cols)] = "white" # removes "Google"
		outputs[[1]] = vals; outputs[[2]] = cols; return(outputs)
	}

map1 = ggmap::get_map(c(left=2.5,bottom=50.5,right=4.5,top=51.5), source="stamen", maptype="watercolor")
map2 = ggmap_rast(map1, F); par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(map2[[1]], col=map2[[2]], legend=F, axes=F, box=F) # ggmap(map1)+labs(x ="",y="")

map1 = ggmap::get_map(c(left=2.5,bottom=50.5,right=4.5,top=51.5), source="google", maptype="terrain")
map2 = ggmap_rast(map1, T); par(oma=c(0,0,0,0), mar=c(0,0,0,0))
plot(map2[[1]], col=map2[[2]], legend=F, axes=F, box=F) # ggmap(map1)+labs(x ="",y="")
