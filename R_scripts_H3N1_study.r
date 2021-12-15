# 1. Loading and preparing the wind data from two meteorological stations
# 2. Preparation of the BEAST inputs for the preliminary phylogenetic analysis
# 3. Investigating the temporal signal using the program TempEst
# 4. Running the preliminary phylogenetic analysis in BEAST
# 5. Subsampling of phylogenetic clusters (i.e. farm clusters)
# 6. Projecting sampling coordinates into Belgian Lambert 2008 
# 7. Preparing the starting tree for the phylogeographic run
# 8. Running the continuous phylogeographic analysis in BEAST
# 9. Extracting the spatio-temporal information embedded in posterior trees
# 10. Plotting the dispersal history of H3N1 lineages
# 11. Estimating the dispersal statistics of viral lineages
# 12. Exploring the correlation between dispersal durations and geographic distances
# 13. Exploring the correlation between patristic and geographic distances
# 14. Generating a null dispersal model (randomisations)
# 15. Investigating the impact of wind direction on dispersal direction of viral lineages
# 16. Investigating the correlation between patristic distances and several covariates
# 17. Analysis of the phylogenetic signal associated with several covariates

library(adephylo)
library(diagram)
library(ecodist)
library(lubridate)
library(maptools)
library(phytools)
library(plotrix)
library(seraphim)
library(yhat)

nberOfExtractionFiles = 1000

# 1. Loading and preparing the wind data from two meteorological stations

wind_data = read.csv("Wind_data_two_stations/Winddata_091121.csv", head=T)
meteorological_stations = matrix(nrow=4, ncol=2)
colnames(meteorological_stations) = c("longitude","latitude")
row.names(meteorological_stations) = c("Beitem","Semmerzake","Koksijde","Zeebrugge")
meteorological_stations[1,] = cbind(3.114230, 50.892830)
meteorological_stations[2,] = cbind(3.664843, 50.943126)
meteorological_stations[3,] = cbind(2.650142, 51.105516)
meteorological_stations[4,] = cbind(3.207273, 51.317804)
coordinates = data.frame(meteorological_stations[,c("longitude","latitude")])
coordinates(coordinates) = c("longitude","latitude")
proj4string(coordinates) = CRS("+init=epsg:4326")
municipalities = shapefile("Shapefile_municipalities/Shapefile_post_codes.shp")
meteorological_stations = spTransform(coordinates, crs(municipalities))@coords[1:2,]
wind_data = wind_data[which(tolower(wind_data[,"name"])%in%tolower(row.names(meteorological_stations))),]
wind_data = wind_data[which(wind_data[,"wind_intensity_average"]!="null"),]

# 2. Preparation of the BEAST inputs for the preliminary phylogenetic analysis

sequences = scan("Analysis_1_141220.fasta", what="", sep="\n", quiet=T)
sequence_IDs = gsub(">","",sequences[which(grepl(">",sequences))])
tab1 = read.csv("Analysis_1_141220.csv", head=T, sep=";") # exported from Excel
tab2 = matrix(nrow=length(sequence_IDs), ncol=4)
colnames(tab2) = c("trait","collection_date","latitude","longitude")
for (i in 1:length(sequence_IDs))
	{
		index = which(tab1[,"FASTA"]==sequence_IDs[i])
		if (length(index) != 1)
			{
				print(sequence_IDs[i])
			}	else	{
				tab2[i,"trait"] = sequence_IDs[i]
				date = unlist(strsplit(tab1[index,"sample.date"],"\\/"))
				tab2[i,"collection_date"] = paste(date[3],date[2],date[1],sep="-")
				tab2[i,"latitude"] = tab1[index,"latitude"]
				tab2[i,"longitude"] = tab1[index,"longitude"]			
			}
	}
write.table(tab2, "Analysis_1_141220.txt", quote=F, row.names=F, sep="\t")

# 3. Investigating the temporal signal using the program TempEst

# 4. Running the preliminary phylogenetic analysis in BEAST

	# - substitution model: GTR+G
	# - molecular clock model: relaxed (lognormal distribution)
	# - coalescent model: skygrid (cut-off of 1 years, 50 grid points)

# 5. Subsampling of phylogenetic clusters (i.e. farm clusters)

tree = readAnnotatedNexus("Analysis_1_141220.tree")
tab = read.csv(paste0("Analysis_1_141220.csv"), head=T, sep=";")
tab[,"beslagnummer.AFSCA"] = gsub("HOK E","",tab[,"beslagnummer.AFSCA"])
tab[,"beslagnummer.AFSCA"] = gsub("\xca","",tab[,"beslagnummer.AFSCA"])
tab[,"beslagnummer.AFSCA"] = gsub(" ","",tab[,"beslagnummer.AFSCA"])
subTrees = ape::subtrees(tree)
c1 = 0; clusters1 = list(); bootstraps1 = list()
c2 = 0; clusters2 = list(); bootstraps2 = list() 
for (i in 2:length(subTrees)) # the 1st subtree is the entire one
	{
		subTree = subTrees[i][[1]]
		labels = subTree$tip.label
		farms = c()
		for (j in 1:length(labels))
			{
				index = which(tab[,"FASTA"]==labels[j])
				farms = c(farms, tab[index,"beslagnummer.AFSCA"])
			}
		if (length(unique(farms)) == 1)
			{
				c1 = c1+1; clusters1[[c1]] = labels
			}
	}
for (i in 1:length(clusters1))
	{
		nested = FALSE
		for (j in 1:length(clusters1))
			{
				if (i != j)
					{
						allSequencesIncluded = TRUE
						for (k in 1:length(clusters1[[i]]))
							{
								if (!clusters1[[i]][k]%in%clusters1[[j]]) allSequencesIncluded = FALSE
							}
						if (allSequencesIncluded == TRUE) nested = TRUE
					}
			}
		if (nested == FALSE)
			{
				c2 = c2 + 1; clusters2[[c2]] = clusters1[[i]]
			}	
	}
sequencesToRemove2 = c()
for (i in 1:length(clusters2))
	{
		sequencesToRemove2 = c(sequencesToRemove2, sample(clusters2[[i]],length(clusters2[[i]])-1,replace=F))
	}
sequences = scan("Analysis_1_141220.fasta", what="", sep="\n", quiet=T)
sequence_IDs = gsub(">","",sequences[which(grepl(">",sequences))])
sequences = sequences[which(!grepl(">",sequences))]; fasta = c()
sequences = sequences[which(!sequence_IDs%in%sequencesToRemove2)]
sequence_IDs = sequence_IDs[which(!sequence_IDs%in%sequencesToRemove2)]
for (i in 1:length(sequence_IDs)) fasta = c(fasta, paste0(">",sequence_IDs[i]), sequences[i])
write(fasta, "Analysis_2_141220.fasta")

# 6. Projecting sampling coordinates into Belgian Lambert 2008 

sequences = scan("Analysis_2_141220.fasta", what="", sep="\n", quiet=T)
sequence_IDs = gsub(">","",sequences[which(grepl(">",sequences))])
tab0 = read.csv("Analysis_1_141220.csv", head=T, sep=";") # exported from Excel
tab1 = read.table("Analysis_1_141220.txt", head=T, sep="\t")
tab1 = tab1[which(tab1[,"trait"]%in%sequence_IDs),]; buffer = tab1
for (i in 1:dim(buffer)[1]) buffer[i,] = tab1[which(tab1[,"trait"]==sequence_IDs[i]),]
tab2 = matrix(nrow=length(sequence_IDs), ncol=4)
colnames(tab2) = c("trait","collection_date","latitude","longitude")
municipalities = shapefile("Shapefile_municipalities/Shapefile_post_codes.shp")
coordinates = tab1[,c("longitude","latitude")]
coordinates(coordinates) = c("longitude","latitude")
proj4string(coordinates) = CRS("+init=epsg:4326")
coordinates = spTransform(coordinates, crs(municipalities))
resamplingIdenticalCoordinatesWithinMunicipalities = FALSE
for (i in 1:length(sequence_IDs))
	{
		index = which(tab1[,"trait"]==sequence_IDs[i])
		tab2[i,"trait"] = tab1[index,"trait"]; tab2[i,"collection_date"] = tab1[index,"collection_date"]
		if (resamplingIdenticalCoordinatesWithinMunicipalities == FALSE)
			{
				tab2[i,"latitude"] = coordinates@coords[index,2]; tab2[i,"longitude"] = coordinates@coords[index,1]
			}	else		{
				indices = which((tab1[,"longitude"]==tab1[i,"longitude"])&(tab1[,"latitude"]==tab1[i,"latitude"]))
				if (length(indices) > 1)	
					{
						original_longitude = coordinates@coords[index,"longitude"]
						original_latitude = coordinates@coords[index,"latitude"]
						zipCode = tab0[which(tab0[,"FASTA"]==sequence_IDs[i]),"Zip.code"]
						index = which(municipalities@data[,"nouveau_PO"]==zipCode)
						if (length(index) == 1)
							{
								pols = municipalities@polygons[[index]]
								areas = rep(NA, length(pols@Polygons))
								for (j in 1:length(areas)) areas[j] = pols@Polygons[[j]]@area
								points = pols@Polygons[[which(areas==max(areas))]]@coords
								p = Polygon(points); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								if (point.in.polygon(original_longitude,original_latitude,points[,1],points[,2]) != 1)
									{
										for (j in 1:length(municipalities@polygons))
											{
												for (k in 1:length(municipalities@polygons[[j]]@Polygons))
													{
														points = municipalities@polygons[[j]]@Polygons[[k]]@coords
														if (point.in.polygon(original_longitude,original_latitude,points[,1],points[,2]) == 1)
															{
																p = Polygon(points); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
															}
													}
											}
									}
								new_point = spsample(sps, 1, type="random")
								tab2[i,"latitude"] = new_point@coords[,2]
								tab2[i,"longitude"] = new_point@coords[,1]			
							}
					}
			}
	}
write.table(tab2, "Analysis_2_141220.txt", quote=F, row.names=F, sep="\t")

# 7. Preparing the starting tree for the phylogeographic run

tree1 = read.nexus("Analysis_1_141220.tree")
sequences = scan("Analysis_2_141220.fasta", what="", sep="\n", quiet=T)
sequence_IDs = gsub(">","",sequences[which(grepl(">",sequences))])
tree2 = drop.tip(tree1, tree1$tip.label[which(!tree1$tip.label%in%sequence_IDs)])
write.tree(tree2, "Starting_tree_RRW.tree")

# 8. Running the continuous phylogeographic analysis in BEAST

	# - substitution model: GTR+G
	# - molecular clock model: relaxed (lognormal distribution)
	# - coalescent model: skygrid (cut-off of 1 years, 50 grid points)
	# - relaxed random walk model: gamma (jitter value: 2000, in Lambert 72)

# 9. Extracting the spatio-temporal information embedded in posterior trees

allTrees = scan(file=paste0("Analysis_2_j_2000.trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
index1 = which(allTrees=="\t\t;")[2]
index2 = which(grepl("tree STATE_10100000 ",allTrees))
index3 = which(grepl("tree STATE_130000000 ",allTrees))
interval = (index3-index2)/nberOfExtractionFiles
selectedTrees = allTrees[c(1:index1,seq(index2+interval,index3,interval))]
write(c(selectedTrees,"End;"), "Analysis_2_selected.trees")

localTreesDirectory = "Analysis_2_extractions"
allTrees = scan(file="Analysis_2_selected.trees", what="", sep="\n", quiet=T)
burnIn = 0; randomSampling = FALSE
nberOfTreesToSample = nberOfExtractionFiles
metadata = read.table("Analysis_2_141220.txt", head=T)
mostRecentSamplingDatum = max(decimal_date(ymd(metadata[,"collection_date"])))
coordinateAttributeName = "location"; nberOfCores = 5
treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)

# 10. Plotting the dispersal history of H3N1 lineages

croppingPolygons = FALSE
colour_scale1 = colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121]
colour_scale1 = c(rep(colour_scale1[1],100), colour_scale1)

source("mccTreeExtractions.r")
tree = readAnnotatedNexus("Analysis_2_selected.tree")
metadata = read.table("Analysis_2_141220.txt", head=T)
mostRecentSamplingDatum = max(decimal_date(ymd(metadata[,"collection_date"])))
mcc = mccTreeExtractions(tree, mostRecentSamplingDatum)

localTreesDirectory = "Analysis_2_extractions"
prob = 0.80; precision = 1/12
rootHeights = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		csv = read.csv(paste(localTreesDirectory,"/TreeExtractions_",i,".csv",sep=""), header=T)
		rootHeights[i] = min(csv[,"startYear"])
	}
startDatum = quantile(rootHeights, probs=c(0.025,0.975))[1]
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))

minYear = startDatum; maxYear = max(mcc[,"endYear"])
mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
tree_node_colours = colour_scale1[((((nodeHeights(tree)[,2]+min(mcc[,"startYear"]))-minYear)/(maxYear-minYear))*200)+1]
endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*200)+1
endYears_colours = colour_scale1[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*200)+1)
		polygons_colours[i] = paste0(colour_scale1[polygon_index],"70")
	}

municipalities = shapefile("Shapefile_municipalities/Shapefile_post_codes.shp")
belgium = spTransform(raster::getData("GADM", country="BEL", level=0), crs(municipalities))
e = extent(530000, 610000, 650000, 730000); municipalities_cropped = crop(municipalities, e); belgium_cropped = crop(belgium, e)

farms1 = read.csv("All_farms_01022020.csv", head=T, sep=";")[,15:16]
colnames(farms1) = c("longitude","latitude"); farms1 = farms1[which(!is.na(farms1[,"longitude"])),]
farms2 = data.frame(lon=as.numeric(farms1[,"longitude"]), lat=as.numeric(farms1[,"latitude"]))
lambert72 = CRS("+init=epsg:31370"); coordinates(farms2) = c("lon", "lat")
proj4string(farms2) = lambert72; farms2 = spTransform(farms2, crs(municipalities))
farms = as.matrix(cbind(farms2@coords)); colnames(farms) = c("longitude","latitude")
pts = point.in.polygon(farms[,"longitude"],farms[,"latitude"],c(e@xmin,e@xmax,e@xmax,e@xmin,e@xmin),c(e@ymin,e@ymin,e@ymax,e@ymax,e@ymin))
farms_cropped = farms[which(pts==1),]

farms1 = read.csv("All_farms_01022020.csv", head=T, sep=";")[,c(2,15:16)]
colnames(farms1) = c("farm","longitude","latitude"); farms1 = farms1[which(!is.na(farms1[,"longitude"])),]
farms2 = data.frame(lon=as.numeric(farms1[,"longitude"]), lat=as.numeric(farms1[,"latitude"]))
lambert72 = CRS("+init=epsg:31370"); coordinates(farms2) = c("lon", "lat")
proj4string(farms2) = lambert72; farms2 = spTransform(farms2, crs(municipalities))
outbreaks1 = read.csv("Outrbreaks_030221.csv", head=T, sep=";")[,c(2,4)]
outbreaks2 = cbind(outbreaks1, matrix(nrow=dim(outbreaks1)[1], ncol=2))
colnames(outbreaks2) = c("farm","onset_of_symptoms","longitude","latitude")
for (i in 1:dim(outbreaks2)[1])
	{
		indices = which(farms1[,1]==outbreaks2[i,"farm"])
		outbreaks2[i,c("longitude","latitude")] = farms2@coords[indices[1],1:2]
	}
outbreaks2 = outbreaks2[which(!is.na(outbreaks2[,"longitude"])),]
outbreaks2[,"onset_of_symptoms"] = gsub(" ","",outbreaks2[,"onset_of_symptoms"])
outbreaks2[,"onset_of_symptoms"] = decimal_date(dmy(gsub("\\/","-",outbreaks2[,"onset_of_symptoms"])))
outbreaks2 = outbreaks2[which(!is.na(outbreaks2[,"onset_of_symptoms"])),]
outbreaks2 = outbreaks2[order(outbreaks2[,"onset_of_symptoms"]),]
outbreak_indices = (((outbreaks2[,"onset_of_symptoms"]-minYear)/(maxYear-minYear))*200)+1
outbreak_colours = colour_scale1[outbreak_indices]
pts = point.in.polygon(outbreaks2[,"longitude"],outbreaks2[,"latitude"],c(e@xmin,e@xmax,e@xmax,e@xmin,e@xmin),c(e@ymin,e@ymin,e@ymax,e@ymax,e@ymin))
outbreaks_cropped = outbreaks2[which(pts==1),]
outbreak_indices = (((outbreaks_cropped[,"onset_of_symptoms"]-minYear)/(maxYear-minYear))*200)+1
outbreak_colours_cropped = colour_scale1[outbreak_indices]

ats = decimal_date(dmy(c("01-01-2019","01-02-2019","01-03-2019","01-04-2019","01-05-2019","01-06-2019","01-07-2019")))
labels = c("01-01-2019","01-02-2019","01-03-2019","01-04-2019","01-05-2019","01-06-2019","01-07-2019")

pdf("Analysis_2_spread1_NEW.pdf", width=11, height=5.6)
par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(tree, show.tip.label=F, show.node.label=F, edge.width=0.2, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
for (i in 1:dim(tree$edge)[1])
	{
		if (i == 1)
			{
				nodelabels(node=tree$edge[i,1], pch=16, cex=0.7, col=colour_scale1[1])
				nodelabels(node=tree$edge[i,1], pch=1, cex=0.7, col="gray30", lwd=0.25)
			}
		if (!tree$edge[i,2]%in%tree$edge[,1])
			{
				nodelabels(node=tree$edge[i,2], pch=15, cex=0.7, col=tree_node_colours[i])
				nodelabels(node=tree$edge[i,2], pch=0, cex=0.7, col="gray30", lwd=0.2)
			}	else		{
				nodelabels(node=tree$edge[i,2], pch=16, cex=0.7, col=tree_node_colours[i])
				nodelabels(node=tree$edge[i,2], pch=1, cex=0.7, col="gray30", lwd=0.2)				
			}
	}
plot(municipalities_cropped, border=NA, col="gray90", axes=F, frame=F) 
plot(municipalities_cropped, border="white", col=NA, lwd=0.2, add=T)
for (i in 1:length(polygons))
	{
		for (j in 1:length(polygons[[i]]@polygons))
			{
				polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
			}
		pol = polygons[[i]]; crs(pol) = crs(municipalities)
		if (croppingPolygons == TRUE) pol = crop(pol, municipalities)
		plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
	}
for (i in 1:dim(mcc)[1])
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		if (!mcc[i,"node2"]%in%mcc[,"node1"])
			{
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.65)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.2, cex=0.65)
			}	else		{
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)				
			}
		if (i == 1)
			{
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale1[1], cex=0.7)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
			}
	}
rect(530000, 650000, 610000, 730000, lwd=0.2, col=NA, border="gray30")
dev.off()

pdf("Analysis_2_spread2_NEW.pdf", width=11, height=4.8)
par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(municipalities, border=NA, col="gray90", axes=F, frame=F) 
plot(municipalities, border="white", col=NA, lwd=0.2, add=T)
plot(belgium, border="gray50", col=NA, lwd=0.5, add=T)
points(farms[,c("longitude","latitude")], cex=0.3, pch=16, col="gray70")
for (i in dim(outbreaks2)[1]:1)
	{
		points(outbreaks2[i,c("longitude","latitude")], cex=0.65, pch=17, col=outbreak_colours[i])
		points(outbreaks2[i,c("longitude","latitude")], cex=0.65, pch=2, col="gray30", lwd=0.2)
	}
plot(municipalities, border=NA, col="gray90", axes=F, frame=F) 
plot(municipalities, border="white", col=NA, lwd=0.2, add=T)
plot(belgium, border="gray50", col=NA, lwd=0.5, add=T)
for (i in 1:length(polygons))
	{
		for (j in 1:length(polygons[[i]]@polygons))
			{
				polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
			}
		pol = polygons[[i]]; crs(pol) = crs(municipalities)
		if (croppingPolygons == TRUE) pol = crop(pol, municipalities)
		plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
	}
for (i in 1:dim(mcc)[1])
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		if (!mcc[i,"node2"]%in%mcc[,"node1"])
			{
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.65)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.2, cex=0.65)
			}	else		{
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)				
			}
		if (i == 1)
			{
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale1[1], cex=0.7)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
			}
	}
rect(530000, 650000, 610000, 730000, lwd=0.2, col=NA, border="gray30")
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
index1 = round((((rast[1]-minYear)/(maxYear-minYear))*200)+1); index2 = round((((rast[2]-minYear)/(maxYear-minYear))*200)+1)
plot(rast, legend.only=T, add=T, col=colour_scale1[index1:index2], legend.width=0.5, legend.shrink=0.3, smallplot=c(0.050,0.550,0.100,0.110),
	 legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
	 axis.args=list(cex.axis=0.7, lwd=0, lwd.tick=0.2, tck=-1.0, col.axis="gray30", line=0, mgp=c(0,0.06,0), at=ats, label=labels))
dev.off()

pdf("Analysis_2_spread3_NEW.pdf", width=11, height=5.6)
par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(municipalities_cropped, border=NA, col="gray90", axes=F, frame=F) 
plot(municipalities_cropped, border="white", col=NA, lwd=0.2, add=T)
points(farms_cropped[,c("longitude","latitude")], cex=0.3, pch=16, col="gray70")
for (i in dim(outbreaks2)[1]:1)
	{
		points(outbreaks_cropped[i,c("longitude","latitude")], cex=0.80, pch=17, col=outbreak_colours_cropped[i])
		points(outbreaks_cropped[i,c("longitude","latitude")], cex=0.80, pch=2, col="gray30", lwd=0.2)
	}
rect(530000, 650000, 610000, 730000, lwd=0.2, col=NA, border="gray30")
plot(municipalities_cropped, border=NA, col="gray90", axes=F, frame=F) 
plot(municipalities_cropped, border="white", col=NA, lwd=0.2, add=T)
for (i in 1:length(polygons))
	{
		for (j in 1:length(polygons[[i]]@polygons))
			{
				polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
			}
		pol = polygons[[i]]; crs(pol) = crs(municipalities)
		if (croppingPolygons == TRUE) pol = crop(pol, municipalities)
		plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
	}
for (i in 1:dim(mcc)[1])
	{
		curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
	}
for (i in dim(mcc)[1]:1)
	{
		if (!mcc[i,"node2"]%in%mcc[,"node1"])
			{
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.65)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.2, cex=0.65)
			}	else		{
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)				
			}
		if (i == 1)
			{
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale1[1], cex=0.7)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
			}
	}
points(meteorological_stations, cex=1, pch=16, col="gray30")
rect(530000, 650000, 610000, 730000, lwd=0.2, col=NA, border="gray30")
dev.off()

cutOffs = c(decimal_date(dmy(c("05-04-2019","26-04-2019","16-05-2019"))), mostRecentSamplingDatum+1)

pdf("Analysis_2_spread4_NEW.pdf", width=11, height=3)
par(mfrow=c(1,4), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (h in 1:length(cutOffs))
	{
		plot(municipalities_cropped, border=NA, col="gray90", axes=F) 
		plot(municipalities_cropped, border="white", col=NA, lwd=0.2, add=T)
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) < cutOffs[h])
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(municipalities)
						if (croppingPolygons == TRUE) pol = crop(pol, municipalities)
						plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] < cutOffs[h])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
								    arr.width=0, lwd=0.2, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
			}
		for (i in dim(mcc)[1]:1)
			{
				if (mcc[i,"endYear"] < cutOffs[h])
					{
						if (!mcc[i,"node2"]%in%mcc[,"node1"])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.65)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.2, cex=0.65)
							}	else		{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.7)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)				
							}
						if (i == 1)
							{
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale1, cex=0.7)
								points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.2, cex=0.7)
							}
					}
			}
		rect(530000, 650000, 610000, 730000, lwd=0.2, col=NA, border="gray30")
	}
dev.off()

pdf("Analysis_2_spread5_NEW.pdf", width=11, height=3)
par(mfrow=c(1,4), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
dates = decimal_date(ymd(wind_data[,"date"]))
for (h in 1:length(cutOffs))
	{
		if (h == 1)
			{
				startingDate = min(mcc[,"startYear"])
			}	else	{
				startingDate = cutOffs[h-1]
			}
		sub = wind_data[which((dates>startingDate)&(dates<cutOffs[h])),]
		intensities = rep(NA, length(unique(sub[,"date"]))); directions = rep(NA, length(unique(sub[,"date"])))
		for (i in 1:length(unique(sub[,"date"])))
			{
				intensities[i] = mean(as.numeric(sub[which(sub[,"date"]==unique(sub[,"date"])[i]),"wind_intensity_average"]))
				directions[i] = mean(as.numeric(sub[which(sub[,"date"]==unique(sub[,"date"])[i]),"wind_direction_average"]))
			}
		directions = directions-180; directions[directions[]<0] = directions[directions[]<0]+360
		directions = c(0, directions); intensities = c(10, intensities)
		polar.plot(intensities, directions, start=90, clockwise=T, labels="")
	}
dev.off()

# 11. Estimating the dispersal statistics of viral lineages

timeSlices = 100; onlyTipBranches = FALSE; showingPlots = FALSE
nberOfCores = 5; slidingWindow = 1/(365/14)
localTreesDirectory = "Analysis_2_extractions"
dir.create(file.path("Analysis_2_disp_stats"), showWarnings=F)
outputName = paste0("Analysis_2_disp_stats/YFV3")

spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow) 
	# mean branch dispersal velocity of 4.69 km/day (95% CI = [2.61-20.39])
	# weighted branch dispersal velocity of 0.74 km/day (95% CI = [0.58-0.91])

# 12. Exploring the correlation between dispersal durations and geographic distances

corrs = rep(NA, nberOfExtractionFiles); lrR2s = rep(NA, nberOfExtractionFiles)
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0("Analysis_2_extractions/TreeExtractions_",i,".csv")); tS = tab[,"length"]
		dS = sqrt(((tab[,"endLon"]-tab[,"startLon"])^2)+((tab[,"endLat"]-tab[,"startLat"])^2))
		corrs[i] = cor(tS, dS, method="spearman"); lrR2s[i] = summary(lm("tS ~ dS"))$r.squared
	}
cat(paste0("correlation (Spearman): ",round(median(corrs),3),", 95% HPD [",round(quantile(corrs,0.025),3),"-",round(quantile(corrs,0.975),3),"]\n"))
cat(paste0("R2 from linear regression: ",round(median(lrR2s),3),", 95% HPD [",round(quantile(lrR2s,0.025),3),"-",round(quantile(lrR2s,0.975),3),"]\n"))
	# correlation (Spearman): 0.375, 95% HPD [0.253-0.493], R2 from linear regression: 0.047, 95% HPD [0.008-0.100]

# 13. Exploring the correlation between patristic and geographic distances

geographic_distances_list = list(); patristic_distances_list = list(); c = 0
tab = read.csv(paste0("Analysis_2_extractions/TreeExtractions_1.csv"))
farms = round(tab[which(!tab[,"node2"]%in%tab[,"node1"]),c("endLon","endLat")])
for (i in 2:dim(farms)[1])
	{
		for (j in 1:(i-1))
			{
				geographic_distances = c(); patristic_distances = c()
				for (k in 1:nberOfExtractionFiles)
					{
						tab = read.csv(paste0("Analysis_2_extractions/TreeExtractions_",k,".csv"), header=T)
						index1 = which((round(tab[,"endLon"])==farms[i,"endLon"])&(round(tab[,"endLat"])==farms[i,"endLat"]))
						index2 = which((round(tab[,"endLon"])==farms[j,"endLon"])&(round(tab[,"endLat"])==farms[j,"endLat"]))
						if ((length(index1)>0)&(length(index2)>0))
							{
								geographic_dis = sqrt(((tab[index1,"endLon"]-tab[index2,"endLon"])^2)+((tab[index1,"endLat"]-tab[index2,"endLat"])^2))
								indices1 = index1; root = FALSE
								while (root == FALSE)
									{	
										if (tab[indices1[length(indices1)],"node1"]%in%tab[,"node2"])
											{
												indices1 = c(indices1, which(tab[,"node2"]==tab[indices1[length(indices1)],"node1"]))
											}	else		{
												root = TRUE
											}
									}
								indices2 = index2; root = FALSE
								while (root == FALSE)
									{	
										if (tab[indices2[length(indices2)],"node1"]%in%tab[,"node2"])
											{
												indices2 = c(indices2, which(tab[,"node2"]==tab[indices2[length(indices2)],"node1"]))
											}	else		{
												root = TRUE
											}
									}
								indices3 = indices1[which(indices1%in%indices2)]
								if (length(indices3) == 0)
									{
										patristic_dis = sum(tab[c(indices1,indices2),"length"])
									}	else		{
										patristic_dis = sum(tab[c(indices1[which(!indices1%in%indices3)],indices2[which(!indices2%in%indices3)]),"length"])
									}
								geographic_distances = c(geographic_distances, geographic_dis)
								patristic_distances = c(patristic_distances, patristic_dis)
							}	else		{
								print(c(i,j,k))
							}
					}
				c = c+1
				geographic_distances_list[[c]] = geographic_distances
				patristic_distances_list[[c]] = patristic_distances
			}
	}
median_geographic_distances = c(); median_patristic_distances = c()
for (i in 1:length(geographic_distances_list))
	{
		median_patristic_distances = c(median_patristic_distances, median(patristic_distances_list[[i]]))
		median_geographic_distances = c(median_geographic_distances, median(geographic_distances_list[[i]]))
	}
cat(paste0("correlation (Spearman): ",round(cor(median_patristic_distances,median_geographic_distances),3)))
cat(paste0("R2 from linear regression: ",round(summary(lm("median_patristic_distances ~ median_geographic_distances"))$r.squared,3)))
	# correlation (Spearman): 0.148, R2 from linear regression: 0.022

# 14. Generating a null dispersal model (randomisations)

localTreesDirectory = "Analysis_2_extractions"
municipalities = shapefile("Shapefile_municipalities/Shapefile_post_codes.shp")
rast = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(municipalities))
rast[!is.na(rast[])] = 0; envVariables = list(rast); randomProcedure = 3; nberOfCores = 1
treesRandomisation(localTreesDirectory, nberOfExtractionFiles, envVariables, randomProcedure, nberOfCores)

# 15. Investigating the impact of wind direction on dispersal direction of viral lineages

angle1 = function(x1, y1, x2, y2)
	{
		a = abs(x2-x1)
		h = sqrt(((x2-x1)^2)+((y2-y1)^2))
		theta = acos(a/h)*(180/pi)
		return(theta)
	}
angle2 = function(x1, y1, x2, y2)
	{
		if (x1 < x2)
			{
				if (y1 < y2) a = 90-angle1(x1,y1,x2,y2)
				if (y1 > y2) a = 90+angle1(x1,y1,x2,y2)
			}
		if (x1 > x2)
			{
				if (y1 > y2) a = 270-angle1(x1,y1,x2,y2)
				if (y1 < y2) a = 270+angle1(x1,y1,x2,y2)
			}
		if (x1 != x2) return(a)
		if (x1 == x2) return(NA)
	}
minStartYears = c()
for (i in 1:nberOfExtractionFiles)
	{
		tab = read.csv(paste0("Analysis_2_extractions/TreeExtractions_",i,".csv"), header=T)
		minStartYears = c(minStartYears, min(tab[,"startYear"]))
	}
metadata = read.table("Analysis_2_141220.txt", head=T)
mostRecentSamplingDatum = max(decimal_date(ymd(metadata[,"collection_date"])))
distance_threshold = 1000; cutOffs = c(mostRecentSamplingDatum+1) # to get the BF for the entire period
cutOffs = c(decimal_date(dmy(c("26-04-2019","16-05-2019"))), mostRecentSamplingDatum+1)
for (h in 1:length(cutOffs))
	{
		if (h == 1)
			{
				startYear = 0
			}	else {
				startYear = cutOffs[h-1]
			}
		differencesM_obs = matrix(nrow=nberOfExtractionFiles, ncol=1)
		differencesM_ran = matrix(nrow=nberOfExtractionFiles, ncol=1)
		dates = decimal_date(ymd(wind_data[,"date"]))
		for (i in 1:nberOfExtractionFiles)
			{
				tab = read.csv(paste0("Analysis_2_extractions/TreeExtractions_",i,".csv"), header=T)
				tab = tab[which(((tab[,"startYear"]+tab[,"endYear"])/2)>=min(dates)),] # to be improved?
				tab = tab[which((tab[,"startYear"]>startYear)&(tab[,"endYear"]<=cutOffs[h])),]
				values1 = matrix(nrow=dim(tab)[1], ncol=2); colnames(values1) = c("angle","angle_intensity")
				values2 = matrix(nrow=dim(tab)[1], ncol=2); colnames(values2) = c("angle","angle_intensity")
				for (j in 1:dim(tab)[1])
					{
						geoDistance = sqrt(((tab[j,"endLon"]-tab[j,"startLon"])^2)+((tab[j,"endLat"]-tab[j,"startLat"])^2))
						if (geoDistance < distance_threshold)
						# if (geoDistance > distance_threshold)
							{
								values1[j,1] = angle2(tab[j,"startLon"],tab[j,"startLat"],tab[j,"endLon"],tab[j,"endLat"])
								values1[j,2] = values1[j,1]*(geoDistance/tab[j,"length"])
								indices = which((dates>tab[j,"startYear"])&(dates<=tab[j,"endYear"]))
								if (length(indices) == 0)
									{
										indices = which((abs(dates-((tab[j,"startYear"]+tab[j,"endYear"])/2)))==min(abs(dates-((tab[j,"startYear"]+tab[j,"endYear"])/2))))
									}
								values2[j,1] = mean(as.numeric(wind_data[indices,"wind_direction_average"]))
								values2[j,1] = values2[j,1]-180
								if (values2[j,1] < 0) values2[j,1] = values2[j,1]+360
								values2[j,2] = values2[j,1]*mean(as.numeric(wind_data[indices,"wind_intensity_average"]))
							}
					}
				differences = abs(values2[,1]-values1[,1])
				differences[(differences[]>180)&(!is.na(differences))] = 360-differences[(differences[]>180)&(!is.na(differences))]
				differencesM_obs[i,1] = median(differences, na.rm=T)
				tab = read.csv(paste0("Analysis_2_extractions/TreeRandomisation_",i,".csv"), header=T)
				tab = tab[which(((tab[,"startYear"]+tab[,"endYear"])/2)>=min(dates)),] # to be improved?
				tab = tab[which((tab[,"startYear"]>startYear)&(tab[,"endYear"]<=cutOffs[h])),]
				values1 = matrix(nrow=dim(tab)[1], ncol=2); colnames(values1) = c("angle","angle_intensity")
				values2 = matrix(nrow=dim(tab)[1], ncol=2); colnames(values2) = c("angle","angle_intensity")
				for (j in 1:dim(tab)[1])
					{
						geoDistance = sqrt(((tab[j,"endLon"]-tab[j,"startLon"])^2)+((tab[j,"endLat"]-tab[j,"startLat"])^2))
						if (geoDistance > distance_threshold)
							{
								values1[j,1] = angle2(tab[j,"startLon"],tab[j,"startLat"],tab[j,"endLon"],tab[j,"endLat"])
								values1[j,2] = values1[j,1]*(geoDistance/tab[j,"length"])
								indices = which((dates>tab[j,"startYear"])&(dates<=tab[j,"endYear"]))
								if (length(indices) == 0)
									{
										indices = which((abs(dates-((tab[j,"startYear"]+tab[j,"endYear"])/2)))==min(abs(dates-((tab[j,"startYear"]+tab[j,"endYear"])/2))))
									}
								values2[j,1] = mean(as.numeric(wind_data[indices,"wind_direction_average"]))
								values2[j,1] = values2[j,1]-180
								if (values2[j,1] < 0) values2[j,1] = values2[j,1]+360
								values2[j,2] = values2[j,1]*mean(as.numeric(wind_data[indices,"wind_intensity_average"]))
							}
					}
				differences = abs(values2[,1]-values1[,1])
				differences[(differences[]>180)&(!is.na(differences))] = 360-differences[(differences[]>180)&(!is.na(differences))]
				differencesM_ran[i,1] = median(differences, na.rm=T)
			}
		p = sum((differencesM_ran-differencesM_obs)>0, na.rm=T)/length(differencesM_obs)
		BF = (p/(1-p))/(0.5/(1-0.5)); cat("Period ",h,"  -  BF = ",round(BF,2),"\n",sep="")
	}
		# BF with cut-off distance < 1000:  1.36 (and for the 3 time periods: 0.86, 1.00, 0.94)
		# BF with cut-off distance > 1000:  2.47 (and for the 3 time periods: 2.09, 1.04, 1.12)
		# BF with cut-off distance > 2000:  2.76 (and for the 3 time periods: 2.62, 0.84, 1.32)
		# BF with cut-off distance > 5000:  3.08 (and for the 3 time periods: 3.33, 0.67, 1.04)
		# BF with cut-off distance > 10000: 3.76 (and for the 3 time periods: 4.05, 0.42, 0.97)

# 16. Investigating the correlation between patristic distances and several covariates

trees = read.nexus("Analysis_2_selected.trees"); tree1 = trees[1] # outbreak "0" was in the same farm as outbreak "1"
labels = as.character(tree1$tip.label$tip.label); labels[which(labels=="0-1-PTL")] = "1-1-PTL"
metadata = read.table("Analysis_2_141220.txt", head=T); metadata[which(metadata[,"trait"]=="0-1-PTL"),"trait"] = "1-1-PTL"
dists = matrix(nrow=length(labels), ncol=length(labels)); row.names(dists) = labels; colnames(dists) = labels
for (i in 2:dim(dists)[1])
	{
		for (j in 1:(i-1))
			{
				index1 = which(metadata[,"trait"]==row.names(dists)[j]); index2 = which(metadata[,"trait"]==row.names(dists)[i])
				dists[i,j] = 1/(1+sqrt(((metadata[index1,"longitude"]-metadata[index2,"longitude"])^2)+((metadata[index1,"latitude"]-metadata[index2,"latitude"])^2)))
				dists[j,i] = 1/(1+sqrt(((metadata[index1,"longitude"]-metadata[index2,"longitude"])^2)+((metadata[index1,"latitude"]-metadata[index2,"latitude"])^2)))
			}
	}
covariates = list(); covariate_names = c()
covariates[[1]] = dists; covariate_names = c(covariate_names, "proximity")
matrice_names = c("families","feed","hatcheries","rendac","veterinaries")
matrice_names = c("human_mvt","transports","SatSQan")
for (i in 1:length(matrice_names))
	{
		dists = matrix(nrow=length(labels), ncol=length(labels))
		row.names(dists) = labels; colnames(dists) = labels
		# mat = read.csv(paste0("Matrices_of_Sciensano/All_previous_versions/Matrix_",matrice_names[i],".csv"), head=T, sep=";")
		mat = read.csv(paste0("Matrices_of_Sciensano/Matrix_",matrice_names[i],".csv"), head=T, sep=";")
		colnames(mat) = gsub("X","",colnames(mat)); farms = gsub("a","",gsub("b","",row.names(mat)))
		for (j in 2:dim(dists)[1])
			{
				for (k in 1:(j-1))
					{
						outbreak1 = unlist(strsplit(row.names(dists)[k],"-"))[1]; index1 = which(farms==outbreak1)
						outbreak2 = unlist(strsplit(row.names(dists)[j],"-"))[1]; index2 = which(farms==outbreak2)
						if ((length(index1)!=0)&(length(index2)!=0))
							{
								dists[j,k] = mean(c(mat[index1[1],index2[1]],mat[index2[1],index1[1]]), na.rm=T)
								dists[k,j] = mean(c(mat[index1[1],index2[1]],mat[index2[1],index1[1]]), na.rm=T)
							}
					}
			}		
		covariates[[length(covariates)+1]] = dists; covariate_names = c(covariate_names, matrice_names[i])
	}
		# - proximity (continuous variable): 1/(1+geographic_distance)
		# - families (symmetric matrix, binary variable): "1" when the two farms are owned by members of the same family, "0" otherwise
		# - hatcheries (symmetric matrix, binary variable): "1" when the two farms belong to the same hatchery network (the birds from the two farms have the same origin), "0" otherwise
		# - rendac (symmetric matrix): number of times both farms were visited by the same Rendac truck
		# - veterinaries (asymmetric matrix*): number of times a farm (row) was visited by a veterinary who also visited another farm (column)
		# - feed (asymmetric matrix**): number of times feed was delivered to a farm (row) by the same feed company delivering another farm (column)
		
		# (*) asymetric because a common vet for two farms A and B can have visited farm A 3 times, while he visited farm B only 1 time. So the value (a,b) and (b,a) can differ
		# (**) same remarks as for the veterinaries: absolute number of feed deliveries counted in a similar way as the visits by the veterinaries
		
		# - human movements (symmetric matrix, binary variable): human movement ("1") or not ("0") between farms 
		# - transports (symmetric matrix, binary variable): transports (feed companies, manure and/or Rendac; "1") or not ("0") between farms
		# - SatSQan (symmetric matrix, binary variable): no connection with the farm in any cluster ("0"), connection with this farm in the same cluster ("1")

	# 17.1. Univariate analyses (Mantel tests)
	
mantel_tests_r = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mantel_tests_r) = covariate_names
mantel_tests_p = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mantel_tests_p) = covariate_names
for (i in 1:length(trees))
	{
		print(i)
		dists = matrix(0, nrow=length(labels), ncol=length(labels))
		row.names(dists) = labels; colnames(dists) = labels
		mat = as.matrix(adephylo::distTips(as(trees[i][[1]],"phylo4"), tips=c(labels[j],labels[k]), method="patristic"))
		row.names(mat)[row.names(mat)=="0-1-PTL"] = "1-1-PTL"; colnames(mat)[colnames(mat)=="0-1-PTL"] = "1-1-PTL"
		for (j in 2:dim(dists)[1])
			{
				for (k in 1:(j-1))
					{
						index1 = which(row.names(mat)==row.names(dists)[k])
						index2 = which(colnames(mat)==row.names(dists)[j])
						dists[j,k] = 1/mat[index1,index2]; dists[k,j] = 1/mat[index2,index1]
					}
			}
		for (j in 1:length(covariates))
			{
				# plot(covariates[[j]][lower.tri(covariates[[j]])], dists[lower.tri(dists)])
				mantel = mantel(covariates[[j]], dists, method="spearman", na.rm=T)
				mantel_tests_r[i,j] = mantel$statistic; mantel_tests_p[i,j] = mantel$signif
			}
	}
write.csv(mantel_tests_r, "Mantel_test_corrs.csv", row.names=F, quote=F)
write.csv(mantel_tests_p, "Mantel_test_pvals.csv", row.names=F, quote=F)
for (i in 1:length(covariates))
	{
		cat(covariate_names[i],": ",round(mean(mantel_tests_r[,i]),2)," 95% HPD [",
		    round(quantile(mantel_tests_r[,i],0.025),2),",",round(quantile(mantel_tests_r[,i],0.975),2),"]\n",sep="")
	}
		# r = 0.19, 95% HPD = [0.14, 0.26]	for proximity
		# r = 0.08, 95% HPD = [0.06, 0.11]	for families
		# r = 0.05, 95% HPD = [0.02, 0.10]	for feed
		# r = 0.03, 95% HPD = [-0.01, 0.07]	for hatcheries
		# r = 0.11, 95% HPD = [0.07, 0.15]	for rendac
		# r = 0.19, 95% HPD = [0.13, 0.24]	for veterinaries
		# r = 0.13, 95% HPD = [0.09, 0.18]  for human movements
		# r = 0.10, 95% HPD = [0.05, 0.13]  for transports
		# r = 0.16, 95% HPD = [0.13, 0.18]  for SatSQan
	
	# 17.2. Multivariate analyses (MRDM+CA)

source("CC4log_R_function.r")
vectorisationAndzTransformation = function(x)
	{ 
		x = data.matrix(x)
		x[which(!is.finite(x))] = NA
		x = x[upper.tri(x, diag=F)]
		x = (x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T))
	}
MRDM_CA = function(dists, covariates, covariate_names, outputName, CA, logistic, NAzero)
	{
		matrices = list(); matrices[[1]] = dists
		for (i in 2:(length(covariates)+1))
			{
				matrices[[i]] = covariates[[i-1]]
				if (NAzero == TRUE) matrices[[i]][is.na(matrices[[i]][])] = 0
			}
		matrices_z = lapply(matrices, vectorisationAndzTransformation)
		matrices_z = as.data.frame(matrices_z)
		names(matrices_z) = c("phylo_dists", covariate_names)
		names(matrices_z) = gsub(",",".",names(matrices_z))
		fm = paste(names(matrices_z)[1]," ~ ",names(matrices_z)[2], sep="")
		if (length(matrices_z) > 2)
			{
				for (i in 3:length(matrices_z)) fm = paste(fm," + ",names(matrices_z)[i],sep="")
			}
		if (logistic == TRUE)
			{
				matrices_z[matrices_z[] >= 0] = 1
				matrices_z[matrices_z[] < 0] = 0
			}
		if (CA == FALSE)
			{
				if (logistic == FALSE)
					{
						lr = lm(fm, matrices_z)
						return(summary(lr))
					}	else	{
						lr = glm(fm, data=matrices_z, family=binomial("logit"))
						return(summary(lr))
					}
			}	else	{
				if (logistic == FALSE)
					{
						lr = lm(fm, matrices_z)
						ca = calc.yhat(lr, prec=5)
					}	else	{
						lr = glm(fm, data=matrices_z, family=binomial("logit"))
						if (lr$converged=="FALSE") stop("Logistic model did not converge!", call.=F)
						ca = CC4log(matrices_z, "phylo_dists", covariate_names, "N")
					}
				return(ca)
			}
	}
mrdm_ca_r = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mrdm_ca_r) = covariate_names
mrdm_ca_b = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mrdm_ca_b) = covariate_names
mrdm_ca_u = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mrdm_ca_u) = covariate_names
mrdm_ca_c = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mrdm_ca_c) = covariate_names
mrdm_ca_2 = matrix(nrow=length(trees), ncol=length(covariates)); colnames(mrdm_ca_2) = covariate_names
for (i in 1:length(trees))
	{
		dists = matrix(0, nrow=length(labels), ncol=length(labels))
		row.names(dists) = labels; colnames(dists) = labels
		mat = as.matrix(adephylo::distTips(as(trees[i][[1]],"phylo4"), tips=c(labels[j],labels[k]), method="patristic"))
		row.names(mat)[row.names(mat)=="0-1-PTL"] = "1-1-PTL"; colnames(mat)[colnames(mat)=="0-1-PTL"] = "1-1-PTL"
		for (j in 2:dim(dists)[1])
			{
				for (k in 1:(j-1))
					{
						index1 = which(row.names(mat)==row.names(dists)[k])
						index2 = which(colnames(mat)==row.names(dists)[j])
						dists[j,k] = 1/mat[index1,index2]; dists[k,j] = 1/mat[index2,index1]
					}
			}
		outputName = paste0("tree_",i); logistic = FALSE; NAzero = TRUE; NAzero = FALSE
		CA = FALSE; lr = MRDM_CA(dists, covariates, covariate_names, outputName, CA, logistic, NAzero)
		CA = TRUE; ca = MRDM_CA(dists, covariates, covariate_names, outputName, CA, logistic, NAzero)		
		ca = ca$PredictorMetrics[1:(dim(ca$PredictorMetrics)[1]-1),c("r","Beta","Unique","Common")]
		colnames(ca) = c("r","beta","U","C"); globalR2 = round(lr$r.squared, 2)
		for (j in 1:length(covariate_names))
			{
				mrdm_ca_r[i,covariate_names[j]] = ca[covariate_names[j],"r"]
				mrdm_ca_b[i,covariate_names[j]] = ca[covariate_names[j],"beta"]
				mrdm_ca_u[i,covariate_names[j]] = ca[covariate_names[j],"U"]
				mrdm_ca_c[i,covariate_names[j]] = ca[covariate_names[j],"C"]
				mrdm_ca_2[i,covariate_names[j]] = globalR2
			}
	}
write.csv(mrdm_ca_r, "MRDM_CA_r_values.csv", row.names=F, quote=F)
write.csv(mrdm_ca_b, "MRDM_CA_b_values.csv", row.names=F, quote=F)
write.csv(mrdm_ca_c, "MRDM_CA_C_values.csv", row.names=F, quote=F)
write.csv(mrdm_ca_u, "MRDM_CA_U_values.csv", row.names=F, quote=F)
write.csv(mrdm_ca_2, "MRDM_CA_R2_values.csv", row.names=F, quote=F)

results = matrix(nrow=length(covariates), ncol=5)
row.names(results) = covariate_names; colnames(results) = c("r","beta","C","U","R2")
for (i in 1:length(covariates))
	{
		results[i,1] = paste0(round(mean(mrdm_ca_r[,i],na.rm=T),3)," [",
		    				  round(quantile(mrdm_ca_r[,i],0.025,na.rm=T),3),", ",round(quantile(mrdm_ca_r[,i],0.975,na.rm=T),3),"]",sep="")
		results[i,2] = paste0(round(mean(mrdm_ca_b[,i],na.rm=T),3)," [",
		    				  round(quantile(mrdm_ca_b[,i],0.025,na.rm=T),3),", ",round(quantile(mrdm_ca_b[,i],0.975,na.rm=T),3),"]",sep="")
		results[i,3] = paste0(round(mean(mrdm_ca_c[,i],na.rm=T),3)," [",
		    				  round(quantile(mrdm_ca_c[,i],0.025,na.rm=T),3),", ",round(quantile(mrdm_ca_c[,i],0.975,na.rm=T),3),"]",sep="")
		results[i,4] = paste0(round(mean(mrdm_ca_u[,i],na.rm=T),3)," [",
		    				  round(quantile(mrdm_ca_u[,i],0.025,na.rm=T),3),", ",round(quantile(mrdm_ca_u[,i],0.975,na.rm=T),3),"]",sep="")
		results[i,5] = paste0(round(mean(mrdm_ca_2[,i],na.rm=T),3)," [",
		    				  round(quantile(mrdm_ca_2[,i],0.025,na.rm=T),3),", ",round(quantile(mrdm_ca_2[,i],0.975,na.rm=T),3),"]",sep="")
	}
write.table(results, "MRDM_CA_results.txt", quote=F, sep="\t")

# 17. Analysis of the phylogenetic signal associated with several covariates

tree = readAnnotatedNexus("Analysis_2_selected.tree")
trees = read.nexus("Analysis_2_selected.trees"); Ks_inf_list = list(); Ks_ran_list = list()
metadata = read.csv("Analysis_1_141220.csv", head=T, sep=";"); cols = list(); variableIDs = list() ; variableNames = list() 
colours = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#d3d3d3")
colours = c("red","deepskyblue2","orange","green3","yellow","pink","purple","brown")
temp = metadata[,"SatSQan"]; temp[temp==0] = NA; cols[[1]] = colours[temp]; variableIDs[[1]] = "SatSQan"; variableNames[[1]] = "Satscan"
temp = metadata[,"transports"]; temp[temp==0] = NA; cols[[2]] = colours[temp]; variableIDs[[2]] = "transports"; variableNames[[2]] = "Transports"
temp = metadata[,"human_movements"]; temp[temp==0] = NA; cols[[3]] = colours[temp]; variableIDs[[3]] = "human_movements"; variableNames[[3]] = "Human movements"

pdf("Analysis_2_selected_NEW.pdf", width=11, height=5.6)
par(mfrow=c(1,3), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
for (i in 1:length(variableNames))
	{
		plot(tree, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, col="gray30", edge.color="gray30")
		for (j in 1:dim(tree$edge)[1])
			{
				if (!tree$edge[j,2]%in%tree$edge[,1])
					{
						tipLabel = tree$tip.label[tree$edge[j,2]]
						index = which(metadata[,"FASTA"]==tipLabel)
						if (length(index) == 0)
							{
								print(c(i,j,tipLabel))
							}	else	{
								if (!is.na(cols[[i]][index]))
									{
										nodelabels(node=tree$edge[j,2], pch=15, cex=0.9, col=cols[[i]][index])
										nodelabels(node=tree$edge[j,2], pch=0, cex=0.9, col="gray30", lwd=0.3)
									}
							}
					}
			}
		mtext(variableNames[[i]], side=3, line=-3, at=0.1, col="gray30", cex=0.8)
	}
dev.off()

for (i in 1:length(variableNames))
	{
		c = 0
		Ks_inf = rep(NA, length(trees)); Ks_inf_pVals = rep(NA, length(trees))
		Ks_ran = rep(NA, length(trees)); Ks_ran_pVals = rep(NA, length(trees))
		for (j in 1:length(trees))
			{
				values = matrix(nrow=length(trees[[j]]$tip.label), ncol=1)
				row.names(values) = trees[[j]]$tip.label
				for (k in 1:dim(values)[1])
					{
						index = which(metadata[,"FASTA"]==trees[[j]]$tip.label[k]) 
						values[k,1] = metadata[index,variableIDs[[i]]]
					}
				values[which(values[,1]==0),1] = NA
				test1 = phytools::phylosig(trees[[j]], values, test=T, method="K")
				Ks_inf[j] = test1$K; Ks_inf_pVals[j] = test1$P; buffer = values
				buffer[which(!is.na(buffer[,1])),1] = sample(buffer[which(!is.na(buffer[,1])),1],length(buffer[which(!is.na(buffer[,1])),1]),replace=F)
				randos = matrix(nrow=length(trees[[j]]$tip.label), ncol=1)
				randos[,1] = as.numeric(buffer); row.names(randos) = row.names(values)
				test2 = phytools::phylosig(trees[[j]], randos, test=T, method="K")
				Ks_ran[j] = test2$K; Ks_ran_pVals[j] = test2$P
				if (Ks_inf[j] > Ks_ran[j]) c = c+1
			}
		Ks_inf_list[[i]] = Ks_inf; Ks_ran_list[[i]] = Ks_ran
	}
for (i in 1:length(variableNames))
	{
		c = 0; Ks_inf = Ks_inf_list[[i]]; Ks_ran = Ks_ran_list[[i]]
		for (j in 1:length(trees))
			{
				if (Ks_inf[j] > Ks_ran[j]) c = c+1
			}
		p = c/length(trees); BF = p/(1-p)
		cat(paste0(variableNames[[i]],": K = ",round(mean(Ks_inf),2),", 95% HPD [",round(quantile(Ks_inf,0.025),2),"-",round(quantile(Ks_inf,0.975),2),"], BF = ",round(BF,1),"\n"))
	}

		# Satscan: K = 1.46, 95% HPD [0.70-2.32], BF = >999
		# Transports: K = 0.81, 95% HPD [0.36-1.37], BF = 199
		# Human movements: K = 0.76, 95% HPD [0.51-1.01], BF = >999

