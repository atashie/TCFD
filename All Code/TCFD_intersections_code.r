library(data.table)

# reading in hazard maps
wfMaps = fread('J:\\Cai_data\\TCFD\\BurntArea\\burntarea_allData_30MAY2022.csv')
gwMaps = fread('J:\\Cai_data\\TCFD\\GWstorage\\groundwater_nonegs_allData_30MAY2022.csv')
wsMaps = fread('J:\\Cai_data\\TCFD\\WaterScarcity\\WaterScarcity_allData_03JUN2022.csv')
#tcMaps = fread('J:\\Cai_data\\TCFD\\TropicalCycloneArea\\tropicalCyclone_allData_30JUN2022.csv')
mzMaps = fread('J:\\Cai_data\\TCFD\\MaizeRainfed\\maizeProductivityUnirrigated_allData_30JUN2022.csv')
rfMaps = fread('J:\\Cai_data\\TCFD\\RiverFloodArea\\RiverFloodArea_allData_01JUL2022.csv')
sweMaps =fread('J:\\Cai_data\\TCFD\\SWE\\SWE_allData_30MAY2022.csv')

	# more is bad
oldRows = which(wfMaps$Decade == '2010s')
theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
theseQuantilesVals = quantile(wfMaps$Q_50_abs[oldRows], theseQuantiles)
wfMaps$Q_25_indx = 0.01
for(i in 2:length(theseQuantilesVals))	{
	wfMaps$Q_25_indx[which(wfMaps$Q_50_abs > theseQuantilesVals[i])] = theseQuantiles[i]
}
wfMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')]	 = round(wfMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] * 100, 6)
plot(wfMaps$Q_50_abs[9:90000], wfMaps$Q_25_indx[9:90000])


	# more is bad
oldRows = which(rfMaps$Decade == '2010s')
theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
theseQuantilesVals = quantile(rfMaps$Q_50_abs[oldRows], theseQuantiles)
rfMaps$Q_25_indx = 0.01
for(i in 1:length(theseQuantilesVals))	{
	rfMaps$Q_25_indx[which(rfMaps$Q_50_abs > theseQuantilesVals[i])] = theseQuantiles[i]
}
rfMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')]	 = round(rfMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')], 6)
plot(rfMaps$Q_50_abs[9:900000], rfMaps$Q_25_indx[9:900000])


#	# more is bad
#oldRows = which(tcMaps$Decade == '2010s')
#theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
#theseQuantilesVals = quantile(tcMaps$Q_50_abs[oldRows], theseQuantiles)
#tcMaps$Q_25_indx = 0.01
#for(i in 2:length(theseQuantilesVals))	{
#	tcMaps$Q_25_indx[which(tcMaps$Q_50_abs > theseQuantilesVals[i])] = theseQuantiles[i]
#}
#tcMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')]	 = round(tcMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] * 100, 6)
#plot(tcMaps$Q_50_abs[9:90000], tcMaps$Q_25_indx[9:90000])


	# less is bad
oldRows = which(mzMaps$Decade == '2010s' & mzMaps$Q_50_abs > 0)
theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
theseQuantilesVals = rev(quantile(mzMaps$Q_50_abs[oldRows], theseQuantiles))
mzMaps$Q_25_indx = 0.01
for(i in 1:length(theseQuantilesVals))	{
	mzMaps$Q_25_indx[which(mzMaps$Q_50_abs < theseQuantilesVals[i])] = theseQuantiles[i]
}
mzMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')]	 = round(mzMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] / 10, 6)
plot(mzMaps$Q_50_abs[1:1500000], mzMaps$Q_25_indx[1:1500000])


	# less is bad
oldRows = which(sweMaps$Decade == '2010s' & sweMaps$Q_50_abs > 0)
theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
theseQuantilesVals = rev(quantile(sweMaps$Q_50_abs[oldRows], theseQuantiles))
sweMaps$Q_25_indx = 0.01
for(i in 1:length(theseQuantilesVals))	{
	sweMaps$Q_25_indx[which(sweMaps$Q_50_abs < theseQuantilesVals[i])] = theseQuantiles[i]
}
sweMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')]	 = round(sweMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] / 10, 6)
plot(sweMaps$Q_50_abs[1:15000], sweMaps$Q_25_indx[1:15000])



	# less is bad
oldRows = which(gwMaps$Decade == '2010s')
theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
theseQuantilesVals = rev(quantile(gwMaps$Q_50_abs[oldRows], theseQuantiles))
gwMaps$Q_25_indx = 0.01
for(i in 1:length(theseQuantilesVals))	{
	gwMaps$Q_25_indx[which(gwMaps$Q_50_abs < theseQuantilesVals[i])] = theseQuantiles[i]
}
gwMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')]	 = round(gwMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] / 10, 6)
plot(gwMaps$Q_50_abs[1:10000], gwMaps$Q_25_indx[1:10000])

	# less is bad
oldRows = which(wsMaps$Decade == '2010s')
wsMaps$Q_25_abs = as.numeric(wsMaps$Q_25_abs)
wsMaps$Q_50_abs = as.numeric(wsMaps$Q_50_abs)
wsMaps$Q_75_abs = as.numeric(wsMaps$Q_75_abs)
theseQuantiles = seq(0.01,0.99,0.01)#c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
theseQuantilesVals = rev(quantile(wsMaps$Q_50_abs[oldRows], theseQuantiles))
wsMaps$Q_25_indx = 0.01
for(i in 1:length(theseQuantilesVals))	{
	wsMaps$Q_25_indx[which(wsMaps$Q_50_abs < theseQuantilesVals[i])] = theseQuantiles[i]
}
wsMaps$Q_25_abs[wsMaps$Q_25_abs < -100]	 = -100
wsMaps$Q_50_abs[wsMaps$Q_50_abs < -100] = -100
wsMaps$Q_75_abs[wsMaps$Q_75_abs < -100] = -100
wsMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] = round(wsMaps[,c('Q_25_abs','Q_50_abs','Q_75_abs')] +0 , 6)
plot(wsMaps$Q_50_abs[1:90000], wsMaps$Q_25_indx[1:90000])


hazardMaps = rbind(wfMaps, gwMaps) ; rm(wfMaps); rm(gwMaps)
hazardMaps = rbind(hazardMaps, wsMaps)	; rm(wsMaps)
hazardMaps = rbind(hazardMaps, mzMaps)	; rm(mzMaps)
hazardMaps = rbind(hazardMaps, rfMaps)	; rm(rfMaps)
hazardMaps = rbind(hazardMaps, sweMaps)	; rm(sweMaps)


# reading in customer database
customerPath = 'J:\\Cai_data\\TCFD\\'
#customerCsv = 'CustomerLocData_Simplot.csv'
#customerTable = fread(paste0(customerPath, customerCsv))
#customerCsv = 'walmartstoreloaction.csv'
#customerTable = fread(paste0(customerPath, customerCsv))[seq(1,4655,10), c('name', 'state', 'city', 'latitude', 'longitude')]
#names(customerTable) = c('Name', 'Region', 'Subregion', 'Lat', 'Lon')
customerCsv = 'Richs-sample-locations.csv'
customerTable = fread(paste0(customerPath, customerCsv))

closeLats = NULL
closeLons = NULL
for(j in 1:nrow(customerTable))	{
	closeLons = c(closeLons, hazardMaps$Lon[which.min(abs(customerTable$Lon[j] - hazardMaps$Lon))])
	closeLats = c(closeLats, hazardMaps$Lat[which.min(abs(customerTable$Lat[j] - hazardMaps$Lat))])
}

closeLonLats = NULL
theseLabels = NULL
theseRegions = NULL
theseSubregions = NULL
for(k in 1:length(closeLons))	{
	newCloseLonLats = which(hazardMaps$Lat == closeLats[k] & hazardMaps$Lon == closeLons[k])
	closeLonLats = c(closeLonLats, newCloseLonLats)
	theseLabels = c(theseLabels, rep(customerTable$Name[k], length(newCloseLonLats)))
	theseRegions = c(theseRegions, rep(customerTable$Region[k], length(newCloseLonLats)))
	theseSubregions = c(theseSubregions, rep(customerTable$Subregion[k], length(newCloseLonLats)))
}


customerHazards = hazardMaps[closeLonLats,]
customerHazards$Regions = theseRegions
customerHazards$Subregions = theseSubregions
#customerHazards$Labels = unlist(strsplit(theseLabels, ','))[seq(1,(length(theseLabels)*2),2)]
customerHazards$Labels = theseLabels

			library(psych)
			customerAggHaz = cbind(as.numeric(subset(customerHazards, Hazard == "Fire (% Area Burned)")$Q_25_indx),
				as.numeric(subset(customerHazards, Hazard == "Groundwater (mm)")$Q_25_indx),
				as.numeric(subset(customerHazards, Hazard == "Water Scarcity (% Storage Remaining)")$Q_25_indx),
				as.numeric(subset(customerHazards, Hazard == 'Maize Productivity (unirrigated)')$Q_25_indx),
				as.numeric(subset(custoemrHazards, Hazard == 'Hundred Year Flood (% Area Impacted)')$Q_25_indx))
			customerAggHaz_df =	subset(customerHazards, Hazard == "Fire (% Area Burned)")
			customerAggHaz_df$Hazard =	'Aggregated Score'
			customerAggHaz_df[, c('Q_25_abs','Q_50_abs','Q_75_abs','Q_25_rel','Q_50_rel','Q_75_rel','Q_75_indx','Likelihood')] = NA
			customerAggHaz_df$Q_25_indx = apply(customerAggHaz, 1, geometric.mean)	

			customerHazards = rbind(customerHazards, customerAggHaz_df)

customerHazards$ID = 1:nrow(customerHazards)

customerHazards$Q_50_indx = 'Low'
customerHazards$Q_50_indx[which(customerHazards$Q_25_indx > 0.50)] = 'Medium'
customerHazards$Q_50_indx[which(customerHazards$Q_25_indx > 0.75)] = 'High'


thisdf = subset(customerHazards, Hazard == 'Maize Productivity (unirrigated)' & Decade == '2010s') 
plot(thisdf$Lon, thisdf$Lat, col=c('blue','orange','red')[as.numeric(thisdf$Q_50_indx)],lwd=3)

fwrite(customerHazards, paste0(customerPath, 'RichsSampleLocs_13JUL2022.csv'))
