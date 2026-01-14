library(jsonlite) # for reading json files

latLonLocator = function(cityName, countryCode, thisHit = 1)	{
	thisURL = paste0('https://nominatim.openstreetmap.org/search.php?city=',
		cityName,
		'&countrycodes=',
		countryCode,
		'&limit=9&format=jsonv2')
	return(fromJSON(thisURL)[thisHit, c('lat','lon','display_name')])	# returning only the first hit; not
}

cityName = 'Saint-Paul'
countryCode = 'US'			# look up country codes here https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2
dataOut = latLonLocator(cityName, countryCode)



library(data.table)	# for cleaner dataframes instead of matrices 
inputFileLocation = 'J:\\Downloads\\'
inputFileName = 'Template for deriving lat lons from city names - Sheet1'
inputFile = read.csv(paste0(inputFileLocation, inputFileName, '.csv'))

outputFile = data.frame(cityName = character(), countryCode = character(), Lat=numeric(), Lon=numeric(), fullLocationName=character())
for(i in 1:nrow(inputFile))	{
	outputFile[i, ] = c(inputFile[i,1], inputFile[i,2],
		latLonLocator(inputFile[i, 1], inputFile[i, 2]))
}

outputFileLocation = 'J:\\Downloads\\'
outputFileName = 'nameOfCustomerAndDate'
fwrite(outputFile, paste0(outputFileLocation, outputFileName, '.csv')) 



nrow <- nrow(test)
counter <- 1
test$lon[counter] <- 0
test$lat[counter] <- 0
while (counter <= nrow){
  CityName <- gsub(' ','%20',test$CityLong[counter]) #remove space for URLs
  CountryCode <- test$Country[counter]
  url <- paste(
    "http://nominatim.openstreetmap.org/search?city="
    , CityName
    , "&countrycodes="
    , CountryCode
    , "&limit=9&format=json"
    , sep="")
  x <- fromJSON(url)
  if(is.vector(x)){
    test$lon[counter] <- x[[1]]$lon
    test$lat[counter] <- x[[1]]$lat    
  }
  counter <- counter + 1
}