########################################################
### library(hydroGOF)		# for nse calculations
library(dataRetrieval)	# for streamflow data (I think)
library(data.table)
library(sf)
	sf::sf_use_s2(FALSE) # for problem with intersecting spherical w flat
library(lubridate)
library(ddplyr) # for left_join() merging by date
library(ncdf4)
library(magrittr)
library(maptools)
#library(FAdist) 	# for the log-pearson type III distribution
#library(moments)	# for skewness
library(lmom)		# for samlmu()



	# function for translating future flood freq according to historic flood freq
newFldFrqFunc = function(oldQ, newQ, probFld)	{
		# generating distribution using old / historic time series data
	volQ = oldQ
	sortedVolQ = sort(volQ, decreasing=TRUE)
	probQ = c(1:length(volQ)) / (length(volQ) + 1)
	trQ = 1 / probQ

	fit = samlmu(volQ)
	para3 = pelpe3(fit)

	lp3Accum = cdfpe3(sortedVolQ, para3)
	fittedTr3 = 1 / (1 - lp3Accum)

	thisProbFld = probFld
	theseFldVols = seq(max(volQ)*5, min(volQ), length.out = 10000)
	thisProbFldVol = theseFldVols[
		which.min(abs(thisProbFld - (1 / (1 -  cdfpe3(theseFldVols, para3)))))
		]

		# creating a new distribution with the future time series data data
	newVolQ = newQ
	newSortedvolQ = sort(newVolQ, decreasing = TRUE)
	newFit = samlmu(newVolQ)
	newPara3 = pelpe3(newFit)

		# identifying the return period of a 'probFld' from the historic data according to the distribution from the future data
	newLp3Accum = cdfpe3(thisProbFldVol, newPara3)
	newFittedTr3 = 1 / (1 - newLp3Accum)
	return(newFittedTr3)
}


#########################################
# reading in climai netcdf data
ncpath = "J:\\Cai_data\\TCFD\\Flash Floods\\"
	# reading ssp126 climate data 
ncname = "clm45_gfdl-esm2m_ewembi_historical_2005soc_co2_maxdis_global_monthly_1861_2005.nc4"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = rev(ncvar_get(ncin, 'lat'))	# lat is given high to low
nc_lon = ncvar_get(ncin, 'lon')
nc_date = as.Date("1661-01-01") + ncvar_get(ncin, 'time') * 30.4375# time is days after jan 1 2015
nc_years = unique(year(nc_date))
nc_maxdis_all = ncvar_get(ncin,"maxdis")

ncname_f = "clm45_gfdl-esm2m_ewembi_rcp26_2005soc_co2_maxdis_global_monthly_2006_2099.nc4"#"clm45_gfdl-esm2m_ewembi_rcp60_2005soc_co2_maxdis_global_monthly_2006_2099.nc4"  
ncfname_f = paste0(ncpath, ncname_f)
ncin_f = nc_open(ncfname_f)
nc_lat_f = ncvar_get(ncin_f, 'lat')
nc_lon_f = ncvar_get(ncin_f, 'lon')
nc_date_f = as.Date("1661-01-01") + ncvar_get(ncin_f, 'time') * 30.4375# time is days after jan 1 2015
nc_years_f = unique(year(nc_date_f))
nc_maxdis_all_f = ncvar_get(ncin_f,"maxdis")	# in lon,lat,time


ftrFldIntrv_df_5_2029 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_10_2029 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_20_2029 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_30_2029 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_50_2029 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_100_2029 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))

ftrFldIntrv_df_5_3039 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_10_3039 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_20_3039 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_30_3039 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_50_3039 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_100_3039 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))

ftrFldIntrv_df_5_4049 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_10_4049 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_20_4049 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_30_4049 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_50_4049 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_100_4049 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))

ftrFldIntrv_df_5_5059 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_10_5059 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_20_5059 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_30_5059 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_50_5059 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_100_5059 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))

ftrFldIntrv_df_5_6069 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_10_6069 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_20_6069 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_30_6069 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_50_6069 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))
ftrFldIntrv_df_100_6069 = data.frame(matrix(NA, nrow=length(nc_lat), ncol=length(nc_lon)))

old_dates_past = which(year(nc_date) == 1920)[1]:length(nc_date)
old_dates_future = 1:which(year(nc_date_f) == 2020)[12]
future_dates2029 = which(year(nc_date_f) == 2010)[1]:which(year(nc_date_f) == 2039)[12]
future_dates3039 = which(year(nc_date_f) == 2020)[1]:which(year(nc_date_f) == 2049)[12]
future_dates4049 = which(year(nc_date_f) == 2030)[1]:which(year(nc_date_f) == 2059)[12]
future_dates5059 = which(year(nc_date_f) == 2040)[1]:which(year(nc_date_f) == 2069)[12]
future_dates6069 = which(year(nc_date_f) == 2050)[1]:length(nc_date_f)

returnFreqs = c(5,10,20,30,50,100)

for(i in 1:length(nc_lat))	{
	for(j in 1:length(nc_lon))	{
		nc_maxdis = c(nc_maxdis_all[j,i,old_dates_past], nc_maxdis_all_f[j,i,old_dates_future])	# in lon,lat,time
		if(any(!is.na(nc_maxdis)))	{
			nc_maxdis_f_2029 = nc_maxdis_all_f[j,i,future_dates2029] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_maxdis_f_3039 = nc_maxdis_all_f[j,i,future_dates3039] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_maxdis_f_4049 = nc_maxdis_all_f[j,i,future_dates4049] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_maxdis_f_5059 = nc_maxdis_all_f[j,i,future_dates5059] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time
			nc_maxdis_f_6069 = nc_maxdis_all_f[j,i,future_dates6069] #ncvar_get(ncin_f,"maxdis")[j,i,]	# in lon,lat,time			
			if(any(nc_maxdis > 0) & any(nc_maxdis_f_2029 > 0) & any(nc_maxdis_f_3039) & any(nc_maxdis_f_4049) & any(nc_maxdis_f_5059) & any(nc_maxdis_f_6069))	{
					
				print(c(i,j))

				randNoise = seq(0.1^3,0.1^2,length.out = 100)
				nc_maxdis_y = NULL
				for(kh in 1:(length(nc_maxdis)/12))	{
					nc_maxdis_y = c(nc_maxdis_y, max(nc_maxdis[(1:12)+(12*(kh-1))])+sample(randNoise,1))
				}

				nc_maxdis_y_f_2029 = NULL
				for(kh in 1:(length(nc_maxdis_f_2029)/12))	{
					nc_maxdis_y_f_2029 = c(nc_maxdis_y_f_2029, max(nc_maxdis_f_2029[(1:12)+(12*(kh-1))])+sample(randNoise,1))
				}

				nc_maxdis_y_f_3039 = NULL
				for(kh in 1:(length(nc_maxdis_f_3039)/12))	{
					nc_maxdis_y_f_3039 = c(nc_maxdis_y_f_3039, max(nc_maxdis_f_3039[(1:12)+(12*(kh-1))])+sample(randNoise,1))
				}

				nc_maxdis_y_f_4049 = NULL
				for(kh in 1:(length(nc_maxdis_f_4049)/12))	{
					nc_maxdis_y_f_4049 = c(nc_maxdis_y_f_4049, max(nc_maxdis_f_4049[(1:12)+(12*(kh-1))])+sample(randNoise,1))
				}

				nc_maxdis_y_f_5059 = NULL
				for(kh in 1:(length(nc_maxdis_f_5059)/12))	{
					nc_maxdis_y_f_5059 = c(nc_maxdis_y_f_5059, max(nc_maxdis_f_5059[(1:12)+(12*(kh-1))])+sample(randNoise,1))
				}

				nc_maxdis_y_f_6069 = NULL
				for(kh in 1:(length(nc_maxdis_f_6069)/12))	{
					nc_maxdis_y_f_6069 = c(nc_maxdis_y_f_6069, max(nc_maxdis_f_6069[(1:12)+(12*(kh-1))])+sample(randNoise,1))
				}

				ftrFldIntrv_df_5_2029[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_2029, probFld = 5)
				ftrFldIntrv_df_5_3039[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_3039, probFld = 5)
				ftrFldIntrv_df_5_4049[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_4049, probFld = 5)
				ftrFldIntrv_df_5_5059[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_5059, probFld = 5)
				ftrFldIntrv_df_5_6069[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_6069, probFld = 5)
				ftrFldIntrv_df_10_2029[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_2029, probFld = 10)
				ftrFldIntrv_df_10_3039[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_3039, probFld = 10)
				ftrFldIntrv_df_10_4049[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_4049, probFld = 10)
				ftrFldIntrv_df_10_5059[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_5059, probFld = 10)
				ftrFldIntrv_df_10_6069[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_6069, probFld = 10)
				ftrFldIntrv_df_20_2029[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_2029, probFld = 20)
				ftrFldIntrv_df_20_3039[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_3039, probFld = 20)
				ftrFldIntrv_df_20_4049[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_4049, probFld = 20)
				ftrFldIntrv_df_20_5059[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_5059, probFld = 20)
				ftrFldIntrv_df_20_6069[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_6069, probFld = 20)
				ftrFldIntrv_df_30_2029[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_2029, probFld = 30)
				ftrFldIntrv_df_30_3039[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_3039, probFld = 30)
				ftrFldIntrv_df_30_4049[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_4049, probFld = 30)
				ftrFldIntrv_df_30_5059[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_5059, probFld = 30)
				ftrFldIntrv_df_30_6069[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_6069, probFld = 30)
				ftrFldIntrv_df_50_2029[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_2029, probFld = 50)
				ftrFldIntrv_df_50_3039[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_3039, probFld = 50)
				ftrFldIntrv_df_50_4049[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_4049, probFld = 50)
				ftrFldIntrv_df_50_5059[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_5059, probFld = 50)
				ftrFldIntrv_df_50_6069[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_6069, probFld = 50)
				ftrFldIntrv_df_100_2029[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_2029, probFld = 100)
				ftrFldIntrv_df_100_3039[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_3039, probFld = 100)
				ftrFldIntrv_df_100_4049[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_4049, probFld = 100)
				ftrFldIntrv_df_100_5059[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_5059, probFld = 100)
				ftrFldIntrv_df_100_6069[i,j] = newFldFrqFunc(oldQ = nc_maxdis_y, newQ = nc_maxdis_y_f_6069, probFld = 100)
			}
		} 
	}
#	print(ftrFldIntrv_df)
}
write.csv(ftrFldIntrv_df_5_2029, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_5yrflod_CLM_2029.csv")
write.csv(ftrFldIntrv_df_5_3039, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_5yrflod_CLM_3039.csv")
write.csv(ftrFldIntrv_df_5_4049, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_5yrflod_CLM_4049.csv")
write.csv(ftrFldIntrv_df_5_5059, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_5yrflod_CLM_5059.csv")
write.csv(ftrFldIntrv_df_5_6069, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_5yrflod_CLM_6069.csv")

write.csv(ftrFldIntrv_df_10_2029, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_10yrflod_CLM_2029.csv")
write.csv(ftrFldIntrv_df_10_3039, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_10yrflod_CLM_3039.csv")
write.csv(ftrFldIntrv_df_10_4049, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_10yrflod_CLM_4049.csv")
write.csv(ftrFldIntrv_df_10_5059, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_10yrflod_CLM_5059.csv")
write.csv(ftrFldIntrv_df_10_6069, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_10yrflod_CLM_6069.csv")

write.csv(ftrFldIntrv_df_20_2029, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_20yrflod_CLM_2029.csv")
write.csv(ftrFldIntrv_df_20_3039, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_20yrflod_CLM_3039.csv")
write.csv(ftrFldIntrv_df_20_4049, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_20yrflod_CLM_4049.csv")
write.csv(ftrFldIntrv_df_20_5059, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_20yrflod_CLM_5059.csv")
write.csv(ftrFldIntrv_df_20_6069, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_20yrflod_CLM_6069.csv")

write.csv(ftrFldIntrv_df_30_2029, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_30yrflod_CLM_2029.csv")
write.csv(ftrFldIntrv_df_30_3039, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_30yrflod_CLM_3039.csv")
write.csv(ftrFldIntrv_df_30_4049, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_30yrflod_CLM_4049.csv")
write.csv(ftrFldIntrv_df_30_5059, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_30yrflod_CLM_5059.csv")
write.csv(ftrFldIntrv_df_30_6069, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_30yrflod_CLM_6069.csv")

write.csv(ftrFldIntrv_df_50_2029, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_50yrflod_CLM_2029.csv")
write.csv(ftrFldIntrv_df_50_3039, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_50yrflod_CLM_3039.csv")
write.csv(ftrFldIntrv_df_50_4049, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_50yrflod_CLM_4049.csv")
write.csv(ftrFldIntrv_df_50_5059, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_50yrflod_CLM_5059.csv")
write.csv(ftrFldIntrv_df_50_6069, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_50yrflod_CLM_6069.csv")

write.csv(ftrFldIntrv_df_100_2029, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_100yrflod_CLM_2029.csv")
write.csv(ftrFldIntrv_df_100_3039, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_100yrflod_CLM_3039.csv")
write.csv(ftrFldIntrv_df_100_4049, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_100yrflod_CLM_4049.csv")
write.csv(ftrFldIntrv_df_100_5059, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_100yrflod_CLM_5059.csv")
write.csv(ftrFldIntrv_df_100_6069, "J:\\Cai_data\\TCFD\\Flash Floods\\rcp26_100yrflod_CLM_6069.csv")



#ftrFldIntrv_df = as.data.frame(read.csv("J:\\Cai_data\\TCFD\\Flash Floods\\test_out.csv"))[,-1]


library(colorRamps)
col5 <- colorRampPalette(c('blue4', 'gray96', 'yellow4'))  #create color ramp starting from blue to red

par(mar=c(2.1,2.1,2.1,1))
image(nc_lon,rev(nc_lat),as.matrix(ftrFldIntrv_df[,ncol(ftrFldIntrv_df):1]), ylim = c(-57,85),
	col=col5(n=11), breaks = seq(0,40,length.out=12),
	ylab='Lat', xlab = 'lon', main = "Changes in flood frequency in 21st C")
data(wrld_simpl)
plot(wrld_simpl,add=TRUE)

if (any(is.na(as.matrix(ftrFldIntrv_df[,ncol(ftrFldIntrv_df):1])))) {
    mmat <- ifelse(is.na(as.matrix(ftrFldIntrv_df[,ncol(ftrFldIntrv_df):1])), 1, NA)
    image(mmat, axes = FALSE, xlab = "", ylab = "", col = 'white', useRaster=TRUE, add = TRUE)
  }

ftrFldIntrv_mt = as.matrix(unlist(ftrFldIntrv_df))

pt1
ftrFldIntrv_sf = st_as_sf(as.vector(ftrFldIntrv_df), coords = 

















ftrFldIntrv_df[is.na(ftrFldIntrv_df)] = 0

# converting df into array 
lon2 = nc_lon	; lat2 = nc_lat
nlon2 = length(nc_lon)	; nlat2 = length(nc_lat)
#time2 = time; tunits2 = tunits	; nt2 = nt

# convert to 2d array
ftrFldIntrv_ar <- array(ftrFldIntrv_df, dim=c(nlon2,nlat2))
dim(ftrFldIntrv_ar)

# plot to check creation of arrays
library(lattice)
library(RColorBrewer)

grid <- expand.grid(lon=nc_lon, lat=nc_lat)

levelplot(ftrFldIntrv_ar ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
  col.regions=(rev(brewer.pal(10,"RdBu"))), main="MTWA (C)")



	#names of path, file, and dname
ncpath = "J:\\Cai_data\\TCFD\\Flash Floods\\test_out.csv"
ncname = "GLM_ftrFlood"
ncfname = paste(ncpath, ncname, '.nc', sep='')
dname = 'fldInt'


	#define dimensions
londim = ncdim_def('lon','degrees',as.double(nc_lon))
latdim = ncdim_def('lat','degrees',as.double(nc_lat))
#timedim <- ncdim_def("time",tunits3,as.double(time3))

# define variables
fillvalue = NA
dlname = "Flood Interval in Years"
fldInt_def = ncvar_def("fldInt","yrs",list(londim,latdim),fillvalue,dlname,prec="single")
#dlname = "mean_temperture_warmest_month"
#mtwa.def <- ncvar_def("mtwa","deg_C",list(londim,latdim,timedim),fillvalue,dlname,prec="single")


# create netCDF file and put arrays
#ncout <- nc_create(ncfname,list(tmp_def,mtco.def,mtwa.def,mat.def),force_v4=TRUE)
ncout <- nc_create(ncfname,list(fldInt_def),force_v4=FALSE)

# put variables
#ncvar_put(ncout,tmp_def,tmp_array3)
#ncvar_put(ncout,mtwa.def,mtwa_array3)
ncvar_put(ncout,fldInt_def,ftrFldIntrv_df)
ncvar_put(ncout,mat.def,mat_array3)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

# add global attributes
ncatt_put(ncout,0,"title",title$value)
ncatt_put(ncout,0,"institution",institution$value)
ncatt_put(ncout,0,"source",datasource$value)
ncatt_put(ncout,0,"references",references$value)
history <- paste("P.J. Bartlein", date(), sep=", ")
ncatt_put(ncout,0,"history",history)
ncatt_put(ncout,0,"Conventions",Conventions$value)

# Get a summary of the created file:
ncout





dim1 <- dim.def.ncdf( name = "lat", units = "degrees", vals = as.double(nc_lat) )
dim2 <- dim.def.ncdf( name = "lon", units = "degrees", vals = as.double(nc_lon) )
ftrFldIntrv_mt = as.matrix(ftrFldIntrv_df)

# define the EMPTY (elevation) netcdf variable

varz <- var.def.ncdf(name = "floodFreq", units =  "years", 
                    dim = list(dim1, dim2), missval = NA, 
                    longname ="Flood Frequency of what was a 20 year flood")		

nc.ex <- create.ncdf( filename = "J:\\Cai_data\\TCFD\\Flash Floods\\test_out.nc", vars = varz )
put.var.ncdf(nc = nc.ex, varid = varz, vals = ftrFldIntrv_mt)




ggplot(aes(x=X1, y=X2, fill=value), data=mpr) + geom_raster() 
		
myClrs =  colorRamp(c("red", "blue"))
plot(x=nc_lon, y=nc_lat, 



lon1 <- ncdim_def( "Lon", "degrees_east", unique(DB$lon))
lat1 <- ncdim_def( "Lat", "degrees_north", unique(DB$lat))
record <- ncdim_def( "record", "files", list_obs$record_uid,unlim=T)
mv <- -999 # missing value to use
tmpmean_ <- ncvar_def( "tmpmean_", "degrees", list(lon1,lat1,record),mv)

time<-list_obs$record_uid

ncnew <- nc_create( "tmpmean_.nc",tmpmean_)

for( i in 1:length(record))
  ncvar_put( ncnew, tmpmean_, DB$tmpmean_, start=c(1,1,i), 
             count=c(-1,-1,1))
		
		
		
		
		

cal_out = data.frame(sfcf=rep(NA,length(HBV_in_ls)))

#"https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23"
for(thisClmMod in c(4,6:15))	{
	HBVoutCol = HBVoutCol + 1

#	ncvar_thismod_ppt =  ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] ugh annoying
	ncvar_thismod_ppt_1 = ncvar_get(ncin_ppt_1, "tp")[ ,thisClmMod , , ]	# is [time,model,lon,lat] and now merging two separate netcdfs
	ncvar_thismod_ppt_2 = ncvar_get(ncin_ppt_2, "tp")[date_range_ppt2,thisClmMod , , ]
	ncvar_thismod_tmin = ncvar_get(ncin, "t2m_min")[, , , thisClmMod]
	ncvar_thismod_tmax = ncvar_get(ncin, "t2m_max")[, , , thisClmMod]



