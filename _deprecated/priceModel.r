library(crypto2)
library(data.table)
library(lubridate)
library(CAST)
library(caret)
library(dplyr)
library(zoo)


library(ncdf4)

#################################################################
## Reading in historical cali FNF data
	# for the Sacramento River
BND = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=BND&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")
BND$Date = ymd(unlist(strsplit(BND$DATE.TIME, " "))[seq(1,nrow(BND)*2,2)])

ORO = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=ORO&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")
ORO$Date = ymd(unlist(strsplit(ORO$DATE.TIME, " "))[seq(1,nrow(ORO)*2,2)])
 
YRS = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=YRS&SensorNums=8&dur_code=D&Start=1906-01-01&End=2020-12-31")
YRS$Date = ymd(unlist(strsplit(YRS$DATE.TIME, " "))[seq(1,nrow(YRS)*2,2)])

FOL = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=FOL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-03-23")
FOL$Date = ymd(unlist(strsplit(FOL$DATE.TIME, " "))[seq(1,nrow(FOL)*2,2)])

	# data for the San Jaoquin Valley
NML = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=NML&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")#Stanislaus River inflow to New Melones Lake
NML$Date = ymd(unlist(strsplit(NML$DATE.TIME, " "))[seq(1,nrow(NML)*2,2)])

DNP = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=DNP&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-03-23")#Tuolumne River inflow to New Don Pedro Reservoir
DNP$Date = ymd(unlist(strsplit(DNP$DATE.TIME, " "))[seq(1,nrow(DNP)*2,2)])
 
EXC = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=EXC&SensorNums=8&dur_code=D&Start=1900-01-01&End=2022-03-23")#Merced River inflow to Lake McClure
EXC$Date = ymd(unlist(strsplit(EXC$DATE.TIME, " "))[seq(1,nrow(EXC)*2,2)])

MIL = read.csv("https://cdec.water.ca.gov/dynamicapp/req/CSVDataServlet?Stations=MIL&SensorNums=76&dur_code=D&Start=1900-01-01&End=2022-03-23")#San Joaquin River inflow to Millerton Lake
MIL$Date = ymd(unlist(strsplit(MIL$DATE.TIME, " "))[seq(1,nrow(MIL)*2,2)])



#########################################
# reading in climai netcdf clim data
ncpath = "C:\\Users\\arik\\Documents\\climateAi\\Wonderful\\"
	# reading in tmin 
ncname = "CA_era5_tmin.nc"  
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_lat = ncvar_get(ncin, 'latitude')
nc_lon = ncvar_get(ncin, 'longitude')
nc_date = as.Date("1990-01-01") + ncvar_get(ncin, 'time') # time is days after jan 1 1990
nc_tmin = ncvar_get(ncin, "t2m_min")

	# reading in tmin 
ncname = "CA_era5_tmax.nc"  		# check to ensure the lat lons / times are the same for all files 
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_tmax = ncvar_get(ncin, "t2m_max")

	# reading in ppt 
ncname = "CA_era5_tp.nc"  			# check to ensure the lat lons / times are the same for all files 
ncfname = paste0(ncpath, ncname)
ncin = nc_open(ncfname)
nc_ppt = ncvar_get(ncin, "tp")







predict_lag = 60 	# how many days in the future do we want to predict
maxlag = 90			# maximum number of days we'll test for correlation

waterPrice = fread("J:\\Cai_data\\WestWater\\NQH2O history.csv")
waterPrice$Date = as.Date(waterPrice$Date, format="%m/%d/%Y")
fullDates = data.frame(seq(as.Date(waterPrice$Date[1]), as.Date(tail(waterPrice$Date,1)),1))
names(fullDates) = "Date"
mp = full_join(fullDates, waterPrice,by='Date')
mp$value = na.approx(mp$N)
mp$BND = left_join(mp, BND[,c('Date','VALUE')], by='Date')$VALUE
mp$ORO = left_join(mp, ORO[,c('Date','VALUE')], by='Date')$VALUE
mp$YRS = left_join(mp, YRS[,c('Date','VALUE')], by='Date')$VALUE
mp$FOL = left_join(mp, FOL[,c('Date','VALUE')], by='Date')$VALUE
mp$NML = left_join(mp, NML[,c('Date','VALUE')], by='Date')$VALUE
mp$DNP = left_join(mp, DNP[,c('Date','VALUE')], by='Date')$VALUE
mp$EXC = left_join(mp, EXC[,c('Date','VALUE')], by='Date')$VALUE
mp$MIL = left_join(mp, MIL[,c('Date','VALUE')], by='Date')$VALUE
mp$month = month(mp$Date)

mp$pct_lag = c(diff(mp$value, predict_lag) / mp$value[-((nrow(mp)-(predict_lag-1)):nrow(mp))],rep(NA, predict_lag))

	for(laggards in 1:maxlag)	{
		rmv_dys = ((nrow(mp)-(laggards-1)):nrow(mp))

		strmnonas = as.numeric(mp$BND)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('BND_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$ORO)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('ORO_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$YRS)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('YRS_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$FOL)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('FOL_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$NML)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('NML_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$DNP)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('DNP_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$EXC)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('EXC_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$MIL)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('MIL_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

#		mp[ , paste0("pct_",laggards)] = 		c(rep(NA, laggards), diff(mp$value, laggards) 				 / mp$value[-rmv_dys])
	}

	
			# correlation of each data
	cor_results = matrix(data = NA, nrow = maxlag, ncol = 8)
	for(laggards in 1:maxlag)	{
		cor_results[laggards,1] = cor.test(unlist(mp[ , paste0("BND_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,2] = cor.test(unlist(mp[ , paste0("ORO_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,3] = cor.test(unlist(mp[ , paste0("YRS_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,4] = cor.test(unlist(mp[ , paste0("FOL_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,5] = cor.test(unlist(mp[ , paste0("NML_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,6] = cor.test(unlist(mp[ , paste0("DNP_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,7] = cor.test(unlist(mp[ , paste0("EXC_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
		cor_results[laggards,8] = cor.test(unlist(mp[ , paste0("MIL_sum",laggards)]), unlist(mp$pct_lag), method="kendall")$p.value
	}
#	print(cor_results, round = 3)

	best_BND_lag = order(cor_results[,1])[1]
	best_ORO_lag = order(cor_results[,2])[1]
	best_YRS_lag = order(cor_results[,3])[1]
	best_FOL_lag = order(cor_results[,4])[1]
	best_NML_lag = order(cor_results[,5])[1]
	best_DNP_lag = order(cor_results[,6])[1]
	best_EXC_lag = order(cor_results[,7])[1]
	best_MIL_lag = order(cor_results[,8])[1]
	
	
	
#	mp$Day = wday(mp$Date, label=FALSE)
#	mp$Month = month(mp$Date)
	mp_subset = mp[,c("Date","value","pct_lag")]
	
	mp_subset$BND1 = mp[ , paste0("BND_sum", best_BND_lag)]
	mp_subset$ORO1 = mp[ , paste0("ORO_sum", best_ORO_lag)]
	mp_subset$YRS1 = mp[ , paste0("YRS_sum", best_YRS_lag)]
	mp_subset$FOL1 = mp[ , paste0("FOL_sum", best_FOL_lag)]
	mp_subset$NML1 = mp[ , paste0("NML_sum", best_NML_lag)]
	mp_subset$DNP1 = mp[ , paste0("DNP_sum", best_DNP_lag)]
	mp_subset$EXC1 = mp[ , paste0("EXC_sum", best_EXC_lag)]
	mp_subset$MIL1 = mp[ , paste0("MIL_sum", best_MIL_lag)]
	mp_subset$month = mp$month
	

	mp_test = subset(mp_subset, Date < as.Date("2020-01-01") & Date > as.Date("2014-01-01"))
	mp_holdout = subset(mp_subset, Date > as.Date("2020-01-01") & Date < as.Date("2021-09-01"))
	
	model_LLO = train(mp_test[,-c(1,2,3)],
						mp_test[,3],
						metric="Rsquared",
						method='rf',
						#tuneLength=?,
						imporance=TRUE,
						ntree=2000,
						trControl = trainControl(
							method='repeatedcv',
							number=5,
							repeats=8))
							#method='cv'))
		#					index = indices$index))
		#					p=.75))

	model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars

	predict_val = predict(model_LLO, newdata = mp_holdout)
	summary(lm(predict_val ~ mp_holdout$pct_lag))$adj.r.squared
	cor.test(mp_holdout$pct_lag, predict_val)$p.value

	summary(lm((mp_holdout$pct_lag * mp_holdout$value + mp_holdout$value) ~ (predict_val *  mp_holdout$value + mp_holdout$value)))

	plot(mp_holdout$Date, mp_holdout$pct_lag * mp_holdout$value + mp_holdout$value, type='l', lwd=1.5, col='black')
	lines(mp_holdout$Date[-c(1:predict_lag)],
		headpredict_val[-c(1:predict_lag)] *  head(mp_holdout$value,nrow(mp_holdout) - predict_lag) + mp_holdout$value[-c(1:predict_lag)], col='red', lwd=2)
	lines(smooth.spline(mp_holdout$Date, predict_val *  mp_holdout$value + mp_holdout$value, spar=0.35), col='red', lwd=2)
	
	predictions = predict(model_LLO, newdata = mp_subset[-c(1:max_lag),])
	


	p_l = length(predictions)
	sb_l = nrow(mp_subset)
	which_col = 5
	for(j in rev(0:predict_lag))	{
		which_col = which_col + 1
		coin_results[which_coin,which_col] = (1 + predictions[p_l - j]) * mp_subset$avg[sb_l - j]
	}
	coin_results[which_coin,which_col+1] = predictions[p_l] * 100
	print(coin_results[order(coin_results[,4], decreasing=TRUE),]) 
	write.csv(coin_results[order(coin_results[,4], decreasing=TRUE),], paste0("C:\\Users\\arik\\Documents\\money_please\\coinbets\\60day",format(last(mp$Date),"%y%m%d"),".csv"))
#	}
	}
}




























#########################################################################################
## predicting / describing annual average price acccording to water year totals
waterPrice = fread("J:\\Cai_data\\WestWater\\NQH2O history.csv")
waterPrice$Date = as.Date(waterPrice$Date, format="%m/%d/%Y")
fullDates = data.frame(seq(as.Date(waterPrice$Date[1]), as.Date(tail(waterPrice$Date,1)),1))
names(fullDates) = "Date"
mp = full_join(fullDates, waterPrice,by='Date')
mp$value = na.approx(mp$N)

BND$VALUE[is.na(as.numeric(BND$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
ORO$VALUE[is.na(as.numeric(ORO$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
YRS$VALUE[is.na(as.numeric(YRS$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
FOL$VALUE[is.na(as.numeric(FOL$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
NML$VALUE[is.na(as.numeric(NML$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
DNP$VALUE[is.na(as.numeric(DNP$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
EXC$VALUE[is.na(as.numeric(EXC$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values
MIL$VALUE[is.na(as.numeric(MIL$VALUE))] = 0# ugh stupid data has NAs and '---' indicating missing values

mp$BND = as.numeric(left_join(mp, BND[,c('Date','VALUE')], by='Date')$VALUE)	# ugh stupid data has NAs and '---' indicating missing values
mp$ORO = as.numeric(left_join(mp, ORO[,c('Date','VALUE')], by='Date')$VALUE)
mp$YRS = as.numeric(left_join(mp, YRS[,c('Date','VALUE')], by='Date')$VALUE) # has substantial missing data
mp$FOL = as.numeric(left_join(mp, FOL[,c('Date','VALUE')], by='Date')$VALUE)
mp$NML = as.numeric(left_join(mp, NML[,c('Date','VALUE')], by='Date')$VALUE)
mp$DNP = as.numeric(left_join(mp, DNP[,c('Date','VALUE')], by='Date')$VALUE)
mp$EXC = as.numeric(left_join(mp, EXC[,c('Date','VALUE')], by='Date')$VALUE)
mp$MIL = as.numeric(left_join(mp, MIL[,c('Date','VALUE')], by='Date')$VALUE)
mp$FNF = apply(mp[,c(4,5,7:11)], 1, sum) / 504166.667 # converting from dumbo cfs to dumbier mil acr ft / day
meanFNFannum = mean(mp$FNF, na.rm=TRUE) * 365.25
mp$month = month(mp$Date)
mp$year = year(mp$Date)

numWYs = length(unique(mp$year))
yearsToTest = unique(mp$year)[-c(1,(numWYs))]


priceMtrx = matrix(data = NA, ncol = 6, nrow = length(yearsToTest))
colnames(priceMtrx) = c('year','mnPrice','mdPrice','lyWtr','tyWtr','nyWtr')
priceDf = as.data.frame(priceMtrx)
iter = 0
for(i in yearsToTest)	{
	iter = iter+1
		#previous water year
	lyLastWY = which(mp$year == (i-2) & mp$month >= 10)
	thLastWY = which(mp$year == (i-1) & mp$month < 10)

		#this water year
	lyThisWY = which(mp$year == (i-1) & mp$month >= 10)
	thThisWY = which(mp$year == i & mp$month < 10)

		#next water year
	lyNextWY = which(mp$year == i & mp$month >= 10)
	thNextWY = which(mp$year == (i+1) & mp$month < 10)

	priceDf$year[iter] = i
	priceDf$mnPrice[iter] = mean(mp$value[c(lyThisWY,thThisWY)])
	priceDf$mdPrice[iter] = median(mp$value[c(lyThisWY,thThisWY)])
	priceDf$lyWtr[iter] = sum(mp$FNF[c(lyLastWY,thLastWY)])
	priceDf$tyWtr[iter] = sum(mp$FNF[c(lyThisWY,thThisWY)])
	priceDf$nyWtr[iter] = sum(mp$FNF[c(lyNextWY,thNextWY)])
}

summary(lm(priceDf$mnPrice ~ priceDf$lyWtr))
summary(lm(priceDf$mnPrice ~ priceDf$tyWtr))
summary(lm(priceDf$mnPrice ~ priceDf$nyWtr))


plot(priceDf$year, (priceDf$mnPrice - mean(priceDf$mnPrice)) / sd(priceDf$mnPrice), type = 'l')
lines(priceDf$year, (priceDf$lyWtr - mean(priceDf$lyWtr)) / sd(priceDf$lyWtr), type = 'l', col='red4')
lines(priceDf$year, (priceDf$tyWtr - mean(priceDf$tyWtr)) / sd(priceDf$tyWtr), type = 'l', col='green3')
lines(priceDf$year, (priceDf$nyWtr - mean(priceDf$nyWtr, na.rm=TRUE)) / sd(priceDf$nyWtr, na.rm=TRUE), type = 'l', col='blue4')


par(cex=1.3)
plot(priceDf$lyWtr, priceDf$mnPrice, lwd=4, col='red2', pch=2, ylab="Annual Average Index Value (USD)", xlab="Total Surface Water Supply (mil-acr-ft)", frame.plot=FALSE)
box(bty='l')
abline(lm(priceDf$mnPrice ~ priceDf$lyWtr), col='red2', lwd=3)
text(35,650,paste0('Prev. Yr r2=',round(summary(lm(priceDf$mnPrice ~ priceDf$lyWtr))$r.squared,2)), col='red2',cex=1.5)
#points(priceDf$tyWtr, priceDf$mnPrice, lwd=3, col='green4', pch=3)
#abline(lm(priceDf$mnPrice ~ priceDf$tyWtr), col='green4', lwd=2)
#text(33,625,paste0('Current Yr r2=',round(summary(lm(priceDf$mnPrice ~ priceDf$tyWtr))$r.squared,2)), col='green4',cex=1.4)
points(priceDf$nyWtr, priceDf$mnPrice, lwd=2, col='blue3', pch=4)
text(35,600,paste0('Next Yr r2=',round(summary(lm(priceDf$mnPrice ~ priceDf$nyWtr))$r.squared,2)), col='blue3',cex=1.3)






mp$rankSeason = (mp$year - mp$year[1]) * 4 + floor(mp$month / 4) + 1 # sequential seasons
numSeasons = length(unique(mp$rankSeason))
seasonsToTest = unique(mp$rankSeason)[-c(1:4,(numSeasons-4):numSeasons)]

priceMtrxS = matrix(data = NA, ncol = 13, nrow = length(seasonsToTest))
colnames(priceMtrxS) = c('year','season','mnPrice','mdPrice','pv4Wtr','pv3Wtr','pv2Wtr','pv1Wtr','thsWtr','ft1Wtr','ft2Wtr','ft3Wtr','ft4Wtr')
priceDfS = as.data.frame(priceMtrxS)
iter = 0
for(thisSeason in seasonsToTest)	{
	iter = iter+1
			#previous seasons
	pv4 = which(mp$rankSeason == thisSeason - 4)
	pv3 = which(mp$rankSeason == thisSeason - 3)
	pv2 = which(mp$rankSeason == thisSeason - 2)
	pv1 = which(mp$rankSeason == thisSeason - 1)
			#this season
	thisS = which(mp$rankSeason == thisSeason)
			#next season + 4
	ft1 = which(mp$rankSeason == thisSeason + 1)
	ft2 = which(mp$rankSeason == thisSeason + 2)
	ft3 = which(mp$rankSeason == thisSeason + 3)
	ft4 = which(mp$rankSeason == thisSeason + 4)

	priceDfS$year[iter] = mp$year[thisS][1]
	priceDfS$season[iter] = thisSeason
	priceDfS$mnPrice[iter] = mean(mp$value[thisS])
	priceDfS$mdPrice[iter] = median(mp$value[thisS])
	priceDfS$pv4Wtr[iter] = sum(mp$FNF[pv4])
	priceDfS$pv3Wtr[iter] = sum(mp$FNF[pv3])
	priceDfS$pv2Wtr[iter] = sum(mp$FNF[pv2])
	priceDfS$pv1Wtr[iter] = sum(mp$FNF[pv1])
	priceDfS$thsWtr[iter] = sum(mp$FNF[thisS])
	priceDfS$ft1Wtr[iter] = sum(mp$FNF[ft1])
	priceDfS$ft2Wtr[iter] = sum(mp$FNF[ft2])
	priceDfS$ft3Wtr[iter] = sum(mp$FNF[ft3])
	priceDfS$ft4Wtr[iter] = sum(mp$FNF[ft4])

}


summary(lm(priceDfS$mnPrice ~ priceDfS$pv4Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$pv3Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$pv2Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$pv1Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$thsWtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$ft1Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$ft2Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$ft3Wtr))$r.squared
summary(lm(priceDfS$mnPrice ~ priceDfS$ft4Wtr))$r.squared


summary(lm(priceDfS$mnPrice ~ priceDfS$pv4Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$pv3Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$pv2Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$pv1Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$thsWtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$ft1Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$ft2Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$ft3Wtr))$adj
summary(lm(priceDfS$mnPrice ~ priceDfS$ft4Wtr))$adj

summary(lm(priceDfS$mdPrice ~ priceDfS$pv4Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$pv3Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$pv2Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$pv1Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$thsWtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$ft1Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$ft2Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$ft3Wtr))
summary(lm(priceDfS$mdPrice ~ priceDfS$ft4Wtr))



predPwr = c(summary(lm(priceDfS$mnPrice ~ priceDfS$pv4Wtr))$r.squared, 
	summary(lm(priceDfS$mnPrice ~ priceDfS$pv3Wtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$pv2Wtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$pv1Wtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$thsWtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$ft1Wtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$ft2Wtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$ft3Wtr))$r.squared,
	summary(lm(priceDfS$mnPrice ~ priceDfS$ft4Wtr))$r.squared)

barplot(predPwr,names.arg=c(-4:4), ylab="Predicative Power (r-squared)", col=c('red4','red3','red2','red1','grey40','blue1','blue2','blue3','blue4'))

















#############################################################################################################
## ML for predicting actual price instead of price swings



predict_lag = 60 	# how many days in the future do we want to predict
maxlag = 90			# maximum number of days we'll test for correlation

waterPrice = fread("J:\\Cai_data\\WestWater\\NQH2O history.csv")
waterPrice$Date = as.Date(waterPrice$Date, format="%m/%d/%Y")
fullDates = data.frame(seq(as.Date(waterPrice$Date[1]), as.Date(tail(waterPrice$Date,1)),1))
names(fullDates) = "Date"
mp = full_join(fullDates, waterPrice,by='Date')
mp$value = na.approx(mp$N)
mp$valueFtr = c(mp$value[-c(1:predict_lag)] ,rep(NA, predict_lag))
mp$BND = left_join(mp, BND[,c('Date','VALUE')], by='Date')$VALUE
mp$ORO = left_join(mp, ORO[,c('Date','VALUE')], by='Date')$VALUE
mp$YRS = left_join(mp, YRS[,c('Date','VALUE')], by='Date')$VALUE
mp$FOL = left_join(mp, FOL[,c('Date','VALUE')], by='Date')$VALUE
mp$NML = left_join(mp, NML[,c('Date','VALUE')], by='Date')$VALUE
mp$DNP = left_join(mp, DNP[,c('Date','VALUE')], by='Date')$VALUE
mp$EXC = left_join(mp, EXC[,c('Date','VALUE')], by='Date')$VALUE
mp$MIL = left_join(mp, MIL[,c('Date','VALUE')], by='Date')$VALUE
mp$month = month(mp$Date)
mp$doy = yday(mp$Date)


	for(laggards in 1:maxlag)	{
		rmv_dys = ((nrow(mp)-(laggards-1)):nrow(mp))

		strmnonas = as.numeric(mp$BND)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('BND_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$ORO)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('ORO_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$YRS)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('YRS_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$FOL)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('FOL_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$NML)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('NML_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$DNP)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('DNP_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$EXC)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('EXC_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

		strmnonas = as.numeric(mp$MIL)
		miss = is.na(strmnonas)
		strmnonas[miss] = 0
		strmCumSum = strmnonas
		mp[ , paste0('MIL_sum',laggards)] = c(rep(NA, laggards), strmCumSum[-c(1:laggards)] - strmCumSum[-rmv_dys])

#		mp[ , paste0("pct_",laggards)] = 		c(rep(NA, laggards), diff(mp$value, laggards) 				 / mp$value[-rmv_dys])
	}

	
			# correlation of each data
	cor_results = matrix(data = NA, nrow = maxlag, ncol = 8)
	for(laggards in 1:maxlag)	{
		cor_results[laggards,1] = cor.test(unlist(mp[ , paste0("BND_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,2] = cor.test(unlist(mp[ , paste0("ORO_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,3] = cor.test(unlist(mp[ , paste0("YRS_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,4] = cor.test(unlist(mp[ , paste0("FOL_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,5] = cor.test(unlist(mp[ , paste0("NML_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,6] = cor.test(unlist(mp[ , paste0("DNP_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,7] = cor.test(unlist(mp[ , paste0("EXC_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
		cor_results[laggards,8] = cor.test(unlist(mp[ , paste0("MIL_sum",laggards)]), unlist(mp$valueFtr), method="kendall")$p.value
	}
#	print(cor_results, round = 3)

	best_BND_lag = order(cor_results[,1])[1]
	best_ORO_lag = order(cor_results[,2])[1]
	best_YRS_lag = order(cor_results[,3])[1]
	best_FOL_lag = order(cor_results[,4])[1]
	best_NML_lag = order(cor_results[,5])[1]
	best_DNP_lag = order(cor_results[,6])[1]
	best_EXC_lag = order(cor_results[,7])[1]
	best_MIL_lag = order(cor_results[,8])[1]
	
	
	
#	mp$Day = wday(mp$Date, label=FALSE)
#	mp$Month = month(mp$Date)
	mp_subset = mp[,c("Date","value","valueFtr")]
	
	mp_subset$BND1 = mp[ , paste0("BND_sum", best_BND_lag)]
	mp_subset$ORO1 = mp[ , paste0("ORO_sum", best_ORO_lag)]
#	mp_subset$YRS1 = mp[ , paste0("YRS_sum", best_YRS_lag)]	# missing a lot of data post 2020
	mp_subset$FOL1 = mp[ , paste0("FOL_sum", best_FOL_lag)]
	mp_subset$NML1 = mp[ , paste0("NML_sum", best_NML_lag)]
	mp_subset$DNP1 = mp[ , paste0("DNP_sum", best_DNP_lag)]
	mp_subset$EXC1 = mp[ , paste0("EXC_sum", best_EXC_lag)]
	mp_subset$MIL1 = mp[ , paste0("MIL_sum", best_MIL_lag)]
	
	mp_subset$todaysPrice = mp$value
	mp_subset$month = mp$month
#	mp_subset$doy = mp$doy
	

	mp_test = subset(mp_subset, Date < as.Date("2020-01-01") & Date > as.Date("2014-01-01"))
	mp_holdout = subset(mp_subset, Date > as.Date("2020-01-01") & Date < as.Date("2021-09-01"))
	
	model_LLO = train(mp_test[,-c(1,2,3)],
						mp_test[,3],
						metric="Rsquared",
						method='rf',
						#tuneLength=?,
						imporance=TRUE,
						ntree=1000,
						trControl = trainControl(
							method='repeatedcv',
							number=5,
							repeats=3))
							#method='cv'))
		#					index = indices$index))
		#					p=.75))

	model_LLO
#	sortVars = varImp(model_LLO)$importance	;	sortVars$tot = apply(sortVars, 1, sum)
#	sortedVars = row.names(sortVars)[rev(order(sortVars$tot))]
#	sortedVars


#model_LLO_14 = model_LLO
#model_LLO_30 = model_LLO
#model_LLO_60 = model_LLO
#model_LLO_90 = model_LLO
#model_LLO = model_LLO_14
modelRMSE = model_LLO$results[which.min(model_LLO$results$RMSE),"RMSE"] * 2

	predict_val = predict(model_LLO, newdata = mp_holdout)
		# raw model results
	summary(lm(predict_val ~ mp_holdout$valueFtr))$adj.r.squared
	cor.test(mp_holdout$valueFtr, predict_val)$p.value
		# smoothed model results
	summary(lm(smooth.spline(mp_holdout$Date, predict_val, spar=0.35)$y ~ 
		mp_holdout$valueFtr))$adj.r.squared



par(cex=1.3)
	plot(mp_holdout$Date, mp_holdout$value, type='l', lwd=1.5, col='black', ylim=c(0,1000), frame.plot=FALSE)
	box(bty='l')
		# high model
#	lines(mp_holdout$Date[-c(1:predict_lag)],
#		predict_val[-c((nrow(mp_holdout)-(predict_lag-1)):nrow(mp_holdout))] + modelRMSE,
#		col=alpha('royalblue1', 0.4), lwd=2)
		# smoothed model high
	lines(smooth.spline(mp_holdout$Date[-c(1:predict_lag)],
		predict_val[-c((nrow(mp_holdout)-(predict_lag-1)):nrow(mp_holdout))]  + modelRMSE, spar=0.05),
		col=alpha('royalblue1', 0.4), lwd=2)
		# low model
#	lines(mp_holdout$Date[-c(1:predict_lag)],
#		predict_val[-c((nrow(mp_holdout)-(predict_lag-1)):nrow(mp_holdout))] - modelRMSE,
#		col=alpha('royalblue1', 0.4), lwd=2)
		# smoothed model low
	lines(smooth.spline(mp_holdout$Date[-c(1:predict_lag)],
		predict_val[-c((nrow(mp_holdout)-(predict_lag-1)):nrow(mp_holdout))]  - modelRMSE, spar=0.05),
		col=alpha('royalblue1', 0.4), lwd=2)
		# raw model
#	lines(mp_holdout$Date[-c(1:predict_lag)],
#		predict_val[-c((nrow(mp_holdout)-(predict_lag-1)):nrow(mp_holdout))],
#		col='red', lwd=2)
		# smoothed model
	lines(smooth.spline(mp_holdout$Date[-c(1:predict_lag)],
		predict_val[-c((nrow(mp_holdout)-(predict_lag-1)):nrow(mp_holdout))], spar=0.35),
		col='royalblue2', lwd=4)
		# real data again 
	lines(mp_holdout$Date, mp_holdout$value, lwd = 3)
	
	headpredict_val[-c(1:predict_lag)] *  head(mp_holdout$value,nrow(mp_holdout) - predict_lag) + mp_holdout$value[-c(1:predict_lag)], col='red', lwd=2)




	
	predictions = predict(model_LLO, newdata = mp_subset[-c(1:max_lag),])
	


	p_l = length(predictions)
	sb_l = nrow(mp_subset)
	which_col = 5
	for(j in rev(0:predict_lag))	{
		which_col = which_col + 1
		coin_results[which_coin,which_col] = (1 + predictions[p_l - j]) * mp_subset$avg[sb_l - j]
	}
	coin_results[which_coin,which_col+1] = predictions[p_l] * 100
	print(coin_results[order(coin_results[,4], decreasing=TRUE),]) 
	write.csv(coin_results[order(coin_results[,4], decreasing=TRUE),], paste0("C:\\Users\\arik\\Documents\\money_please\\coinbets\\60day",format(last(mp$Date),"%y%m%d"),".csv"))
#	}
	}
}

























			# data smoothed to a monthly mean
#	mp$avg = apply(mp[,c("open","close","high","low")],1,median)
#	mp$mnth_avg = c(rep(NA, maxlag/2 - 1), rollmean(rollmean(mp$avg,27),3), rep(NA, maxlag/2 - 1))
#	mp$pct_avg =  c(NA,(diff(mp$mnth_avg) / mp$mnth_avg[-1]))

#	acf_results = acf(subset(mp, !is.na(pct_avg))$pct_avg,
#		lag.max=maxlag, ylim=c(-.5,.5))
#	strng_acf = order(abs(acf_results$acf), decreasing=TRUE)[2:6] - 1
#	sum(abs(acf_results$acf[strng_acf+1]))

#	ccf_results = ccf(mp$vlm_pctChng, mp$pct_avg,
#		lag.max=maxlag, ylim=c(-.5,.5))
#	strng_ccf = (maxlag+1) - order(abs(ccf_results$acf[1:21]), decreasing=TRUE)[1:3]



	