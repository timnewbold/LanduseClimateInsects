# LanduseClimateInsects

There are 5 scripts that should be run in the following order:

1. convenience_functions.r
  This contains required functions for the rest of the scripts.

2.climate_percNH.r
  This will obtain climate and percNH data for all sites in the predicts database

3.models.r
  Create models with the data

4.plot.r
  Create plots of model results

5.maps.r 
  Create maps of standardised climate anomaly in 1970, 2005, 2018
  This script can be run at anytime.
 

Descriptions of variables created by the climate_percNH script + convience_functions are as follows:

LandUse is Landcover excluding urban sites

log_abundance is the log(x+1) transformation of Total_abundance

Scaled_abundance scales abundance within studies in predicts, giving the maximum abundance in a study a value of 1

and all other abundances with that study relative to that value

log_scaled_abundance is the log(x+0.01) transformation of Scaled_abundance

avg_temp is the Mean annual temp of the site in year preceding the end sample date

tmax is the mean maximum temperature of the warmest month in the year preceding the end sample date

avg_temp_sd is the standard deviation of mean monthly temperatures in the year preceding the end sample date

historic_sd is the standard deviation of mean monthly temperatures between 1901 and 1905

climate_anomaly is avg_temp - the mean temperature for all months between 1901 and 1905

tmax_anomaly is tmax - mean of the hottest month of each year from 1901 to 1905

percNH is the percentage natural habitat in the by 5km grid cell in which the site is found

percNH.s is the rescaled and centered version of percNH

std_climate anomaly is climate_anomaly / historic_sd

std_tmax_anomaly is tmax_anomaly / historic_sd
