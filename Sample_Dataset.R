# This is a sample code for sample dataset
# Annas Hummingbird - Temperate
# Tufted Coquet - Tropical
# Steps:
# 1. Pull GBIF Data
# 2. Get biogeographic realm shapefile
# 3. Load in climate data and clip it by shapefile
# 4. Generate available points within shapefile
# 5. Run maxent
# 6. 

## Bioclim variables: BIO1=Ann Mean Temp; BIO2=Mean Diurnal Range (Mean of monthly 
# (max temp - min temp)); BIO3=Isothermality (BIO2/BIO7) (* 100); BIO4=Temp Seasonality 
# (SD *100);BIO5 = Max Temp Warmest Mon;BIO6=Min Temp Coldest Mon; BIO7=Temp Ann Range;
# BIO8=Mean Temp Wettest Q;BIO9=Mean Temp Driest Q; BIO10=Mean Temp Warmest Q;
# BIO11=Mean Temp Coldest Q; BIO12 = Annual Precip; BIO13=Precip Wettest Month;
# BIO14=Precip Driest Month; BIO15=Precip Seasonality (CV); BIO16=Precip Wettest Q;
# BIO17=Precip Driest Q; BIO18=Precip Warmest Q; BIO19 = Precip Coldest Q

