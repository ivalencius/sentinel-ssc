library(dataRetrieval)
library(dplyr)
library(USAboundaries)
library(data.table)

us_states_abbr <-state.abb

ssc_code <- c("80154", "00076") # SSC, turbidity
#wqp_stations <- data.table(whatWQPdata(parameterCd = ssc_code))
#                                           [,":=" (
#                                             agency_cd = "USGS",
#                                             site_no = unlist(strsplit(MonitoringLocationIdentifier, split='-'))[c(FALSE, TRUE)],
#                                             sample_dt = ActivityStartDate,
#                                             SSC_mgL = ResultMeasureValue,
#                                             sample_depth_m = ActivityDepthHeightMeasure.MeasureValue * 0.3048,
#                                             sample_method = SampleCollectionMethod.MethodName
#                                             )][,.(agency_cd, site_no, sample_dt, SSC_mgL, sample_depth_m, sample_method)]

# Extract all NWIS stations
nwis_stations <- data.frame(matrix(ncol=4,nrow=0, 
    dimnames=list(NULL, c("site_no", "begin_date", "end_date", "count_nu"))))

for (state in us_states_abbr){
  print(state)
  nwis_data <- data.table(whatNWISdata(stateCd=state, parameterCd = ssc_code))
  
  if(nrow(nwis_data)==0){} # Do nothing if no data returned
  else{
    nwis_data_clean <- nwis_data[,.(site_no, begin_date, end_date, count_nu)]
    # If station has both, combine across site_no
    # begin_date = min, end_date = max, count_nu = sum
    nwis_data_clean <- nwis_data_clean[, .(
      begin_date = min(begin_date),
      end_date   = max(end_date),
      count_nu   = sum(count_nu)
    ), by = site_no]
    nwis_stations <- rbind(nwis_stations, nwis_data_clean)
  }
}

nwis_stations$start_year <- as.numeric(format(as.Date(
    nwis_stations$begin_date, format="%Y-%m-%d"),"%Y"))
nwis_stations$end_year <- as.numeric(format(as.Date(
    nwis_stations$end_date, format="%Y-%m-%d"),"%Y"))
stations_year <- nwis_stations %>% mutate(
    years_of_operation = end_year - start_year + 1
    # average_samples = count_nu/years_of_operation
)


# Get years of operation
monitoring <- data.frame(
    year = min(nwis_stations$start_year):2023,
    n_stations = rep(0, each=length(year)),
    n_samples = rep(0, each=length(year))
)

for(i in seq_along(monitoring[,1])){
    year <- monitoring$year[i]
    stations_i <- stations_year %>% filter((start_year<=year) & (end_year>=year))
    monitoring$n_stations[i] <- nrow(stations_i)
    monitoring$n_samples[i] <- sum(stations_i$count_nu)
}

write.csv(monitoring, "figure-data/stations-over-time.csv", row.names=FALSE)
