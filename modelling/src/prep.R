



##-----------------------------------------------------------------
## Read data
dataFull = read.csv2('./data/high_tc_20190913.csv', fileEncoding="latin1")
N = dim(dataFull)[1]
set.seed(144257)


# Subsampling
nSkip = 100
Idx = seq(1,N,nSkip)
data = dataFull #dataFull[Idx,]
n = dim(data)[1]


# Define provided status changes
changeTimes <- c(as.POSIXct("2019-09-17 13:00:00"),
                 as.POSIXct("2019-09-17 13:30:00"),
                 as.POSIXct("2019-09-17 14:00:00"),
                 as.POSIXct("2019-09-17 14:45:00"),
                 as.POSIXct("2019-09-17 15:30:00"),
                 as.POSIXct("2019-09-17 16:30:00"),
                 as.POSIXct("2019-09-17 17:30:00"),
                 as.POSIXct("2019-09-17 19:00:00"))


# Data
time = as.POSIXct(data$Time)
hour_index = which(c(1,diff(as.numeric(strftime(time, format="%H"))))==1)

# Convert to step function
status = createMeltSteps(changeTimes,time)
ice_data_full = data.frame("status" = status, "CompCap" = data$A40_A_CompCap/100)
ice_data_full$t = seq(0,n-1)*6/60

ice_data = ice_data_full[Idx,]