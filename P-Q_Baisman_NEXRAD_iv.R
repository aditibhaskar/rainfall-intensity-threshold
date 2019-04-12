# aditi bhaskar, aditi.bhaskar@colostate.edu
# 6 july 2017
# precipitation - discharge threshold analysis and uses the NEXRAD output.

# 2019 03 18 aditi is making modifications looking at what stephanie sent for code. 
# 2019 04 10 aditi is trying out iv's and this is a modification from P-Q_DRKR_NEXRAD_justbaseflow_exclude_P_lt_1mm_2019 03 18_20190405


library(EcoHydRology)
library(dplyr)
library(zoo)
library(dygraphs)
library(ggplot2)
library(ggthemes)
library(xts)

# import RIST file. 
# RIST <- read.csv("./DRKR NEX_RIST_2005-2009.csv", sep="", header=TRUE, skip=5, col.names = c("date", "precip.mm", "duration.hrs", "max_5.mmhr", "max_10.mmhr", "max_15.mmhr", "max_30.mmhr", "max_60.mmhr", "energy.mjha", "ei30.mjmmhahr"))
# RIST$date <- as.Date(RIST$date, format("%m/%d/%Y"), tz="EST")

RIST <- read.csv("G:/My Drive/PROJECTS/2017/Kampf Green P threshold analysis/analysis files/Baltimore/Baisman/RIST_Baisman_rain_start_times.csv", header=TRUE, skip=5, col.names = c("datetime", "date", "precip.mm", "duration.hrs", "max_5.mmhr", "max_10.mmhr", "max_15.mmhr", "max_30.mmhr", "max_60.mmhr", "energy.mjha", "ei30.mjmmhahr"))
RIST$date <- as.Date(RIST$date, format("%m/%d/%Y"), tz="EST")
RIST$datetime <- as.POSIXct(RIST$datetime, format = "%m/%d/%Y %H:%M", tz="EST")
RIST$precip.mm <- as.numeric(RIST$precip.mm)
RIST$duration.hrs <- as.numeric(RIST$duration.hrs)
RIST$max_60.mmhr <- as.numeric(RIST$max_60.mmhr)
RIST$max_30.mmhr <- as.numeric(RIST$max_30.mmhr)
RIST$max_15.mmhr <- as.numeric(RIST$max_15.mmhr)


# suggestion from Stephanie on 9/19/2018: How about trying to adjust RIST dates based on duration of event? We could use the RIST output to compute event start, event end, then assign the event to the day with the majority of time?
# so I am going to work with an assumption (arbitrary) that the event starts at 12 noon on the RIST date in question, because RIST itself doesn't give a time. so therefore, if the duration is > 12 hours, I move the RIST$date to the following day.  If the RIST$date > 24+12 or 36 hours, I move the RIST$date + 2 days, etc..  
# RIST$old.date <- RIST$date
# RIST$date <- RIST$date + ceiling((RIST$duration.hrs-12)/24)
# update on 9/19/2018 Stephanie sent me a file with the actual dates.  
#https://stackoverflow.com/questions/11922181/adding-time-to-posixct-object-in-r
#POSIXct uses seconds.  so go from hours to 60 minutes/hour to 60 seconds/min
# RIST$date <- as.Date(RIST$datetime + (RIST$duration.hrs*60*60)/2)

# SO ASB is realizing on 2019 04 05 that when I got the RIST date/times of the storms from Stephanie - I never seem to actually have used them!?!?! 
# RIST$date <- as.Date(RIST$datetime + (RIST$duration.hrs*60*60)/2)

RIST <- subset(RIST, RIST$date >= "2005-01-01" & RIST$date <= "2009-12-31")


# I DONT THINK i NEED TO DO THIS ANYMORE NOW THAT i HAVE THE DATE/TIME OF THE START OF THE RAIN EVENT. 
# # subset RIST to remove rows that are on the same date (two events on the same date, the one that will be removed has lower MI60)
# ordered.RIST <- RIST[order(RIST$max_60.mmhr, decreasing = TRUE),]  # order RIST with highest MI60 up top
# duplicated.return <- duplicated(ordered.RIST$date) # find duplicate dates with second date indicated as TRUE, with second date also having lower MI60
# RIST.not.duplicated <- ordered.RIST[!duplicated.return,] # remove duplicated dates
# RIST <- RIST.not.duplicated[order(RIST.not.duplicated$date),] #return in original date order

## remove storm events that were not in April - October to avoid snow time periods.  do this by deleting these lines in rist. 
library("lubridate")
RIST <- subset(RIST,month(RIST$date) > 3 & month(RIST$date) < 11)

library("dataRetrieval") # USGS-R package to get NWIS data
startDateTime <- as.POSIXct("2005-01-01 00:00:00 EST", tz="EST")
endDateTime <- as.POSIXct("2009-12-31 23:59:59 EST", tz="EST")

######### DOWNLOAD IV'S FOR DRKR
library(hydrostats)

# read from dataRetrieval streamflow - only do this in the beginning, otherwise read it in. 
# iv.streamflow <- readNWISdata(sites="01583580", service="iv", parameterCd="00060",startDate="2005-01-01T00:00", endDate="2009-12-13T23:59", tz="America/Jamaica")
# write.csv(iv.streamflow, "./Baisman_Q_iv_2005-2009.csv")
iv.streamflow <- read.csv("./Baisman_Q_iv_2005-2009.csv")

######## COPIED FROM KRISSY HOPKINS' CODE 3_IDENTIFY_HIGH_FLOW_EVENTS.R: 

final <- iv.streamflow
final$dateTime <- as.POSIXct(strptime(final$dateTime, "%Y-%m-%d %H:%M:%S"), tz = "America/Jamaica")
final <- data.frame(final$dateTime, final$X_00060_00000)
colnames(final) <- c("Date", "Q_cfs")
final <- final[order(as.POSIXct(final$Date, format="%Y-%m-%d %H:%M:%S")),]

# Run baseflow separation filter
# Must not have an NAs for streamflow (filter out NAs with na.omit())
# Use 0.99 for the filter parameter 
Q = final %>% na.omit()

Q_Base <- BaseflowSeparation(Q$Q_cfs, filter_parameter = 0.99, passes = 3)
Q_Base <- cbind(Q, Q_Base)  
colnames(Q_Base) <- c("dateTime", "Q_cfs", "Qbase_cfs", "Qquick_cfs")

# Create dataframe with timesteps every 15-min
# Join with the flow data
Q_Com <- data.frame(seq.POSIXt(startDateTime, endDateTime, by="15 mins"))
colnames(Q_Com) <- "dateTime"
final <- left_join(Q_Com, Q_Base, by="dateTime")

rm(Q, Q_Base, Q_Com, iv.streamflow)

##############################
# Calculate Rolling minimums #
##############################
# Values are the number of time intervals
# I'm using 5-min intervals - ASB MODIFIED TO 15 MINUTE INTERVALS
final$Qquick_leadrollmin_24hr <- rollapply(final$Qquick_cfs, 96, fill = NA, align = "left", min) # lead, falling limb 
# 288 is how many 5 minute intervals there are in 24 hours, 144 is how many 5 minute intervals there are in 12 hours. 
final$Qquick_lagrollmin_12hr <- rollapply(final$Qquick_cfs, 48, fill = NA, align = "right", min) # lag rising limb
final$Q_leadrollmin_12hr <- rollapply(final$Q_cfs, 48, fill = NA, align = "left", min) # lead, falling limb

# Calculate discharge minus rolling mins
final <- final %>% 
  mutate(Qquick_minLead24hr = Qquick_cfs - Qquick_leadrollmin_24hr) %>%  # 24 hour window
  mutate(Qquick_minlag12hr = Qquick_cfs - Qquick_lagrollmin_12hr) %>% 
  mutate(Q_cfs_minLead12hr = Q_cfs - Q_leadrollmin_12hr) # 12 hours window


#########################
# Create indicator code #
#########################
# # Use for TR104 and TR109
# final$EventInd <- 0
# final$EventInd[final$Qquick_cfs > 0.25] <- 1
# final$EventInd[final$Q_cfs > 2] <- 1
# final$EventInd[final$Qquick_minLead24hr > 0.2] <- 1
# final$EventInd[final$Qquick_minlag12hr > 0.2] <- 1

# For Soper and CR
# see notes word document from 2019 04 10 with Aditi's notes as to why these values are being changed. 
final$EventInd <- 0
# final$EventInd[final$Qquick_cfs > 0.25] <- 1
final$EventInd[final$Qquick_cfs > 0.6] <- 1
# final$EventInd[final$Q_cfs > 3] <- 1
final$EventInd[final$Q_cfs > 5] <- 1
# final$EventInd[final$Qquick_minLead24hr > 0.2] <- 1
final$EventInd[final$Qquick_minLead24hr > 1] <- 1
# final$EventInd[final$Qquick_minlag12hr > 0.2] <- 1
final$EventInd[final$Qquick_minlag12hr > 0.5] <- 1

# RISING LIMBS ONLY - this Aditi added, not in Krissy's IDs
# final$diff <- 0 # set first row difference from previous to be 0
# final$diff[2:nrow(final)] <- diff(final$Q_cfs) # all other rows, take difference in Q from previous row's Q
# # if the flow was smaller than the one on the previous timestep, make it not an event. 
# final$EventInd[final$diff < 0] <- 0 
# final$EventInd[final$Qquick_minLead24hr == final$Qquick_cfs] <- 0
# THE ABOVE isn't working right.  entire storms are being missed. 

####################
# check indicators #
####################
# final_Check = final %>% filter(final$Qquick_minLead24hr > 0.2)
# sum(final_Check$EventInd)

# Plot indicator
ggplot(final  %>% 
         filter(dateTime >= as.POSIXct("2007-05-01") & 
                  dateTime <= as.POSIXct("2007-05-20")), 
       aes(dateTime, Q_cfs, color = factor(EventInd))) + 
  geom_point() + 
  scale_colour_brewer(palette = "Set1")+
  ylim(0,10)+
  theme_few()+
  labs(colour = "Variable", y = "Q (cfs)", 
       title = "Baisman")+
  geom_hline(yintercept = 0.25, color = "black", size=1)


write.csv(final, paste0("./Baisman_Processed_", format(Sys.time(), "%d-%b-%Y %H.%M"), ".csv"), row.names = FALSE)

####### DYGRAPHS COMPARISON PLOTS FOR IV #############

not.baseflow <- subset(final, final$EventInd == 1)
xts.flow <- xts(final, order.by=final$dateTime)
xts.not.baseflow <- xts(not.baseflow, order.by=not.baseflow$dateTime)
combined.flow <- merge.xts(xts.flow$Q_cfs, xts.not.baseflow$Q_cfs)
xts.RIST <- xts(RIST, order.by=RIST$date)
dygraph(cbind(combined.flow, xts.RIST$max_60.mmhr), main = "Baisman") %>% 
  dyAxis("y", label = "Flow (cfs)") %>%
  dyAxis("y2", label = "MI60 (mm/hr)", independentTicks = TRUE) %>%
  dyRangeSelector()  %>%
  dyOptions(drawPoints = TRUE, pointSize = 2) %>%
  dySeries("max_60.mmhr", axis = 'y2')

###########################  END PASTING FROM KRISSY HOPKINS' CODE
int <- interval(ymd("2001-01-01", tz="EST"), ymd("2002-01-01", tz="EST"))
int_start(int)
int_start(int) <- RIST$datetime
int_end(int) <- RIST$datetime + RIST$duration.hrs*60*60

precip.events <- int

RIST$not.baseflow.cutoff <- FALSE

for (precip.events.index in 1:length(precip.events)){
  RIST$not.baseflow.cutoff[precip.events.index] <- any(not.baseflow$dateTime %within% precip.events[precip.events.index])
}

# ASB modifying on 2019 04 05 to also create a storm end date matching
# RIST$not.baseflow.cutoff.start <- as.Date(as.character(RIST$date), tz="EST") %in% as.Date(as.character(not.baseflow$Date), tz="EST")
# RIST$date.end <- as.Date(RIST$datetime + (RIST$duration.hrs*60*60))
# RIST$not.baseflow.cutoff.end <- as.Date(as.character(RIST$date.end), tz="EST") %in% as.Date(as.character(not.baseflow$Date), tz="EST")
# 
# RIST$not.baseflow.cutoff.end.2 <- as.Date(as.character(RIST$date.end-1), tz="EST") %in% as.Date(as.character(not.baseflow$Date), tz="EST")
# 
# RIST$not.baseflow.cutoff.end.3 <- as.Date(as.character(RIST$date.end-2), tz="EST") %in% as.Date(as.character(not.baseflow$Date), tz="EST")
# 
# # plot(daily.streamflow$dateTime, daily.streamflow$Flow, "l", log="y")
# RIST$not.baseflow.cutoff <- RIST$not.baseflow.cutoff.end | RIST$not.baseflow.cutoff.start
# 
# RIST$not.baseflow.cutoff[RIST$date.end - RIST$date == 2] <- RIST$not.baseflow.cutoff[RIST$date.end - RIST$date == 2] | RIST$not.baseflow.cutoff.end.2[RIST$date.end - RIST$date == 2]

# for (i in nrow(RIST)){
#   if (RIST$date.end[i] - RIST$date[i] >= 2){
#     RIST$not.baseflow.cutoff[i] <- RIST$not.baseflow.cutoff.end[i] | RIST$not.baseflow.cutoff.start[i] | RIST.not.baseflow.cutoff.end.2[i] 
#   }
#   if (RIST$date.end[i] - RIST$date[i] >= 3){
#     RIST$not.baseflow.cutoff[i] <- RIST$not.baseflow.cutoff.end[i] | RIST$not.baseflow.cutoff.start[i] | RIST.not.baseflow.cutoff.end.2[i] | RIST.not.baseflow.cutoff.end.3[i]
#   }
# # }
# 
# 
# plot(subset(daily.streamflow$dateTime, daily.streamflow$dateTime > as.Date("2009-06-01")), subset(daily.streamflow$Flow, daily.streamflow$dateTime > as.Date("2009-06-01")), 'l', log="y", xlab="", ylab=expression(paste("Daily Stream Discharge (m"^"3","/s)")))
# points(daily.streamflow$dateTime, daily.streamflow$Flow)
# points(not.baseflow$Date,not.baseflow$Q, col="red")
# 
# 
# 
# # non dygraphs
# par(mfrow=c(2,1))
# plot(final$dateTime, final$Q_cfs, 'l', log="y")
# points(final$dateTime, final$Q_cfs)
# points(not.baseflow$dateTime,not.baseflow$Q_cfs, col="blue", pch=20)
# legend('bottomleft', c("Flow event (not baseflow)", "Baseflow"), col=c("blue", "black"), pch=c(20,1))
# par(new=TRUE)
# plot(RIST$date, RIST$max_60.mmhr, pch=4, xlim=range(as.Date(daily.streamflow$dateTime)))
# axis(side = 4, ylab="RIST$MI60 mm/hr")
# points(RIST[RIST$not.baseflow.cutoff,]$date, RIST[RIST$not.baseflow.cutoff,]$max_60.mmhr, pch=4,col="red")
# png(filename=paste0("./for_paper_05-09_baseflow.separation.alpha",alpha,"_cfsabove",cfs.above.baseflow,".png"), width=7, height=4, units="in", res=300)
# par(mar=c(3, 6, 3, 3),xpd=TRUE)
# plot(subset(daily.streamflow$dateTime, daily.streamflow$dateTime > as.Date("2009-09-01")), subset(daily.streamflow$Flow*0.0283168, daily.streamflow$dateTime > as.Date("2009-09-01")), 'l', log="y", xlab="", ylab=expression(paste("Daily Stream Discharge (m"^"3","/s)")), col="blue", ylim=c(0.02,6.5))
# points(subset(daily.streamflow$dateTime, daily.streamflow$dateTime > as.Date("2009-09-01")), subset(daily.streamflow$Flow*0.0283168, daily.streamflow$dateTime > as.Date("2009-09-01")))
# points(subset(not.baseflow$Date, not.baseflow$Date > "2009-09-01"),subset(not.baseflow$Q*0.0283168, not.baseflow$Date > "2009-09-01"), col="red")
# # legend('topleft', c("Flow event", "Baseflow", "Daily Discharge"), col=c("red", "black", "blue"), pch=c(1,1,NA), lty=c(NA, NA, 1), horiz=TRUE)
# title('b. Urban', adj=0)
# # plot(RIST$date, RIST$precip.mm, 'h')
# dev.off()


# ASB modification 2019 03 11, to remove storms with <1 mm/hr. 
# RIST <- subset(RIST, RIST$max_60.mmhr > 1) # ASB modification 2019 04 05: > instead of >= based on Stephanie's email from April 3, 2019
# MI60over1=filter(RIST,RIST$max_60.mmhr>1)
# other notes from Stephanie in the same email: 
# #here need to calculate F based on all predictions, not just where flow is correctly predicted. in my code this would look like
# #if predict = flow, correct = 1, else correct = 0
# #p0 = sum(correct)/length(correct)
# #but I am not sure how to update your code to do the same
# for (threshold.index in thresholds$tested.value) {
#   thresholds$prediction[thresholds$tested.value==threshold.index] <- 
#     sum((MI60over1$max_60.mmhr > threshold.index) == 
#           MI60over1$not.baseflow.cutoff)/length(MI60over1$not.baseflow.cutoff)
# }

## simplified plot
png(filename=paste0("./05-09_MI60not.baseflow.cutoff", format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"))#, width=4, height=3, units="in", res=300)
# thresholds <- data.frame(tested.value = seq(0, max(RIST$max_60.mmhr), length.out=4000), prediction = 0)
thresholds <- data.frame(tested.value = seq(from=0, to=50, by=0.1), prediction = 0) # modifying to use what stephanie is doing. 
for (threshold.index in thresholds$tested.value) {
  thresholds$prediction[thresholds$tested.value==threshold.index] <- sum((RIST$max_60.mmhr > threshold.index) == RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff)
}
threshold <- thresholds$tested.value[which.max(thresholds$prediction)]
p0 <- max(thresholds$prediction) # this part matches fine

# it is the pe that isn't matching.\
#pe <- sum(RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff)*sum(RIST$max_60.mmhr>threshold)/length(RIST$not.baseflow.cutoff)  + (1-sum(RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff))*sum(RIST$max_60.mmhr<= threshold)/length(RIST$not.baseflow.cutoff)

pobs_flow <- sum(RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff) # matches

# THIS IS THE PROBLEM CALCULATION. 
pt_flow <- sum(RIST$max_60.mmhr>threshold)/length(RIST$not.baseflow.cutoff) 

pobs_noflow <- (1-sum(RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff))
pt_noflow <- sum(RIST$max_60.mmhr<= threshold)/length(RIST$not.baseflow.cutoff)

pe <- pt_flow*pobs_flow + pt_noflow*pobs_noflow

kappa <- (p0-pe)/(1-pe)
plot(subset(RIST, !RIST$not.baseflow.cutoff)$max_60.mmhr, subset(RIST, !RIST$not.baseflow.cutoff)$not.baseflow.cutoff, pch=4, col="red", xlim=c(0, max(RIST$max_60.mmhr)), ylim=c(0,1), yaxt="n", xlab="RIST MI60 mm/hr", ylab="not baseflow based on separation"
     , main=paste0("Baisman threshold of ", round(thresholds$tested.value[which.max(thresholds$prediction)],2), " mm/hr. kappa=",round(kappa,2)))
axis(side=2, at =c(0, 0.25, 0.5, 0.75, 1))
points(subset(RIST, RIST$not.baseflow.cutoff)$max_60.mmhr, subset(RIST, RIST$not.baseflow.cutoff)$not.baseflow.cutoff, pch=4, col="black")
############## THRESHOLDS value based on procedure from Stephanie's email on 7/14/17
# test out threshold values in a iteration loop
points(thresholds$tested.value, thresholds$prediction, col="green", pch=1)
legend(10,0.4, c("storm responses", "not storm responses","# events correctly predicted at this MI60/total # events"), pch=c(4,4,1), col=c("black","red","green"))
# text(10, 1.1, paste0(round(thresholds$tested.value[which.max(thresholds$prediction)],2), "=threshold"))
dev.off()


write.csv(RIST, file=paste0("./05-09_Baisman_iv",format(Sys.time(), "%d-%b-%Y %H.%M"), ".csv"))

### COMPARISON PLOTS$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# 
# 
# xts.flow <- xts(daily.streamflow, order.by=daily.streamflow$dateTime)
# xts.baseflow <- xts(not.baseflow, order.by=not.baseflow$Date)
# combined.flow <- merge.xts(xts.flow$Flow, xts.baseflow$Q)
# xts.RIST <- xts(RIST, order.by=RIST$date)
# rainfall <- merge.xts(xts.RIST$precip.mm, xts.RIST$max_60.mmhr)
# dygraph(cbind(combined.flow, rainfall$precip.mm), main = "Baisman") %>% 
#   dyAxis("y", label = "Flow (cfs)") %>%
#   dyAxis("y2", label = "Rainfall (mm)", independentTicks = TRUE) %>%
#   dyRangeSelector()  %>%
#   dyOptions(drawPoints = TRUE, pointSize = 2) %>%
#   dySeries("precip.mm", axis = 'y2')
# dygraph(cbind(combined.flow, rainfall$max_60.mmhr), main = "Baisman") %>% 
#   dyAxis("y", label = "Flow (cfs)") %>%
#   dyAxis("y2", label = "MI60 (mm/hr)", independentTicks = TRUE) %>%
#   dyRangeSelector()  %>%
#   dyOptions(drawPoints = TRUE, pointSize = 2) %>%
#   dySeries("max_60.mmhr", axis = 'y2')
# 
# # non dygraphs
# # par(mfrow=c(2,1))
# plot(daily.streamflow$dateTime, daily.streamflow$Flow, 'l', log="y")
# points(daily.streamflow$dateTime, daily.streamflow$Flow)
# points(not.baseflow$Date,not.baseflow$Q, col="blue", pch=20)
# legend('bottomleft', c("Flow event (not baseflow)", "Baseflow"), col=c("blue", "black"), pch=c(20,1))
# par(new=TRUE)
# plot(RIST$date, RIST$max_60.mmhr, pch=4, xlim=range(as.Date(daily.streamflow$dateTime)))
# axis(side = 4, ylab="RIST$MI60 mm/hr")
# points(RIST[RIST$not.baseflow.cutoff,]$date, RIST[RIST$not.baseflow.cutoff,]$max_60.mmhr, pch=4,col="red")
# 
