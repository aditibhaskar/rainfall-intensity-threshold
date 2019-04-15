# aditi bhaskar, aditi.bhaskar@colostate.edu
# precipitation - discharge threshold analysis and uses the NEXRAD output.
# 2019 04 10 aditi is trying out iv's and this is a modification from P-Q_DRKR_NEXRAD_justbaseflow_exclude_P_lt_1mm_2019 03 18_20190405

library(EcoHydRology)
library(dplyr)
library(zoo)
library(dygraphs)
library(ggplot2)
library(ggthemes)
library(xts)

#DRKR drainage area km2 
drainage_area <- 14

RIST <- read.csv("./RIST_DRKR_rain_start_times.csv", header=TRUE, skip=5, col.names = c("datetime", "date", "precip.mm", "duration.hrs", "max_5.mmhr", "max_10.mmhr", "max_15.mmhr", "max_30.mmhr", "max_60.mmhr", "energy.mjha", "ei30.mjmmhahr"))
RIST$date <- as.Date(RIST$date, format("%m/%d/%Y"), tz="EST")
RIST$datetime <- as.POSIXct(RIST$datetime, format = "%m/%d/%Y %H:%M", tz="EST")

RIST <- subset(RIST, RIST$date >= "2005-01-01" & RIST$date <= "2009-12-31")

# ASB commenting this out on 2019-04-11 because now we can consider multiple events on the same day. 
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

# ASB modification 2019 03 11, to remove storms with <1 mm/hr. 
# RIST <- subset(RIST, RIST$max_60.mmhr >= 1)

######### DOWNLOAD IV'S FOR DRKR
library(hydrostats)

# read from dataRetrieval streamflow - only do this in the beginning, otherwise read it in. 
# iv.streamflow <- readNWISdata(sites="01589330", service="iv", parameterCd="00060",startDate="2005-01-01T00:00", endDate="2009-12-13T23:59", tz="America/Jamaica")
# write.csv(iv.streamflow, "./DRKR_Q_iv_2005-2009.csv")
iv.streamflow <- read.csv("./DRKR_Q_iv_2005-2009.csv")

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

# Create dataframe with timesteps every 5-min
# Join with the flow data
Q_Com <- data.frame(seq.POSIXt(startDateTime, endDateTime, by="5 mins"))
colnames(Q_Com) <- "dateTime"
final <- left_join(Q_Com, Q_Base, by="dateTime")

rm(Q, Q_Base, Q_Com, iv.streamflow)

##############################
# Calculate Rolling minimums #
##############################
# Values are the number of time intervals
# I'm using 5-min intervals
final$Qquick_leadrollmin_24hr <- rollapply(final$Qquick_cfs, 288, fill = NA, align = "left", min) # lead, falling limb 
# 288 is how many 5 minute intervals there are in 24 hours, 144 is how many 5 minute intervals there are in 12 hours. 
final$Qquick_lagrollmin_12hr <- rollapply(final$Qquick_cfs, 144, fill = NA, align = "right", min) # lag rising limb
final$Q_leadrollmin_12hr <- rollapply(final$Q_cfs, 144, fill = NA, align = "left", min) # lead, falling limb

# Calculate discharge minus rolling mins
final <- final %>% 
  mutate(Qquick_minLead24hr = Qquick_cfs - Qquick_leadrollmin_24hr) %>%  # 24 hour window
  mutate(Qquick_minlag12hr = Qquick_cfs - Qquick_lagrollmin_12hr) %>% 
  mutate(Q_cfs_minLead12hr = Q_cfs - Q_leadrollmin_12hr) # 12 hours window


#########################
# Create indicator code #
#########################
# see notes word document from 2019 04 10 with Aditi's notes as to why these values are being changed. 
final$EventInd <- 0
# final$EventInd[final$Qquick_cfs > 1.3] <- 1
# final$EventInd[final$Q_cfs > 4] <- 1
# final$EventInd[final$Qquick_minLead24hr > 10] <- 1
# final$EventInd[final$Qquick_minlag12hr > 10] <- 1

# based on drainage areas: 
final$EventInd[final$Qquick_cfs > 0.2*drainage_area] <- 1
final$EventInd[final$Q_cfs > 1.5*drainage_area] <- 1
final$EventInd[final$Qquick_minLead24hr > 0.4*drainage_area] <- 1
final$EventInd[final$Qquick_minlag12hr > 0.4*drainage_area] <- 1

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
         filter(dateTime >= as.POSIXct("2006-07-25") & 
                  dateTime <= as.POSIXct("2006-07-30")), 
       aes(dateTime, Q_cfs, color = factor(EventInd))) + 
  geom_point() + 
  scale_colour_brewer(palette = "Set1")+
  ylim(0,10)+
  theme_few()+
  labs(colour = "Variable", y = "Q (cfs)", 
       title = "DRKR")+
  geom_hline(yintercept = 0.25, color = "black", size=1)

write.csv(final, paste0("./DRKR_Processed_", format(Sys.time(), "%d-%b-%Y %H.%M"), ".csv"), row.names = FALSE)

####### DYGRAPHS COMPARISON PLOTS FOR IV #############

not.baseflow <- subset(final, final$EventInd == 1)
xts.flow <- xts(final, order.by=final$dateTime)
xts.not.baseflow <- xts(not.baseflow, order.by=not.baseflow$dateTime)
combined.flow <- merge.xts(xts.flow$Q_cfs, xts.not.baseflow$Q_cfs)
xts.RIST <- xts(RIST, order.by=RIST$date)
dygraph(cbind(combined.flow, xts.RIST$max_60.mmhr), main = "DRKR") %>% 
  dyAxis("y", label = "Flow (cfs)") %>%
  dyAxis("y2", label = "MI60 (mm/hr)", independentTicks = TRUE) %>%
  dyRangeSelector()  %>%
  dyOptions(drawPoints = TRUE, pointSize = 2) %>%
  dySeries("max_60.mmhr", axis = 'y2')

###########################  END PASTING FROM KRISSY HOPKINS' CODE
int <- interval(ymd("2001-01-01", tz="EST"), ymd("2002-01-01", tz="EST"))
int_start(int)
int_start(int) <- RIST$datetime
# update on 9/19/2018 Stephanie sent me a file with the actual dates.  therefore, allow rainfall and flow event matching if the flow event occurs anytime within the duration of the rainfall event. 
#https://stackoverflow.com/questions/11922181/adding-time-to-posixct-object-in-r
#POSIXct uses seconds.  so go from hours to 60 minutes/hour to 60 seconds/min
# RIST$date <- as.Date(RIST$datetime + (RIST$duration.hrs*60*60)/2)
int_end(int) <- RIST$datetime + RIST$duration.hrs*60*60

precip.events <- int

RIST$not.baseflow.cutoff <- FALSE

for (precip.events.index in 1:length(precip.events)){
  RIST$not.baseflow.cutoff[precip.events.index] <- any(not.baseflow$dateTime %within% precip.events[precip.events.index])
}


png(filename=paste0("./for_paper_DRKR_iv_", format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"), width=7, height=4, units="in", res=300)
final.forplot <- subset(final, final$dateTime > as.Date("2009-09-01"))
final.forplot <- subset(final.forplot, final.forplot$dateTime < as.Date("2009-12-01"))
not.baseflow.forplot <- subset(not.baseflow, not.baseflow$dateTime > as.Date("2009-09-01"))
not.baseflow.forplot <- subset(not.baseflow.forplot, not.baseflow.forplot$dateTime < as.Date("2009-12-01"))
par(mar=c(3, 6, 3, 3),xpd=TRUE)
plot(final.forplot$dateTime, final.forplot$Q_cfs*0.0283168, 'l', log="y", xlab="", ylab=expression(paste("Stream Discharge (m"^"3","/s)")), col="blue", ylim=c(0.02,50))
points(final.forplot$dateTime, final.forplot$Q_cfs*0.0283168, cex=0.3)
points(not.baseflow.forplot$dateTime, not.baseflow.forplot$Q_cfs*0.0283168, col="red", cex=0.3)
#legend('topleft', c("Flow event", "Baseflow", "Daily Discharge"), col=c("red", "black", "blue"), pch=c(1,1,NA), lty=c(NA, NA, 1), horiz=TRUE)
title('b. Urban', adj=0)
# plot(RIST$date, RIST$precip.mm, 'h')
dev.off()

png(filename=paste0("./05-09_DRKR_MI60not.baseflow.cutoff", format(Sys.time(), "%d-%b-%Y %H.%M"), ".png"))#, width=4, height=3, units="in", res=300)
thresholds <- data.frame(tested.value = seq(from=0, to=50, by=0.1), prediction = 0)
for (threshold.index in thresholds$tested.value) {
  thresholds$prediction[thresholds$tested.value==threshold.index] <- sum((RIST$max_60.mmhr > threshold.index) == RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff)
}
threshold <- thresholds$tested.value[which.max(thresholds$prediction)]
p0 <- max(thresholds$prediction) 
pobs_flow <- sum(RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff)
pt_flow <- sum(RIST$max_60.mmhr>threshold)/length(RIST$not.baseflow.cutoff) 
pobs_noflow <- (1-sum(RIST$not.baseflow.cutoff)/length(RIST$not.baseflow.cutoff))
pt_noflow <- sum(RIST$max_60.mmhr<= threshold)/length(RIST$not.baseflow.cutoff)
pe <- pt_flow*pobs_flow + pt_noflow*pobs_noflow
kappa <- (p0-pe)/(1-pe)
plot(subset(RIST, !RIST$not.baseflow.cutoff)$max_60.mmhr, subset(RIST, !RIST$not.baseflow.cutoff)$not.baseflow.cutoff, pch=4, col="red", xlim=c(0, max(RIST$max_60.mmhr)), ylim=c(0,1), yaxt="n", xlab="RIST MI60 mm/hr", ylab="not baseflow based on separation"
     , main=paste0("DRKR threshold of ", round(thresholds$tested.value[which.max(thresholds$prediction)],2), " mm/hr. kappa=",round(kappa,2)))
axis(side=2, at =c(0, 0.25, 0.5, 0.75, 1))
points(subset(RIST, RIST$not.baseflow.cutoff)$max_60.mmhr, subset(RIST, RIST$not.baseflow.cutoff)$not.baseflow.cutoff, pch=4, col="black")
points(thresholds$tested.value, thresholds$prediction, col="green", pch=1)
legend(10,0.4, c("storm responses", "not storm responses","# events correctly predicted at this MI60/total # events"), pch=c(4,4,1), col=c("black","red","green"))
dev.off()

write.csv(RIST, file=paste0("./05-09_DRKR_iv",format(Sys.time(), "%d-%b-%Y %H.%M"), ".csv"))