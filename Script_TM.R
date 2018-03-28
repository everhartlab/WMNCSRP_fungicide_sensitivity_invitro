#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
TMdata <- read.csv("data/TM_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
TMdata <- subset(TMdata, select = - ( X: X.6 ))
TMdata$repeats <- rep_len(1:4,length.out = 2184)
TMdata <-
  TMdata %>% mutate(growth = ((ecuatorial + polar) / 2))

get_range <- function(mynumber){
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}

TM_filtered <- TMdata %>%
  group_by(ID, experimental_replicate, dose) %>%
  mutate(growth_range = list(get_range(growth))) %>%
  unnest() %>%
  filter(growth <= upper & growth >= lower) %>%
  rename(response = growth) %>% 
  ungroup()
#Plotting data
p <- ggplot(TM_filtered, aes(x = ecuatorial, y = polar, fill= ID))
p +  geom_boxplot(notch = TRUE, aes(group =dose, experimental_replicate), outlier.colour = "#000000") + facet_wrap(~ ID)  +labs(title="Plot of dose ~ response by experimental_replicate", x ="Dose (ppm)", y = "Response (cm)")
#Getting EC50 by ECtable function
xx <- EC_table(TM_filtered, form = response ~ dose)
# Filterring by baseline
baseline <- xx %>% filter(sample=="1"|sample=="118"|sample=="123"|sample=="12B"|sample=="20"|sample=="21"|sample=="449"|sample=="461"|sample=="467"|sample=="475"|sample=="558"|sample=="564"|sample=="568"|sample=="581"|sample=="645"|sample=="667"|sample=="74SS1"|sample=="8"|sample=="87")
summary(baseline)
#Filterring by baseline
survey<- xx %>% filter(sample=="1025"|sample=="1026"|sample=="1027"|sample=="1029"|sample=="1032"|sample=="1033"|sample=="318"|sample=="413"|sample=="62-02"|sample=="62-03"|sample=="62-04"|sample=="78-01"|sample=="78-02"|sample=="78-05"|sample=="H-01"|sample=="H-03"|sample=="H-04")
summary(survey)

# Getting the Relative growth (RG) at each dose
RG <- TM_filtered %>% group_by(ID,dose) %>% summarise(mean_response=mean(response,na.rm=TRUE)) %>% spread(dose,mean_response) %>%  mutate(RG0.75=((`0.75`/`0`)*100)) %>% mutate(RG1=((`1`/`0`)*100)) %>% mutate(RG1.5=((`1.5`/`0`)*100)) %>% mutate(RG2=((`2`/`0`)*100)) %>% mutate(RG2.5=((`2.5`/`0`)*100)) %>%  mutate(RG10=((`10`/`0`)*100))



TM_xx <-  xx %>% 
  select(-Estimate.10,-SE.10,-Estimate.90,-SE.90)

names(TM_xx)[names(TM_xx) == "sample"] <- "ID"
final_TM<-left_join (RG, TM_xx) %>% mutate(logEC50=(log(Estimate.50)))

##
finalRG0.75 <- lm(logEC50~RG0.75,final_TM) 
summary(finalRG0.75)
check_assumptions(finalRG0.75)

ggplot(finalRG0.75,aes(x=logEC50, y=RG0.75))+  geom_point() + geom_smooth()
##
finalRG1 <- lm(logEC50~RG1,final_TM) 
summary(finalRG1)
check_assumptions(finalRG1)
ggplot(finalRG1,aes(x=logEC50, y=RG1))+  geom_point() + geom_smooth()

###
finalRG1.5 <- lm(logEC50~RG1.5,final_TM) 
summary(finalRG1.5)
check_assumptions(finalRG1.5)

ggplot(finalRG1.5,aes(x=logEC50, y=RG1.5))+  geom_point() + geom_smooth()


###
finalRG2 <- lm(logEC50~RG2,final_TM) 
summary(finalRG2)
check_assumptions(finalRG2)

ggplot(finalRG2,aes(x=logEC50, y=RG2))+  geom_point() + geom_smooth()

###

finalRG2.5<- lm(logEC50~RG2.5,final_TM) 
summary(finalRG2.5)
check_assumptions(finalRG2.5)

ggplot(finalRG2.5,aes(x=logEC50, y=RG2.5))+  geom_point() + geom_smooth()

##
finalRG10<- lm(logEC50~RG10,final_TM) 
summary(finalRG10)
check_assumptions(finalRG2.5)

ggplot(finalRG10,aes(x=logEC50, y=RG10))+  geom_point() + geom_smooth()
