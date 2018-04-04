#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
TMdata <- read.csv("data/TM_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
# Deleting columns
TMdata <- subset(TMdata, select = - ( X: X.6 ))
# Setting 4 repeats for each treatment
TMdata$repeats <- rep_len(1:4,length.out = 2184)
# Getting growth as an average between ecuatorial and polar
TMdata <-
  TMdata %>% mutate(growth = ((ecuatorial + polar) / 2))
# Creating a function for using the outputs of boxplot.stats function and notice the outlayers
get_range <- function(mynumber){
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}
# Using the function in order to take out the outlayers
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

#Setting the groups of isolates
baseline_isolates <- c("1","118", "123", "12B", "129", "20", "21", "449", "461", "467", "475", "558", "564", "568", "581", "645", "800", "667", "74SS1", "8", "87")
survey_isolates <-  c("318", "413", "419", "62-02", "62-03", "62-04", "78-01", "78-02", "78-05", "H-01", "H-03", "H-04", "I-20", "S-01", "W212")
treatmentyear2016_isolates <- c("1025", "1026", "1027", "1029", "1032","1032")
# Filtering by baseline group
baseline <-  xx %>% filter(sample %in% baseline_isolates)
summary(baseline)
# Filtering by survey group
survey <-  xx %>% filter(sample %in% survey_isolates)
summary(survey)

# Filtering by treatmentyear2016 group
treatmentyear2016 <- xx %>% filter(sample %in% treatmentyear2016_isolates)
summary(treatmentyear2016)
# Getting relative growth column at each dose where RG= Relative Growth, eg. RG0.75= Relative Growth at 0.75 ppm
RG <- TM_filtered %>% group_by(ID,dose) %>% summarise(mean_response=mean(response,na.rm=TRUE)) %>% spread(dose,mean_response) %>%  mutate(RG0.75=((`0.75`/`0`)*100)) %>% mutate(RG1=((`1`/`0`)*100)) %>% mutate(RG1.5=((`1.5`/`0`)*100)) %>% mutate(RG2=((`2`/`0`)*100)) %>% mutate(RG2.5=((`2.5`/`0`)*100)) %>%  mutate(RG10=((`10`/`0`)*100))


# Taking out some columns
TM_xx <-  xx %>% 
  select(-Estimate.10,-SE.10,-Estimate.90,-SE.90)
# Replacing "sample" for "ID"
names(TM_xx)[names(TM_xx) == "sample"] <- "ID"
# Getting the log of the EC50
final_TM<-left_join (RG, TM_xx) %>% mutate(logEC50=(log(Estimate.50)))
pdf("TM_assumptions_linearmodel_each_dose.pdf")
# Linear model of log EC50 and relative growth at 0.75ppm, check normality and homogeneity of variancesfinalRG0.75 <- lm(logEC50~RG0.75,final_TM) 
summary(finalRG0.75)
check_assumptions(finalRG0.75)
    # Linear Regression graph oflog EC50 and relative growth at 0.75 ppm  
ggplot(finalRG0.75,aes(x=logEC50, y=RG0.75))+  geom_point() + geom_smooth(method = "lm")
# Linear model of log EC50 and relative growth at 1 ppm, check normality and homogeneity of variances
finalRG1 <- lm(logEC50~RG1,final_TM) 
summary(finalRG1)
check_assumptions(finalRG1)
    # Linear Regression graph of log EC50 and relative growth at 1 ppm 
ggplot(finalRG1,aes(x=logEC50, y=RG1))+  geom_point() + geom_smooth(method = "lm")

# Linear model of log EC50 and relative growth at 1.5 ppm, check normality and homogeneity of variances
finalRG1.5 <- lm(logEC50~RG1.5,final_TM) 
summary(finalRG1.5)
check_assumptions(finalRG1.5)
    # Linear Regression graph of log EC50 and relative growth at 1.5 ppm
ggplot(finalRG1.5,aes(x=logEC50, y=RG1.5))+  geom_point() + geom_smooth(method = "lm")


# Linear model of log EC50 and relative growth at 2 ppm, check normality and homogeneity of variances
finalRG2 <- lm(logEC50~RG2,final_TM) 
summary(finalRG2)
check_assumptions(finalRG2)
      # Linear Regression graph of log EC50 and relative growth at 2 ppm
ggplot(finalRG2,aes(x=logEC50, y=RG2))+  geom_point() + geom_smooth(method = "lm")

# Linear model of log EC50 and relative growth at 2.5 ppm, check normality and homogeneity of variances

finalRG2.5<- lm(logEC50~RG2.5,final_TM) 
summary(finalRG2.5)
check_assumptions(finalRG2.5)
      # Linear Regression graph of log EC50 and relative growth at 2.5 ppm
ggplot(finalRG2.5,aes(x=logEC50, y=RG2.5))+  geom_point() + geom_smooth()

#Dose chosen as DD
# Linear model of log EC50 and relative growth at 10 ppm, check normality and homogeneity of variances
finalRG10<- lm(logEC50~RG10,final_TM) 
summary(finalRG10)
check_assumptions(finalRG2.5)
    # Linear Regression graph of log EC50 and relative growth at 10 ppm
ggplot(finalRG10,aes(x=logEC50, y=RG10))+  geom_point() + geom_smooth()

# Getting the EC50DD according to the model
final_TM_DD <- final_TM %>% mutate(Estimate.50DD = exp(1.19249 - (0.02646*RG10)))
# Linear model of EC50  and EC50DD, check normality and homogeneity of variances
final_TM_DD_10 <- lm (Estimate.50DD ~ Estimate.50, final_TM_DD)
summary(final_TM_DD_10)
check_assumptions(final_TM_DD_10)
pdf("Linear regression of thiophanate methyl 10 ppm.pdf")   
# Linear Regression graph of EC50 and EC50DD
ggplot(final_TM_DD_10, aes(x = Estimate.50, y = Estimate.50DD)) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 1)) + labs(title = "Linear regression of thiophanate methyl 10 ppm", x ="EC50 (ppm)", y = "EC50DD (ppm)" ) + geom_point() + geom_smooth(method = "lm")

dev.off()