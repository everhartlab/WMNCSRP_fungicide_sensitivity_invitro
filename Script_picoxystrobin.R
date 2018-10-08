#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
library(broom)
#READING AND SUBSETTING DATA
picoxystrobin_data <-
  read.csv("data/picoxystrobin_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
# Deleting columns
picoxystrobin_data <- subset(picoxystrobin_data, select = -(X:X.9))
# Setting 4 repeats for each treatment
picoxystrobin_data$repeats <- rep_len(1:4, length.out = 2016)
# Getting growth as an average between ecuatorial and polar

picoxystrobin_data <- picoxystrobin_data %>%  
    group_by(ID, experimental_replicate, dose, ecuatorial, polar, repeats) %>%
  filter( !ID == "78-05" | !experimental_replicate ==1 | !dose== 0 | !repeats== 2 , !ID == "62-04" | !experimental_replicate ==2 | !dose== 0 | !repeats== 3, !ID == "H-04" | !experimental_replicate ==2 | !dose== 0 | !repeats== 2 )
picoxystrobin_data <-
  picoxystrobin_data %>% mutate(growth = ((ecuatorial + polar) / 2))
# Creating a function for using the outputs of boxplot.stats function and notice the outlayers
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}
# Using the function in order to take out the outlayers
picoxystrobin_filtered <- picoxystrobin_data %>%
  group_by(ID, dose) %>%
  mutate(growth_range = list(get_range(growth))) %>%
  unnest() %>%
  filter(growth <= upper & growth >= lower) %>%
  rename(response = growth) %>%
  ungroup()%>% 
  select(c(ID, experimental_replicate, repeats, dose, response))

#FIRST NORMALITY

##Shapiro_test

shapiro.test_picoxystrobin <- picoxystrobin_filtered %>% 
  group_by(ID) %>% 
  do(tidy(shapiro.test(.$response))) 

picoxystrobinfinal <- left_join (shapiro.test_picoxystrobin, picoxystrobin_filtered) %>% 
  mutate(normality = ifelse(p.value >0.05, "normal", "nonormal")) %>% 
  spread(experimental_replicate, response)
View(picoxystrobinfinal)

##Wilcox test
wilcox_picoxystrobin<- picoxystrobinfinal%>%
  group_by(ID, dose) %>%
  do(tidy(wilcox.test(.$`1`, .$`2`, paired=TRUE)))
View(wilcox_picoxystrobin)


# Plotting data
p <-
  ggplot(picoxystrobin_filtered, aes(x = ecuatorial, y = polar, fill = ID))
p +  geom_boxplot(
  notch = TRUE,
  aes(group = dose, experimental_replicate),
  outlier.colour = "#000000"
) + facet_wrap( ~ ID)  + labs(title = "Plot of dose ~ response by experimental_replicate", x =
                                "Dose (ppm)", y = "Response (cm)")
# Getting EC50 by EC_table function 
xx <- EC_table(picoxystrobin_filtered, form = response ~ dose)
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
#Getting relative growth column at each dose where RG= Relative Growth, eg. RG0.01= Relative Growth at 0.01 ppm
RG <-
  picoxystrobin_filtered %>% group_by(ID, dose) %>% summarise(mean_response =
                                                                mean(response, na.rm = TRUE)) %>% spread(dose, mean_response) %>%  mutate(RG0.01 =
                                                                                                                                            ((`0.01` / `0`) * 100)) %>% mutate(RG0.02 = ((`0.02` / `0`) * 100)) %>% mutate(RG0.04 =
                                                                                                                                                                                                                             ((`0.04` / `0`) * 100)) %>% mutate(RG0.06 = ((`0.06` / `0`) * 100)) %>% mutate(RG0.1 = ((`0.1` / `0`) * 100))

# Taking out some columns                                                                                                                                                                                                                                                                                                              
picoxystrobin_xx <-  xx %>%
  select(-Estimate.10, -SE.10, -Estimate.90, -SE.90)
# Replacing "sample" for "ID"
names(picoxystrobin_xx)[names(picoxystrobin_xx) == "sample"] <- "ID"

# Getting the log of the EC50
final_picoxystrobin <-
  left_join (RG, picoxystrobin_xx) %>% mutate(logEC50 = (log(Estimate.50)))
pdf("picoxystrobin_assumptions_linearmodel_each_dose.pdf")

#Dose chosen as DD
# Linear model of log EC50 and relative growth at 0.01 ppm, check normality and homogeneity of variances

finalRG0.01 <- lm(logEC50 ~ RG0.01, final_picoxystrobin)
summary(finalRG0.01)
plot(finalRG0.01)
check_assumptions(finalRG0.01)
    # Linear Regression graph of log EC50 and relative growth at 0.01 ppm 

ggplot(finalRG0.01, aes(x = RG0.01, y = logEC50)) +  geom_point() + geom_smooth(method = "lm")

    # Getting the EC50DD according to the model
final_picoxystrobin_DD <- final_picoxystrobin %>% mutate(Estimate.50DD = exp(-6.479912 + (0.036911*RG0.01)))
    # Linear model of EC50  and EC50DD, check normality and homogeneity of variances
final_picoxystrobin_DD_0.01 <- lm (Estimate.50DD ~ Estimate.50, final_picoxystrobin_DD)
summary(final_picoxystrobin_DD_0.01 )
plot(final_picoxystrobin_DD_0.01)
check_assumptions(final_picoxystrobin_DD_0.01 )
pdf("Linear regression of picoxystrobin 0.01 ppm.pdf")  
 # Linear Regression graph of EC50 and EC50DD
ggplot(final_picoxystrobin_DD_0.01, aes(x = Estimate.50, y = Estimate.50DD)) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 1)) + labs(title = "Linear regression of picoxystrobin 0.01 ppm", x ="EC50 (ppm)", y = "EC50DD (ppm)" ) + geom_point() + geom_smooth(method = "lm")
dev.off()


# Linear model of log EC50 and relative growth at 0.02 ppm, check normality and homogeneity of variances

finalRG0.02 <- lm(logEC50 ~ RG0.02, final_picoxystrobin)
summary(finalRG0.02)
plot(finalRG0.02)
check_assumptions(finalRG0.02)
    # Linear Regression graph of log EC50 and relative growth at 0.02 ppm 
ggplot(finalRG0.02, aes(x = RG0.02, y = logEC50)) +  geom_point() + geom_smooth(method = "lm")

# Linear model of log EC50 and relative growth at 0.04 ppm, check normality and homogeneity of variances

finalRG0.04 <- lm(logEC50 ~ RG0.04, final_picoxystrobin)
summary(finalRG0.04)
plot(finalRG0.04)
check_assumptions(finalRG0.04)
    # Linear Regression graph of log EC50 and relative growth at 0.04 ppm 

ggplot(finalRG0.04, aes(x = RG0.04, y = logEC50)) +  geom_point() + geom_smooth(method = "lm")


# Linear model of log EC50 and relative growth at 0.06 ppm, check normality and homogeneity of variances
finalRG0.06 <- lm(logEC50 ~ RG0.06, final_picoxystrobin)
summary(finalRG0.06)
plot(finalRG0.06)
check_assumptions(finalRG0.06)
    # Linear Regression graph of log EC50 and relative growth at 0.06 ppm 

ggplot(finalRG0.06, aes(x = RG0.06, y = logEC50)) +  geom_point() + geom_smooth(method = "lm")


# Linear model of log EC50 and relative growth at 0.1 ppm, check normality and homogeneity of variances


finalRG0.1<- lm(logEC50 ~ RG0.1, final_picoxystrobin)
summary(finalRG0.1)
plot(RG0.1)
check_assumptions(finalRG0.1)
    # Linear Regression graph of log EC50 and relative growth at 0.1 ppm 

ggplot(finalRG0.1, aes(x = RG0.1, y = logEC50)) +  geom_point() + geom_smooth(method = "lm")
dev.off()