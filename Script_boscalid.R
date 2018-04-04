

#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
boscaliddata <-
  read.csv("data/boscalid_serialdilution_baseline(validationdata).csv")
# Deleting column
boscaliddata <- subset(boscaliddata, select = -X)
# Setting 4 repeats for each treatment
boscaliddata$repeats <- rep_len(1:4, length.out = 2186)
# Discarding dose at 0.4 ppm and isolate 800 because it has just one experimental replication
boscaliddata <- boscaliddata %>%  filter(!dose ==  0.4 & !ID == 800)
# Getting growth as an average between ecuatorial and polar
boscaliddata <-
  boscaliddata %>% mutate(growth = ((ecuatorial + polar) / 2))
# Creating a function for using the outputs of boxplot.stats function and notice the outlayers
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}
# Using the function in order to take out the outlayers
boscalid_filtered <- boscaliddata %>%
  group_by(ID, experimental_replicate, dose) %>%
  mutate(growth_range = list(get_range(growth))) %>%
  unnest() %>%
  filter(growth <= upper & growth >= lower) %>%
  rename(response = growth) %>%
  ungroup()
# Plotting data
p <-
  ggplot(boscalid_filtered, aes(x = ecuatorial, y = polar, fill = ID))
p +  geom_boxplot(
  notch = TRUE,
  aes(group = dose, experimental_replicate),
  outlier.colour = "#000000"
) + facet_wrap(~ ID)  + labs(title = "Plot of dose ~ response by experimental_replicate", x =
                               "Dose (ppm)", y = "Response (cm)")
# Getting EC50 by EC_table function
xx <- EC_table(boscalid_filtered, form = response ~ dose)
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
# Getting relative growth column at each dose where RG = Relative Growth, eg. RG0.025 = Relative Growth at 0.025 ppm
RG <-
  boscalid_filtered %>% group_by(ID, dose) %>% summarise(mean_response =
                                                           mean(response, na.rm = TRUE)) %>% spread(dose, mean_response) %>%  mutate(RG0.025 =
                                                                                                                                       ((`0.025` / `0`) * 100)) %>% mutate(RG0.05 = ((`0.05` / `0`) * 100)) %>% mutate(RG0.05 =
                                                                                                                                                                                                                         ((`0.05` / `0`) * 100)) %>% mutate(RG0.1 = ((`0.1` / `0`) * 100)) %>% mutate(RG0.2 =
                                                                                                                                                                                                                                                                                                        ((`0.2` / `0`) * 100)) %>%  mutate(RG0.8 = ((`0.8` / `0`) * 100))


# Taking out some columns
boscalid_xx <-  xx %>%
  select(-Estimate.10,-SE.10,-Estimate.90,-SE.90)
# Replacing "sample" for "ID"
names(boscalid_xx)[names(boscalid_xx) == "sample"] <- "ID"
# Getting the log of the EC50
final_boscalid <-
  left_join (RG, boscalid_xx) %>% mutate(logEC50 = (log(Estimate.50)))
pdf("boscalid_assumptions_linearmodel_each_dose.pdf")
# Linear model of log EC50 and relative growth at 0.025 ppm, check normality and homogeneity of variances
finalRG0.025 <- lm(logEC50 ~ RG0.025, final_boscalid)
summary(finalRG0.025)
check_assumptions(finalRG0.025)
    # Linear Regression graph oflog EC50 and relative growth at 0.025 ppm
ggplot(finalRG0.025, aes(x = logEC50, y = RG0.025)) +  geom_point() + geom_smooth(method = "lm")
# Linear model of log EC50 and relative growth at 0.05 ppm, check normality and homogeneity of variances
finalRG0.05 <- lm(logEC50 ~ RG0.05, final_boscalid)
summary(finalRG0.05)
check_assumptions(finalRG0.05)
    # Linear Regression graph of log EC50 and relative growth at 0.05 ppm
ggplot(finalRG0.05, aes(x = RG0.05, y = logEC50)) +  geom_point() + geom_smooth(method = "lm")

# Linear model of log EC50 and relative growth at 0.1 ppm, check normality and homogeneity of variances
finalRG0.1 <- lm(logEC50 ~ RG0.1, final_boscalid)
summary(finalRG0.1)
check_assumptions(finalRG0.1)
    # Linear Regression graph of log EC50 and relative growth at 0.1 ppm

ggplot(finalRG0.1, aes(x = logEC50, y = RG0.1)) +  geom_point() + geom_smooth(method = "lm")

#Dose chosen as DD
# Linear model of log EC50 and relative growth at 0.2 ppm, check normality and homogeneity of variances

finalRG0.2 <- lm(logEC50 ~ RG0.2, final_boscalid)
summary(finalRG0.2)
check_assumptions(finalRG0.2)
    # Linear Regression graph of log EC50 and relative growth at 0.2 ppm
ggplot(finalRG0.2, aes(x = logEC50, y = RG0.2)) +  geom_point() + geom_smooth(method = "lm")
    # Getting the EC50DD according to the model
final_boscalid_DD <- final_boscalid %>% mutate(Estimate.50DD = exp(-4.328466 + (0.054537*RG0.2)))
    # Linear model of EC50  and EC50DD, check normality and homogeneity of variances
final_boscalid_DD_0.2 <- lm (Estimate.50DD ~ Estimate.50, final_boscalid_DD)
summary(final_boscalid_DD_0.2)
check_assumptions(final_boscalid_DD_0.2)
    # Linear Regression graph of EC50 and EC50DD
pdf("Linear regression of boscalid 0.2 ppm.pdf")  
ggplot(final_boscalid_DD_0.2, aes(x = Estimate.50, y = Estimate.50DD)) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 1)) + labs(title = "Linear regression of boscalid 0.2 ppm", x ="EC50 (ppm)", y = "EC50DD (ppm)" ) + geom_point() + geom_smooth(method = "lm")
# Linear model of log EC50 and relative growth at 0.8 ppm, check normality and homogeneity of variances
dev.off()
finalRG0.8 <- lm(logEC50 ~ RG0.8, final_boscalid)
summary(finalRG0.8)
check_assumptions(finalRG0.8)
    # Linear Regression graph of log EC50 and relative growth at 0.8 ppm
ggplot(finalRG0.8, aes(x = logEC50, y = RG0.8)) +  geom_point() + geom_smooth(method = "lm")
dev.off()