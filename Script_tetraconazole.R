#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
tetraconazole_data <-
  read.csv("data/tetraconazole_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
# tetraconazole_data <- subset(tetraconazole_data, select = -(X:X.9))

# Setting 4 repeats for each treatment
tetraconazole_data$repeats <- rep_len(1:4, length.out = 1680)
# tetraconazole_data <-
# tetraconazole_data %>% mutate(growth = ((ecuatorial + polar) / 2))

# Creating a function for using the outputs of boxplot.stats function and notice the outlayers
tetraconazole_data <-  tetraconazole_data%>% select(ID:dose,repeats,response)
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}

# Using the function in order to take out the outlayers
tetraconazole_filtered <- tetraconazole_data %>%
  group_by(ID, experimental_replicate, dose) %>%
  mutate(response_range = list(get_range(response))) %>%
  unnest() %>%
  filter(response <= upper & response >= lower) %>%
  # rename(response = growth) %>%
  ungroup()

# p <-
#   ggplot(tetraconazole_filtered, aes(x = ecuatorial, y = polar, fill = ID))
# p +  geom_boxplot(
#   notch = TRUE,
#   aes(group = dose, experimental_replicate),
#   outlier.colour = "#000000"
# ) + facet_wrap( ~ ID)  + labs(title = "Plot of dose ~ response by experimental_replicate", x =
#                                 "Dose (ppm)", y = "Response (cm)")

# Getting EC50 by EC_table function
xx <- EC_table(tetraconazole_filtered, form = response ~ dose)
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


# Getting relative growth column at each dose where RG= Relative Growth, eg. RG0.5= Relative Growth at 0.5 ppm
RG <-
  tetraconazole_filtered %>% group_by(ID, dose) %>% summarise(mean_response =
                                                                mean(response, na.rm = TRUE)) %>% spread(dose, mean_response) %>%  mutate(RG0.5 =
                                                                                                                                            ((`0.5` / `0`) * 100)) %>% mutate(RG1 = ((`1` / `0`) * 100)) %>% mutate(RG2 =
                                                                                                                                                                                                                      ((`2` / `0`) * 100)) %>% mutate(RG3 = ((`3` / `0`) * 100)) %>% mutate(RG5 =
                                                                                                                                                                                                                                                                                              ((`5` / `0`) * 100)) 

# Taking out some columns
tetraconazole_xx <-  xx %>%
  select(-Estimate.10, -SE.10, -Estimate.90, -SE.90)
# Replacing "sample" for "ID"
names(tetraconazole_xx)[names(tetraconazole_xx) == "sample"] <- "ID"

# Getting the log of the EC50
final_tetraconazole <-
  left_join (RG, tetraconazole_xx) %>% mutate(logEC50 = (log(Estimate.50)))
pdf("tetraconazole_assumptions_linearmodel_each_dose.pdf")
# Linear model of log EC50 and relative growth at 0.5 ppm, check normality and homogeneity of variancesfinalRG0.5 <- lm(logEC50 ~ RG0.5, final_tetraconazole)
summary(finalRG0.5)
check_assumptions(finalRG0.5)
    # Linear Regression graph oflog EC50 and relative growth at 0.5 ppm
ggplot(finalRG0.5, aes(x = logEC50, y = RG0.5)) +  geom_point() + geom_smooth(method = "lm")
# Linear model of log EC50 and relative growth at 1 ppm, check normality and homogeneity of variances
finalRG1 <- lm(logEC50 ~ RG1, final_tetraconazole)
summary(finalRG1)
check_assumptions(finalRG1)
    # Linear Regression graph of log EC50 and relative growth at 1 ppm
ggplot(finalRG1, aes(x = logEC50, y = RG1)) +  geom_point() + geom_smooth(method = "lm")

#Dose chosen as DD
# Linear model of log EC50 and relative growth at 2 ppm, check normality and homogeneity of variances
finalRG2 <- lm(logEC50 ~ RG2, final_tetraconazole)
summary(finalRG2)
check_assumptions(finalRG2)
    # Linear Regression graph of log EC50 and relative growth at 2 ppm
ggplot(finalRG2, aes(x = logEC50, y = RG2)) +  geom_point() + geom_smooth(method = "lm")

    # Getting the EC50DD according to the model
final_tetraconazole_DD <- final_tetraconazole %>% mutate(Estimate.50DD = exp(-1.336487 + (0.036323*RG2)))
    # Linear model of EC50  and EC50DD, check normality and homogeneity of variances
final_tetraconazole_DD_2 <- lm (Estimate.50DD ~ Estimate.50, final_tetraconazole_DD)
summary(final_tetraconazole_DD_2)
check_assumptions(final_tetraconazole_DD_2)
pdf("Linear regression of tetraconazole 2 ppm.pdf")   
 # Linear Regression graph of EC50 and EC50DD
ggplot(final_tetraconazole_DD_2, aes(x = Estimate.50, y = Estimate.50DD)) +  theme(plot.title = element_text(size = 16, face = "bold", hjust = 1)) + labs(title = "Linear regression of tetraconazole 2 ppm", x ="EC50 (ppm)", y = "EC50DD (ppm)" ) + geom_point() + geom_smooth(method = "lm")
dev.off()
# Linear model of log EC50 and relative growth at 3 ppm, check normality and homogeneity of variances
finalRG3 <- lm(logEC50 ~ RG3, final_tetraconazole)
summary(finalRG3)
check_assumptions(finalRG3)
    # Linear Regression graph of log EC50 and relative growth at 3 ppm
ggplot(finalRG3, aes(x = logEC50, y = RG3)) +  geom_point() + geom_smooth(method = "lm")

# Linear model of log EC50 and relative growth at 5 ppm, check normality and homogeneity of variances
finalRG5 <- lm(logEC50 ~ RG5, final_tetraconazole)
summary(finalRG5)
check_assumptions(finalRG5)
      # Linear Regression graph of log EC50 and relative growth at 5 ppmggplot(finalRG5, aes(x = logEC50, y = RG5)) +  geom_point() + geom_smooth(method = "lm")
ggplot(finalRG5, aes(x = logEC50, y = RG5)) +  geom_point() + geom_smooth(method = "lm")

dev.off()