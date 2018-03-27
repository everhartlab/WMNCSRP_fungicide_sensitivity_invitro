#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
picoxystrobin_data <-
  read.csv("data/picoxystrobin_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
picoxystrobin_data <- subset(picoxystrobin_data, select = -(X:X.9))
picoxystrobin_data$repeats <- rep_len(1:4, length.out = 2017)
picoxystrobin_data <-
  picoxystrobin_data %>% mutate(growth = ((ecuatorial + polar) / 2))
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}

picoxystrobin_filtered <- picoxystrobin_data %>%
  group_by(ID, experimental_replicate, dose) %>%
  mutate(growth_range = list(get_range(growth))) %>%
  unnest() %>%
  filter(growth <= upper & growth >= lower) %>%
  rename(response = growth) %>%
  ungroup()
#Plotting data
p <-
  ggplot(picoxystrobin_filtered, aes(x = ecuatorial, y = polar, fill = ID))
p +  geom_boxplot(
  notch = TRUE,
  aes(group = dose, experimental_replicate),
  outlier.colour = "#000000"
) + facet_wrap( ~ ID)  + labs(title = "Plot of dose ~ response by experimental_replicate", x =
                                "Dose (ppm)", y = "Response (cm)")
#getting EC50 by EC_table
xx <- EC_table(picoxystrobin_filtered, form = response ~ dose)
#Filterring by baseline
baseline <-
  xx %>% filter(
    sample == "1" |
      sample == "118" |
      sample == "123" |
      sample == "12B" |
      sample == "20" |
      sample == "21" |
      sample == "449" |
      sample == "461" |
      sample == "467" |
      sample == "475" |
      sample == "558" |
      sample == "564" |
      sample == "568" |
      sample == "581" |
      sample == "645" |
      sample == "667" | sample == "74SS1" | sample == "8" | sample == "87"
  )
summary(baseline)
#Filterring by baseline
survey <-
  xx %>% filter(
    sample == "1025" |
      sample == "1026" |
      sample == "1027" |
      sample == "1029" |
      sample == "1032" |
      sample == "1033" |
      sample == "318" |
      sample == "413" |
      sample == "62-02" |
      sample == "62-03" |
      sample == "62-04" |
      sample == "78-01" |
      sample == "78-02" |
      sample == "78-05" | sample == "H-01" | sample == "H-03" |
      sample == "H-04"
  )
summary(survey)


RG <-
  picoxystrobin_filtered %>% group_by(ID, dose) %>% summarise(mean_response =
                                                                mean(response, na.rm = TRUE)) %>% spread(dose, mean_response) %>%  mutate(RG0.01 =
                                                                                                                                            ((`0.01` / `0`) * 100)) %>% mutate(RG0.02 = ((`0.02` / `0`) * 100)) %>% mutate(RG0.04 =
                                                                                                                                                                                                                             ((`0.04` / `0`) * 100)) %>% mutate(RG0.06 = ((`0.06` / `0`) * 100)) %>% mutate(RG0.1 =
                                                                                                                                                                                                                                                                                                              ((`0.1` / `0`) * 100)) 
picoxystrobin_xx <-  xx %>%
  select(-Estimate.10, -SE.10, -Estimate.90, -SE.90)

names(picoxystrobin_xx)[names(picoxystrobin_xx) == "sample"] <- "ID"
final_picoxystrobin <-
  left_join (RG, picoxystrobin_xx) %>% mutate(logEC50 = (log(Estimate.50)))
pdf("sample.pdf")
##
finalRG0.01 <- lm(logEC50 ~ RG0.01, final_picoxystrobin)
summary(finalRG0.01)
check_assumptions(finalRG0.01)

ggplot(finalRG0.01, aes(x = logEC50, y = RG0.01)) +  geom_point() + geom_smooth()
##
finalRG0.02 <- lm(logEC50 ~ RG0.02, final_picoxystrobin)
summary(finalRG0.02)
check_assumptions(finalRG0.02)
ggplot(finalRG0.02, aes(x = logEC50, y = RG0.02)) +  geom_point() + geom_smooth()

###
finalRG0.04 <- lm(logEC50 ~ RG0.04, final_picoxystrobin)
summary(finalRG0.04)
check_assumptions(finalRG0.04)

ggplot(finalRG0.04, aes(x = logEC50, y = RG0.04)) +  geom_point() + geom_smooth()


###
finalRG0.06 <- lm(logEC50 ~ RG0.06, final_picoxystrobin)
summary(finalRG0.06)
check_assumptions(finalRG0.06)

ggplot(finalRG0.06, aes(x = logEC50, y = RG0.06)) +  geom_point() + geom_smooth()
dev.off()
#

# finalRG0.8 <- lm(logEC50 ~ RG0.8, final_picoxystrobin)
# summary(finalRG0.8)
# check_assumptions(finalRG0.8)
# 
# ggplot(finalRG0.8, aes(x = logEC50, y = RG0.8)) +  geom_point() + geom_smooth()