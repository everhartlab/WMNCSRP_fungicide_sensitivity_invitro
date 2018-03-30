#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
tetraconazole_data <-
  read.csv("data/tetraconazole_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
# tetraconazole_data <- subset(tetraconazole_data, select = -(X:X.9))
tetraconazole_data$repeats <- rep_len(1:4, length.out = 1680)
# tetraconazole_data <-
# tetraconazole_data %>% mutate(growth = ((ecuatorial + polar) / 2))

tetraconazole_data <-  tetraconazole_data%>% select(ID:dose,repeats,response)
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}

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
#Getting EC50 by ECtable
xx <- EC_table(tetraconazole_filtered, form = response ~ dose)
#filterring by baseline
baseline <-
  xx %>% filter(
    sample == "1" |
      sample == "118" |
      sample == "123" |
      sample == "12B" |
      sample == "129" |
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
      sample == "800" |
      sample == "667" | sample == "74SS1" | sample == "8" | sample == "87"
  )
summary(baseline)
#filterring by survey
survey <-
  xx %>% filter(
    sample == "318" |
      sample == "413" |
      sample == "419" |
      sample == "62-02" |
      sample == "62-03" |
      sample == "62-04" |
      sample == "78-01" |
      sample == "78-02" |
      sample == "78-05" | 
      sample == "H-01" |
      sample == "H-03" |
      sample == "H-04"|
      sample == "I-20" |sample == "S-01"| sample == "W212"
  )
summary(survey)

#filterring by treatmentyear2016
treatmentyear2016 <-
  xx %>% filter(
    sample == "1025" |
      sample == "1026" |
      sample == "1027" |
      sample == "1029" |
      sample == "1032" |
      sample == "1033" 
  )
summary(treatmentyear2016)

RG <-
  tetraconazole_filtered %>% group_by(ID, dose) %>% summarise(mean_response =
                                                                mean(response, na.rm = TRUE)) %>% spread(dose, mean_response) %>%  mutate(RG0.5 =
                                                                                                                                            ((`0.5` / `0`) * 100)) %>% mutate(RG1 = ((`1` / `0`) * 100)) %>% mutate(RG2 =
                                                                                                                                                                                                                      ((`2` / `0`) * 100)) %>% mutate(RG3 = ((`3` / `0`) * 100)) %>% mutate(RG5 =
                                                                                                                                                                                                                                                                                              ((`5` / `0`) * 100)) 
tetraconazole_xx <-  xx %>%
  select(-Estimate.10, -SE.10, -Estimate.90, -SE.90)

names(tetraconazole_xx)[names(tetraconazole_xx) == "sample"] <- "ID"
final_tetraconazole <-
  left_join (RG, tetraconazole_xx) %>% mutate(logEC50 = (log(Estimate.50)))
pdf("tetraconazole_assumptions_linearmodel_each_dose.pdf")
##
finalRG0.5 <- lm(logEC50 ~ RG0.5, final_tetraconazole)
summary(finalRG0.5)
check_assumptions(finalRG0.5)

ggplot(finalRG0.5, aes(x = logEC50, y = RG0.5)) +  geom_point() + geom_smooth()
##
finalRG1 <- lm(logEC50 ~ RG1, final_tetraconazole)
summary(finalRG1)
check_assumptions(finalRG1)
ggplot(finalRG1, aes(x = logEC50, y = RG1)) +  geom_point() + geom_smooth()

###
finalRG2 <- lm(logEC50 ~ RG2, final_tetraconazole)
summary(finalRG2)
check_assumptions(finalRG2)

ggplot(finalRG2, aes(x = logEC50, y = RG2)) +  geom_point() + geom_smooth()


###
finalRG3 <- lm(logEC50 ~ RG3, final_tetraconazole)
summary(finalRG3)
check_assumptions(finalRG3)

ggplot(finalRG3, aes(x = logEC50, y = RG3)) +  geom_point() + geom_smooth()

##

finalRG5 <- lm(logEC50 ~ RG5, final_tetraconazole)
summary(finalRG5)
check_assumptions(finalRG5)

ggplot(finalRG5, aes(x = logEC50, y = RG5)) +  geom_point() + geom_smooth()

dev.off()