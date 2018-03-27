#running data
library(tidyverse)
library(ggplot2)
library(ezec)
library(NRES803)
#READING AND SUBSETTING DATA
boscaliddata <-
  read.csv("data/boscalid_serialdilution_baseline(validationdata).csv")
boscaliddata <- subset(boscaliddata, select = -X)
boscaliddata$repeats <- rep_len(1:4, length.out = 2186)
boscaliddata <- boscaliddata %>%  filter(!dose ==  0.4 & !ID == 800)
boscaliddata <-
  boscaliddata %>% mutate(growth = ((ecuatorial + polar) / 2))
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}

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
) + facet_wrap( ~ ID)  + labs(title = "Plot of dose ~ response by experimental_replicate", x =
                                "Dose (ppm)", y = "Response (cm)")
#Getting EC50 by ECtable
xx <- EC_table(boscalid_filtered, form = response ~ dose)
#filterring by baseline
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
#filterring by survey
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
  boscalid_filtered %>% group_by(ID, dose) %>% summarise(mean_response =
                                                           mean(response, na.rm = TRUE)) %>% spread(dose, mean_response) %>%  mutate(RG0.025 =
                                                                                                                                       ((`0.025` / `0`) * 100)) %>% mutate(RG0.05 = ((`0.05` / `0`) * 100)) %>% mutate(RG0.05 =
                                                                                                                                                                                                                         ((`0.05` / `0`) * 100)) %>% mutate(RG0.1 = ((`0.1` / `0`) * 100)) %>% mutate(RG0.2 =
                                                                                                                                                                                                                                                                                                        ((`0.2` / `0`) * 100)) %>%  mutate(RG0.8 = ((`0.8` / `0`) * 100))



boscalid_xx <-  xx %>%
  select(-Estimate.10, -SE.10, -Estimate.90, -SE.90)

names(boscalid_xx)[names(boscalid_xx) == "sample"] <- "ID"
final_boscalid <-
  left_join (RG, boscalid_xx) %>% mutate(logEC50 = (log(Estimate.50)))
pdf("sample.pdf")
##
finalRG0.025 <- lm(logEC50 ~ RG0.025, final_boscalid)
summary(finalRG0.025)
check_assumptions(finalRG0.025)

ggplot(finalRG0.025, aes(x = logEC50, y = RG0.025)) +  geom_point() + geom_smooth()
##
finalRG0.05 <- lm(logEC50 ~ RG0.05, final_boscalid)
summary(finalRG0.05)
check_assumptions(finalRG0.05)
ggplot(finalRG0.05, aes(x = logEC50, y = RG0.05)) +  geom_point() + geom_smooth()

###
finalRG0.1 <- lm(logEC50 ~ RG0.1, final_boscalid)
summary(finalRG0.1)
check_assumptions(finalRG0.1)

ggplot(finalRG0.1, aes(x = logEC50, y = RG0.1)) +  geom_point() + geom_smooth()


###
finalRG0.2 <- lm(logEC50 ~ RG0.2, final_boscalid)
summary(finalRG0.2)
check_assumptions(finalRG0.2)

ggplot(finalRG0.2, aes(x = logEC50, y = RG0.2)) +  geom_point() + geom_smooth()

#

finalRG0.8 <- lm(logEC50 ~ RG0.8, final_boscalid)
summary(finalRG0.8)
check_assumptions(finalRG0.8)

ggplot(finalRG0.8, aes(x = logEC50, y = RG0.8)) +  geom_point() + geom_smooth()
dev.off()