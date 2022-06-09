WMNCSRP_fungicide_sensitivity_invitro
================
2022-06-09

# About

This is an example of fungicide sensitivity test *in vitro* for
*Sclerotinia sclerotiorum* by discriminatory concentrations approach
**WMNCSRP_fungicide_sensitivity_invitro**:

1.  Serial dilution method with 6 concentrations for each fungicide from
    22 baseline isolates (never exposed to fungicides) + 22 from betwenn
    farmer fields and fungicide field trials

2.  EC50 estimated from dose-response curve by ezec-package R (Kamvar,
    2016)

3.  Identify a single concentration that can be a discriminatory
    concentration

4.  Use DC to assess sensitivity of large collection

5.  Data analyses using reproducible methods

The 23 own created functions were compelled in a package FastconcR in
<https://github.com/EdgarNietoFungi/FastconcR>

# Background

Knowledge of baseline sensitivity of a pathogen population to fungicides
is important for detection of resistance and for achieving effective
disease management strategies (Brent and Holloman 2007). QoI, DMI, MBC
and SDHI fungicides are used to control *S. sclerotiorum* in an IPM
approach in soybean and dry bean in the U.S. DMI fungicides are
considered to medium risk of resistance emergence, whereas QoI, SDHI and
MBC are high risk (Fungicide Resistance Action Committee 2016). MBC and
dicarboximide fungicides are used for white mold in Brazil, wherein
resistance to MBC fungicide, Thiophanate methyl, has been reported
(Lehner et al. 2015). Fungicide resistance has also been reported in
China to the MBC, Carbendazim and dicarboximide, Dimethachlon (Ma et
al. 2009).

## **Commonly Used Fungicides For WM**

|        Fungicide         |  Efficacy of Control   | Launched |     |
|:------------------------:|:----------------------:|:--------:|:---:|
| Thiophanate methyl (MBC) |      fair control      |   1970   |     |
|   Tetraconazole (DMI)    |      fair control      |   1988   |     |
|   Picoxystrobin (QoI)    | good-very good control |   2000   |     |
|     Boscalid (SDHI)      |   very good control    |   2003   |     |

Examples of these include: **Topsin, Domark, Aproach & Endura**. Several
other products are mixtures

MBC: Methyl Benzimidazole Carbamates DMI: DeMethylation Inhibitors QoI:
Quinone outside Inhibitors SDHI: Succinate dehydrogenase inhibitors

Field applications of fungicides in the United States on soybean show
control of *S. sclerotiorum* is low to fair for QoI, DMI and MBC
fungicides, and good for SDHI fungicides (The North Central Regional
Committee on Soybean Diseases 2015). However, it is unknown whether
efficacy is related to reduce sensitivity within populations under
study.

# Hypotheses

1.  Isolates *S. sclerotiorum* collected from producer fields and
    fungicide research trials, which already have been in contact with
    **MBC (FRAC group 1), DMI (FRAC group 3), SdhI (FRAC group 7), and
    QoI (FRAC group 11) fungicides** will have less sensitivity
    **(except FRAC group 3)** compared to the baseline, because these
    fungicides have a high risk of resistance and resistance to mostly
    these groups have been reported previously in *S. sclerotiorum*;
2.  Isolates *S. sclerotiorum* collected from dry bean will have lower
    sensitivity compared to the isolates from soybean for the three
    group fungicides **(except FRAC group 3)**, which have fewer
    fungicide applications;
3.  Isolates *S. sclerotiorum* from Brazil and Mexico will have less
    sensitivity compared to the USA isolates for the three group
    fungicides **(except FRAC group 3)** because of different
    populations and management practices used in tropical and
    sub-tropical environment.

``` r
baseline.isolates <- c(1, 8, 12, 20, 21, 74, 87, 118, 123, 129, 449, 461, 467, 475, 558, 564, 568, 581, 645, 667, 800)

survey.isolates <-  c( 2385, 2386, 2098, 2099, 2100, 2139, 2140, 2143, 2220, 2222, 2223, 2320, 2362, 2388, 2390 )

Farmer.fields <- c(136:442, 455, 456, 466, 468, 470, 471, 478:554, 602:610, 612:635, 671, 672, 682:690, 695:797, 810:823, 834, 835, 848:1024, 1255: 1326, 1491: 1500, 1661:1670, 1831: 2246, 2303: 2342,2362: 2381, 2393: 2573)

baseline.isolates.2 <- c(1: 136, 444: 454, 457: 465, 467, 469, 472:477, 555:601, 611, 636:670, 673: 680, 691:694 , 798: 809,  824: 833, 836: 847)

drybean <- c(1, 5, 12,13:118, 123:128, 133:135, 145, 146, 152, 155:160, 182: 185, 194:200, 205, 220, 223, 248, 253:255, 274, 279, 280, 290, 294, 304: 309, 323, 358, 393, 395, 396, 397, 400:402, 405, 408, 409, 434, 443:465, 467, 469: 493, 495:505, 555:613, 615: 760, 762: 819, 824:833, 835:855, 858:921, 966: 971, 980, 981, 985: 996, 998: 1009, 1220:1225 , 1857: 1940, 2442: 2569)

soybean <- c(143, 147, 181, 187, 188, 189,202, 257: 259, 264:268, 276, 281, 289, 293, 295, 310 , 399, 412:417, 419, 425:427, 439, 440, 494, 506:554, 834, 972:979, 982, 997, 1010: 1022, 1025:1219, 1229: 1856, 1941:2441 )

SI.Production.Field <- c( 698:741, 743:744, 746:760, 762: 778, 
786: 797, 810:819, 848: 855, 858: 914)

SI.Screening.Nursery.Field <- c(444:454, 457:465, 467, 469, 472:477, 555: 582, 584: 601, 611, 636:670, 673:681,691:694, 798:809, 824: 833 , 836:847)# also taking into account WM Monitor 
mexican.isolates <- c(242, 248, 615:632, 779, 1857: 1940)
brazilian.isolates <-  c(400:402, 972:979, 1010:1022, 1227:1254, 2417:2569)


whitemold.inventory <- read_csv("data/whitemold.inventory.csv") %>% mutate(Field_Year= paste (Field, Year,sep = "_year_"))
```

    ## New names:
    ## Rows: 2326 Columns: 40
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (30): Host, State, County, Country, Field, Form.ID, Serial_agar_dilution... dbl
    ## (8): ...1, ID, Year, lat, long, Applications, Rate, Number.of.Years.1 lgl (2):
    ## company_current_season, company_previous_seasons
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
for (variable in unique(whitemold.inventory$State)) {
    assign(variable, whitemold.inventory %>% filter (State == variable) %>% pull(ID), envir = .GlobalEnv)
}

for (variable in unique(whitemold.inventory$Field_Year)) {
    assign(variable, whitemold.inventory %>% filter (Field_Year == variable)%>% pull(ID), envir = .GlobalEnv)
}
IDA <- ID
#MAKE THE SAME BY COUNTY
for (variable in unique(whitemold.inventory$County)) {
    assign(variable, whitemold.inventory %>% filter (County == variable)%>% pull(ID), envir = .GlobalEnv)
}
```

# Functions

``` r
#Get range function
get_range <- function(mynumber) {
  bb <- boxplot.stats(mynumber)
  cc <- bb$stats
  dd <- max(cc)
  ee <- min(cc)
  return(data.frame(upper = dd, lower = ee))
}

##Function.Reading data for serial.dilution
reading_data_serial_dilution <- function(filename){
  data <-
  read.csv(filename)
 # data <- subset(dareading_data("data/Benomyl.Colletotrichum.csv")ta, select = -X)
  data$repeats <- rep_len(1:4, length.out = nrow(data))
  
   # renaming ID for the current names
 data$ID[ data$ID =="12B"]  <- "12" 
 #Although 74SS1 belong to # 74, the one what I used was specifically this:74SS1
 data$ID[ data$ID =="74SS1"]  <- "74" 
 data$ID[ data$ID =="62-02"]  <- "2098" 
 data$ID[ data$ID =="62-03"]  <- "2099" 
 data$ID[ data$ID =="62-04"]  <- "2100" 
 data$ID[ data$ID =="78-01"]  <- "2139" 
 data$ID[ data$ID =="78-02"]  <- "2140" 
 data$ID[ data$ID =="78-05"]  <- "2143" 
 data$ID[ data$ID =="S-01"]  <- "2320" 
 data$ID[ data$ID =="H-01"]  <- "2220" 
 data$ID[ data$ID =="H-03"]  <- "2222" 
 data$ID[ data$ID =="H-04"]  <- "2223"
 data$ID[ data$ID =="I-20"]  <- "2362"
 data$ID[ data$ID =="419"]  <- "2385"
 data$ID[ data$ID =="413"]  <- "2386"
 data$ID[ data$ID =="W212"]  <- "2388"
 data$ID[ data$ID =="318"]  <- "2390"
 data$ID[data$ID == "54C"]  <- "2453"
 data$ID[data$ID == "65B"]  <- "2528"
 data$ID[data$ID == "51C"]  <- "2449"
 data$ID[data$ID == "71B"]  <- "2505"
 data$ID[data$ID == "64D"]  <- "2518"
 data$ID[data$ID == "53B"]  <- "2536"
 data$ID[data$ID == "60A"]  <- "2559"
 data$ID[data$ID == "H-01"]  <- "2220"
 data$ID[data$ID == "H-03"]  <- "2222"
 data$ID[data$ID == "419"]  <- "2385"
 data$ID[data$ID == "78-02"]  <- "2140"

 
 data$ID <- as.numeric(data$ID)
  data <- data %>%
         mutate(polar= replace(polar, polar == 0, 0.6)) %>%  #replacing 0 cm growth for the size of plug that is 0.6
    mutate(ecuatorial= replace(ecuatorial, ecuatorial == 0, 0.6)) %>% #replacing 0 cm growth for the size of plug that is 0.6
    # group_by(ID, experimental_replicate, dose, ecuatorial, polar, repeats) %>%
    mutate(growth = ((ecuatorial + polar) / 2)) %>% 
  group_by(ID, dose) %>%
  mutate(growth_range = list(get_range(growth))) %>%
  unnest(cols = c(growth_range)) %>% 
  filter(growth <= upper & growth >= lower) %>%
  dplyr::rename(response = growth) %>%
  ungroup() %>%
  select(c(ID, experimental_replicate, repeats, dose, response))
    }


##Function.Reading data for Discriminatory concentration
reading_data_discriminatory_concentration <- function(filename) {
  data <- filename
  # data <- subset(data, select = -X)
  #data$repeats <- rep_len(1:3, length.out = nrow(data))
  # renaming ID for the current names
  data$ID[data$ID == "12B"]  <- "12"
  #Although 74SS1 belong to # 74, the one what I used was specifically this:74SS1
  data$ID[data$ID == "74SS1"]  <- "74"
  data$ID[data$ID == "62-02"]  <- "2098"
  data$ID[data$ID == "62-03"]  <- "2099"
  data$ID[data$ID == "62-04"]  <- "2100"
  data$ID[data$ID == "78-01"]  <- "2139"
  data$ID[data$ID == "78-02"]  <- "2140"
  data$ID[data$ID == "78-05"]  <- "2143"
  data$ID[data$ID == "S-01"]  <- "2320"
  data$ID[data$ID == "H-01"]  <- "2220"
  data$ID[data$ID == "H-03"]  <- "2222"
  data$ID[data$ID == "H-04"]  <- "2223"
  data$ID[data$ID == "I-20"]  <- "2362"
  data$ID[data$ID == "419"]  <- "2385"
  data$ID[data$ID == "413"]  <- "2386"
  data$ID[data$ID == "W212"]  <- "2388"
  data$ID[data$ID == "318"]  <- "2390"
  data$ID[data$ID == "54C"]  <- "2453"
  data$ID[data$ID == "65B"]  <- "2528"
  data$ID[data$ID == "51C"]  <- "2449"
  data$ID[data$ID == "71B"]  <- "2505"
  data$ID[data$ID == "64D"]  <- "2518"
  data$ID[data$ID == "53B"]  <- "2536"
  data$ID[data$ID == "60A"]  <- "2559"
  data$ID[data$ID == "H-01"]  <- "2220"
  data$ID[data$ID == "H-03"]  <- "2222"
  data$ID[data$ID == "419"]  <- "2385"
  data$ID[data$ID == "78-02"]  <- "2140"
  data$ID <- as.numeric(data$ID)
  data <- data %>%
    mutate(response = replace(response, response == 0, 0.6)) %>%
    #replacing 0 cm growth for the size of plug that is 0.6
    # group_by(ID, experimental_replicate, concentration, growth, repeats) %>%
    group_by(ID) %>%
    mutate(response_range = list(get_range(response))) %>%
    unnest() %>%
    mutate(control_range = list(get_range(control))) %>%
    unnest() %>%
    filter(response <= upper &
             response >= lower,
           control <= upper1 & control >= lower1) %>%
    ungroup() 
}

#getting EC50
getting_EC50 <- function(filename){
  
  getting.EC50 <- EC_table(filename, form = response ~ dose)
getting.EC50 <- getting.EC50 %>% 
  dplyr::rename(ID = sample ) %>% # renaming
 # mutate(ID = as.factor(ID)) %>% 
  dplyr::rename(EC50 = Estimate.50) %>%
   group_by(ID, EC50) %>%
  mutate(source = ifelse(
    ID %in% baseline.isolates.2,# if they are in this group
    "baseline",
    ifelse(
      ID %in% survey.isolates, # if they are in this group, that is named farmer fiels also ahead
      "survey_isolates",
      "treatmentyear2016_isolates" # if they are in this group, that is named fungicide field trials also ahead
      )), Host = ifelse(
    ID %in% drybean,
    "Drybean", ifelse(
      ID %in% soybean,
      "Soybean","different_host"))) %>% ungroup() %>% 
  mutate(
  source = as.factor(source),
  Host = as.factor(Host)) %>%
#taking out "different host" becuse does not belong to drybean nor soybean
  filter(!ID == 8, !ID == 129) %>% 
      mutate(source = recode(source, survey_isolates = "Farmer fields", treatmentyear2016_isolates = "Fungicide field trials", baseline = "Baseline")) %>% mutate(source = as.factor(source)) %>% 
  ungroup()
}
##
getting_DC <- function(filename){
   data <-filename
  
  data1<- data %>% group_by(ID, dose) %>% summarise(mean_response = mean(response, na.rm = TRUE)) %>%
    ungroup() %>% 
  spread(dose, mean_response) %>% 
    dplyr::rename(control = `0`) %>% #taking out "different host" becuse does not belong to drybean nor soybean
  filter(!ID == 8, !ID == 129) 
    
    
  data2 <- as_data_frame( lapply(data1,function(x) {x[1] ; x}))
    colnames(data2)[-c(1:2)] <- paste("RG_", colnames(data2[,-c(1:2)]), sep = "")
  data3 <- data2 %>% mutate_each(funs((./control)*100), starts_with("RG_")) %>% 
    select(-c(control)) 
} 
##
  getting_DC_2 <- function(filename){
   data <-filename
  
data1 <-data %>% #joining
  # normalizing data
  mutate(logEC50 = (log(EC50))) #%>% 
  }

  lineal_model <- function (vardep, varindep, data) {
    lm(paste(vardep, "~", varindep), data = data)
  }

  wrangling_DC<- function(filename){
   data <-filename

data1 <-data %>%
 mutate(source = ifelse(
    ID %in% baseline.isolates.2,
    "Baseline",
    ifelse(
      ID %in% Farmer.fields,
      "Producer Fields","Fungicide Field Trials"
    )), Host = ifelse(
    ID %in% drybean,
    "Drybean", ifelse(
      ID %in% soybean,
      "Soybean","different_host")))  %>%
  mutate(
  source = as.factor(source),
  Host = as.factor(Host)) %>%
#taking out "different host" becuse does not belong to drybean nor soybean
  filter(!ID == 8, !ID == 129) %>%
  mutate(
    country = ifelse(
      ID %in% mexican.isolates,
      "Mexico",
      ifelse(ID %in% brazilian.isolates,
             "Brazil", "USA")))%>% 
 mutate(Field_Year= case_when(
ID %in% NA_year_1973~"NA_year_1973",
ID %in% NA_year_1974~"NA_year_1974",
ID %in% NA_year_1963~"NA_year_1963",
ID %in% NA_year_1975~"NA_year_1975",
ID %in% NA_year_1977~"NA_year_1977",
ID %in% NA_year_1976~"NA_year_1976",
ID %in% NA_year_1978~"NA_year_1978",
ID %in% NA_year_1979~"NA_year_1979",
ID %in% NA_year_1980~"NA_year_1980",
ID %in% NA_year_2003~"NA_year_2003",
ID %in% NA_year_1982~"NA_year_1982",
ID %in% NA_year_1983~"NA_year_1983",
ID %in% NA_year_1984~"NA_year_1984",
ID %in% NA_year_1985~"NA_year_1985",
ID %in% NA_year_1987~"NA_year_1987",
ID %in% NA_year_1990~"NA_year_1990",
ID %in% NA_year_1965~"NA_year_1965",
ID %in% NA_year_1989~"NA_year_1989",
ID %in% NA_year_1991~"NA_year_1991",
ID %in% NA_year_1992~"NA_year_1992",
ID %in% NA_year_1969~"NA_year_1969",
ID %in% NA_year_1981~"NA_year_1981",
ID %in% NA_year_1993~"NA_year_1993",
ID %in% NA_year_1994~"NA_year_1994",
ID %in% NA_year_1995~"NA_year_1995",
ID %in% NA_year_1996~"NA_year_1996",
ID %in% NA_year_1997~"NA_year_1997",
ID %in% WM_production_31_year_1997~"WM_production_31_year_1997",
ID %in% WM_production_32_year_1997~"WM_production_32_year_1997",
ID %in% NA_year_1988~"NA_year_1988",
ID %in% WM_production_33_year_1997~"WM_production_33_year_1997",
ID %in% NA_year_NA~"NA_year_NA",
ID %in% NA_year_1998~"NA_year_1998",
ID %in% NA_year_1999~"NA_year_1999",
ID %in% NA_year_2000~"NA_year_2000",
ID %in% NA_year_2002~"NA_year_2002",
ID %in% WM_nursery_1_year_2003~"WM_nursery_1_year_2003",
ID %in% WM_nursery_2_year_2003~"WM_nursery_2_year_2003",
ID %in% WM_nursery_3_year_2003~"WM_nursery_3_year_2003",
ID %in% WM_nursery_4_year_2003~"WM_nursery_4_year_2003",
ID %in% WM_nursery_5_year_2003~"WM_nursery_5_year_2003",
ID %in% WM_production_35_year_2003~"WM_production_35_year_2003",
ID %in% NA_year_2004~"NA_year_2004",
ID %in% NA_year_2005~"NA_year_2005",
ID %in% NA_year_2006~"NA_year_2006",
ID %in% WM_production_36_year_2004~"WM_production_36_year_2004",
ID %in% WM_production_37_year_2004~"WM_production_37_year_2004",
ID %in% WM_production_38_year_2004~"WM_production_38_year_2004",
ID %in% WM_production_39_year_2004~"WM_production_39_year_2004",
ID %in% WM_production_40_year_2004~"WM_production_40_year_2004",
ID %in% WM_production_41_year_2004~"WM_production_41_year_2004",
ID %in% WM_production_42_year_2004~"WM_production_42_year_2004",
ID %in% WM_production_43_year_2004~"WM_production_43_year_2004",
ID %in% WM_production_44_year_2004~"WM_production_44_year_2004",
ID %in% WM_nursery_6_year_2004~"WM_nursery_6_year_2004",
ID %in% WM_nursery_7_year_2004~"WM_nursery_7_year_2004",
ID %in% WM_nursery_8_year_2004~"WM_nursery_8_year_2004",
ID %in% WM_nursery_9_year_2004~"WM_nursery_9_year_2004",
ID %in% WM_nursery_10_year_2004~"WM_nursery_10_year_2004",
ID %in% WM_nursery_11_year_2004~"WM_nursery_11_year_2004",
ID %in% WM_nursery_12_year_2004~"WM_nursery_12_year_2004",
ID %in% WM_nursery_24_year_2005~"WM_nursery_24_year_2005",
ID %in% WM_nursery_13_year_2005~"WM_nursery_13_year_2005",
ID %in% WM_nursery_14_year_2005~"WM_nursery_14_year_2005",
ID %in% WM_nursery_15_year_2005~"WM_nursery_15_year_2005",
ID %in% WM_nursery_16_year_2005~"WM_nursery_16_year_2005",
ID %in% WM_nursery_17_year_2005~"WM_nursery_17_year_2005",
ID %in% WM_nursery_18_year_2005~"WM_nursery_18_year_2005",
ID %in% WM_nursery_19_year_2006~"WM_nursery_19_year_2006",
ID %in% NA_year_2007~"NA_year_2007",
ID %in% WM_production_1_year_2007~"WM_production_1_year_2007",
ID %in% WM_production_2_year_2007~"WM_production_2_year_2007",
ID %in% WM_production_3_year_2007~"WM_production_3_year_2007",
ID %in% WM_production_4_year_2007~"WM_production_4_year_2007",
ID %in% WM_production_5_year_2007~"WM_production_5_year_2007",
ID %in% WM_production_6_year_2007~"WM_production_6_year_2007",
ID %in% WM_production_7_year_2007~"WM_production_7_year_2007",
ID %in% WM_production_8_year_2007~"WM_production_8_year_2007",
ID %in% WM_production_9_year_2007~"WM_production_9_year_2007",
ID %in% NA_year_2008~"NA_year_2008",
ID %in% WM_production_10_year_2008~"WM_production_10_year_2008",
ID %in% WM_nursery_20_year_2008~"WM_nursery_20_year_2008",
ID %in% WM_production_11_year_2008~"WM_production_11_year_2008",
ID %in% WM_production_12_year_2008~"WM_production_12_year_2008",
ID %in% WM_production_13_year_2008~"WM_production_13_year_2008",
ID %in% NA_year_2009~"NA_year_2009",
ID %in% WM_nursery_21_year_2008~"WM_nursery_21_year_2008",
ID %in% WM_production_38_year_2009~"WM_production_38_year_2009",
ID %in% WM_nursery_22_year_2009~"WM_nursery_22_year_2009",
ID %in% WM_nursery_23_year_2009~"WM_nursery_23_year_2009",
ID %in% WM_production_14_year_2009~"WM_production_14_year_2009",
ID %in% WM_production_15_year_2009~"WM_production_15_year_2009",
ID %in% WM_production_16_year_2009~"WM_production_16_year_2009",
ID %in% WM_production_17_year_2009~"WM_production_17_year_2009",
ID %in% WM_production_18_year_2009~"WM_production_18_year_2009",
ID %in% WM_production_19_year_2009~"WM_production_19_year_2009",
ID %in% WM_production_20_year_2009~"WM_production_20_year_2009",
ID %in% WM_production_21_year_2010~"WM_production_21_year_2010",
ID %in% WM_production_22_year_2010~"WM_production_22_year_2010",
ID %in% WM_production_23_year_2010~"WM_production_23_year_2010",
ID %in% WM_production_24_year_2010~"WM_production_24_year_2010",
ID %in% WM_production_25_year_2010~"WM_production_25_year_2010",
ID %in% WM_production_26_year_2010~"WM_production_26_year_2010",
ID %in% WM_production_27_year_2010~"WM_production_27_year_2010",
ID %in% WM_production_28_year_2010~"WM_production_28_year_2010",
ID %in% WM_production_29_year_2010~"WM_production_29_year_2010",
ID %in% WM_production_30_year_2010~"WM_production_30_year_2010",
ID %in% NA_year_2011~"NA_year_2011",
ID %in% NA_year_2012~"NA_year_2012",
ID %in% NA_year_2013~"NA_year_2013",
ID %in% NA_year_2014~"NA_year_2014",
ID %in% NA_year_2016~"NA_year_2014",
ID %in% po_year_2016~"NA_year_2016",
ID %in% flo_year_2016~"flo_year_2016",
ID %in% lan_year_2016~"lan_year_2016",
ID %in% ho_year_2017~"ho_year_2017",
ID %in% an_year_2017~"an_year_2017",
ID %in% dod_year_2017~"dod_year_2017",
ID %in% va_year_2017~"va_year_2017",
ID %in% na_year_2017~"na_year_2017",
ID %in% sa_year_2017~"sa_year_2017",
ID %in% han_year_2017~"han_year_2017",
ID %in% cu_year_2017~"cu_year_2017",
ID %in% mon_year_2017~"mon_year_2017",
ID %in% ant_year_2017~"ant_year_2017",
ID %in% holt_year_2017~"holt_year_2017",
ID %in% arr_year_2017~"arr_year_2017",
ID %in% ini_year_2017~"ini_year_2017",
ID %in% ley_year_2017~"ley_year_2017",
ID %in% pre_year_2017~"pre_year_2017",
ID %in% ru_year_2017~"ru_year_2017",
ID %in% ca_year_2017~"ca_year_2017",
ID %in% st_year_2017~"st_year_2017",
ID %in% Al_year_2015~"Al_year_2015",
ID %in% Hi_year_2015~"Hi_year_2015",
ID %in% In_year_2015~"In_year_2015",
ID %in% mont_year_2015~"mont_year_2015",
ID %in% Sand_year_2015~"Sand_year_2015",
ID %in% Sani_year_2015~"Sani_year_2015",
ID %in% St_year_2015~"St_year_2015",
ID %in% wau_year_2015~"wau_year_2015",
ID %in% la_year_2015~"la_year_2015",
ID %in% co_year_2015~"co_year_2015",
ID %in% cern_year_2015~"cern_year_2015",
ID %in% cer_year_2015~"cer_year_2015",
ID %in% NA_year_2010~"NA_year_2010",
ID %in% NA_year_2017~"NA_year_2017",
TRUE ~ "no specific field"),
  State= case_when(
ID %in% NE~"NE",
ID %in% CO~"CO",
ID %in% IDA~"ID",
ID %in% CA~"CA",
ID %in% NY~"NY",
ID %in% NA~"NA",
ID %in% WA~"WA",
ID %in% NC~"NC",
ID %in% AZ~"AZ",
ID %in% MT~"MT",
ID %in% WI~"WI",
ID %in% MI~"MI",
ID %in% ON~"ON",
ID %in% AB~"AB",
ID %in% KS~"KS",
ID %in% ND~"ND",
ID %in% MD~"MD",
ID %in% Arusha~"Arusha",
ID %in% SK~"SK",
ID %in% OK~"OK",
ID %in% DF~"DF",
ID %in% OH~"OH",
ID %in% SIN~"SIN",
ID %in% IA~"IA",
ID %in% PA~"PA",
ID %in% Herts~"Herts",
ID %in% Unknown~"Unknown",
ID %in% FL~"FL",
ID %in% GA~"GA",
ID %in% MS~"MS",
ID %in% LA~"LA",
ID %in% KY~"KY",
ID %in% OR~"OR",
ID %in% MG~"MG",
ID %in% GO~"GO",
ID %in% Herverlee~"Herverlee",
# ID %in% Mooi Rive~"Mooi Rive",
ID %in% Marathon~"Marathon",
ID %in% SC~"SC",
ID %in% MN~"MN",
ID %in% Blenheim~"Blenheim",
ID %in% Barranquita~"Barranquita",
ID %in% Victoria~"Victoria",
ID %in% Tasmania~"Tasmania",
ID %in% DE~"DE",
ID %in% Guerbigny~"Guerbigny",
ID %in% QUE~"QUE",
ID %in% Kigali~"Kigali",
ID %in% birdseed~"birdseed",
ID %in% Kessenich~"Kessenich",
ID %in% Elen~"Elen",
ID %in% Moshi~"Moshi",
ID %in% Peronne~"Peronne",
ID %in% Loudéac~"Loudéac",
# ID %in% Bordeaux-Hourtin~"Bordeaux-Hourtin",
ID %in% BA~"BA",
ID %in% RS~"RS",
ID %in% MO~"MO",
ID %in% PR~"PR",
TRUE ~ "no specific state"
), 
County= case_when(
ID %in% Lambari~"Lambari",
ID %in% `Porto Firme`~"Porto Firme",
ID %in% `Unaí`~"Unai",
TRUE ~ "no specific county"
)) %>% 
  mutate(
    Field_Year =  ifelse(
      County == "Lambari",
      "Lambari",
      ifelse(
        County == "Porto Firme",
        "Porto Firme",
        ifelse(County == "Unai", "Unai", Field_Year
        )))) %>% mutate(
Field_Year = recode(
      Field_Year,
      WM_production_21_year_2010 = "NE.10.fld21.DB",
      WM_production_1_year_2007 = "ND.07.fld01.DB",
      WM_production_9_year_2007 = "WA.07.fld09.DB",
      lan_year_2016 = "MI.16.East.Lansing.SB",
      WM_production_10_year_2008 = "WA.08.fld10.DB",
      WM_production_25_year_2010 = "ND.10.fld25.DB",
      mont_year_2015= "MI.15.Montcalm.SB",
      na_year_2017 = "IA.17.Nashua.SB",
      mon_year_2017 = "MI.17.Montcalm.SB",
      Baseline = "Baseline",
      han_year_2017 = "WI.17.Hancock.SB",
      NA_year_2010 = "BR.2010.fld01.DB",
      pre_year_2017 = "MX.2017.pre.DB",
      ca_year_2017 = "MX.2017.ca.DB",
      ley_year_2017 = "MX.2017.ley.DB",
      cer_year_2015 = "WI.2015.cer.SB",
      `Porto Firme` = "BR.2010.por.DB",
      `Unai` = "BR.2010.unai.DB",
      `Lambari` = "BR.2010.lam.DB"
        )) 
  
  }

  #NEED TO MODIFI HERE WHEN BRAZIL, MAKE BT COUNTY NOR EVEN STATE
##wrangling survey to get the EC50DC
  wrangling_survey <- function(filename) {
    data <- filename
    data1 <- data %>%
      group_by(ID) %>% 
      summarise(
        mean_response = mean(response, na.rm = TRUE),
        mean_control = mean(control, na.rm = TRUE)
      ) %>% ungroup() %>%
      mutate(RG = (mean_response / mean_control) * 100)
  }

   wrangling_DC_2<- function(filename){
   data <-filename
   data1 <-data %>%
    group_by(ID) %>% #keeping the highest value of each ID value
  top_n(1, EC50DC) %>% 
  ungroup()
   }
  
  wrangling_DC_3<- function(filename){
   data <-filename
   data1 <-data %>%
   mutate(Field_Year =  ifelse(source == "Baseline", "Baseline", Field_Year),
    Field_Year = as.factor(Field_Year)) %>% 
   
  group_by(Field_Year) %>% filter(n()>=9) %>% 
     ungroup() %>% 
     filter(!Field_Year =="no specific field" ) %>% 
     mutate(Field_Year = fct_reorder(Field_Year, EC50DC, mean))
}

   #test source
  
  test_kruskal_source <- function(filename, vardep, varindep) {
    data <- filename
    vardep <- data$EC50DC
    varindep <- data$source
    data1 <-  data %>%
      do(tidy(shapiro.test(vardep)))
    if (data1$p.value > 0.05)  {
      print ("normal")# I have to add the ANOVA
    } else {
      kruskal.test(as.numeric(vardep) ~ as.factor(varindep), data = data1)
         }}
  
  test_Dunn_source <- function(filename, vardep, varindep) {
    data <- filename
    vardep <- data$EC50DC
    varindep <- data$source
    data1 <-  data %>%
      do(tidy(shapiro.test(vardep)))
    if (data1$p.value > 0.05)  {
      print ("normal")# I have to add the ANOVA
    } else {
      kruskal.test(as.numeric(vardep) ~ as.factor(varindep), data = data1)
     if (data1$p.value > 0.05) {
        print ("No stats difference")
      } else{
        DunnTest(as.numeric(vardep) ~ as.factor(varindep), data = data1, method="bonferroni")}}}
  
   #test host

   test_kruskal_host <- function(filename, vardep, varindep) {
    data <- filename
    vardep <- data$EC50DC
    varindep <- data$Host
    data1 <-  data %>%
      do(tidy(shapiro.test(vardep))) 
      if (data1$p.value > 0.05)  {
      print ("normal")# I have to add the ANOVA
    } else {
      kruskal.test(as.numeric(vardep) ~ as.factor(varindep), data = data1)
      }}    
  
  test_Dunn_host <- function(filename, vardep, varindep) {
    data <- filename
    vardep <- data$EC50DC
    varindep <- data$Host
    data1 <-  data %>%
      do(tidy(shapiro.test(vardep)))
    if (data1$p.value > 0.05)  {
      print ("normal")# I have to add the ANOVA
    } else {
      kruskal.test(as.numeric(vardep) ~ as.factor(varindep), data = data1)
     if (data1$p.value > 0.05) {
        print ("No stats difference")
      } else{
        DunnTest(as.numeric(vardep) ~ as.factor(varindep), data = data1, method="bonferroni")}}}
  
  test_field <- function(filename, vardep, varindep) {
    data <- filename
    vardep <- data$EC50DC
    varindep <- data$Field_Year
    data1 <-  data %>%
      do(tidy(shapiro.test(vardep)))
    if (data1$p.value > 0.05)  {
      print ("normal")# I have to add the ANOVA
    } else {
      kruskal.test(as.numeric(vardep) ~ as.factor(varindep), data = data1)
     if (data1$p.value > 0.05) {
        print ("No stats difference")
     } 
      
      else{
      DunnTest(as.numeric(vardep) ~ as.factor(varindep), data = data1, method="bonferroni")}
    
    
    }}
  ###
      myplot_model_1_boscalid <- function(filename, RG, logEC50, lm)  { 
    data <- filename
    RG <- data$RG
    logEC50 <- data$logEC50
    Host <- data$Host
    source <- data$source
ggplot(data = filename,  aes(x = RG, y = logEC50))  +
  geom_point(aes(
    shape = Host,
    fill= source
  ),
  size = 2, stroke= 1, colour= "black")  + scale_shape_manual(values = c(21, 24))+ scale_fill_grey(start = 0, end = 1) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        colour = "black",
        size = 0.3
      ) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          family = "Arial"
        ),
        axis.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
        axis.text = element_text(
          face = "bold",
          size = 18,
          family = "Arial"
        ),
        panel.background = element_rect(fill = "white", colour = "grey50")
      ) + labs(x =
                 "Relative growth (%)", y = expression(bold(Log  (EC[bold("50")])))) + geom_abline(
                   intercept = 0,
                   slope = 1,
                   color = "black",
                   lty = "dashed",
                   size = 0.3
                  ) + expand_limits(x = c(15, 65)) + scale_x_continuous(name = waiver(), breaks = c(15, 25, 35, 45, 55, 65))  + geom_label(aes(x = 53,
                                                                                                                                                                                        y = -2.75)
                                                                                                                                                                                    ,
                                                                                                                                                                                    label = c(paste  (
                                                                                                                                                                                      " Y =",
                                                                                                                                                                                      paste0 (round(summary(lm)[[4]][2], 4), "x"),
                                                                                                                                                                                      "-",
                                                                                                                                                                                      round(summary(lm)[[4]][1], 3) *
                                                                                                                                                                                        -1,
                                                                                                                                                                                      paste  ("\n R =",
                                                                                                                                                                                              paste0 (round(
                                                                                                                                                                                                summary(lm)[[9]][1], 4
                                                                                                                                                                                              ), "\n p < 0.001"))
                                                                                                                                                                                    )),size = 3.5,
                                                                                                                                                                                    fontface = "bold")

 }
    ### tetraconazole
  
   myplot_model_1_tetraconazole <- function(filename, RG, logEC50, lm)  { 
    data <- filename
    RG <- data$RG
    logEC50 <- data$logEC50
    Host <- data$Host
    source <- data$source
ggplot(data = filename,  aes(x = RG, y = logEC50))  +
  geom_point(aes(
    shape = Host,
    fill= source
  ),
  size = 2, stroke= 1, colour= "black")  + scale_shape_manual(values = c(21, 24))+ scale_fill_grey(start = 0, end = 1) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        colour = "black",
        size = 0.3
      ) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          family = "Arial"
        ),
        axis.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
        axis.text = element_text(
          face = "bold",
          size = 18,
          family = "Arial"
        ),
        panel.background = element_rect(fill = "white", colour = "grey50")
      ) + labs(x =
                 "Relative growth (%)", y = expression(bold(Log  (EC[bold("50")])))) + geom_abline(
                   intercept = 0,
                   slope = 1,
                   color = "black",
                   lty = "dashed",
                   size = 0.3
                  ) + expand_limits(x = c(15, 65)) + scale_x_continuous(name = waiver(), breaks = c(15, 25, 35, 45, 55, 65))  + geom_label(aes(x = 53,
                                                                                                                                                                                        y = -0.75)
                                                                                                                                                                                    ,
                                                                                                                                                                                    label = c(paste  (
                                                                                                                                                                                      " Y =",
                                                                                                                                                                                      paste0 (round(summary(lm)[[4]][2], 4), "x"),
                                                                                                                                                                                      "-",
                                                                                                                                                                                      round(summary(lm)[[4]][1], 3) *
                                                                                                                                                                                        -1,
                                                                                                                                                                                      paste  ("\n R =",
                                                                                                                                                                                              paste0 (round(
                                                                                                                                                                                                summary(lm)[[9]][1], 4
                                                                                                                                                                                              ), "\n p < 0.001"))
                                                                                                                                                                                    )),size = 3.5,
                                                                                                                                                                                    fontface = "bold")

 }
  
   #### picoxystrobin
  
   myplot_model_1_picoxystrobin <- function(filename, RG, logEC50, lm)  { 
    data <- filename
    RG <- data$RG
    logEC50 <- data$logEC50
    Host <- data$Host
    source <- data$source
ggplot(data = filename,  aes(x = RG, y = logEC50))  +
  geom_point(aes(
    shape = Host,
    fill= source
  ),
  size = 2, stroke= 1, colour= "black")  + scale_shape_manual(values = c(21, 24))+ scale_fill_grey(start = 0, end = 1) +
      geom_smooth(
        method = "lm",
        se = FALSE,
        colour = "black",
        size = 0.3
      ) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          family = "Arial"
        ),
        axis.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5
        ),
        axis.text = element_text(
          face = "bold",
          size = 18,
          family = "Arial"
        ),
        axis.text.y = element_text(
                                                                                                                                                                                                                                angle = 20,
                                                                                                                                       hjust = 1),
        panel.background = element_rect(fill = "white", colour = "grey50")
      ) + labs(x =
                 "Relative growth (%)", y = expression(bold(Log  (EC[bold("50")])))) + geom_abline(
                   intercept = 0,
                   slope = 1,
                   color = "black",
                   lty = "dashed",
                   size = 0.3
                  ) + expand_limits(x = c(15, 65)) + scale_x_continuous(name = waiver(), breaks = c(15, 25, 35, 45, 55, 65))  + geom_label(aes(x = 28,
                                                                                                                                                                                        y = -4.30)
                                                                                                                                                                                    ,
                                                                                                                                                                                    label = c(paste  (
                                                                                                                                                                                      " Y =",
                                                                                                                                                                                      paste0 (round(summary(lm)[[4]][2], 4), "x"),
                                                                                                                                                                                      "-",
                                                                                                                                                                                      round(summary(lm)[[4]][1], 3) *
                                                                                                                                                                                        -1,
                                                                                                                                                                                      paste  ("\n R =",
                                                                                                                                                                                              paste0 (round(
                                                                                                                                                                                                summary(lm)[[9]][1], 4
                                                                                                                                                                                              ), "\n p < 0.001"))
                                                                                                                                                                                    )),size = 3.5,
                                                                                                                                                                                    fontface = "bold")

 }
  
  
    ####
   myplot_model_2 <- function(filename, EC50, Estimate.50DC, Host, source) {
     data <- filename
    EC50 <- data$EC50
    Estimate.50DC <- data$Estimate.50DC
    Host <- data$Host
    source <- data$source
    ggplot(data = filename,  aes(x = EC50, y = Estimate.50DC))  +
  geom_point(aes(
    shape = Host,
    fill= source
  ),
  size = 2, stroke= 1, colour= "black")  + scale_shape_manual(values = c(21, 24))+ scale_fill_grey(start = 0, end = 1)  +
   geom_smooth(method = "lm",
               se = FALSE,
               colour = "black",
               size = 0.3) +
   theme(
     plot.title = element_text(
       size = 18,
       face = "bold",
       hjust = 0.5,
       family = "Arial"
     ),
     axis.title = element_text(size = 18, face = "bold", hjust = 0.5),
     axis.text = element_text(
       face = "bold",
       size = 12,
       family = "Arial"
     ),
     panel.background = element_rect(fill = "white", colour = "grey50")
   ) +
   labs(x =
          expression(bold(EC[bold("50")]) ~ (bold(mg/L ~ bold(
            "a.i."
          )))), y = expression(bold(EC[bold("50") ~ (bold("D"))]) ~ (bold(mg/L ~ bold(
            "a.i."
          ))))) + geom_abline(
            intercept = 0,
            slope = 1,
            color = "black",
            lty = "dashed",
            size = 0.3)
   }
#label 
  generate_label <- function(HSD){
    # Extract labels and factor levels from Tukey post-hoc 
    Dunn.levels <- HSD[[1]][,2]
    Dunn.labels <- multcompLetters(Dunn.levels)['Letters']
    plot.labels <- names(Dunn.labels[['Letters']])
    plot.levels <- data.frame(plot.labels, Dunn.levels, labels = Dunn.labels[['Letters']],
                              stringsAsFactors = FALSE)
        return(plot.levels)
}
   ###
generate_label_2 <- function(HSD, flev, surveyy, highlabel){
  # Extract labels and factor levels from Tukey post-hoc 
  Dunn.levels <- HSD[[1]][,2]
  Dunn.labels <- multcompLetters(Dunn.levels)['Letters']
  plot.labels <- names(Dunn.labels[['Letters']])
  # Get highest quantile for Dunn's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- plyr:: ddply(surveyy, flev, function (x) max(fivenum(x$EC50DC)) + highlabel)
    # Create a data frame out of the factor levels and Dunn's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Dunn.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}
###
myplot_comparison <- function(filename, Field_Year,EC50DC, highlabel) {
  data <- filename
  Field_Year <- data$Field_Year
  EC50DC<- data$EC50DC
  
  ggplot(data = filename,  aes(x = Field_Year, y = EC50DC)) + geom_boxplot() +
    geom_text(data = generate_label_2(test_field(filename), 'Field_Year', filename,highlabel),
              aes(x = plot.labels, y = V1, label = labels), size= 3) +     labs(x = "Field", y = expression(bold(EC[bold("50") ~ (bold("D"))]) ~
                                                                                                     (bold(mg/L ~ bold(
                                                                                                       "a.i."
                                                                                                     ))))) + theme(
                                                                                                       panel.border = element_rect(
                                                                                                         colour = "black",
                                                                                                         fill = NA,
                                                                                                         size = 0.5
                                                                                                       ),
                                                                                                       axis.title = element_text(size = 11, face = "bold", hjust = 0.5),
                                                                                                       axis.text.x = element_text(
                                                                                                         face = "bold",
                                                                                                         size = 8,
                                                                                                         family = "Arial",
                                                                                                         angle = 30,
                                                                                                         hjust = 1
                                                                                                       ),
                                                                                                       axis.text.y = element_text(
                                                                                                         face = "bold",
                                                                                                         size = 8,
                                                                                                         family = "Arial"
                                                                                                       ),
                                                                                                       panel.background = element_rect(fill = "white", colour = "grey50")
                                                                                                     ) 
  
  
}

 myplot_baseline <- function(filename, EC50DC, source) {
  data <- filename
  EC50DC <- data$EC50DC
  source <- data$source
  
  data %>% filter(source == "Baseline") %>% mutate(
    Baseline = ifelse(
      ID %in% SI.Screening.Nursery.Field,
      # if they are in this group
      "Contemporary isolates*",
      "Historic isolates**"
    )
  ) %>%   ggplot(aes(x = Baseline, y = EC50DC, )) + geom_violin() + geom_jitter(aes(shape = Baseline), size = 6) +     labs(x = "Baseline", y = expression(bold(EC[bold("50") ~ (bold("D"))]) ~
                                                                                                                                                             (bold(mg/L ~ bold(
                                                                                                                                                               "a.i."
                                                                                                                                                             ))))) + theme(
                                                                                                                                                               panel.border = element_rect(
                                                                                                                                                                 colour = "black",
                                                                                                                                                                 fill = NA,
                                                                                                                                                                 size = 1
                                                                                                                                                               ),
                                                                                                                                                               axis.title = element_text(
                                                                                                                                                                 size = 18,
                                                                                                                                                                 face = "bold",
                                                                                                                                                                 hjust = 0.5
                                                                                                                                                               ),
                                                                                                                                                               axis.text.x = element_text(
                                                                                                                                                                 face = "bold",
                                                                                                                                                                 size = 10,
                                                                                                                                                                 family = "Arial",
                                                                                                                                                                 hjust = 0.5
                                                                                                                                                               ),
                                                                                                                                                               axis.text.y = element_text(
                                                                                                                                                                 face = "bold",
                                                                                                                                                                 size = 10,
                                                                                                                                                                 family = "Arial"
                                                                                                                                                               ),
                                                                                                                                                               panel.background = element_rect(fill = "white", colour = "grey50")
                                                                                                                                                             )
}
```

# Fungicides:

**1. Thiophanate methyl**  
**2. Tetraconazole**  
**3. Boscalid**  
**4. Picoxystrobin**

``` r
#TM
hey <-
  reading_data_serial_dilution("data/TM_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv") %>%
  filter(!experimental_replicate == 1 |
           !dose == 1.5 | !repeats == 3)  #something weird happened with the repeats 3 so Iam revoming those
#Check here because there are more repeats for control, in that case create a new variable or column
## First try with postdoc Thomas
survey.TM.thomas <- read.csv("data/DD_survey-WMNCSRP.csv") %>%
  select(
    ID,
    repeats,
    experimental_replicate,
    TM.ecuatorial,
    TM.polar,
    control.ecuatorial,
    control.polar
  ) %>%
  mutate(TM.growth = ((TM.ecuatorial + TM.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = TM.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


## Second try with Edgar
survey.TM.USA <- 
      read.csv("data/DD_survey-WMNCSRP_TM.csv") %>% 
mutate(TM.growth = ((TM.ecuatorial + TM.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = TM.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()

## Third try with Cristian
survey.TM.more.USA <- 
 read.csv("data/TM.Cristian.csv") %>% 
mutate(TM.growth = ((TM.ecuatorial + TM.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = TM.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()

#Tetraconazole
hey2 <- reading_data_serial_dilution("data/tetraconazole_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
## First try with postdoc Thomas
survey.tetraconazole.thomas <- read.csv("data/DD_survey-WMNCSRP.csv") %>%
  select(
    ID,
    repeats,
    experimental_replicate,
    tetraconazole.ecuatorial,
    tetraconazole.polar,
    control.ecuatorial,
    control.polar
  ) %>%
  mutate(tetraconazole.growth = ((tetraconazole.ecuatorial + tetraconazole.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = tetraconazole.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()

## Second try with Edgar
survey.tetraconazole.USA <- 
      read.csv("data/DD_survey-WMNCSRP_tetraconazole.csv") %>% 
mutate(tetraconazole.growth = ((tetraconazole.ecuatorial + tetraconazole.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = tetraconazole.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


## Third try with Cristian
survey.tetraconazole.more.USA <- 
 read.csv("data/tetraconazole.Cristian.csv") %>% 
mutate(tetraconazole.growth = ((tetraconazole.ecuatorial + tetraconazole.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = tetraconazole.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>%  
reading_data_discriminatory_concentration()

## Fourth try with Edgar 
survey.tetraconazole.MX.BZ <- 
 read.csv("data/DD_survey-WMNCSRP_tetraconazole.Mexican.Brazilian.csv") %>% mutate(tetraconazole.growth = ((tetraconazole.ecuatorial + tetraconazole.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = tetraconazole.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


#Boscalid
hey3 <- reading_data_serial_dilution("data/boscalid_serialdilution_baseline(validationdata).csv") %>% filter(!dose ==  0.4 & !ID == 800)
## First try with postdoc Thomas
survey.boscalid.thomas <- read.csv("data/DD_survey-WMNCSRP.csv") %>%
  select(
    ID,
    repeats,
    experimental_replicate,
    boscalid.ecuatorial,
    boscalid.polar,
    control.ecuatorial,
    control.polar
  ) %>%
  mutate(boscalid.growth = ((boscalid.ecuatorial + boscalid.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = boscalid.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


   ## Second try with Edgar
survey.boscalid.USA <- 
      read.csv("data/DD_survey-WMNCSRP_boscalid.csv") %>% 
mutate(boscalid.growth = ((boscalid.ecuatorial + boscalid.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = boscalid.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


 ## Third try with Cristian
survey.boscalid.more.USA <- 
 read.csv("data/Boscalid.Cristian.csv") %>% 
mutate(boscalid.growth = ((boscalid.ecuatorial + boscalid.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = boscalid.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


## Fourth try with Edgar 
survey.boscalid.MX.BZ <- 
 read.csv("data/DD_survey-WMNCSRP_boscalid.Brazilian.Mexican.csv") %>% mutate(boscalid.growth = ((boscalid.ecuatorial + boscalid.polar) / 2),
         control.growth = ((control.ecuatorial + control.polar) / 2)) %>% 
dplyr::rename(response = boscalid.growth, control = control.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


# #Picoxystrobin
hey4 <- reading_data_serial_dilution("data/picoxystrobin_serialdi_BL-WMN_survey-WMNCSRP(validationdata).csv")
## First try with postdoc Thomas
survey.picoxystrobin.thomas <- read.csv("data/DD_survey-WMNCSRP.csv") %>%
  select(
    ID,
    repeats,
    experimental_replicate,
    picoxystrobin.ecuatorial,
    picoxystrobin.polar,
     sham.ecuatorial,
    sham.polar
  ) %>%
  mutate(picoxystrobin.growth = ((picoxystrobin.ecuatorial + picoxystrobin.polar) / 2),
         sham.growth = ((sham.ecuatorial + sham.polar) / 2)) %>% 
dplyr::rename(response = picoxystrobin.growth, control = sham.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
reading_data_discriminatory_concentration()
 
## Second try with Edgar
survey.picoxystrobin.USA <- 
      read.csv("data/DD_survey-WMNCSRP_picoxystrobin.csv") %>% 
mutate(picoxystrobin.growth = ((picoxystrobin.ecuatorial + picoxystrobin.polar) / 2),
         sham.growth = ((sham.ecuatorial + sham.polar) / 2)) %>% 
dplyr::rename(response = picoxystrobin.growth, control = sham.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
  reading_data_discriminatory_concentration()


## Third try with Cristian
survey.picoxystrobin.more.USA <- 
 read.csv("data/picoxystrobin.Cristian.csv") %>% 
mutate(picoxystrobin.growth = ((picoxystrobin.ecuatorial + picoxystrobin.polar) / 2),
         sham.growth = ((sham.ecuatorial + sham.polar) / 2)) %>% 
dplyr::rename(response = picoxystrobin.growth, control = sham.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
 reading_data_discriminatory_concentration()

## Fourth try with Edgar 
survey.picoxystrobin.MX.BZ <- 
 read.csv("data/DD_survey-WMNCSRP_picoxystrobin.Mexican.Brazilian.csv") %>% mutate(picoxystrobin.growth = ((picoxystrobin.ecuatorial + picoxystrobin.polar) / 2),
         sham.growth = ((sham.ecuatorial + sham.polar) / 2)) %>% 
dplyr::rename(response = picoxystrobin.growth, control = sham.growth) %>% 
  select(c(ID, experimental_replicate, repeats, response, control)) %>% 
 reading_data_discriminatory_concentration()
```

# Fungicide sensitivity by fungicide

## Serial dilution

Fungicide sensitivity was determined from 22 baseline isolates plus 21
additional isolates from FF and FFT randomly selected. Serial dilution
by fungicide eas carried out from 5-6 concentrations with 4 repetitions
per isolate and 2 experimental replications. Plates were incubated in
darkness at 23 ± 2 °C for 42 h and diameter measured with digital
calipers. A dose-response curve was fit to estimate the EC50 for each
fungicide.

``` r
# TM
TM.EC50 <-  getting_EC50(hey)
```

![](READM_files/figure-gfm/EC50%20serial%20dilution-1.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-2.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-3.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-4.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-5.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-6.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-7.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-8.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-9.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-10.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-11.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-12.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-13.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-14.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-15.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-16.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-17.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-18.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-19.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-20.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-21.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-22.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-23.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-24.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-25.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-26.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-27.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-28.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-29.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-30.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-31.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-32.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-33.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-34.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-35.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-36.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-37.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-38.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-39.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-40.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-41.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-42.png)<!-- -->

``` r
# Tetraconazole
tetraconazole.EC50 <- getting_EC50(hey2)
```

![](READM_files/figure-gfm/EC50%20serial%20dilution-43.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-44.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-45.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-46.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-47.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-48.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-49.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-50.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-51.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-52.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-53.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-54.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-55.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-56.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-57.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-58.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-59.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-60.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-61.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-62.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-63.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-64.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-65.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-66.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-67.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-68.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-69.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-70.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-71.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-72.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-73.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-74.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-75.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-76.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-77.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-78.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-79.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-80.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-81.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-82.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-83.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-84.png)<!-- -->

``` r
# Boscalid
boscalid.EC50 <- getting_EC50(hey3)
```

![](READM_files/figure-gfm/EC50%20serial%20dilution-85.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-86.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-87.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-88.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-89.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-90.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-91.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-92.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-93.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-94.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-95.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-96.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-97.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-98.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-99.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-100.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-101.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-102.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-103.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-104.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-105.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-106.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-107.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-108.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-109.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-110.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-111.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-112.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-113.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-114.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-115.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-116.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-117.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-118.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-119.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-120.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-121.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-122.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-123.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-124.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-125.png)<!-- -->

``` r
#Picoxystrobin
picoxystrobin.EC50 <- getting_EC50(hey4)
```

![](READM_files/figure-gfm/EC50%20serial%20dilution-126.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-127.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-128.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-129.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-130.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-131.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-132.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-133.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-134.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-135.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-136.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-137.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-138.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-139.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-140.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-141.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-142.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-143.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-144.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-145.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-146.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-147.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-148.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-149.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-150.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-151.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-152.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-153.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-154.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-155.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-156.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-157.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-158.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-159.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-160.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-161.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-162.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-163.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-164.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-165.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-166.png)<!-- -->![](READM_files/figure-gfm/EC50%20serial%20dilution-167.png)<!-- -->

## Procedure to get the Discriminatory Concentration (DC)

We identified the concentration with the best prediction of EC50 for
each fungicide, known as discriminatory concentration (DC), by linear
regression of **% mycelial growth vs. log(EC50)**. The concentration
yielding the **highest coefficient of determination (r2)** was selected
as the DC for each fungicide. For the remaining isolates among the 207,
growth on the DC was used to estimate EC50, which is termed EC50(D).

``` r
#TM
  #getting_EC50 releases ID with factor
TM <- getting_DC(hey) 
```

    ## `summarise()` has grouped output by 'ID'. You can override using the `.groups`
    ## argument.

``` r
  TM.joined <-  TM %>% left_join(TM.EC50)
```

    ## Joining, by = "ID"

``` r
TM.joined.2 <-  getting_DC_2(TM.joined) 
  
#Boscalid
boscalid <- getting_DC(hey3) 
```

    ## `summarise()` has grouped output by 'ID'. You can override using the `.groups`
    ## argument.

``` r
  boscalid.joined <-  boscalid%>% left_join(boscalid.EC50)
```

    ## Joining, by = "ID"

``` r
 boscalid.joined.2 <-  getting_DC_2(boscalid.joined) 
 ###Dose chosen as DC. Although there is one dose with higher r2, there is no difference between and this was selected
# Linear model of log EC50 and relative growth at 0.2 ppm, check normality and homogeneity of variances
finalRG_0.2 <- lm(logEC50 ~ RG_0.2,  boscalid.joined.2)
summary(finalRG_0.2)
```

    ## 
    ## Call:
    ## lm(formula = logEC50 ~ RG_0.2, data = boscalid.joined.2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.27096 -0.02645  0.01943  0.05137  0.25614 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.338863   0.128376  -33.80   <2e-16 ***
    ## RG_0.2       0.054589   0.003376   16.17   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1056 on 37 degrees of freedom
    ## Multiple R-squared:  0.8761, Adjusted R-squared:  0.8727 
    ## F-statistic: 261.5 on 1 and 37 DF,  p-value: < 2.2e-16

``` r
#plot(finalRG0.2)

##taking out outliers from Q-Q PLOT $ LEVERAGE TEST, ISOLATES # 1032 from treatment, 2222 & 2223 from Farmer fields##
boscalid.joined.2.clean <-
  boscalid.joined.2 %>% filter(!ID == "1032", !ID == "2222",!ID == "2223") %>% 
  select(ID, logEC50, RG_0.2, EC50) %>% 
  dplyr::rename(RG = RG_0.2) %>% 
   wrangling_DC()

summary(boscalid.lm <- lineal_model("logEC50","RG",data =  boscalid.joined.2.clean))
```

    ## 
    ## Call:
    ## lm(formula = paste(vardep, "~", varindep), data = data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.22927 -0.03254  0.01107  0.03859  0.16920 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -4.412096   0.105765  -41.72   <2e-16 ***
    ## RG           0.056687   0.002757   20.56   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08019 on 34 degrees of freedom
    ## Multiple R-squared:  0.9256, Adjusted R-squared:  0.9234 
    ## F-statistic: 422.8 on 1 and 34 DF,  p-value: < 2.2e-16

``` r
   # Getting the EC50DC according to the model
boscalid.EC50DC <- boscalid.joined.2.clean %>% mutate(Estimate.50DC =
                                     #using the model the intercept and the coefficient of the model
                                    exp(boscalid.lm[[1]][[1]] + boscalid.lm[[1]][[2]] * RG)) %>%# here I need to do this dynamic
  select(ID, Estimate.50DC, EC50)  %>%  
  wrangling_DC()

# Linear Regression graph of log EC50 and relative growth at 0.2 ppm
###

myplot_model_1_boscalid(boscalid.joined.2.clean, lm = boscalid.lm) 
```

    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Getting%20Discriminatory%20dose-1.png)<!-- -->

``` r
######### Linear model of EC50  and EC50DC just to corroborate is a GOOD RELATIONSHIP, check normality and homogeneity of variances

summary(boscalid.lm.2 <- lineal_model("Estimate.50DC","EC50",data =  boscalid.EC50DC))
```

    ## 
    ## Call:
    ## lm(formula = paste(vardep, "~", varindep), data = data)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -0.0189959 -0.0045486 -0.0005029  0.0039102  0.0234536 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 0.014242   0.005271   2.702   0.0107 *  
    ## EC50        0.865738   0.046429  18.646   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.008695 on 34 degrees of freedom
    ## Multiple R-squared:  0.9109, Adjusted R-squared:  0.9083 
    ## F-statistic: 347.7 on 1 and 34 DF,  p-value: < 2.2e-16

``` r
round(boscalid.lm.2[[1]][[2]],4)
```

    ## [1] 0.8657

``` r
round(boscalid.lm.2[[1]][[1]],3) 
```

    ## [1] 0.014

``` r
myplot_model_2(boscalid.EC50DC)+ expand_limits(x = c(0.04, 0.25), y = c(0.04, 0.25))  + scale_x_continuous(name = waiver(),
                                                                             breaks = c(0.05, 0.1, 0.15, 0.20, 0.25)) + scale_y_continuous(breaks = c(0.05, 0.1, 0.15, 0.20, 0.25)) + geom_label(aes(x = 0.20,
                                                    y = 0.10)
                                                ,
                                                label = c(paste  (
                                                  " Y =",
                                                  paste0 (round(summary(boscalid.lm.2)[[4]][2], 4), "x"),
                                                  
                                                  round(summary(boscalid.lm.2)[[4]][1], 3) *
                                                    -1,
                                                  paste  ("\n R =",
                                                          paste0 (round(
                                                           summary(boscalid.lm.2)[[9]][1], 4
                                                          ), "\n p < 0.001"))
                                                )),size = 5,
                                                fontface = "bold")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Getting%20Discriminatory%20dose-2.png)<!-- -->

``` r
round(boscalid.lm[[1]][[2]],4)
```

    ## [1] 0.0567

``` r
round(boscalid.lm[[1]][[1]],3) 
```

    ## [1] -4.412

``` r
#Tetraconazole
  
tetraconazole <- getting_DC(hey2) 
```

    ## `summarise()` has grouped output by 'ID'. You can override using the `.groups`
    ## argument.

``` r
  tetraconazole.joined <-  tetraconazole%>% left_join(tetraconazole.EC50)
```

    ## Joining, by = "ID"

``` r
tetraconazole.joined.2 <-  getting_DC_2(tetraconazole.joined) 
 ###Dose chosen as DC. Although there is one dose with higher r2, there is no difference between and this was selected
# Linear model of log EC50 and relative growth at 2 ppm, check normality and homogeneity of variances

finalRG_2 <- lm(logEC50 ~ RG_2, tetraconazole.joined.2)
summary(finalRG_2)
```

    ## 
    ## Call:
    ## lm(formula = logEC50 ~ RG_2, data = tetraconazole.joined.2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15348 -0.04122  0.01195  0.04652  0.19668 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -1.429298   0.052190  -27.39   <2e-16 ***
    ## RG_2         0.039056   0.001389   28.11   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08135 on 38 degrees of freedom
    ## Multiple R-squared:  0.9541, Adjusted R-squared:  0.9529 
    ## F-statistic: 790.1 on 1 and 38 DF,  p-value: < 2.2e-16

``` r
#plot(finalRG0.2)

##taking out outliers from Q-Q PLOT $ LEVERAGE TEST, ISOLATES # 1032 from treatment, 2222 & 2223 from Farmer fields##
tetraconazole.joined.2.clean <-
  tetraconazole.joined.2 %>% filter(!ID == "21", !ID == "558",!ID == "667") %>% 
   select(ID, logEC50, RG_2, EC50) %>% 
  dplyr::rename(RG = RG_2) %>% 
   wrangling_DC()

summary(tetraconazole.lm <- lineal_model("logEC50","RG",data =  tetraconazole.joined.2.clean))
```

    ## 
    ## Call:
    ## lm(formula = paste(vardep, "~", varindep), data = data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.14413 -0.04247  0.02069  0.03689  0.11121 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -1.669458   0.076294  -21.88   <2e-16 ***
    ## RG           0.045156   0.001958   23.06   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06204 on 35 degrees of freedom
    ## Multiple R-squared:  0.9382, Adjusted R-squared:  0.9365 
    ## F-statistic: 531.8 on 1 and 35 DF,  p-value: < 2.2e-16

``` r
   # Getting the EC50DC according to the model
tetraconazole.EC50DC <- tetraconazole.joined.2.clean %>% mutate(Estimate.50DC =
                                     #using the model the intercept and the coefficient of the model
                                    exp(tetraconazole.lm[[1]][[1]] + tetraconazole.lm[[1]][[2]] * RG)) %>%# here I need to do this dynamic
  select(ID, Estimate.50DC, EC50)  %>%  wrangling_DC()

# Linear Regression graph of log EC50 and relative growth at 0.2 ppm
###

myplot_model_1_tetraconazole(tetraconazole.joined.2.clean, lm = tetraconazole.lm) 
```

    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Getting%20Discriminatory%20dose-3.png)<!-- -->

``` r
######### Linear model of EC50  and EC50DC just to corroborate is a GOOD RELATIONSHIP, check normality and homogeneity of variances

summary(tetraconazole.lm.2 <- lineal_model("Estimate.50DC","EC50",data =  tetraconazole.EC50DC))
```

    ## 
    ## Call:
    ## lm(formula = paste(vardep, "~", varindep), data = data)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.10710 -0.03721 -0.01551  0.05445  0.15171 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.14852    0.05618   2.644   0.0122 *  
    ## EC50         0.86314    0.05000  17.264   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06418 on 35 degrees of freedom
    ## Multiple R-squared:  0.8949, Adjusted R-squared:  0.8919 
    ## F-statistic:   298 on 1 and 35 DF,  p-value: < 2.2e-16

``` r
round(tetraconazole.lm.2[[1]][[2]],4)
```

    ## [1] 0.8631

``` r
round(tetraconazole.lm.2[[1]][[1]],3)
```

    ## [1] 0.149

``` r
myplot_model_2(tetraconazole.EC50DC)+ expand_limits(x = c(0.27, 1.8), y = c(0.27, 1.8))  + scale_x_continuous(name = waiver(),
                                                                             breaks = c(0.3, 0.6, 0.9, 1.2, 1.5,1.8)) + scale_y_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5,1.8)) + geom_label(aes(x = 1.5,
                  y = 0.6)
              ,
              label = c(paste  (
                " Y =",
                paste0 (round(summary(tetraconazole.lm.2)[[4]][2], 4), "x"),
                round(summary(tetraconazole.lm.2)[[4]][1], 3) *
                  -1,
                paste  ("\n R =",
                        paste0 (round(
                          summary(tetraconazole.lm.2)[[9]][1], 4
                        ), "\n p < 0.001"))
              )),size = 5,
              fontface = "bold")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Getting%20Discriminatory%20dose-4.png)<!-- -->

``` r
 ## Picoxystrobin

picoxystrobin<- getting_DC(hey4) 
```

    ## `summarise()` has grouped output by 'ID'. You can override using the `.groups`
    ## argument.

``` r
  picoxystrobin.joined <-  picoxystrobin%>% left_join(picoxystrobin.EC50)
```

    ## Joining, by = "ID"

``` r
 picoxystrobin.joined.2 <-  getting_DC_2(picoxystrobin.joined) 
 ###Dose chosen as DC. Although there is one dose with higher r2, there is no difference between and this was selected
# Linear model of log EC50 and relative growth at 0.01 ppm, check normality and homogeneity of variances
finalRG_0.01 <- lm(logEC50 ~ RG_0.01,  picoxystrobin.joined.2)
summary(finalRG_0.01)
```

    ## 
    ## Call:
    ## lm(formula = logEC50 ~ RG_0.01, data = picoxystrobin.joined.2)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.249660 -0.053292  0.004241  0.055571  0.188057 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -6.685437   0.158498  -42.18   <2e-16 ***
    ## RG_0.01      0.040251   0.002815   14.30   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09447 on 38 degrees of freedom
    ## Multiple R-squared:  0.8433, Adjusted R-squared:  0.8392 
    ## F-statistic: 204.5 on 1 and 38 DF,  p-value: < 2.2e-16

``` r
##taking out outliers from Q-Q PLOT $ LEVERAGE TEST, ISOLATES # 1032 from treatment, 2222 & 2223 from Farmer fields##
picoxystrobin.joined.2.clean <-
  picoxystrobin.joined.2 %>% filter(!ID == "1032", !ID == "558",!ID == "581") %>% 
  select(ID, logEC50, RG_0.01, EC50) %>% 
  dplyr::rename(RG = RG_0.01) %>% 
   wrangling_DC()


summary(picoxystrobin.lm <- lineal_model("logEC50","RG",data =  picoxystrobin.joined.2.clean))
```

    ## 
    ## Call:
    ## lm(formula = paste(vardep, "~", varindep), data = data)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.180702 -0.055184 -0.000729  0.049534  0.147093 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) -6.571213   0.139542  -47.09   <2e-16 ***
    ## RG           0.038131   0.002507   15.21   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07633 on 35 degrees of freedom
    ## Multiple R-squared:  0.8686, Adjusted R-squared:  0.8648 
    ## F-statistic: 231.3 on 1 and 35 DF,  p-value: < 2.2e-16

``` r
   # Getting the EC50DC according to the model
picoxystrobin.EC50DC <- picoxystrobin.joined.2.clean %>% mutate(Estimate.50DC =
                                     #using the model the intercept and the coefficient of the model
                                    exp(picoxystrobin.lm[[1]][[1]] + picoxystrobin.lm[[1]][[2]] * RG)) %>%# here I need to do this dynamic
  select(ID, Estimate.50DC, EC50)  %>%  wrangling_DC()

# Linear Regression graph of log EC50 and relative growth at 0.2 ppm
###
summary(picoxystrobin.lm.2 <- lineal_model("Estimate.50DC","EC50",data =  picoxystrobin.EC50DC))
```

    ## 
    ## Call:
    ## lm(formula = paste(vardep, "~", varindep), data = data)
    ## 
    ## Residuals:
    ##        Min         1Q     Median         3Q        Max 
    ## -0.0015409 -0.0006240 -0.0002075  0.0007106  0.0026937 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 0.001884   0.000798   2.361   0.0239 *  
    ## EC50        0.838025   0.066339  12.632 1.35e-14 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0008963 on 35 degrees of freedom
    ## Multiple R-squared:  0.8201, Adjusted R-squared:  0.815 
    ## F-statistic: 159.6 on 1 and 35 DF,  p-value: 1.347e-14

``` r
round(picoxystrobin.lm.2[[1]][[2]],4)
```

    ## [1] 0.838

``` r
round(picoxystrobin.lm.2[[1]][[1]],3)
```

    ## [1] 0.002

``` r
myplot_model_1_picoxystrobin(picoxystrobin.joined.2.clean, lm = picoxystrobin.lm) 
```

    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Getting%20Discriminatory%20dose-5.png)<!-- -->

``` r
######### Linear model of EC50  and EC50DC just to corroborate is a GOOD RELATIONSHIP, check normality and homogeneity of variances
myplot_model_2(picoxystrobin.EC50DC)+ expand_limits(x = c(0.005, 0.018), y = c(0.005, 0.018))  + scale_x_continuous(name = waiver(),
                                                                             breaks = c(0.0075, 0.01, 0.0125, 0.0150, 0.018)) + scale_y_continuous(breaks = c(0.0075, 0.01, 0.0125, 0.0150, 0.018)) + geom_label(aes(x = 0.015,
                  y = 0.009)
              ,
              label = c(paste  (
                " Y =",
                paste0 (round(summary(picoxystrobin.lm.2)[[4]][2], 4), "x"),
                round(summary(picoxystrobin.lm.2)[[4]][1], 3) *
                  -1,
                paste  ("\n R =",
                        paste0 (round(
                          summary(picoxystrobin.lm.2)[[9]][1], 4
                        ), "\n p < 0.001"))
              )),size = 5,
              fontface = "bold")
```

    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Getting%20Discriminatory%20dose-6.png)<!-- -->

``` r
#Function to clean
picoxystrobin.joined.2 %>% filter(!EC50 %in% Outlier(EC50))
```

    ## # A tibble: 36 × 15
    ##       ID RG_0.01 RG_0.02 RG_0.04 RG_0.06 RG_0.1 Estimate.10    SE.10    EC50
    ##    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>       <dbl>    <dbl>   <dbl>
    ##  1     1    55.8    37.4    19.6   15.4   10.0      0.00169 0.000264 0.0123 
    ##  2    12    55.8    39.7    25.3   19.5   12.2      0.00124 0.000198 0.0128 
    ##  3    20    56.7    37.5    21.4   15.0    8.74     0.00181 0.000291 0.0127 
    ##  4    21    56.9    23.3    14.0    8.90   8.90     0.00268 0.000564 0.0111 
    ##  5    74    49.9    24.1    14.5   12.0    9.84     0.00135 0.000497 0.00906
    ##  6    87    49.8    25.0    15.0   11.7    9.10     0.00136 0.000557 0.00916
    ##  7   118    55.0    36.1    24.0   18.9    9.89     0.00126 0.000226 0.0120 
    ##  8   123    54.2    38.4    20.7   11.6    8.17     0.00179 0.000337 0.0121 
    ##  9   449    56.0    30.8    19.0   18.0   11.5      0.00135 0.000494 0.0111 
    ## 10   461    53.9    34.8    20.7   18.3   11.8      0.00115 0.000458 0.0111 
    ## # … with 26 more rows, and 6 more variables: SE.50 <dbl>, Estimate.90 <dbl>,
    ## #   SE.90 <dbl>, source <fct>, Host <fct>, logEC50 <dbl>

### Reading from White Mold Inventory; needed downstream

``` r
#based on object table1 which has all the complete isolates used in this study,  

 # I have to keep Field adn Year separate tahts why reading again the main table from data
 
 whitemold.inventory.2 <- read_csv("data/whitemold.inventory.csv")
```

    ## New names:
    ## • `` -> `...1`

    ## Warning: One or more parsing issues, see `problems()` for details

    ## Rows: 2326 Columns: 40
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (30): Host, State, County, Country, Field, Form.ID, Serial_agar_dilution...
    ## dbl  (8): ...1, ID, Year, lat, long, Applications, Rate, Number.of.Years.1
    ## lgl  (2): company_current_season, company_previous_seasons
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
 Colletotrichum <-
    whitemold.inventory.2 %>% 
   # filter(ID %in% lasiodiplodia) %>% 
   mutate(country = ifelse(
      ID %in% mexican.isolates,
      "Mexico",
      ifelse(ID %in% brazilian.isolates,
             "Brazil", "USA")),
      source = ifelse(
    ID %in% baseline.isolates.2,
    "Baseline",
    ifelse(
      ID %in% Farmer.fields,
      "Producer Fields","Fungicide Field Trials"
    )), Host = ifelse(
    ID %in% drybean,
    "Drybean", ifelse(
      ID %in% soybean,
      "Soybean","different_host")),
 source = as.factor(source),
  Host = as.factor(Host),
 country = as.factor(country),
 Field = as.factor(Field),
 Year = as.factor(Year),
 State = as.factor(State),
 County = as.factor(County),
 Fungicide_current_season = as.factor(Fungicide_current_season),
 molecule.s._current_season = as.factor(molecule.s._current_season),
 Group_MOA_current_season = as.factor(Group_MOA_current_season),
 Applications = as.factor(Applications),
 Rate = as.factor(Rate)
      ) %>%    select(
     ID,
     source,#we will change later to the real variable "source"
     Host,
     country, 
     Field,
     Year,
     State,
     County,
     Fungicide_current_season,
     molecule.s._current_season,
     Group_MOA_current_season,
     Applications,
     Rate
   ) 
 
 controls.in.the.plots <- Colletotrichum %>% filter( Fungicide_current_season == "control") %>% select(ID) %>% pull()

 controls.in.the.plots.2 <- c(1175, 1195, 1327, 1443, 1446, 1461, 1541, 1710)
```

## Use of the Discriminatory Concentration DC

We randomly selected 154 isolates: 86 from dry beans (from 9 states with
1-17 isolates/state, and 2 states from Brazil) and 66 from soybeans (6
NE, 18 IA, 13 WI & 18 MI and 4 from Mexico) from the collection
(baseline, FF, and FF) plus 5 same isolates used in the serial dilution
approach

``` r
TM.joined.2 <- TM.joined.2 %>% select(ID, EC50)

survey.TM <- bind_rows(survey.TM.thomas,
survey.TM.USA,
survey.TM.more.USA,
#survey.TM.MX.BZ, there are the same the other fungicides
TM.joined.2) 

#Have to analyze the controls

#wranglind_DC_2

survey.TM.2 <- survey.TM%>% 
  wrangling_DC() %>% 
  mutate(Field_Year =  ifelse(source == "Baseline", "Baseline", Field_Year),
    Field_Year = as.factor(Field_Year)) %>% 
  dplyr::rename(EC50DC = EC50) %>% #change here to be able to join to other datframes to get the total used
  select(ID, EC50DC) %>% 
   # is here where now I have to remove the isolates from controls from fields and potato host isolates provided by Dr.Chilvers at the end
  dplyr::filter(! ID %in% controls.in.the.plots,  ! ID %in% 2570:2574)
  
survey.boscalid.thomas.2 <- wrangling_survey(survey.boscalid.thomas) %>%
  mutate(EC50DC = exp (boscalid.lm [[1]][[1]] + boscalid.lm[[1]][[2]] * RG))# exponential is the opposite of log
        #using the model the intercept and the coefficient of the model

survey.boscalid.USA.2 <- wrangling_survey(survey.boscalid.USA) %>%
  mutate(EC50DC = exp (boscalid.lm [[1]][[1]] + boscalid.lm[[1]][[2]] * RG))#

##
survey.boscalid.more.USA <- survey.boscalid.more.USA %>% mutate(response = (response/10), control = (control/10))




survey.boscalid.more.USA.2 <- wrangling_survey(survey.boscalid.more.USA) %>%
  mutate(EC50DC = exp (boscalid.lm [[1]][[1]] + boscalid.lm[[1]][[2]] * RG))#
##

survey.boscalid.MX.BZ.2 <- wrangling_survey(survey.boscalid.MX.BZ) %>%
  mutate(EC50DC = exp (boscalid.lm [[1]][[1]] + boscalid.lm[[1]][[2]] * RG))#
       
#Gathering 
boscalid.EC50.2 <- boscalid.EC50 %>% select(ID, EC50) %>% 
  dplyr::rename(EC50DC = EC50)#just renaming in purpose to join the following dataframe

survey.boscalid <- bind_rows(survey.boscalid.thomas.2,
                             survey.boscalid.USA.2,
                             survey.boscalid.more.USA.2,
                             survey.boscalid.MX.BZ.2, 
                             boscalid.EC50.2) %>% select(ID,EC50DC)%>% 
   # is here where now I have to remove the isolates from controls from fields and potato host isolates provided by Dr. Chilvers at the end
  dplyr::filter(! ID %in% controls.in.the.plots, ! ID %in% 2570:2574 )

#Have to analyze the controls
#fix the function here
#wranglind_DC_2
#Check this
IDA <- ID
survey.boscalid.complete <- survey.boscalid%>% 
  wrangling_DC() %>%  
 wrangling_DC_2() 

survey.boscalid.fields <- survey.boscalid%>% 
  wrangling_DC() %>% 
   wrangling_DC_2() %>% 
  wrangling_DC_3()
  
###For comparisons between hosts I had to take out the ones coming from baseline otherwise would not be fair

survey.boscalid.complete.host <- survey.boscalid.complete %>% 
  filter(!source=="Baseline")

# Tetraconazole
#any particular case check so far by this way
survey.tetraconazole.thomas.2 <- wrangling_survey(survey.tetraconazole.thomas) %>%
  mutate(EC50DC = exp (tetraconazole.lm [[1]][[1]] + tetraconazole.lm[[1]][[2]] * RG))# exponential is the opposite of log
        #using the model the intercept and the coefficient of the model
##
survey.tetraconazole.USA.2 <- wrangling_survey(survey.tetraconazole.USA) %>%
  mutate(EC50DC = exp (tetraconazole.lm [[1]][[1]] + tetraconazole.lm[[1]][[2]] * RG))#

##
survey.tetraconazole.more.USA <- survey.tetraconazole.more.USA %>% mutate(response = (response/10), control = (control/10))

survey.tetraconazole.more.USA.2 <- wrangling_survey(survey.tetraconazole.more.USA) %>%
  mutate(EC50DC = exp (tetraconazole.lm [[1]][[1]] + tetraconazole.lm[[1]][[2]] * RG))#

survey.tetraconazole.MX.BZ.2 <- wrangling_survey(survey.tetraconazole.MX.BZ) %>%
  mutate(EC50DC = exp (tetraconazole.lm [[1]][[1]] + tetraconazole.lm[[1]][[2]] * RG))#
       
#Gathering 
tetraconazole.EC50.2 <- tetraconazole.EC50 %>% select(ID, EC50) %>% 
  dplyr::rename(EC50DC = EC50)#just rebaming in purpose to join the following dataframe

survey.tetraconazole <- bind_rows(survey.tetraconazole.thomas.2,
                             survey.tetraconazole.USA.2,
                             survey.tetraconazole.more.USA.2,
                             survey.tetraconazole.MX.BZ.2, 
                             tetraconazole.EC50.2) %>% select(ID,  EC50DC) %>% 
# is here where now I have to remove the isolates from controls from fields and potato host isolates provided by Dr. Chilvers at the end
  dplyr::filter(! ID %in% controls.in.the.plots, ! ID %in% 2570:2574 )
#Have to analyze the controls

#wranglind_DC_2

survey.tetraconazole.complete <- survey.tetraconazole%>% 
  wrangling_DC() %>% 
  wrangling_DC_2()

survey.tetraconazole.fields <- survey.tetraconazole%>% 
  wrangling_DC() %>%
  wrangling_DC_2() %>% 
  wrangling_DC_3()
  # count(Field_Year) %>% 
  #  filter(n >= 9)




###For comparisons between hosts I had to take out the ones coming from baseline otherwise would not be fair

survey.tetraconazole.complete.host <- survey.tetraconazole.complete %>% 
  filter(!source=="Baseline")



# picoxystrobin

#any particular case check so far by this way
survey.picoxystrobin.thomas.2 <- wrangling_survey(survey.picoxystrobin.thomas) %>%
  mutate(EC50DC = exp (picoxystrobin.lm [[1]][[1]] + picoxystrobin.lm[[1]][[2]] * RG))# exponential is the opposite of log
        #using the model the intercept and the coefficient of the model
##
# 
survey.picoxystrobin.USA.2 <- wrangling_survey(survey.picoxystrobin.USA) %>%
  mutate(EC50DC = exp (picoxystrobin.lm [[1]][[1]] + picoxystrobin.lm[[1]][[2]] * RG))#

##
survey.picoxystrobin.more.USA <- survey.picoxystrobin.more.USA %>% mutate(response = (response/10), control = (control/10))

survey.picoxystrobin.more.USA.2 <- wrangling_survey(survey.picoxystrobin.more.USA) %>%
  mutate(EC50DC = exp (picoxystrobin.lm [[1]][[1]] + picoxystrobin.lm[[1]][[2]] * RG))#
##

survey.picoxystrobin.MX.BZ.2 <- wrangling_survey(survey.picoxystrobin.MX.BZ) %>%
  mutate(EC50DC = exp (picoxystrobin.lm [[1]][[1]] + picoxystrobin.lm[[1]][[2]] * RG))#
       
#Gathering 
picoxystrobin.EC50.2 <- picoxystrobin.EC50 %>% select(ID, EC50) %>% 
  dplyr::rename(EC50DC = EC50)#just rebaming in purpose to join the following dataframe

survey.picoxystrobin <- bind_rows(survey.picoxystrobin.thomas.2,
                             survey.picoxystrobin.USA.2,
                             survey.picoxystrobin.more.USA.2,
                             survey.picoxystrobin.MX.BZ.2, 
                             picoxystrobin.EC50.2) %>% select(ID,  EC50DC) %>% 
# is here where now I have to remove the isolates from controls from fields and potato host isolates provided by Dr. Chilvers at the end
  dplyr::filter(! ID %in% controls.in.the.plots, ! ID %in% 2570:2574 )
#Have to analyze the controls

survey.picoxystrobin.complete <- survey.picoxystrobin%>% 
  wrangling_DC() %>% 
  wrangling_DC_2()

survey.picoxystrobin.fields <- survey.picoxystrobin%>% 
  wrangling_DC() %>%
  wrangling_DC_2() %>% 
  wrangling_DC_3()
  
###For comparisons between hosts I had to take out the ones coming from baseline otherwise would not be fair

survey.picoxystrobin.complete.host <- survey.picoxystrobin.complete %>% 
  filter(!source=="Baseline")
```

# Statistic Analysis

## Boscalid

### By Source

``` r
test_kruskal_source(survey.boscalid.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() )
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  as.numeric(vardep) by as.factor(varindep)
    ## Kruskal-Wallis chi-squared = 18.701, df = 2, p-value = 8.694e-05

``` r
test_Dunn_source(survey.boscalid.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() )
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                                        mean.rank.diff    pval    
    ## Fungicide Field Trials-Baseline              4.187608 1.00000    
    ## Producer Fields-Baseline                   -51.975684 0.11835    
    ## Producer Fields-Fungicide Field Trials     -56.163291 0.00017 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
generate_label( test_Dunn_source(survey.boscalid.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() ))
```

    ##                                                   plot.labels  Dunn.levels
    ## Fungicide Field Trials-Baseline        Fungicide Field Trials 1.0000000000
    ## Producer Fields-Baseline                      Producer Fields 0.1183489091
    ## Producer Fields-Fungicide Field Trials               Baseline 0.0001666306
    ##                                        labels
    ## Fungicide Field Trials-Baseline             a
    ## Producer Fields-Baseline                    b
    ## Producer Fields-Fungicide Field Trials     ab

``` r
survey.boscalid.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup()  %>% 
  group_by(source) %>% summarize(round(range(EC50DC),3), average= round(mean(EC50DC), 2), N= n(),  sd = sd(EC50DC), se = round(sd/sqrt(N),3))
```

    ## `summarise()` has grouped output by 'source'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 6 × 6
    ## # Groups:   source [3]
    ##   source                 `round(range(EC50DC), 3)` average     N     sd    se
    ##   <fct>                                      <dbl>   <dbl> <int>  <dbl> <dbl>
    ## 1 Baseline                                   0.066    0.11    21 0.0353 0.008
    ## 2 Baseline                                   0.222    0.11    21 0.0353 0.008
    ## 3 Fungicide Field Trials                     0.062    0.11    83 0.0219 0.002
    ## 4 Fungicide Field Trials                     0.163    0.11    83 0.0219 0.002
    ## 5 Producer Fields                            0.042    0.1    282 0.0290 0.002
    ## 6 Producer Fields                            0.18     0.1    282 0.0290 0.002

### By Host

``` r
#generate_label(
  test_kruskal_host(survey.boscalid.complete.host %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() )
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  as.numeric(vardep) by as.factor(varindep)
    ## Kruskal-Wallis chi-squared = 1.6187, df = 1, p-value = 0.2033

``` r
#generate_label
  test_Dunn_host(survey.boscalid.complete.host %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() )
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                 mean.rank.diff   pval    
    ## Soybean-Drybean       14.53249 0.2033    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
survey.boscalid.complete.host %>%
  group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() %>% 
  group_by(Host) %>% summarize(round(range(EC50DC),3), average= round(mean(EC50DC),2), N= n(),  sd = sd(EC50DC), se = round(sd/sqrt(N),3))
```

    ## `summarise()` has grouped output by 'Host'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 4 × 6
    ## # Groups:   Host [2]
    ##   Host    `round(range(EC50DC), 3)` average     N     sd    se
    ##   <fct>                       <dbl>   <dbl> <int>  <dbl> <dbl>
    ## 1 Drybean                     0.049     0.1   229 0.0298 0.002
    ## 2 Drybean                     0.18      0.1   229 0.0298 0.002
    ## 3 Soybean                     0.042     0.1   136 0.0246 0.002
    ## 4 Soybean                     0.163     0.1   136 0.0246 0.002

### By Field

``` r
generate_label(  test_field (survey.boscalid.fields))
```

    ##                                                   plot.labels  Dunn.levels
    ## ND.07.fld01.DB-NE.10.fld21.DB                  ND.07.fld01.DB 1.000000e+00
    ## WA.07.fld09.DB-NE.10.fld21.DB                  WA.07.fld09.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-NE.10.fld21.DB    MI.16.East.Lansing.SB 1.000000e+00
    ## WA.08.fld10.DB-NE.10.fld21.DB                  WA.08.fld10.DB 1.000000e+00
    ## ND.10.fld25.DB-NE.10.fld21.DB                  ND.10.fld25.DB 2.477204e-01
    ## MI.15.Montcalm.SB-NE.10.fld21.DB            MI.15.Montcalm.SB 4.244188e-02
    ## IA.17.Nashua.SB-NE.10.fld21.DB                IA.17.Nashua.SB 2.054563e-02
    ## Baseline-NE.10.fld21.DB                              Baseline 2.087706e-02
    ## WI.17.Hancock.SB-NE.10.fld21.DB              WI.17.Hancock.SB 2.713616e-02
    ## MI.17.Montcalm.SB-NE.10.fld21.DB            MI.17.Montcalm.SB 3.951429e-03
    ## BR.2010.por.DB-NE.10.fld21.DB                  BR.2010.por.DB 9.775535e-04
    ## BR.2010.unai.DB-NE.10.fld21.DB                BR.2010.unai.DB 2.807856e-04
    ## MX.2017.pre.DB-NE.10.fld21.DB                  MX.2017.pre.DB 2.305319e-04
    ## MX.2017.ca.DB-NE.10.fld21.DB                    MX.2017.ca.DB 9.975110e-05
    ## BR.2010.lam.DB-NE.10.fld21.DB                  BR.2010.lam.DB 9.373366e-06
    ## MX.2017.ley.DB-NE.10.fld21.DB                  MX.2017.ley.DB 4.057426e-06
    ## WA.07.fld09.DB-ND.07.fld01.DB                  NE.10.fld21.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-ND.07.fld01.DB           ND.07.fld01.DB 1.000000e+00
    ## WA.08.fld10.DB-ND.07.fld01.DB                  WA.07.fld09.DB 1.000000e+00
    ## ND.10.fld25.DB-ND.07.fld01.DB           MI.16.East.Lansing.SB 1.000000e+00
    ## MI.15.Montcalm.SB-ND.07.fld01.DB               WA.08.fld10.DB 3.562009e-01
    ## IA.17.Nashua.SB-ND.07.fld01.DB                 ND.10.fld25.DB 2.152783e-01
    ## Baseline-ND.07.fld01.DB                     MI.15.Montcalm.SB 2.312318e-01
    ## WI.17.Hancock.SB-ND.07.fld01.DB               IA.17.Nashua.SB 2.586017e-01
    ## MI.17.Montcalm.SB-ND.07.fld01.DB                     Baseline 4.618881e-02
    ## BR.2010.por.DB-ND.07.fld01.DB                WI.17.Hancock.SB 1.080346e-02
    ## BR.2010.unai.DB-ND.07.fld01.DB              MI.17.Montcalm.SB 3.178176e-03
    ## MX.2017.pre.DB-ND.07.fld01.DB                  BR.2010.por.DB 2.699022e-03
    ## MX.2017.ca.DB-ND.07.fld01.DB                  BR.2010.unai.DB 1.199445e-03
    ## BR.2010.lam.DB-ND.07.fld01.DB                  MX.2017.pre.DB 1.167211e-04
    ## MX.2017.ley.DB-ND.07.fld01.DB                   MX.2017.ca.DB 5.256504e-05
    ## MI.16.East.Lansing.SB-WA.07.fld09.DB           BR.2010.lam.DB 1.000000e+00
    ## WA.08.fld10.DB-WA.07.fld09.DB                  MX.2017.ley.DB 1.000000e+00
    ## ND.10.fld25.DB-WA.07.fld09.DB                  NE.10.fld21.DB 1.000000e+00
    ## MI.15.Montcalm.SB-WA.07.fld09.DB               ND.07.fld01.DB 1.000000e+00
    ## IA.17.Nashua.SB-WA.07.fld09.DB                 WA.07.fld09.DB 1.000000e+00
    ## Baseline-WA.07.fld09.DB                 MI.16.East.Lansing.SB 1.000000e+00
    ## WI.17.Hancock.SB-WA.07.fld09.DB                WA.08.fld10.DB 1.000000e+00
    ## MI.17.Montcalm.SB-WA.07.fld09.DB               ND.10.fld25.DB 6.327472e-01
    ## BR.2010.por.DB-WA.07.fld09.DB               MI.15.Montcalm.SB 1.384048e-01
    ## BR.2010.unai.DB-WA.07.fld09.DB                IA.17.Nashua.SB 4.323582e-02
    ## MX.2017.pre.DB-WA.07.fld09.DB                        Baseline 3.913987e-02
    ## MX.2017.ca.DB-WA.07.fld09.DB                 WI.17.Hancock.SB 1.854095e-02
    ## BR.2010.lam.DB-WA.07.fld09.DB               MI.17.Montcalm.SB 2.007667e-03
    ## MX.2017.ley.DB-WA.07.fld09.DB                  BR.2010.por.DB 1.030144e-03
    ## WA.08.fld10.DB-MI.16.East.Lansing.SB          BR.2010.unai.DB 1.000000e+00
    ## ND.10.fld25.DB-MI.16.East.Lansing.SB           MX.2017.pre.DB 1.000000e+00
    ## MI.15.Montcalm.SB-MI.16.East.Lansing.SB         MX.2017.ca.DB 1.000000e+00
    ## IA.17.Nashua.SB-MI.16.East.Lansing.SB          BR.2010.lam.DB 1.000000e+00
    ## Baseline-MI.16.East.Lansing.SB                 MX.2017.ley.DB 1.000000e+00
    ## WI.17.Hancock.SB-MI.16.East.Lansing.SB         NE.10.fld21.DB 1.000000e+00
    ## MI.17.Montcalm.SB-MI.16.East.Lansing.SB        ND.07.fld01.DB 1.000000e+00
    ## BR.2010.por.DB-MI.16.East.Lansing.SB           WA.07.fld09.DB 1.000000e+00
    ## BR.2010.unai.DB-MI.16.East.Lansing.SB   MI.16.East.Lansing.SB 6.220960e-01
    ## MX.2017.pre.DB-MI.16.East.Lansing.SB           WA.08.fld10.DB 6.198408e-01
    ## MX.2017.ca.DB-MI.16.East.Lansing.SB            ND.10.fld25.DB 3.352475e-01
    ## BR.2010.lam.DB-MI.16.East.Lansing.SB        MI.15.Montcalm.SB 4.755599e-02
    ## MX.2017.ley.DB-MI.16.East.Lansing.SB          IA.17.Nashua.SB 3.227069e-02
    ## ND.10.fld25.DB-WA.08.fld10.DB                        Baseline 1.000000e+00
    ## MI.15.Montcalm.SB-WA.08.fld10.DB             WI.17.Hancock.SB 1.000000e+00
    ## IA.17.Nashua.SB-WA.08.fld10.DB              MI.17.Montcalm.SB 1.000000e+00
    ## Baseline-WA.08.fld10.DB                        BR.2010.por.DB 1.000000e+00
    ## WI.17.Hancock.SB-WA.08.fld10.DB               BR.2010.unai.DB 1.000000e+00
    ## MI.17.Montcalm.SB-WA.08.fld10.DB               MX.2017.pre.DB 1.000000e+00
    ## BR.2010.por.DB-WA.08.fld10.DB                   MX.2017.ca.DB 1.000000e+00
    ## BR.2010.unai.DB-WA.08.fld10.DB                 BR.2010.lam.DB 1.000000e+00
    ## MX.2017.pre.DB-WA.08.fld10.DB                  MX.2017.ley.DB 1.000000e+00
    ## MX.2017.ca.DB-WA.08.fld10.DB                   NE.10.fld21.DB 1.000000e+00
    ## BR.2010.lam.DB-WA.08.fld10.DB                  ND.07.fld01.DB 1.952766e-01
    ## MX.2017.ley.DB-WA.08.fld10.DB                  WA.07.fld09.DB 1.606252e-01
    ## MI.15.Montcalm.SB-ND.10.fld25.DB        MI.16.East.Lansing.SB 1.000000e+00
    ## IA.17.Nashua.SB-ND.10.fld25.DB                 WA.08.fld10.DB 1.000000e+00
    ## Baseline-ND.10.fld25.DB                        ND.10.fld25.DB 1.000000e+00
    ## WI.17.Hancock.SB-ND.10.fld25.DB             MI.15.Montcalm.SB 1.000000e+00
    ## MI.17.Montcalm.SB-ND.10.fld25.DB              IA.17.Nashua.SB 1.000000e+00
    ## BR.2010.por.DB-ND.10.fld25.DB                        Baseline 1.000000e+00
    ## BR.2010.unai.DB-ND.10.fld25.DB               WI.17.Hancock.SB 1.000000e+00
    ## MX.2017.pre.DB-ND.10.fld25.DB               MI.17.Montcalm.SB 1.000000e+00
    ## MX.2017.ca.DB-ND.10.fld25.DB                   BR.2010.por.DB 1.000000e+00
    ## BR.2010.lam.DB-ND.10.fld25.DB                 BR.2010.unai.DB 1.000000e+00
    ## MX.2017.ley.DB-ND.10.fld25.DB                  MX.2017.pre.DB 1.000000e+00
    ## IA.17.Nashua.SB-MI.15.Montcalm.SB               MX.2017.ca.DB 1.000000e+00
    ## Baseline-MI.15.Montcalm.SB                     BR.2010.lam.DB 1.000000e+00
    ## WI.17.Hancock.SB-MI.15.Montcalm.SB             MX.2017.ley.DB 1.000000e+00
    ## MI.17.Montcalm.SB-MI.15.Montcalm.SB            NE.10.fld21.DB 1.000000e+00
    ## BR.2010.por.DB-MI.15.Montcalm.SB               ND.07.fld01.DB 1.000000e+00
    ## BR.2010.unai.DB-MI.15.Montcalm.SB              WA.07.fld09.DB 1.000000e+00
    ## MX.2017.pre.DB-MI.15.Montcalm.SB        MI.16.East.Lansing.SB 1.000000e+00
    ## MX.2017.ca.DB-MI.15.Montcalm.SB                WA.08.fld10.DB 1.000000e+00
    ## BR.2010.lam.DB-MI.15.Montcalm.SB               ND.10.fld25.DB 1.000000e+00
    ## MX.2017.ley.DB-MI.15.Montcalm.SB            MI.15.Montcalm.SB 1.000000e+00
    ## Baseline-IA.17.Nashua.SB                      IA.17.Nashua.SB 1.000000e+00
    ## WI.17.Hancock.SB-IA.17.Nashua.SB                     Baseline 1.000000e+00
    ## MI.17.Montcalm.SB-IA.17.Nashua.SB            WI.17.Hancock.SB 1.000000e+00
    ## BR.2010.por.DB-IA.17.Nashua.SB              MI.17.Montcalm.SB 1.000000e+00
    ## BR.2010.unai.DB-IA.17.Nashua.SB                BR.2010.por.DB 1.000000e+00
    ## MX.2017.pre.DB-IA.17.Nashua.SB                BR.2010.unai.DB 1.000000e+00
    ## MX.2017.ca.DB-IA.17.Nashua.SB                  MX.2017.pre.DB 1.000000e+00
    ## BR.2010.lam.DB-IA.17.Nashua.SB                  MX.2017.ca.DB 1.000000e+00
    ## MX.2017.ley.DB-IA.17.Nashua.SB                 BR.2010.lam.DB 1.000000e+00
    ## WI.17.Hancock.SB-Baseline                      MX.2017.ley.DB 1.000000e+00
    ## MI.17.Montcalm.SB-Baseline                     NE.10.fld21.DB 1.000000e+00
    ## BR.2010.por.DB-Baseline                        ND.07.fld01.DB 1.000000e+00
    ## BR.2010.unai.DB-Baseline                       WA.07.fld09.DB 1.000000e+00
    ## MX.2017.pre.DB-Baseline                 MI.16.East.Lansing.SB 1.000000e+00
    ## MX.2017.ca.DB-Baseline                         WA.08.fld10.DB 1.000000e+00
    ## BR.2010.lam.DB-Baseline                        ND.10.fld25.DB 1.000000e+00
    ## MX.2017.ley.DB-Baseline                     MI.15.Montcalm.SB 1.000000e+00
    ## MI.17.Montcalm.SB-WI.17.Hancock.SB            IA.17.Nashua.SB 1.000000e+00
    ## BR.2010.por.DB-WI.17.Hancock.SB                      Baseline 1.000000e+00
    ## BR.2010.unai.DB-WI.17.Hancock.SB             WI.17.Hancock.SB 1.000000e+00
    ## MX.2017.pre.DB-WI.17.Hancock.SB             MI.17.Montcalm.SB 1.000000e+00
    ## MX.2017.ca.DB-WI.17.Hancock.SB                 BR.2010.por.DB 1.000000e+00
    ## BR.2010.lam.DB-WI.17.Hancock.SB               BR.2010.unai.DB 1.000000e+00
    ## MX.2017.ley.DB-WI.17.Hancock.SB                MX.2017.pre.DB 1.000000e+00
    ## BR.2010.por.DB-MI.17.Montcalm.SB                MX.2017.ca.DB 1.000000e+00
    ## BR.2010.unai.DB-MI.17.Montcalm.SB              BR.2010.lam.DB 1.000000e+00
    ## MX.2017.pre.DB-MI.17.Montcalm.SB               MX.2017.ley.DB 1.000000e+00
    ## MX.2017.ca.DB-MI.17.Montcalm.SB                NE.10.fld21.DB 1.000000e+00
    ## BR.2010.lam.DB-MI.17.Montcalm.SB               ND.07.fld01.DB 1.000000e+00
    ## MX.2017.ley.DB-MI.17.Montcalm.SB               WA.07.fld09.DB 1.000000e+00
    ## BR.2010.unai.DB-BR.2010.por.DB          MI.16.East.Lansing.SB 1.000000e+00
    ## MX.2017.pre.DB-BR.2010.por.DB                  WA.08.fld10.DB 1.000000e+00
    ## MX.2017.ca.DB-BR.2010.por.DB                   ND.10.fld25.DB 1.000000e+00
    ## BR.2010.lam.DB-BR.2010.por.DB               MI.15.Montcalm.SB 1.000000e+00
    ## MX.2017.ley.DB-BR.2010.por.DB                 IA.17.Nashua.SB 1.000000e+00
    ## MX.2017.pre.DB-BR.2010.unai.DB                       Baseline 1.000000e+00
    ## MX.2017.ca.DB-BR.2010.unai.DB                WI.17.Hancock.SB 1.000000e+00
    ## BR.2010.lam.DB-BR.2010.unai.DB              MI.17.Montcalm.SB 1.000000e+00
    ## MX.2017.ley.DB-BR.2010.unai.DB                 BR.2010.por.DB 1.000000e+00
    ## MX.2017.ca.DB-MX.2017.pre.DB                  BR.2010.unai.DB 1.000000e+00
    ## BR.2010.lam.DB-MX.2017.pre.DB                  MX.2017.pre.DB 1.000000e+00
    ## MX.2017.ley.DB-MX.2017.pre.DB                   MX.2017.ca.DB 1.000000e+00
    ## BR.2010.lam.DB-MX.2017.ca.DB                   BR.2010.lam.DB 1.000000e+00
    ## MX.2017.ley.DB-MX.2017.ca.DB                   MX.2017.ley.DB 1.000000e+00
    ## MX.2017.ley.DB-BR.2010.lam.DB                  NE.10.fld21.DB 1.000000e+00
    ##                                         labels
    ## ND.07.fld01.DB-NE.10.fld21.DB               ab
    ## WA.07.fld09.DB-NE.10.fld21.DB              abc
    ## MI.16.East.Lansing.SB-NE.10.fld21.DB      abcd
    ## WA.08.fld10.DB-NE.10.fld21.DB            abcde
    ## ND.10.fld25.DB-NE.10.fld21.DB            abcde
    ## MI.15.Montcalm.SB-NE.10.fld21.DB          acde
    ## IA.17.Nashua.SB-NE.10.fld21.DB            acde
    ## Baseline-NE.10.fld21.DB                   acde
    ## WI.17.Hancock.SB-NE.10.fld21.DB           acde
    ## MI.17.Montcalm.SB-NE.10.fld21.DB           cde
    ## BR.2010.por.DB-NE.10.fld21.DB              cde
    ## BR.2010.unai.DB-NE.10.fld21.DB              de
    ## MX.2017.pre.DB-NE.10.fld21.DB               de
    ## MX.2017.ca.DB-NE.10.fld21.DB                de
    ## BR.2010.lam.DB-NE.10.fld21.DB                e
    ## MX.2017.ley.DB-NE.10.fld21.DB                e
    ## WA.07.fld09.DB-ND.07.fld01.DB                b
    ## MI.16.East.Lansing.SB-ND.07.fld01.DB        ab
    ## WA.08.fld10.DB-ND.07.fld01.DB              abc
    ## ND.10.fld25.DB-ND.07.fld01.DB             abcd
    ## MI.15.Montcalm.SB-ND.07.fld01.DB         abcde
    ## IA.17.Nashua.SB-ND.07.fld01.DB           abcde
    ## Baseline-ND.07.fld01.DB                   acde
    ## WI.17.Hancock.SB-ND.07.fld01.DB           acde
    ## MI.17.Montcalm.SB-ND.07.fld01.DB          acde
    ## BR.2010.por.DB-ND.07.fld01.DB             acde
    ## BR.2010.unai.DB-ND.07.fld01.DB             cde
    ## MX.2017.pre.DB-ND.07.fld01.DB              cde
    ## MX.2017.ca.DB-ND.07.fld01.DB                de
    ## BR.2010.lam.DB-ND.07.fld01.DB               de
    ## MX.2017.ley.DB-ND.07.fld01.DB               de
    ## MI.16.East.Lansing.SB-WA.07.fld09.DB         e
    ## WA.08.fld10.DB-WA.07.fld09.DB                e
    ## ND.10.fld25.DB-WA.07.fld09.DB                b
    ## MI.15.Montcalm.SB-WA.07.fld09.DB            ab
    ## IA.17.Nashua.SB-WA.07.fld09.DB             abc
    ## Baseline-WA.07.fld09.DB                   abcd
    ## WI.17.Hancock.SB-WA.07.fld09.DB          abcde
    ## MI.17.Montcalm.SB-WA.07.fld09.DB         abcde
    ## BR.2010.por.DB-WA.07.fld09.DB             acde
    ## BR.2010.unai.DB-WA.07.fld09.DB            acde
    ## MX.2017.pre.DB-WA.07.fld09.DB             acde
    ## MX.2017.ca.DB-WA.07.fld09.DB              acde
    ## BR.2010.lam.DB-WA.07.fld09.DB              cde
    ## MX.2017.ley.DB-WA.07.fld09.DB              cde
    ## WA.08.fld10.DB-MI.16.East.Lansing.SB        de
    ## ND.10.fld25.DB-MI.16.East.Lansing.SB        de
    ## MI.15.Montcalm.SB-MI.16.East.Lansing.SB     de
    ## IA.17.Nashua.SB-MI.16.East.Lansing.SB        e
    ## Baseline-MI.16.East.Lansing.SB               e
    ## WI.17.Hancock.SB-MI.16.East.Lansing.SB       b
    ## MI.17.Montcalm.SB-MI.16.East.Lansing.SB     ab
    ## BR.2010.por.DB-MI.16.East.Lansing.SB       abc
    ## BR.2010.unai.DB-MI.16.East.Lansing.SB     abcd
    ## MX.2017.pre.DB-MI.16.East.Lansing.SB     abcde
    ## MX.2017.ca.DB-MI.16.East.Lansing.SB      abcde
    ## BR.2010.lam.DB-MI.16.East.Lansing.SB      acde
    ## MX.2017.ley.DB-MI.16.East.Lansing.SB      acde
    ## ND.10.fld25.DB-WA.08.fld10.DB             acde
    ## MI.15.Montcalm.SB-WA.08.fld10.DB          acde
    ## IA.17.Nashua.SB-WA.08.fld10.DB             cde
    ## Baseline-WA.08.fld10.DB                    cde
    ## WI.17.Hancock.SB-WA.08.fld10.DB             de
    ## MI.17.Montcalm.SB-WA.08.fld10.DB            de
    ## BR.2010.por.DB-WA.08.fld10.DB               de
    ## BR.2010.unai.DB-WA.08.fld10.DB               e
    ## MX.2017.pre.DB-WA.08.fld10.DB                e
    ## MX.2017.ca.DB-WA.08.fld10.DB                 b
    ## BR.2010.lam.DB-WA.08.fld10.DB               ab
    ## MX.2017.ley.DB-WA.08.fld10.DB              abc
    ## MI.15.Montcalm.SB-ND.10.fld25.DB          abcd
    ## IA.17.Nashua.SB-ND.10.fld25.DB           abcde
    ## Baseline-ND.10.fld25.DB                  abcde
    ## WI.17.Hancock.SB-ND.10.fld25.DB           acde
    ## MI.17.Montcalm.SB-ND.10.fld25.DB          acde
    ## BR.2010.por.DB-ND.10.fld25.DB             acde
    ## BR.2010.unai.DB-ND.10.fld25.DB            acde
    ## MX.2017.pre.DB-ND.10.fld25.DB              cde
    ## MX.2017.ca.DB-ND.10.fld25.DB               cde
    ## BR.2010.lam.DB-ND.10.fld25.DB               de
    ## MX.2017.ley.DB-ND.10.fld25.DB               de
    ## IA.17.Nashua.SB-MI.15.Montcalm.SB           de
    ## Baseline-MI.15.Montcalm.SB                   e
    ## WI.17.Hancock.SB-MI.15.Montcalm.SB           e
    ## MI.17.Montcalm.SB-MI.15.Montcalm.SB          b
    ## BR.2010.por.DB-MI.15.Montcalm.SB            ab
    ## BR.2010.unai.DB-MI.15.Montcalm.SB          abc
    ## MX.2017.pre.DB-MI.15.Montcalm.SB          abcd
    ## MX.2017.ca.DB-MI.15.Montcalm.SB          abcde
    ## BR.2010.lam.DB-MI.15.Montcalm.SB         abcde
    ## MX.2017.ley.DB-MI.15.Montcalm.SB          acde
    ## Baseline-IA.17.Nashua.SB                  acde
    ## WI.17.Hancock.SB-IA.17.Nashua.SB          acde
    ## MI.17.Montcalm.SB-IA.17.Nashua.SB         acde
    ## BR.2010.por.DB-IA.17.Nashua.SB             cde
    ## BR.2010.unai.DB-IA.17.Nashua.SB            cde
    ## MX.2017.pre.DB-IA.17.Nashua.SB              de
    ## MX.2017.ca.DB-IA.17.Nashua.SB               de
    ## BR.2010.lam.DB-IA.17.Nashua.SB              de
    ## MX.2017.ley.DB-IA.17.Nashua.SB               e
    ## WI.17.Hancock.SB-Baseline                    e
    ## MI.17.Montcalm.SB-Baseline                   b
    ## BR.2010.por.DB-Baseline                     ab
    ## BR.2010.unai.DB-Baseline                   abc
    ## MX.2017.pre.DB-Baseline                   abcd
    ## MX.2017.ca.DB-Baseline                   abcde
    ## BR.2010.lam.DB-Baseline                  abcde
    ## MX.2017.ley.DB-Baseline                   acde
    ## MI.17.Montcalm.SB-WI.17.Hancock.SB        acde
    ## BR.2010.por.DB-WI.17.Hancock.SB           acde
    ## BR.2010.unai.DB-WI.17.Hancock.SB          acde
    ## MX.2017.pre.DB-WI.17.Hancock.SB            cde
    ## MX.2017.ca.DB-WI.17.Hancock.SB             cde
    ## BR.2010.lam.DB-WI.17.Hancock.SB             de
    ## MX.2017.ley.DB-WI.17.Hancock.SB             de
    ## BR.2010.por.DB-MI.17.Montcalm.SB            de
    ## BR.2010.unai.DB-MI.17.Montcalm.SB            e
    ## MX.2017.pre.DB-MI.17.Montcalm.SB             e
    ## MX.2017.ca.DB-MI.17.Montcalm.SB              b
    ## BR.2010.lam.DB-MI.17.Montcalm.SB            ab
    ## MX.2017.ley.DB-MI.17.Montcalm.SB           abc
    ## BR.2010.unai.DB-BR.2010.por.DB            abcd
    ## MX.2017.pre.DB-BR.2010.por.DB            abcde
    ## MX.2017.ca.DB-BR.2010.por.DB             abcde
    ## BR.2010.lam.DB-BR.2010.por.DB             acde
    ## MX.2017.ley.DB-BR.2010.por.DB             acde
    ## MX.2017.pre.DB-BR.2010.unai.DB            acde
    ## MX.2017.ca.DB-BR.2010.unai.DB             acde
    ## BR.2010.lam.DB-BR.2010.unai.DB             cde
    ## MX.2017.ley.DB-BR.2010.unai.DB             cde
    ## MX.2017.ca.DB-MX.2017.pre.DB                de
    ## BR.2010.lam.DB-MX.2017.pre.DB               de
    ## MX.2017.ley.DB-MX.2017.pre.DB               de
    ## BR.2010.lam.DB-MX.2017.ca.DB                 e
    ## MX.2017.ley.DB-MX.2017.ca.DB                 e
    ## MX.2017.ley.DB-BR.2010.lam.DB                b

``` r
survey.boscalid.fields.2 <- survey.boscalid.fields %>% 
mutate(Field_difference = ifelse(
    Field_Year == "Baseline"| Field_Year == "NE.10.fld21.DB"|
      Field_Year == "ND.07.fld01.DB"|Field_Year == "WA.07.fld09.DB",
    "Yes","No")) %>% 
  mutate(Country_difference = ifelse(
    country == "USA",
    "Yes","No"))
  
####
myplot_comparison(survey.boscalid.fields.2,highlabel = 0.01) + expand_limits(y = c(0.06, 0.18)) +
  scale_y_continuous(breaks = c(0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18)) + geom_hline(yintercept =
                                                                                         0.1606226, linetype = "dashed")# make dynamic the value of intercept, create a new fucntion
```

![](READM_files/figure-gfm/Stat%20Anal%20Boscalid%20By%20Field-1.png)<!-- -->

## Tetraconazole

### By Source

``` r
test_kruskal_source(survey.tetraconazole.complete  %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() )
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  as.numeric(vardep) by as.factor(varindep)
    ## Kruskal-Wallis chi-squared = 2.6865, df = 2, p-value = 0.261

``` r
test_Dunn_source(survey.tetraconazole.complete  %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() )
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                                        mean.rank.diff   pval    
    ## Fungicide Field Trials-Baseline              37.89652 0.5794    
    ## Producer Fields-Baseline                     10.90838 1.0000    
    ## Producer Fields-Fungicide Field Trials      -26.98814 0.3759    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
generate_label(test_Dunn_source(survey.tetraconazole.complete  %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() ))
```

    ##                                                   plot.labels Dunn.levels
    ## Fungicide Field Trials-Baseline        Fungicide Field Trials   0.5793826
    ## Producer Fields-Baseline                      Producer Fields   1.0000000
    ## Producer Fields-Fungicide Field Trials               Baseline   0.3758844
    ##                                        labels
    ## Fungicide Field Trials-Baseline             a
    ## Producer Fields-Baseline                    a
    ## Producer Fields-Fungicide Field Trials      a

``` r
survey.tetraconazole.complete %>% 
 group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() %>% 
  group_by(source) %>% summarize(range(EC50DC), average= mean(EC50DC), N= n(),  sd = sd(EC50DC), se = round(sd/sqrt(N),3))
```

    ## `summarise()` has grouped output by 'source'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 6 × 6
    ## # Groups:   source [3]
    ##   source                 `range(EC50DC)` average     N    sd    se
    ##   <fct>                            <dbl>   <dbl> <int> <dbl> <dbl>
    ## 1 Baseline                         0.372    1.03    22 0.368 0.079
    ## 2 Baseline                         1.58     1.03    22 0.368 0.079
    ## 3 Fungicide Field Trials           0.690    1.19    47 0.183 0.027
    ## 4 Fungicide Field Trials           1.83     1.19    47 0.183 0.027
    ## 5 Producer Fields                  0.197    1.04   321 0.416 0.023
    ## 6 Producer Fields                  2.27     1.04   321 0.416 0.023

### By Host

``` r
#generate_label
  test_kruskal_host(survey.tetraconazole.complete.host %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup())
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  as.numeric(vardep) by as.factor(varindep)
    ## Kruskal-Wallis chi-squared = 2.93, df = 1, p-value = 0.08694

``` r
#generate_label
  test_Dunn_host(survey.tetraconazole.complete.host %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup())
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                 mean.rank.diff   pval    
    ## Soybean-Drybean       19.57865 0.0869 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
  #)

survey.tetraconazole.complete.host %>%
group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() %>% 
  group_by(Host) %>% summarize(range(EC50DC), average= mean(EC50DC), N= n(),  sd = sd(EC50DC), se = round(sd/sqrt(N),3))
```

    ## `summarise()` has grouped output by 'Host'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 4 × 6
    ## # Groups:   Host [2]
    ##   Host    `range(EC50DC)` average     N    sd    se
    ##   <fct>             <dbl>   <dbl> <int> <dbl> <dbl>
    ## 1 Drybean           0.197    1.03   229 0.442 0.029
    ## 2 Drybean           2.27     1.03   229 0.442 0.029
    ## 3 Soybean           0.361    1.12   139 0.302 0.026
    ## 4 Soybean           1.83     1.12   139 0.302 0.026

### By Field

``` r
generate_label( test_field (survey.tetraconazole.fields))
```

    ##                                                   plot.labels  Dunn.levels
    ## ND.07.fld01.DB-NE.10.fld21.DB                  ND.07.fld01.DB 1.000000e+00
    ## MX.2017.pre.DB-NE.10.fld21.DB                  MX.2017.pre.DB 1.000000e+00
    ## MX.2017.ca.DB-NE.10.fld21.DB                    MX.2017.ca.DB 1.000000e+00
    ## Baseline-NE.10.fld21.DB                              Baseline 7.288913e-02
    ## MX.2017.ley.DB-NE.10.fld21.DB                  MX.2017.ley.DB 1.450813e-01
    ## BR.2010.por.DB-NE.10.fld21.DB                  BR.2010.por.DB 6.266112e-01
    ## MI.17.Montcalm.SB-NE.10.fld21.DB            MI.17.Montcalm.SB 5.673756e-02
    ## WA.07.fld09.DB-NE.10.fld21.DB                  WA.07.fld09.DB 2.160034e-02
    ## WA.08.fld10.DB-NE.10.fld21.DB                  WA.08.fld10.DB 6.496539e-03
    ## IA.17.Nashua.SB-NE.10.fld21.DB                IA.17.Nashua.SB 1.536112e-03
    ## MI.16.East.Lansing.SB-NE.10.fld21.DB    MI.16.East.Lansing.SB 1.357638e-02
    ## BR.2010.lam.DB-NE.10.fld21.DB                  BR.2010.lam.DB 2.558938e-04
    ## ND.10.fld25.DB-NE.10.fld21.DB                  ND.10.fld25.DB 1.804037e-05
    ## BR.2010.unai.DB-NE.10.fld21.DB                BR.2010.unai.DB 9.497434e-07
    ## MX.2017.pre.DB-ND.07.fld01.DB                  NE.10.fld21.DB 1.000000e+00
    ## MX.2017.ca.DB-ND.07.fld01.DB                   ND.07.fld01.DB 1.000000e+00
    ## Baseline-ND.07.fld01.DB                        MX.2017.pre.DB 1.000000e+00
    ## MX.2017.ley.DB-ND.07.fld01.DB                   MX.2017.ca.DB 1.000000e+00
    ## BR.2010.por.DB-ND.07.fld01.DB                        Baseline 1.000000e+00
    ## MI.17.Montcalm.SB-ND.07.fld01.DB               MX.2017.ley.DB 1.000000e+00
    ## WA.07.fld09.DB-ND.07.fld01.DB                  BR.2010.por.DB 1.000000e+00
    ## WA.08.fld10.DB-ND.07.fld01.DB               MI.17.Montcalm.SB 1.000000e+00
    ## IA.17.Nashua.SB-ND.07.fld01.DB                 WA.07.fld09.DB 3.232033e-01
    ## MI.16.East.Lansing.SB-ND.07.fld01.DB           WA.08.fld10.DB 1.000000e+00
    ## BR.2010.lam.DB-ND.07.fld01.DB                 IA.17.Nashua.SB 6.508085e-02
    ## ND.10.fld25.DB-ND.07.fld01.DB           MI.16.East.Lansing.SB 7.468676e-03
    ## BR.2010.unai.DB-ND.07.fld01.DB                 BR.2010.lam.DB 7.015396e-04
    ## MX.2017.ca.DB-MX.2017.pre.DB                   ND.10.fld25.DB 1.000000e+00
    ## Baseline-MX.2017.pre.DB                       BR.2010.unai.DB 1.000000e+00
    ## MX.2017.ley.DB-MX.2017.pre.DB                  NE.10.fld21.DB 1.000000e+00
    ## BR.2010.por.DB-MX.2017.pre.DB                  ND.07.fld01.DB 1.000000e+00
    ## MI.17.Montcalm.SB-MX.2017.pre.DB               MX.2017.pre.DB 1.000000e+00
    ## WA.07.fld09.DB-MX.2017.pre.DB                   MX.2017.ca.DB 1.000000e+00
    ## WA.08.fld10.DB-MX.2017.pre.DB                        Baseline 9.410956e-01
    ## IA.17.Nashua.SB-MX.2017.pre.DB                 MX.2017.ley.DB 3.054874e-01
    ## MI.16.East.Lansing.SB-MX.2017.pre.DB           BR.2010.por.DB 1.000000e+00
    ## BR.2010.lam.DB-MX.2017.pre.DB               MI.17.Montcalm.SB 6.406962e-02
    ## ND.10.fld25.DB-MX.2017.pre.DB                  WA.07.fld09.DB 8.050138e-03
    ## BR.2010.unai.DB-MX.2017.pre.DB                 WA.08.fld10.DB 8.690408e-04
    ## Baseline-MX.2017.ca.DB                        IA.17.Nashua.SB 1.000000e+00
    ## MX.2017.ley.DB-MX.2017.ca.DB            MI.16.East.Lansing.SB 1.000000e+00
    ## BR.2010.por.DB-MX.2017.ca.DB                   BR.2010.lam.DB 1.000000e+00
    ## MI.17.Montcalm.SB-MX.2017.ca.DB                ND.10.fld25.DB 1.000000e+00
    ## WA.07.fld09.DB-MX.2017.ca.DB                  BR.2010.unai.DB 1.000000e+00
    ## WA.08.fld10.DB-MX.2017.ca.DB                   NE.10.fld21.DB 1.000000e+00
    ## IA.17.Nashua.SB-MX.2017.ca.DB                  ND.07.fld01.DB 8.595254e-01
    ## MI.16.East.Lansing.SB-MX.2017.ca.DB            MX.2017.pre.DB 1.000000e+00
    ## BR.2010.lam.DB-MX.2017.ca.DB                    MX.2017.ca.DB 1.960945e-01
    ## ND.10.fld25.DB-MX.2017.ca.DB                         Baseline 2.884335e-02
    ## BR.2010.unai.DB-MX.2017.ca.DB                  MX.2017.ley.DB 3.764441e-03
    ## MX.2017.ley.DB-Baseline                        BR.2010.por.DB 1.000000e+00
    ## BR.2010.por.DB-Baseline                     MI.17.Montcalm.SB 1.000000e+00
    ## MI.17.Montcalm.SB-Baseline                     WA.07.fld09.DB 1.000000e+00
    ## WA.07.fld09.DB-Baseline                        WA.08.fld10.DB 1.000000e+00
    ## WA.08.fld10.DB-Baseline                       IA.17.Nashua.SB 1.000000e+00
    ## IA.17.Nashua.SB-Baseline                MI.16.East.Lansing.SB 1.000000e+00
    ## MI.16.East.Lansing.SB-Baseline                 BR.2010.lam.DB 1.000000e+00
    ## BR.2010.lam.DB-Baseline                        ND.10.fld25.DB 1.000000e+00
    ## ND.10.fld25.DB-Baseline                       BR.2010.unai.DB 4.796523e-01
    ## BR.2010.unai.DB-Baseline                       NE.10.fld21.DB 6.935555e-02
    ## BR.2010.por.DB-MX.2017.ley.DB                  ND.07.fld01.DB 1.000000e+00
    ## MI.17.Montcalm.SB-MX.2017.ley.DB               MX.2017.pre.DB 1.000000e+00
    ## WA.07.fld09.DB-MX.2017.ley.DB                   MX.2017.ca.DB 1.000000e+00
    ## WA.08.fld10.DB-MX.2017.ley.DB                        Baseline 1.000000e+00
    ## IA.17.Nashua.SB-MX.2017.ley.DB                 MX.2017.ley.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-MX.2017.ley.DB           BR.2010.por.DB 1.000000e+00
    ## BR.2010.lam.DB-MX.2017.ley.DB               MI.17.Montcalm.SB 1.000000e+00
    ## ND.10.fld25.DB-MX.2017.ley.DB                  WA.07.fld09.DB 1.000000e+00
    ## BR.2010.unai.DB-MX.2017.ley.DB                 WA.08.fld10.DB 6.181705e-01
    ## MI.17.Montcalm.SB-BR.2010.por.DB              IA.17.Nashua.SB 1.000000e+00
    ## WA.07.fld09.DB-BR.2010.por.DB           MI.16.East.Lansing.SB 1.000000e+00
    ## WA.08.fld10.DB-BR.2010.por.DB                  BR.2010.lam.DB 1.000000e+00
    ## IA.17.Nashua.SB-BR.2010.por.DB                 ND.10.fld25.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-BR.2010.por.DB          BR.2010.unai.DB 1.000000e+00
    ## BR.2010.lam.DB-BR.2010.por.DB                  NE.10.fld21.DB 1.000000e+00
    ## ND.10.fld25.DB-BR.2010.por.DB                  ND.07.fld01.DB 9.403029e-01
    ## BR.2010.unai.DB-BR.2010.por.DB                 MX.2017.pre.DB 2.171075e-01
    ## WA.07.fld09.DB-MI.17.Montcalm.SB                MX.2017.ca.DB 1.000000e+00
    ## WA.08.fld10.DB-MI.17.Montcalm.SB                     Baseline 1.000000e+00
    ## IA.17.Nashua.SB-MI.17.Montcalm.SB              MX.2017.ley.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-MI.17.Montcalm.SB        BR.2010.por.DB 1.000000e+00
    ## BR.2010.lam.DB-MI.17.Montcalm.SB            MI.17.Montcalm.SB 1.000000e+00
    ## ND.10.fld25.DB-MI.17.Montcalm.SB               WA.07.fld09.DB 1.000000e+00
    ## BR.2010.unai.DB-MI.17.Montcalm.SB              WA.08.fld10.DB 1.000000e+00
    ## WA.08.fld10.DB-WA.07.fld09.DB                 IA.17.Nashua.SB 1.000000e+00
    ## IA.17.Nashua.SB-WA.07.fld09.DB          MI.16.East.Lansing.SB 1.000000e+00
    ## MI.16.East.Lansing.SB-WA.07.fld09.DB           BR.2010.lam.DB 1.000000e+00
    ## BR.2010.lam.DB-WA.07.fld09.DB                  ND.10.fld25.DB 1.000000e+00
    ## ND.10.fld25.DB-WA.07.fld09.DB                 BR.2010.unai.DB 1.000000e+00
    ## BR.2010.unai.DB-WA.07.fld09.DB                 NE.10.fld21.DB 6.968058e-01
    ## IA.17.Nashua.SB-WA.08.fld10.DB                 ND.07.fld01.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-WA.08.fld10.DB           MX.2017.pre.DB 1.000000e+00
    ## BR.2010.lam.DB-WA.08.fld10.DB                   MX.2017.ca.DB 1.000000e+00
    ## ND.10.fld25.DB-WA.08.fld10.DB                        Baseline 1.000000e+00
    ## BR.2010.unai.DB-WA.08.fld10.DB                 MX.2017.ley.DB 1.000000e+00
    ## MI.16.East.Lansing.SB-IA.17.Nashua.SB          BR.2010.por.DB 1.000000e+00
    ## BR.2010.lam.DB-IA.17.Nashua.SB              MI.17.Montcalm.SB 1.000000e+00
    ## ND.10.fld25.DB-IA.17.Nashua.SB                 WA.07.fld09.DB 1.000000e+00
    ## BR.2010.unai.DB-IA.17.Nashua.SB                WA.08.fld10.DB 1.000000e+00
    ## BR.2010.lam.DB-MI.16.East.Lansing.SB          IA.17.Nashua.SB 1.000000e+00
    ## ND.10.fld25.DB-MI.16.East.Lansing.SB    MI.16.East.Lansing.SB 1.000000e+00
    ## BR.2010.unai.DB-MI.16.East.Lansing.SB          BR.2010.lam.DB 1.000000e+00
    ## ND.10.fld25.DB-BR.2010.lam.DB                  ND.10.fld25.DB 1.000000e+00
    ## BR.2010.unai.DB-BR.2010.lam.DB                BR.2010.unai.DB 1.000000e+00
    ## BR.2010.unai.DB-ND.10.fld25.DB                 NE.10.fld21.DB 1.000000e+00
    ##                                         labels
    ## ND.07.fld01.DB-NE.10.fld21.DB               ab
    ## MX.2017.pre.DB-NE.10.fld21.DB               ab
    ## MX.2017.ca.DB-NE.10.fld21.DB                ab
    ## Baseline-NE.10.fld21.DB                    abc
    ## MX.2017.ley.DB-NE.10.fld21.DB              abc
    ## BR.2010.por.DB-NE.10.fld21.DB              abc
    ## MI.17.Montcalm.SB-NE.10.fld21.DB           abc
    ## WA.07.fld09.DB-NE.10.fld21.DB               ac
    ## WA.08.fld10.DB-NE.10.fld21.DB               ac
    ## IA.17.Nashua.SB-NE.10.fld21.DB              ac
    ## MI.16.East.Lansing.SB-NE.10.fld21.DB        ac
    ## BR.2010.lam.DB-NE.10.fld21.DB               ac
    ## ND.10.fld25.DB-NE.10.fld21.DB                c
    ## BR.2010.unai.DB-NE.10.fld21.DB               c
    ## MX.2017.pre.DB-ND.07.fld01.DB                b
    ## MX.2017.ca.DB-ND.07.fld01.DB                ab
    ## Baseline-ND.07.fld01.DB                     ab
    ## MX.2017.ley.DB-ND.07.fld01.DB               ab
    ## BR.2010.por.DB-ND.07.fld01.DB              abc
    ## MI.17.Montcalm.SB-ND.07.fld01.DB           abc
    ## WA.07.fld09.DB-ND.07.fld01.DB              abc
    ## WA.08.fld10.DB-ND.07.fld01.DB              abc
    ## IA.17.Nashua.SB-ND.07.fld01.DB              ac
    ## MI.16.East.Lansing.SB-ND.07.fld01.DB        ac
    ## BR.2010.lam.DB-ND.07.fld01.DB               ac
    ## ND.10.fld25.DB-ND.07.fld01.DB               ac
    ## BR.2010.unai.DB-ND.07.fld01.DB              ac
    ## MX.2017.ca.DB-MX.2017.pre.DB                 c
    ## Baseline-MX.2017.pre.DB                      c
    ## MX.2017.ley.DB-MX.2017.pre.DB                b
    ## BR.2010.por.DB-MX.2017.pre.DB               ab
    ## MI.17.Montcalm.SB-MX.2017.pre.DB            ab
    ## WA.07.fld09.DB-MX.2017.pre.DB               ab
    ## WA.08.fld10.DB-MX.2017.pre.DB              abc
    ## IA.17.Nashua.SB-MX.2017.pre.DB             abc
    ## MI.16.East.Lansing.SB-MX.2017.pre.DB       abc
    ## BR.2010.lam.DB-MX.2017.pre.DB              abc
    ## ND.10.fld25.DB-MX.2017.pre.DB               ac
    ## BR.2010.unai.DB-MX.2017.pre.DB              ac
    ## Baseline-MX.2017.ca.DB                      ac
    ## MX.2017.ley.DB-MX.2017.ca.DB                ac
    ## BR.2010.por.DB-MX.2017.ca.DB                ac
    ## MI.17.Montcalm.SB-MX.2017.ca.DB              c
    ## WA.07.fld09.DB-MX.2017.ca.DB                 c
    ## WA.08.fld10.DB-MX.2017.ca.DB                 b
    ## IA.17.Nashua.SB-MX.2017.ca.DB               ab
    ## MI.16.East.Lansing.SB-MX.2017.ca.DB         ab
    ## BR.2010.lam.DB-MX.2017.ca.DB                ab
    ## ND.10.fld25.DB-MX.2017.ca.DB               abc
    ## BR.2010.unai.DB-MX.2017.ca.DB              abc
    ## MX.2017.ley.DB-Baseline                    abc
    ## BR.2010.por.DB-Baseline                    abc
    ## MI.17.Montcalm.SB-Baseline                  ac
    ## WA.07.fld09.DB-Baseline                     ac
    ## WA.08.fld10.DB-Baseline                     ac
    ## IA.17.Nashua.SB-Baseline                    ac
    ## MI.16.East.Lansing.SB-Baseline              ac
    ## BR.2010.lam.DB-Baseline                      c
    ## ND.10.fld25.DB-Baseline                      c
    ## BR.2010.unai.DB-Baseline                     b
    ## BR.2010.por.DB-MX.2017.ley.DB               ab
    ## MI.17.Montcalm.SB-MX.2017.ley.DB            ab
    ## WA.07.fld09.DB-MX.2017.ley.DB               ab
    ## WA.08.fld10.DB-MX.2017.ley.DB              abc
    ## IA.17.Nashua.SB-MX.2017.ley.DB             abc
    ## MI.16.East.Lansing.SB-MX.2017.ley.DB       abc
    ## BR.2010.lam.DB-MX.2017.ley.DB              abc
    ## ND.10.fld25.DB-MX.2017.ley.DB               ac
    ## BR.2010.unai.DB-MX.2017.ley.DB              ac
    ## MI.17.Montcalm.SB-BR.2010.por.DB            ac
    ## WA.07.fld09.DB-BR.2010.por.DB               ac
    ## WA.08.fld10.DB-BR.2010.por.DB               ac
    ## IA.17.Nashua.SB-BR.2010.por.DB               c
    ## MI.16.East.Lansing.SB-BR.2010.por.DB         c
    ## BR.2010.lam.DB-BR.2010.por.DB                b
    ## ND.10.fld25.DB-BR.2010.por.DB               ab
    ## BR.2010.unai.DB-BR.2010.por.DB              ab
    ## WA.07.fld09.DB-MI.17.Montcalm.SB            ab
    ## WA.08.fld10.DB-MI.17.Montcalm.SB           abc
    ## IA.17.Nashua.SB-MI.17.Montcalm.SB          abc
    ## MI.16.East.Lansing.SB-MI.17.Montcalm.SB    abc
    ## BR.2010.lam.DB-MI.17.Montcalm.SB           abc
    ## ND.10.fld25.DB-MI.17.Montcalm.SB            ac
    ## BR.2010.unai.DB-MI.17.Montcalm.SB           ac
    ## WA.08.fld10.DB-WA.07.fld09.DB               ac
    ## IA.17.Nashua.SB-WA.07.fld09.DB              ac
    ## MI.16.East.Lansing.SB-WA.07.fld09.DB        ac
    ## BR.2010.lam.DB-WA.07.fld09.DB                c
    ## ND.10.fld25.DB-WA.07.fld09.DB                c
    ## BR.2010.unai.DB-WA.07.fld09.DB               b
    ## IA.17.Nashua.SB-WA.08.fld10.DB              ab
    ## MI.16.East.Lansing.SB-WA.08.fld10.DB        ab
    ## BR.2010.lam.DB-WA.08.fld10.DB               ab
    ## ND.10.fld25.DB-WA.08.fld10.DB              abc
    ## BR.2010.unai.DB-WA.08.fld10.DB             abc
    ## MI.16.East.Lansing.SB-IA.17.Nashua.SB      abc
    ## BR.2010.lam.DB-IA.17.Nashua.SB             abc
    ## ND.10.fld25.DB-IA.17.Nashua.SB              ac
    ## BR.2010.unai.DB-IA.17.Nashua.SB             ac
    ## BR.2010.lam.DB-MI.16.East.Lansing.SB        ac
    ## ND.10.fld25.DB-MI.16.East.Lansing.SB        ac
    ## BR.2010.unai.DB-MI.16.East.Lansing.SB       ac
    ## ND.10.fld25.DB-BR.2010.lam.DB                c
    ## BR.2010.unai.DB-BR.2010.lam.DB               c
    ## BR.2010.unai.DB-ND.10.fld25.DB               b

``` r
survey.tetraconazole.fields.2 <- survey.tetraconazole.fields %>% 
mutate(Field_difference = ifelse(
    Field_Year == "Baseline"| Field_Year == "NE.10.fld21.DB"|
      Field_Year == "BR.2010.unai.DB",
    "Yes","No")) %>% 
  mutate(Country_difference = ifelse(
    country == "USA",
    "Yes","No"))
####
myplot_comparison(survey.tetraconazole.fields.2, highlabel = 0.1)  + expand_limits(y = c(0.5, 1.6)) + scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)) 
```

![](READM_files/figure-gfm/Stat%20Anal%20Tetraconazole%20By%20Field-1.png)<!-- -->

## Picoxystrobin

### By Source

``` r
 test_kruskal_source(survey.picoxystrobin.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup())
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  as.numeric(vardep) by as.factor(varindep)
    ## Kruskal-Wallis chi-squared = 13.434, df = 2, p-value = 0.00121

``` r
 ##
test_Dunn_source(survey.picoxystrobin.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup())
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                                        mean.rank.diff   pval    
    ## Fungicide Field Trials-Baseline              50.17540 0.2006    
    ## Producer Fields-Baseline                     80.97263 0.0041 ** 
    ## Producer Fields-Fungicide Field Trials       30.79723 0.0876 .  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# generate the label
generate_label( test_Dunn_source(survey.picoxystrobin.complete %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup()))
```

    ##                                                   plot.labels Dunn.levels
    ## Fungicide Field Trials-Baseline        Fungicide Field Trials 0.200586495
    ## Producer Fields-Baseline                      Producer Fields 0.004142246
    ## Producer Fields-Fungicide Field Trials               Baseline 0.087632317
    ##                                        labels
    ## Fungicide Field Trials-Baseline            ab
    ## Producer Fields-Baseline                    a
    ## Producer Fields-Fungicide Field Trials      b

``` r
survey.picoxystrobin.complete %>% 
  group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup()%>% 
  group_by(source) %>% summarize(round(range(EC50DC),3), average= round(mean(EC50DC), 2), N= n(),  sd = sd(EC50DC), se = round(sd/sqrt(N),3))
```

    ## `summarise()` has grouped output by 'source'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 6 × 6
    ## # Groups:   source [3]
    ##   source                 `round(range(EC50DC), 3)` average     N      sd    se
    ##   <fct>                                      <dbl>   <dbl> <int>   <dbl> <dbl>
    ## 1 Baseline                                   0.011    0.01    22 0.00279 0.001
    ## 2 Baseline                                   0.021    0.01    22 0.00279 0.001
    ## 3 Fungicide Field Trials                     0.008    0.01    85 0.00241 0    
    ## 4 Fungicide Field Trials                     0.021    0.01    85 0.00241 0    
    ## 5 Producer Fields                            0.006    0.02   289 0.00322 0    
    ## 6 Producer Fields                            0.027    0.02   289 0.00322 0

### By Host

``` r
#generate_label
  test_kruskal_host(survey.picoxystrobin.complete.host  %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup())
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  as.numeric(vardep) by as.factor(varindep)
    ## Kruskal-Wallis chi-squared = 38.07, df = 1, p-value = 6.827e-10

``` r
#generate_label
  test_Dunn_host(survey.picoxystrobin.complete.host  %>% group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup())
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                 mean.rank.diff    pval    
    ## Soybean-Drybean      -70.70362 6.8e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
survey.picoxystrobin.complete.host %>%
   group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() %>% 
  group_by(Host) %>% summarize(round(range(EC50DC),3), average= round(mean(EC50DC),2), N= n(),  sd = sd(EC50DC), se = round(sd/sqrt(N),3))
```

    ## `summarise()` has grouped output by 'Host'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 4 × 6
    ## # Groups:   Host [2]
    ##   Host    `round(range(EC50DC), 3)` average     N      sd    se
    ##   <fct>                       <dbl>   <dbl> <int>   <dbl> <dbl>
    ## 1 Drybean                     0.009    0.02   228 0.00326     0
    ## 2 Drybean                     0.027    0.02   228 0.00326     0
    ## 3 Soybean                     0.006    0.01   146 0.00228     0
    ## 4 Soybean                     0.021    0.01   146 0.00228     0

### By Field

``` r
test_field (survey.picoxystrobin.fields)
```

    ## 
    ##  Dunn's test of multiple comparisons using rank sums : bonferroni  
    ## 
    ##                                         mean.rank.diff    pval    
    ## Baseline-WI.2015.cer.SB                    -0.36363636 1.00000    
    ## IA.17.Nashua.SB-WI.2015.cer.SB             10.97058824 1.00000    
    ## MI.17.Montcalm.SB-WI.2015.cer.SB           16.85294118 1.00000    
    ## ND.10.fld25.DB-WI.2015.cer.SB              30.27777778 1.00000    
    ## NE.10.fld21.DB-WI.2015.cer.SB              24.16666667 1.00000    
    ## WA.08.fld10.DB-WI.2015.cer.SB              42.25000000 1.00000    
    ## MI.16.East.Lansing.SB-WI.2015.cer.SB       36.93750000 1.00000    
    ## BR.2010.lam.DB-WI.2015.cer.SB              48.83333333 1.00000    
    ## WI.17.Hancock.SB-WI.2015.cer.SB            48.80769231 1.00000    
    ## MI.15.Montcalm.SB-WI.2015.cer.SB           55.30000000 1.00000    
    ## WA.07.fld09.DB-WI.2015.cer.SB              60.56250000 1.00000    
    ## MX.2017.pre.DB-WI.2015.cer.SB              63.30000000 1.00000    
    ## BR.2010.por.DB-WI.2015.cer.SB              73.00000000 1.00000    
    ## ND.07.fld01.DB-WI.2015.cer.SB             100.75000000 0.04059 *  
    ## BR.2010.unai.DB-WI.2015.cer.SB             84.00000000 0.55079    
    ## MX.2017.ca.DB-WI.2015.cer.SB              113.70000000 0.01244 *  
    ## MX.2017.ley.DB-WI.2015.cer.SB             126.77272727 0.00105 ** 
    ## IA.17.Nashua.SB-Baseline                   11.33422460 1.00000    
    ## MI.17.Montcalm.SB-Baseline                 17.21657754 1.00000    
    ## ND.10.fld25.DB-Baseline                    30.64141414 1.00000    
    ## NE.10.fld21.DB-Baseline                    24.53030303 1.00000    
    ## WA.08.fld10.DB-Baseline                    42.61363636 1.00000    
    ## MI.16.East.Lansing.SB-Baseline             37.30113636 1.00000    
    ## BR.2010.lam.DB-Baseline                    49.19696970 1.00000    
    ## WI.17.Hancock.SB-Baseline                  49.17132867 1.00000    
    ## MI.15.Montcalm.SB-Baseline                 55.66363636 1.00000    
    ## WA.07.fld09.DB-Baseline                    60.92613636 0.61999    
    ## MX.2017.pre.DB-Baseline                    63.66363636 1.00000    
    ## BR.2010.por.DB-Baseline                    73.36363636 0.43890    
    ## ND.07.fld01.DB-Baseline                   101.11363636 0.00193 ** 
    ## BR.2010.unai.DB-Baseline                   84.36363636 0.09286 .  
    ## MX.2017.ca.DB-Baseline                    114.06363636 0.00054 ***
    ## MX.2017.ley.DB-Baseline                   127.13636364 1.5e-05 ***
    ## MI.17.Montcalm.SB-IA.17.Nashua.SB           5.88235294 1.00000    
    ## ND.10.fld25.DB-IA.17.Nashua.SB             19.30718954 1.00000    
    ## NE.10.fld21.DB-IA.17.Nashua.SB             13.19607843 1.00000    
    ## WA.08.fld10.DB-IA.17.Nashua.SB             31.27941176 1.00000    
    ## MI.16.East.Lansing.SB-IA.17.Nashua.SB      25.96691176 1.00000    
    ## BR.2010.lam.DB-IA.17.Nashua.SB             37.86274510 1.00000    
    ## WI.17.Hancock.SB-IA.17.Nashua.SB           37.83710407 1.00000    
    ## MI.15.Montcalm.SB-IA.17.Nashua.SB          44.32941176 1.00000    
    ## WA.07.fld09.DB-IA.17.Nashua.SB             49.59191176 1.00000    
    ## MX.2017.pre.DB-IA.17.Nashua.SB             52.32941176 1.00000    
    ## BR.2010.por.DB-IA.17.Nashua.SB             62.02941176 1.00000    
    ## ND.07.fld01.DB-IA.17.Nashua.SB             89.77941176 0.03422 *  
    ## BR.2010.unai.DB-IA.17.Nashua.SB            73.02941176 0.68972    
    ## MX.2017.ca.DB-IA.17.Nashua.SB             102.72941176 0.00988 ** 
    ## MX.2017.ley.DB-IA.17.Nashua.SB            115.80213904 0.00054 ***
    ## ND.10.fld25.DB-MI.17.Montcalm.SB           13.42483660 1.00000    
    ## NE.10.fld21.DB-MI.17.Montcalm.SB            7.31372549 1.00000    
    ## WA.08.fld10.DB-MI.17.Montcalm.SB           25.39705882 1.00000    
    ## MI.16.East.Lansing.SB-MI.17.Montcalm.SB    20.08455882 1.00000    
    ## BR.2010.lam.DB-MI.17.Montcalm.SB           31.98039216 1.00000    
    ## WI.17.Hancock.SB-MI.17.Montcalm.SB         31.95475113 1.00000    
    ## MI.15.Montcalm.SB-MI.17.Montcalm.SB        38.44705882 1.00000    
    ## WA.07.fld09.DB-MI.17.Montcalm.SB           43.70955882 1.00000    
    ## MX.2017.pre.DB-MI.17.Montcalm.SB           46.44705882 1.00000    
    ## BR.2010.por.DB-MI.17.Montcalm.SB           56.14705882 1.00000    
    ## ND.07.fld01.DB-MI.17.Montcalm.SB           83.89705882 0.08613 .  
    ## BR.2010.unai.DB-MI.17.Montcalm.SB          67.14705882 1.00000    
    ## MX.2017.ca.DB-MI.17.Montcalm.SB            96.84705882 0.02532 *  
    ## MX.2017.ley.DB-MI.17.Montcalm.SB          109.91978610 0.00163 ** 
    ## NE.10.fld21.DB-ND.10.fld25.DB              -6.11111111 1.00000    
    ## WA.08.fld10.DB-ND.10.fld25.DB              11.97222222 1.00000    
    ## MI.16.East.Lansing.SB-ND.10.fld25.DB        6.65972222 1.00000    
    ## BR.2010.lam.DB-ND.10.fld25.DB              18.55555556 1.00000    
    ## WI.17.Hancock.SB-ND.10.fld25.DB            18.52991453 1.00000    
    ## MI.15.Montcalm.SB-ND.10.fld25.DB           25.02222222 1.00000    
    ## WA.07.fld09.DB-ND.10.fld25.DB              30.28472222 1.00000    
    ## MX.2017.pre.DB-ND.10.fld25.DB              33.02222222 1.00000    
    ## BR.2010.por.DB-ND.10.fld25.DB              42.72222222 1.00000    
    ## ND.07.fld01.DB-ND.10.fld25.DB              70.47222222 1.00000    
    ## BR.2010.unai.DB-ND.10.fld25.DB             53.72222222 1.00000    
    ## MX.2017.ca.DB-ND.10.fld25.DB               83.42222222 0.74837    
    ## MX.2017.ley.DB-ND.10.fld25.DB              96.49494949 0.13406    
    ## WA.08.fld10.DB-NE.10.fld21.DB              18.08333333 1.00000    
    ## MI.16.East.Lansing.SB-NE.10.fld21.DB       12.77083333 1.00000    
    ## BR.2010.lam.DB-NE.10.fld21.DB              24.66666667 1.00000    
    ## WI.17.Hancock.SB-NE.10.fld21.DB            24.64102564 1.00000    
    ## MI.15.Montcalm.SB-NE.10.fld21.DB           31.13333333 1.00000    
    ## WA.07.fld09.DB-NE.10.fld21.DB              36.39583333 1.00000    
    ## MX.2017.pre.DB-NE.10.fld21.DB              39.13333333 1.00000    
    ## BR.2010.por.DB-NE.10.fld21.DB              48.83333333 1.00000    
    ## ND.07.fld01.DB-NE.10.fld21.DB              76.58333333 1.00000    
    ## BR.2010.unai.DB-NE.10.fld21.DB             59.83333333 1.00000    
    ## MX.2017.ca.DB-NE.10.fld21.DB               89.53333333 0.38642    
    ## MX.2017.ley.DB-NE.10.fld21.DB             102.60606061 0.06162 .  
    ## MI.16.East.Lansing.SB-WA.08.fld10.DB       -5.31250000 1.00000    
    ## BR.2010.lam.DB-WA.08.fld10.DB               6.58333333 1.00000    
    ## WI.17.Hancock.SB-WA.08.fld10.DB             6.55769231 1.00000    
    ## MI.15.Montcalm.SB-WA.08.fld10.DB           13.05000000 1.00000    
    ## WA.07.fld09.DB-WA.08.fld10.DB              18.31250000 1.00000    
    ## MX.2017.pre.DB-WA.08.fld10.DB              21.05000000 1.00000    
    ## BR.2010.por.DB-WA.08.fld10.DB              30.75000000 1.00000    
    ## ND.07.fld01.DB-WA.08.fld10.DB              58.50000000 1.00000    
    ## BR.2010.unai.DB-WA.08.fld10.DB             41.75000000 1.00000    
    ## MX.2017.ca.DB-WA.08.fld10.DB               71.45000000 1.00000    
    ## MX.2017.ley.DB-WA.08.fld10.DB              84.52272727 0.25989    
    ## BR.2010.lam.DB-MI.16.East.Lansing.SB       11.89583333 1.00000    
    ## WI.17.Hancock.SB-MI.16.East.Lansing.SB     11.87019231 1.00000    
    ## MI.15.Montcalm.SB-MI.16.East.Lansing.SB    18.36250000 1.00000    
    ## WA.07.fld09.DB-MI.16.East.Lansing.SB       23.62500000 1.00000    
    ## MX.2017.pre.DB-MI.16.East.Lansing.SB       26.36250000 1.00000    
    ## BR.2010.por.DB-MI.16.East.Lansing.SB       36.06250000 1.00000    
    ## ND.07.fld01.DB-MI.16.East.Lansing.SB       63.81250000 1.00000    
    ## BR.2010.unai.DB-MI.16.East.Lansing.SB      47.06250000 1.00000    
    ## MX.2017.ca.DB-MI.16.East.Lansing.SB        76.76250000 0.48390    
    ## MX.2017.ley.DB-MI.16.East.Lansing.SB       89.83522727 0.05784 .  
    ## WI.17.Hancock.SB-BR.2010.lam.DB            -0.02564103 1.00000    
    ## MI.15.Montcalm.SB-BR.2010.lam.DB            6.46666667 1.00000    
    ## WA.07.fld09.DB-BR.2010.lam.DB              11.72916667 1.00000    
    ## MX.2017.pre.DB-BR.2010.lam.DB              14.46666667 1.00000    
    ## BR.2010.por.DB-BR.2010.lam.DB              24.16666667 1.00000    
    ## ND.07.fld01.DB-BR.2010.lam.DB              51.91666667 1.00000    
    ## BR.2010.unai.DB-BR.2010.lam.DB             35.16666667 1.00000    
    ## MX.2017.ca.DB-BR.2010.lam.DB               64.86666667 1.00000    
    ## MX.2017.ley.DB-BR.2010.lam.DB              77.93939394 1.00000    
    ## MI.15.Montcalm.SB-WI.17.Hancock.SB          6.49230769 1.00000    
    ## WA.07.fld09.DB-WI.17.Hancock.SB            11.75480769 1.00000    
    ## MX.2017.pre.DB-WI.17.Hancock.SB            14.49230769 1.00000    
    ## BR.2010.por.DB-WI.17.Hancock.SB            24.19230769 1.00000    
    ## ND.07.fld01.DB-WI.17.Hancock.SB            51.94230769 1.00000    
    ## BR.2010.unai.DB-WI.17.Hancock.SB           35.19230769 1.00000    
    ## MX.2017.ca.DB-WI.17.Hancock.SB             64.89230769 1.00000    
    ## MX.2017.ley.DB-WI.17.Hancock.SB            77.96503497 0.48668    
    ## WA.07.fld09.DB-MI.15.Montcalm.SB            5.26250000 1.00000    
    ## MX.2017.pre.DB-MI.15.Montcalm.SB            8.00000000 1.00000    
    ## BR.2010.por.DB-MI.15.Montcalm.SB           17.70000000 1.00000    
    ## ND.07.fld01.DB-MI.15.Montcalm.SB           45.45000000 1.00000    
    ## BR.2010.unai.DB-MI.15.Montcalm.SB          28.70000000 1.00000    
    ## MX.2017.ca.DB-MI.15.Montcalm.SB            58.40000000 1.00000    
    ## MX.2017.ley.DB-MI.15.Montcalm.SB           71.47272727 1.00000    
    ## MX.2017.pre.DB-WA.07.fld09.DB               2.73750000 1.00000    
    ## BR.2010.por.DB-WA.07.fld09.DB              12.43750000 1.00000    
    ## ND.07.fld01.DB-WA.07.fld09.DB              40.18750000 1.00000    
    ## BR.2010.unai.DB-WA.07.fld09.DB             23.43750000 1.00000    
    ## MX.2017.ca.DB-WA.07.fld09.DB               53.13750000 1.00000    
    ## MX.2017.ley.DB-WA.07.fld09.DB              66.21022727 1.00000    
    ## BR.2010.por.DB-MX.2017.pre.DB               9.70000000 1.00000    
    ## ND.07.fld01.DB-MX.2017.pre.DB              37.45000000 1.00000    
    ## BR.2010.unai.DB-MX.2017.pre.DB             20.70000000 1.00000    
    ## MX.2017.ca.DB-MX.2017.pre.DB               50.40000000 1.00000    
    ## MX.2017.ley.DB-MX.2017.pre.DB              63.47272727 1.00000    
    ## ND.07.fld01.DB-BR.2010.por.DB              27.75000000 1.00000    
    ## BR.2010.unai.DB-BR.2010.por.DB             11.00000000 1.00000    
    ## MX.2017.ca.DB-BR.2010.por.DB               40.70000000 1.00000    
    ## MX.2017.ley.DB-BR.2010.por.DB              53.77272727 1.00000    
    ## BR.2010.unai.DB-ND.07.fld01.DB            -16.75000000 1.00000    
    ## MX.2017.ca.DB-ND.07.fld01.DB               12.95000000 1.00000    
    ## MX.2017.ley.DB-ND.07.fld01.DB              26.02272727 1.00000    
    ## MX.2017.ca.DB-BR.2010.unai.DB              29.70000000 1.00000    
    ## MX.2017.ley.DB-BR.2010.unai.DB             42.77272727 1.00000    
    ## MX.2017.ley.DB-MX.2017.ca.DB               13.07272727 1.00000    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
survey.picoxystrobin.fields.2  <- survey.picoxystrobin.fields %>% 
mutate(Field_difference = ifelse(
    Field_Year == "Baseline"| Field_Year == "ND.07.fld01.DB"|
      Field_Year == "BR.2010.unai.DB"|Field_Year == "MX.2017.ca.DB"|Field_Year == "MX.2017.ley.DB" ,
    "Yes","No"))%>% 
  mutate(Country_difference = ifelse(
    country == "USA",
    "Yes","No"))
####
myplot_comparison(survey.picoxystrobin.fields.2,highlabel = 0.005)  + scale_y_continuous(breaks = c(0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024)) + geom_hline(yintercept =
                                                                                                                                                       0.02148932, linetype = "dashed")
```

![](READM_files/figure-gfm/Stat%20Anal%20Picoxystrobin%20By%20Field-1.png)<!-- -->

# Table Publication

``` r
##
table1 <- 
survey.boscalid.complete %>% 
 bind_rows(wrangling_DC(survey.TM.2)) %>% 
  bind_rows(survey.tetraconazole.complete) %>% 
  bind_rows(survey.picoxystrobin.complete) %>% 
select(-c(EC50DC)) %>% 
  distinct(ID, .keep_all = TRUE) %>% 
group_by(source, Host, country, Field_Year, State ) %>%
  summarize(N= n()) %>%
    ungroup() %>%
        add_row(source = "Total", N  = sum(.$N))
```

    ## `summarise()` has grouped output by 'source', 'Host', 'country', 'Field_Year'.
    ## You can override using the `.groups` argument.

``` r
hola.2 <- survey.TM.2 %>% select(-c(EC50DC))  %>% 
  pull(ID) %>% unique()
## Corroborating
survey.boscalid.complete %>% 
 bind_rows(survey.TM.2) %>% 
  bind_rows(survey.tetraconazole.complete) %>% 
  bind_rows(survey.picoxystrobin.complete) %>% 
select(-c(EC50DC)) %>% 
  pull(ID) %>% 
  unique()
```

    ##   [1]   26  274  307  456  478  479  482  504  505  698  699  710  711  724  725
    ##  [16]  731  732  738  739  745  746  751  752  755  756  757  764  765  771  772
    ##  [31]  786  787  811  812  813  814  817  818  824  838  851  852  853  854  855
    ##  [46]  858  859  860  861  862  867  870  871  877  878  884  885  891  892  896
    ##  [61]  897  901  902  905  906  908  909  911  912  914 1058 1081 1087 1109 1127
    ##  [76] 1128 1134 1135 1139 1328 1329 1330 1331 1332 1340 1345 1365 1366 1392 1393
    ##  [91] 1502 1582 1620 1622 1671 1672 1691 1692 1712 1713 1721 1722 1731 1791 1832
    ## [106] 1844 1849 1870 1872 1885 1941 1942 1962 1963 1965 1972 2001 2002 2012 2013
    ## [121] 2014 2015 2027 2029 2030 2031 2033 2097 2108 2110 2111 2119 2120 2121 2122
    ## [136] 2129 2130 2131 2132 2149 2150 2160 2169 2170 2189 2190 2199 2383 2384 2388
    ## [151] 2407 2408 2449 2453 2518 2528 2536 2559 1082 1084 1092 1096 1098 1099 1108
    ## [166] 1110 1400 1403 1427 1434 1436 1587 1590 1593 1595 1612 1614 1618 1628 1725
    ## [181] 1727 1732 1747 1767 1771 1793 2293 2294 2295 2296 2297 2298 2299 2301 2302
    ## [196]  700  701  702  703  704  705  706  707  708  709  713  714  715  716  717
    ## [211]  718  719  720  721  726  727  728  729  730  733  734  735  736  737  740
    ## [226]  741  743  744  747  748  749  750  753  754  758  759  760  762  763  766
    ## [241]  767  768  769  770  788  789  790  791  792  793  794  795  796  797  810
    ## [256]  815  816  819  863  864  865  866  868  869  872  873  874  875  876  879
    ## [271]  880  881  882  883  886  887  888  889  890  893  894  895  898  899  900
    ## [286]  903  904  907  910 1884 1887 1888 1889 1890 1891 1892 1893 1894 1895 1896
    ## [301] 1898 1899 1901 1902 1903 1904 1905 1906 1907 1920 1921 1922 1923 1924 1926
    ## [316] 1928 1929 1930 1931 2457 2481 2483 2484 2485 2486 2487 2488 2489 2490 2504
    ## [331] 2505 2506 2507 2508 2509 2510 2511 2512 2513 2542 2545 2546 2548 2550 2551
    ## [346] 2552 2554 2557    1   12   20   21   74   87  118  123  449  461  467  475
    ## [361]  558  564  568  581  645  667 1025 1026 1027 1029 1032 1033 2098 2099 2100
    ## [376] 2139 2140 2143 2220 2222 2223 2320 2362 2385 2386 2390 1059 1064 1066 1070
    ## [391] 1071 1072 1113 1115 1122 1123 1124 1126 1382 1384 1385 1387 1388 1389 1390
    ## [406] 1391 1501 1503 1506 1520 1522 1530 1533 1537 1680 1739 1746 1751 1763 1769
    ## [421] 1775 1961 1967 1968 1969 2382  800 1842 1843 1845 1847 1848 1852 2035 2036
    ## [436] 2109 2112 2113 2114 2115 2117 2133 2134 2136 2151 2152 2153 2154 2155 2156
    ## [451] 2157 2158 2161 2162 2163 2164 2165 2167 2168 2191 2192 2193 2194 2195 2196
    ## [466] 2198 2547 1136 1140 1145 1146 1150 1163 1167 1349 1379 1401 1423 1451 1487
    ## [481] 1633 1634 1639 1643 1644 1652 1653 1656 1660 1720 1736 1738 1759 1807 2251
    ## [496] 2258 2261 2265 2266 2274 2277 2278 2281 2289 2409 2410 2411 2412 2413 2414
    ## [511] 2415 2416

## Less sensive isolates from boscalid and picoxystrobin

``` r
boscalid.less.sensive <-survey.boscalid.fields.2 %>% arrange(desc(EC50DC))

boscalid.less.sensive[1:11,] %>% arrange(ID)
```

    ## # A tibble: 11 × 10
    ##       ID EC50DC source    Host  country Field_Year State County Field_difference
    ##    <dbl>  <dbl> <fct>     <fct> <chr>   <fct>      <chr> <chr>  <chr>           
    ##  1   449  0.175 Baseline  Dryb… USA     Baseline   MI    no sp… Yes             
    ##  2   564  0.222 Baseline  Dryb… USA     Baseline   WA    no sp… Yes             
    ##  3   790  0.172 Producer… Dryb… USA     WA.08.fld… WA    no sp… No              
    ##  4  1436  0.163 Fungicid… Soyb… USA     IA.17.Nas… IA    no sp… No              
    ##  5  1612  0.161 Fungicid… Soyb… USA     WI.17.Han… WI    no sp… No              
    ##  6  1884  0.160 Producer… Dryb… Mexico  MX.2017.l… SIN   no sp… No              
    ##  7  1889  0.164 Producer… Dryb… Mexico  MX.2017.l… SIN   no sp… No              
    ##  8  1902  0.172 Producer… Dryb… Mexico  MX.2017.p… SIN   no sp… No              
    ##  9  1922  0.180 Producer… Dryb… Mexico  MX.2017.c… SIN   no sp… No              
    ## 10  2488  0.161 Producer… Dryb… Brazil  BR.2010.l… MG    Lamba… No              
    ## 11  2557  0.161 Producer… Dryb… Brazil  BR.2010.u… MG    Unai   No              
    ## # … with 1 more variable: Country_difference <chr>

``` r
boscalid.less.sensive.2 <-boscalid.less.sensive %>% select(ID) %>% pull()

boscalid.less.sensive.3 <- sort(boscalid.less.sensive.2[1:10])
boscalid.less.sensive.3
```

    ##  [1]  449  564  790 1436 1612 1889 1902 1922 2488 2557

``` r
picoxystrobin.less.sensive <-survey.picoxystrobin.fields.2 %>% arrange(desc(EC50DC))

picoxystrobin.less.sensive[1:11,] %>% arrange(ID)
```

    ## # A tibble: 11 × 10
    ##       ID EC50DC source    Host  country Field_Year State County Field_difference
    ##    <dbl>  <dbl> <fct>     <fct> <chr>   <fct>      <chr> <chr>  <chr>           
    ##  1   706 0.0211 Producer… Dryb… USA     ND.07.fld… ND    no sp… Yes             
    ##  2  1884 0.0272 Producer… Dryb… Mexico  MX.2017.l… SIN   no sp… Yes             
    ##  3  1887 0.0251 Producer… Dryb… Mexico  MX.2017.l… SIN   no sp… Yes             
    ##  4  1890 0.0240 Producer… Dryb… Mexico  MX.2017.l… SIN   no sp… Yes             
    ##  5  1895 0.0230 Producer… Dryb… Mexico  MX.2017.l… SIN   no sp… Yes             
    ##  6  1921 0.0231 Producer… Dryb… Mexico  MX.2017.c… SIN   no sp… Yes             
    ##  7  1922 0.0235 Producer… Dryb… Mexico  MX.2017.c… SIN   no sp… Yes             
    ##  8  1926 0.0215 Producer… Dryb… Mexico  MX.2017.c… SIN   no sp… Yes             
    ##  9  1928 0.0238 Producer… Dryb… Mexico  MX.2017.c… SIN   no sp… Yes             
    ## 10  2548 0.0244 Producer… Dryb… Brazil  BR.2010.u… MG    Unai   Yes             
    ## 11  2550 0.0253 Producer… Dryb… Brazil  BR.2010.u… MG    Unai   Yes             
    ## # … with 1 more variable: Country_difference <chr>

``` r
picoxystrobin.less.sensive.2 <-picoxystrobin.less.sensive %>% select(ID) %>% pull()
picoxystrobin.less.sensive.3 <- sort(picoxystrobin.less.sensive.2[1:10])
picoxystrobin.less.sensive.3
```

    ##  [1] 1884 1887 1890 1895 1921 1922 1926 1928 2548 2550

## Corroboration complete isolates

``` r
survey.picoxystrobin.complete %>% group_by( Host) %>% summarise(ave=mean(EC50DC))
```

    ## # A tibble: 2 × 2
    ##   Host       ave
    ##   <fct>    <dbl>
    ## 1 Drybean 0.0160
    ## 2 Soybean 0.0142

``` r
 survey.boscalid.fields %>% select(-c(EC50DC)) %>% 
 distinct(ID, .keep_all = TRUE) %>% 
group_by(Field_Year ) %>%
summarize(N= n()) 
```

    ## # A tibble: 17 × 2
    ##    Field_Year                N
    ##    <fct>                 <int>
    ##  1 NE.10.fld21.DB            9
    ##  2 ND.07.fld01.DB           12
    ##  3 WA.07.fld09.DB           16
    ##  4 MI.16.East.Lansing.SB    17
    ##  5 WA.08.fld10.DB           12
    ##  6 ND.10.fld25.DB            9
    ##  7 MI.15.Montcalm.SB         9
    ##  8 IA.17.Nashua.SB          16
    ##  9 Baseline                 21
    ## 10 WI.17.Hancock.SB         12
    ## 11 MI.17.Montcalm.SB        17
    ## 12 BR.2010.por.DB           10
    ## 13 BR.2010.unai.DB           9
    ## 14 MX.2017.pre.DB           10
    ## 15 MX.2017.ca.DB            10
    ## 16 BR.2010.lam.DB            9
    ## 17 MX.2017.ley.DB           11

# Supplementary Table: Fungicides applied in the field

``` r
#based on object table1 which is correct, just 
lasiodiplodia <- survey.boscalid.complete %>% 
 bind_rows(wrangling_DC(survey.TM.2)) %>% 
  bind_rows(survey.tetraconazole.complete) %>% 
  bind_rows(survey.picoxystrobin.complete) %>% 
select(-c(EC50DC)) %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  select(ID) %>% 
  pull()

 colnames(whitemold.inventory)
```

    ##  [1] "...1"                                               
    ##  [2] "ID"                                                 
    ##  [3] "Year"                                               
    ##  [4] "Host"                                               
    ##  [5] "State"                                              
    ##  [6] "County"                                             
    ##  [7] "Country"                                            
    ##  [8] "Field"                                              
    ##  [9] "Form.ID"                                            
    ## [10] "Serial_agar_dilution_for_creating_the_model_..N.20."
    ## [11] "Origin"                                             
    ## [12] "hyphal_tip..ht."                                    
    ## [13] "DNA.Extraction"                                     
    ## [14] "lat"                                                
    ## [15] "long"                                               
    ## [16] "Size.of.field..acres."                              
    ## [17] "Percentage.Affected"                                
    ## [18] "History"                                            
    ## [19] "Plot"                                               
    ## [20] "Number.of.Years"                                    
    ## [21] "Fungicide.This.Year"                                
    ## [22] "Chemigation"                                        
    ## [23] "Applied.for.White.Mold"                             
    ## [24] "Fungicide_current_season"                           
    ## [25] "molecule.s._current_season"                         
    ## [26] "company_current_season"                             
    ## [27] "Timing"                                             
    ## [28] "Applications"                                       
    ## [29] "Rate"                                               
    ## [30] "Pre.2015.Fungicides"                                
    ## [31] "Number.of.Years.1"                                  
    ## [32] "Fungicide_previous_seasons"                         
    ## [33] "molecule.s._previous_seasons"                       
    ## [34] "company_previous_seasons"                           
    ## [35] "Timing..Soybean"                                    
    ## [36] "Notes"                                              
    ## [37] "Apparently.low.sensitiviy"                          
    ## [38] "Group_MOA_current_season"                           
    ## [39] "Group_MOA_previous_seasons"                         
    ## [40] "Field_Names"                                        
    ## [41] "Field_Year"

``` r
 # I have to keep Field and Year separate tahts why reading again the main table from data
 
tableS1 <-  wrangling_DC(filename = whitemold.inventory) %>% mutate(source = as.factor(source),
                        Host = as.factor(Host),
                        country = as.factor(country),
                        Field = as.factor(Field),
                        Year = as.factor(Year),
                        State = as.factor(State),
                        County = as.factor(County),
                        Fungicide_current_season = as.factor(Fungicide_current_season),
                        molecule.s._current_season = as.factor(molecule.s._current_season),
                        Group_MOA_current_season = as.factor(Group_MOA_current_season),
                        Applications = as.factor(Applications),
                        Rate = as.factor(Rate))%>% select(
                          ID,
                          source,#we will change later to the real variable "source"
                          Host,
                          country, 
                          Field,
                          Year,
                          State,
                          County,
                          Fungicide_current_season,
                          molecule.s._current_season,
                          Group_MOA_current_season,
                          Applications,
                          Rate, Field_Year) %>% filter(ID %in% lasiodiplodia )
```

# Supplementary Figure Baseline Distribution

``` r
a.tet.baseline <-
  myplot_baseline(survey.tetraconazole.complete)  + expand_limits(y = c(0.2, 1.8)) + scale_y_continuous(breaks = c(0.2,0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6,1.8 )) +   theme(
      legend.position = "none"
    ) + labs(tag = "A")  +     theme(plot.tag = element_text(face = "bold")) 
  ##  
b.bos.baseline <- myplot_baseline(survey.boscalid.complete)+ expand_limits(y = c(0.04, 0.24)) +  scale_y_continuous(breaks = c(0.04,0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24))+   theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none"
    ) + labs(tag = "B")  +     theme(plot.tag = element_text(face = "bold"))
####
c.pico.baseline <- myplot_baseline(survey.picoxystrobin.complete)  + expand_limits(y = c(0.008, 0.024)) +  scale_y_continuous(breaks = c(0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.020, 0.022, 0.024 ))+   theme(
      axis.title.y = element_blank(),
      legend.position = "none"
    ) + labs(tag = "C")  +     theme(plot.tag = element_text(face = "bold"))
#Graph together
a.tet.baseline| (b.bos.baseline/c.pico.baseline)
```

![](READM_files/figure-gfm/Supplementary%20Figure%20Baseline-1.png)<!-- -->

# Figures 1 and 2 article

**A. Tetraconazole**  
**B. Boscalid**  
**C. Picoxystrobin**

``` r
a.tet <-
  myplot_model_1_tetraconazole (tetraconazole.joined.2.clean, lm = tetraconazole.lm)+   theme(
     legend.position = "none"
  ) + labs(tag = "A")  +     theme(  plot.tag = element_text(face = "bold")) 
#
b.bos <-
  myplot_model_1_boscalid (boscalid.joined.2.clean, lm = boscalid.lm)  +   theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) + labs(tag = "B")  +     theme(  plot.tag = element_text(face = "bold")) 
###
c.pico <-
  myplot_model_1_picoxystrobin (picoxystrobin.joined.2.clean, lm = picoxystrobin.lm)  +   theme(axis.title.y = element_blank(), legend.position = "none")+ labs(tag = "C")  +     theme(  plot.tag = element_text(face = "bold")) 

#Graph together
a.tet | (b.bos / c.pico)
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Figures%201%20and%202%20Article-1.png)<!-- -->

``` r
###
aa.tet <-
  myplot_model_2(tetraconazole.EC50DC) + expand_limits(x = c(0.27, 1.8), y = c(0.27, 1.8))  + scale_x_continuous(name = waiver(),
                                                                                                                 breaks = c(0.3, 0.6, 0.9, 1.2, 1.5, 1.8)) + scale_y_continuous(breaks = c(0.3, 0.6, 0.9, 1.2, 1.5, 1.8)) +   theme(
                                                                                                                                                                                                       legend.position = "none"
                                                                                                                 ) + geom_label(
                                                                                                                   aes(x = 0.6,
                                                                                                                       y = 1.5)
                                                                                                                   ,
                                                                                                                   label = c(paste  (
                                                                                                                     " Y =",
                                                                                                                     paste0 (round(summary(tetraconazole.lm.2)[[4]][2], 4), "x"),
                                                                                                                     round(summary(tetraconazole.lm.2)[[4]][1], 3) *
                                                                                                                       -1,
                                                                                                                     paste  ("\n R =",
                                                                                                                             paste0 (round(
                                                                                                                               summary(tetraconazole.lm.2)[[9]][1], 4
                                                                                                                             ), "\n p < 0.001"))
                                                                                                                   )),
                                                                                                                   size = 3,
                                                                                                                   fontface = "bold"
                                                                                                                 ) + labs(tag = "A")  +     theme(  plot.tag = element_text(face = "bold")) 

bb.bos <-
  myplot_model_2(boscalid.EC50DC) + expand_limits(x = c(0.04, 0.25), y = c(0.04, 0.25))  + scale_x_continuous(name = waiver(),
                                                                                                              breaks = c(0.05, 0.1, 0.15, 0.20, 0.25)) + scale_y_continuous(breaks = c(0.05, 0.1, 0.15, 0.20, 0.25)) +   theme(axis.title.x   = element_blank(), axis.title.y = element_blank(), legend.position = "none") + geom_label(
                                                                                                                aes(x = 0.08,
                                                                                                                    y = 0.2)
                                                                                                                ,
                                                                                                                label = c(paste  (
                                                                                                                  " Y =",
                                                                                                                  paste0 (round(summary(boscalid.lm.2)[[4]][2], 4), "x"),
                                                                                                                  
                                                                                                                  round(summary(boscalid.lm.2)[[4]][1], 3) *
                                                                                                                    -1,
                                                                                                                  paste  ("\n R =",
                                                                                                                          paste0 (round(
                                                                                                                            summary(boscalid.lm.2)[[9]][1], 4
                                                                                                                          ), "\n p < 0.001"))
                                                                                                                )),
                                                                                                                size = 3,
                                                                                                                fontface = "bold"
                                                                                                              ) + labs(tag = "B")  +     theme(  plot.tag = element_text(face = "bold")) 
###
 cc.pico <-
   myplot_model_2(picoxystrobin.EC50DC) + expand_limits(x = c(0.005, 0.018), y = c(0.005, 0.018))  + scale_x_continuous(name = waiver(),
                                                                                                                        breaks = c(0.0075, 0.01, 0.0125, 0.0150, 0.018)) + scale_y_continuous(breaks = c(0.0075, 0.01, 0.0125, 0.0150, 0.018))   + theme(
                                                                                                                          axis.title.y = element_blank(),
                                                                                                                          legend.position = "none",
                                                                                                                          axis.text.y = element_text(angle = 20,
                                                                                                                                                     hjust = 1), axis.text.x = element_text(angle = 20,
                                                                                                                                                     hjust = 1)
                                                                                                                        ) + geom_label(
                                                                                                                          aes(x = 0.0075,
                                                                                                                              y = 0.0148)
                                                                                                                          ,
                                                                                                                          label = c(paste  (
                                                                                                                            " Y =",
                                                                                                                            paste0 (round(summary(picoxystrobin.lm.2)[[4]][2], 4), "x"),
                                                                                                                            round(summary(picoxystrobin.lm.2)[[4]][1], 3) *
                                                                                                                              -1,
                                                                                                                            paste  ("\n R =",
                                                                                                                                    paste0 (round(
                                                                                                                                      summary(picoxystrobin.lm.2)[[9]][1], 4
                                                                                                                                    ), "\n p < 0.001"))
                                                                                                                          )),
                                                                                                                          size = 3,
                                                                                                                          fontface = "bold"
                                                                                                                        )  + labs(tag = "C")  +     theme(plot.tag = element_text(face = "bold"))
 ## together
aa.tet | (bb.bos / cc.pico)
```

    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'
    ## `geom_smooth()` using formula 'y ~ x'

![](READM_files/figure-gfm/Figures%201%20and%202%20Article-2.png)<!-- -->

### Figure 3 By Fields Article

**A. Tetraconazole**  
**B. Boscalid**  
**C. Picoxystrobin**

``` r
aaa.tet <-
  myplot_comparison(survey.tetraconazole.fields.2, highlabel = 0.1)  + expand_limits(y = c(0.4, 2.5)) + scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4)) + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.25), "cm")) +   
   labs(tag = "A")  +     theme(plot.tag = element_text(face = "bold"))

bbb.bos <-
  myplot_comparison(survey.boscalid.fields.2, highlabel = 0.017) + expand_limits(y = c(0.05, 0.24)) +
  scale_y_continuous(breaks = c(0.06, 0.1, 0.14,  0.18,  0.22)) + geom_hline(yintercept =  boscalid.less.sensive$EC50DC[9]
                                                                                                  , linetype = "dashed") +  theme(plot.margin = unit(c(0, 0, 0, 1), "cm")) +
  theme(
    axis.title.x   = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  ) +  #theme(plot.margin = unit(c(-1, 1, 1, 1), "cm"))
  labs(tag = "B")  +     theme(plot.tag = element_text(face = "bold")) # make dynamic the value of intercept, create a new fucntion

ccc.pico <-
  myplot_comparison(survey.picoxystrobin.fields.2, highlabel = 0.007)  + expand_limits(y = c(0.008, 0.032)) +
  scale_y_continuous(
    breaks = c(
      0.008,
      0.012,
      0.016,
      0.020,
      0.024,
      0.028,
      0.032
    )
  ) + geom_hline(yintercept =
                   picoxystrobin.less.sensive$EC50DC[10] , linetype = "dashed") +# +   theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))  
theme(axis.title.y = element_blank(), legend.position = "none") +   theme(plot.margin = unit(c(0, 0, -1, 0), "cm"))  + labs(tag = "C")  +     theme(plot.tag = element_text(face = "bold"))

#Together
aaa.tet| (bbb.bos/ccc.pico)
```

![](READM_files/figure-gfm/Figure%203%20Article-1.png)<!-- -->

Shield: [![CC BY-SA
4.0](https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by-sa/4.0/)

This work is licensed under a [Creative Commons Attribution-ShareAlike
4.0 International
License](http://creativecommons.org/licenses/by-sa/4.0/).

[![CC BY-SA
4.0](https://licensebuttons.net/l/by-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-sa/4.0/)
