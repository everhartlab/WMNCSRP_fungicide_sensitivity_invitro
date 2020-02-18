boscalid.data.cristian <-
  read.csv("data/boscalid_cristian_01.15.20.csv")
 glimpse(boscalid.data.cristian)
 
 
 
 boscalid.data.cristian.2 <- boscalid.data.cristian %>%
 rename(ecuatorial= Radial_growth1..mm., polar= Radial_growth2..mm.) %>% 
 mutate(polar= replace(polar, polar == 0, 6)) %>%  #replacing 0 cm growth for the size of plug that is 0.6 cm or 6 mm
   mutate(ecuatorial= replace(ecuatorial, ecuatorial == 0, 6)) %>% #replacing 0 cm growth for the size of plug that is 0.6 or 6 mm
      group_by(Isolate_ID, Treatment, Rep) %>%
   mutate(growth = ((
     ecuatorial + polar
   ) / 2) / 10) %>% # getting average and transforming in cm
   select(-c(ecuatorial, polar, Notes)) %>%
   ungroup() %>%
   mutate(experimental_replicate = 1) %>%
   rename(repeats = Rep, ID = Isolate_ID) %>%
   select(ID, repeats,  experimental_replicate, growth, Treatment) %>%
   group_by(Treatment) %>% 
   mutate(grouped_id = row_number()) %>%
   spread(Treatment, growth) %>%
   select(-grouped_id) %>% 
   ungroup() %>% 
   rename(control= Control, response= '0.2 ppm')
 
   ##
 isolates.cristian <- boscalid.data.cristian.2 %>% 
   pull(ID) %>% 
   unique() %>%
   as.numeric()
   
 isolates.cristian <- isolates.cristian[!isolates.cristian == 461]
 # I cannot use the function beacuse I need values of the second experimental replicate
 
   
 #  Using the model to estimate the EC50DC
 boscalid.cristian.filtered. <- boscalid.data.cristian.2 %>%
   group_by(ID) %>%
   summarise(
     mean_response = mean(response, na.rm = TRUE),
     mean_control = mean(control, na.rm = TRUE)
   ) %>%
   mutate(RG = (mean_response / mean_control) * 100) %>%
   
   mutate(Estimate.50DC = exp # exponential is the opposite of log
          #using the model the intercept and the coefficient of the model
          (finalRG0.2[[1]][[1]] + finalRG0.2[[1]][[2]] * RG)) %>% ungroup() %>%
   rename(EC50DC = Estimate.50DC) %>%
   mutate(ID = as.numeric(ID)) %>%
   select(c(ID, EC50DC)) %>%
   mutate(new = ifelse(ID %in% boscalid.filtered.survey.yes,
                       "Yes", "No")) #%>%
 #mutate(new= as.factor(new))
  
 #Selecting the EC50D from the serial dilution
 boscalid.complete.cristian <- final.boscalid.DC %>% 
   ####WATCH OUT, Based on EC50DC
   rename(EC50DC = Estimate.50DC) %>% 
   select(ID, EC50DC) %>% 
   bind_rows(boscalid.filtered.2 ) %>% 
   mutate(new= as.character(new)) %>% 
   bind_rows(boscalid.cristian.filtered. ) %>% 
   mutate(new= as.factor(new)) %>% 
   ungroup()
 
 ## Testing controls between Serial dilution and Duscriminatory concentration
 
 boscalid.controls.cristian <- boscalid.complete.cristian$ID[duplicated(boscalid.complete.cristian$ID)]
 # 24 repeats
 
 boscalid.complete.controls.cristian <- boscalid.complete.cristian %>% 
   mutate(controls = ifelse(
     ID %in% boscalid.controls.cristian,
     "Yes", "No")) %>% 
   filter(controls== "Yes")
 
 ## Normality test
 
 shapiro.test.boscal.controls.cristian <- boscalid.complete.controls.cristian  %>%
   do(tidy(shapiro.test(.$EC50DC)))
 
 shapiro.test.boscal.controls.cristian
 
 
 # All normals, t- test since the sample size for each sampe is 2, in other words grouping by ID
 
 tested.bocalid.controls.cristian <- boscalid.complete.controls.cristian %>% 
   group_by(ID) %>% 
   do(tidy(t.test(.$EC50DC))) %>% 
   filter(p.value<= 0.05)
 
 tested.bocalid.controls.cristian
 
 ##From the 21, there are  that are statistically differenet
 
 # keeping juts the highest value of the repeated ID value
 boscalid.complete.cristian<- boscalid.complete.cristian%>%  
   mutate( ID= as.numeric(ID)) %>% 
   group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() 
 #Filterring specific fields by field scale
 boscalid.complete.cris <- boscalid.complete.cristian%>% 
   mutate(source = ifelse(
     ID %in% baseline_isolates_2,
     "Baseline",
     ifelse(
       ID %in% Farmer_fields,
       "Producer Fields","Fungicide Field Trials"
     )), Host = ifelse(
       ID %in% drybean,
       "Drybean", ifelse(
         ID %in% soybean,
         "Soybean","different_host")),
     Field = ifelse(
       ID %in% MI_East_Lansing,
       "MI:East Lansing",
       ifelse(
         ID %in% IA_Nashua,
         "IA:Nashua",
         ifelse(
           ID %in% WI_Hancock,
           "WI:Hancock",
           ifelse(
             ID %in% MI_Montcalm_2017,
             "MI:Montcalm 2",
             ifelse(
               ID %in% MI_Montcalm_2015,
               "MI: Montcalm 1",
               ifelse(ID %in% potato, 
                      "potato", ifelse(ID %in% baseline_isolates_2, 
                                       "Baseline","other_fields"
                                       
                      )
               )
             )
           )
         )
       )
     )
   )%>%
   mutate(
     source = as.factor(source),
     Host = as.factor(Host),
     Field = as.factor(Field)
   ) %>%    
   #taking out "different host" becuse does not belong to drybean nor soybean
   filter(!ID == 8, !ID == 129)%>%
   mutate(Country = ifelse(
     ID %in% mexican.isolates,
     "Mexico",
     ifelse(ID %in% brazilian.isolates,
            "Brazil", "USA")
   )) 
 
 boscalid.complete.cris <-  boscalid.complete.cris %>% 
   mutate(Year_Field= case_when(
     ID %in%WM_1963~"1963",
     ID %in%WM_1965~"1965",
     ID %in%WM_1969~"1969",
     ID %in%WM_1973~"1973",
     ID %in%WM_1974~"1974",
     ID %in%WM_1975~"1975",
     ID %in%WM_1976~"1976",
     ID %in%WM_1977~"1977",
     ID %in%WM_1978~"1978",
     ID %in%WM_1979~"1979",
     ID %in%WM_1980~"1980",
     ID %in%WM_1981~"1981",
     ID %in%WM_1982~"1982",
     ID %in%WM_1983~"1983",
     ID %in%WM_1984~"1984",
     ID %in%WM_1985~"1985",
     ID %in%WM_1987~"1987",
     ID %in%WM_1989~"1989",
     ID %in%WM_1990~"1990",
     ID %in%WM_1991~"1991",
     ID %in%WM_1992~"1992",
     ID %in%WM_1993~"1993",
     ID %in%WM_1994~"1994",
     ID %in%WM_1995~"1995",
     ID %in%WM_1996~"1996",
     ID %in%WM_1997~"1997",
     ID %in%WM_1997_WM_production_31~"1997_production_31",
     ID %in%WM_1997_WM_production_32~"1997_production_32",
     ID %in%WM_1988~"1988",
     ID %in%WM_1998~"1998",
     ID %in%WM_1999~"1999",
     ID %in%WM_2000~"2000",
     ID %in%WM_2002~"2002",
     ID %in%WM_2003~"2003",
     ID %in%WM_1997_WM_production_33~"1997_production_33",
     ID %in%WM_2003_WM_nursery_1~"2003_nursery_1",
     ID %in%WM_2003_WM_nursery_2~"2003_nursery_2",
     ID %in%WM_2003_WM_nursery_3~"2003_nursery_3",
     ID %in%WM_1996_WM_nursery_4~"1996_nursery_4",
     ID %in%WM_2003_WM_nursery_5~"2003_nursery_5",
     ID %in%WM_2004~"2004",
     ID %in%WM_2005~"2005",
     ID %in%WM_2006~"2006",
     ID %in%WM_2004_WM_nursery_6~"2004_nursery_6",
     ID %in%WM_2004_WM_nursery_7~"2004_nursery_7",
     ID %in%WM_2004_WM_nursery_8~"2004_nursery_8",
     ID %in%WM_2004_WM_nursery_9~"2004_nursery_9",
     ID %in%WM_2004_WM_nursery_10~"2004_nursery_10",
     ID %in%WM_2004_WM_nursery_11~"2004_nursery_11",
     ID %in%WM_2004_WM_nursery_12~"2004_nursery_12",
     ID %in%WM_2005_WM_nursery_13~"2005_nursery_13",
     ID %in%WM_2005_WM_nursery_14~"2005_nursery_14",
     ID %in%WM_2005_WM_nursery_15~"2005_nursery_15",
     ID %in%WM_2005_WM_nursery_16~"2005_nursery_16",
     ID %in%WM_2005_WM_nursery_17~"2005_nursery_17",
     ID %in%WM_2005_WM_nursery_18~"2005_nursery_18",
     ID %in%WM_2006_WM_nursery_19~"2006_nursery_19",
     ID %in%WM_2007~"2007",
     ID %in%WM_2007_WM_production_1~"2007_production_1",
     ID %in%WM_2007_WM_production_2~"2007_production_2",
     ID %in%WM_2007_WM_production_3~"2007_production_3",
     ID %in%WM_2007_WM_production_4~"2007_production_4",
     ID %in%WM_2007_WM_production_5~"2007_production_5",
     ID %in%WM_2007_WM_production_6~"2007_production_6",
     ID %in%WM_2007_WM_production_7~"2007_production_7",
     ID %in%WM_2007_WM_production_8~"2007_production_8",
     ID %in%WM_2007_WM_production_9~"2007_production_9",
     ID %in%WM_2008~"2008",
     ID %in%WM_2008_WM_production_10~"2008_production_10",
     ID %in%WM_2008_WM_nursery_20~"2008_nursery_20",
     ID %in%WM_2008_WM_production_11~"2008_production_11",
     ID %in%WM_2008_WM_production_12~"2008_production_12",
     ID %in%WM_2008_WM_production_13~"2008_production_13",
     ID %in%WM_2009~"2009",
     ID %in%WM_2008_WM_nursery_21~"2008_nursery_21",
     ID %in%WM_2009_WM_nursery_22~"2009_nursery_22",
     ID %in%WM_2009_WM_nursery_23~"2009_nursery_23",
     ID %in%WM_2009_WM_production_14~"2009_production_14",
     ID %in%WM_2009_WM_production_15~"2009_production_15",
     ID %in%WM_2009_WM_production_16~"2009_production_16",
     ID %in%WM_2009_WM_production_17~"2009_production_17",
     ID %in%WM_2009_WM_production_18~"2009_production_18",
     ID %in%WM_2009_WM_production_19~"2009_production_19",
     ID %in%WM_2009_WM_production_20~"2009_production_20",
     ID %in%WM_2010_WM_production_21~"2010_production_21",
     ID %in%WM_2010_WM_production_22~"2010_production_22",
     ID %in%WM_2010_WM_production_23~"2010_production_23",
     ID %in%WM_2010_WM_production_24~"2010_production_24",
     ID %in%WM_2010_WM_production_25~"2010_production_25",
     ID %in%WM_2010_WM_production_26~"2010_production_26",
     ID %in%WM_2010_WM_production_27~"2010_production_27",
     ID %in%WM_2010_WM_production_28~"2010_production_28",
     ID %in%WM_2010_WM_production_29~"2010_production_29",
     ID %in%WM_2010_WM_production_30~"2010_production_30",
     ID %in%WM_2011~"2011",
     ID %in%WM_2013~"2013",
     ID %in%WM_2014~"2014",
     ID %in%WM_2016~"2016",
     ID %in%WM_2016_po~"2016_po",
     ID %in%WM_2016_flo~"2016_flo",
     ID %in%WM_2016_lan~"2016_lan",
     ID %in%WM_2017_ho~"2017_ho",
     ID %in%WM_2017_an~"2017_an",
     ID %in%WM_2017_dod~"2017_dod",
     ID %in%WM_2017_va~"2017_va",
     ID %in%WM_2017_na~"2017_na",
     ID %in%WM_2017_sa~"2017_sa",
     ID %in%WM_2017_ha~"2017_ha",
     ID %in%WM_2017_cu~"2017_cu",
     ID %in%WM_2017_mon~"2017_mon",
     ID %in%WM_2017_ant~"2017_ant",
     ID %in%WM_2017_holt~"2017_holt",
     ID %in%WM_2017_arr~"2017_arr",
     ID %in%WM_2017_ini~"2017_ini",
     ID %in%WM_2017_ley~"2017_ley",
     ID %in%WM_2017_pre~"2017_pre",
     ID %in%WM_2017_ru~"2017_ru",
     ID %in%WM_2017_ca~"2017_ca",
     ID %in%WM_2017_st~"2017_st",
     ID %in%WM_2015_6~"2015_6",
     ID %in%WM_2015_8~"2015_8",
     ID %in%WM_2015_9~"2015_9",
     ID %in%WM_2015_10~"2015_10",
     ID %in%WM_2015_11~"2015_11",
     ID %in%WM_2015_18~"2015_18",
     ID %in%WM_2015_19~"2015_19",
     ID %in%WM_2015_20~"2015_20",
     ID %in%WM_2015_21~"2015_21",
     ID %in%WM_2015_22~"2015_22",
     ID %in%WM_2015_23~"2015_23",
     ID %in%WM_2015_24~"2015_24",
     ID %in%WM_2015_25~"2015_25",
     ID %in%WM_2015_34~"2015_34",
     ID %in%WM_2015_35~"2015_35",
     ID %in%WM_2015_45~"2015_45",
     ID %in%WM_2015_46~"2015_46",
     ID %in%WM_2015_47~"2015_47",
     ID %in%WM_2015_48~"2015_48",
     ID %in%WM_2015_60~"2015_60",
     ID %in%WM_2015_62~"2015_62",
     ID %in%WM_2015_66~"2015_66",
     ID %in%WM_2015_67~"2015_67",
     ID %in%WM_2015_68~"2015_68",
     ID %in%WM_2015_69~"2015_69",
     ID %in%WM_2015_75~"2015_75",
     ID %in%WM_2015_78~"2015_78",
     ID %in%WM_2015_97~"2015_97",
     ID %in%WM_2015_98~"2015_98",
     ID %in%WM_2015_99~"2015_99",
     ID %in%WM_2015_100~"2015_100",
     ID %in%WM_2015_101~"2015_101",
     ID %in%WM_2015_102~"2015_102",
     ID %in%WM_2015_103~"2015_103",
     ID %in%WM_2015_104~"2015_104",
     ID %in%WM_2015_Al~"2015_Al",
     ID %in%WM_2015_Hi~"2015_Hi",
     ID %in%WM_2015_In~"2015_In",
     ID %in%WM_2015_mont~"2015_mont",
     ID %in%WM_2015_Sand~"2015_Sand",
     ID %in%WM_2015_Sani~"2015_Sani",
     ID %in%WM_2015_St~"2015_St",
     ID %in%WM_2015_wau~"2015_wau",
     ID %in%WM_2015_la~"2015_la",
     ID %in%WM_2015_co~"2015_co",
     ID %in%WM_2015_cern~"2015_cern",
     ID %in%WM_2015_cer~"2015_cer",
     ID %in%WM_2010~"2010",
     ID %in%WM_2017~"2017",
     TRUE ~ "no specific field"
   ),
   State= case_when(
     ID %in%WM_AB~"AB",
     ID %in%WM_AZ~"AZ",
     ID %in%WM_BA~"BA",
     ID %in%WM_CA~"CA",
     ID %in%WM_CO~"CO",
     ID %in%WM_DE~"DE",
     ID %in%WM_DF~"DF",
     ID %in%WM_FL~"FL",
     ID %in%WM_GA~"GA",
     ID %in%WM_GO~"GO",
     
     ID %in%WM_IA~"IA",
     ID %in%WM_ID~"ID",
     ID %in%WM_KS~"KS",
     ID %in%WM_KY~"KY",
     ID %in%WM_LA~"LA",
     ID %in%WM_MD~"MD",
     ID %in%WM_MG~"MG",
     ID %in%WM_MI~"MI",
     ID %in%WM_MN~"MN",
     ID %in%WM_MO~"MO",
     
     ID %in%WM_MS~"MS",
     ID %in%WM_MT~"MT",
     ID %in%WM_NC~"NC",
     ID %in%WM_ND~"ND",
     ID %in%WM_NE~"NE",
     ID %in%WM_NY~"NY",
     ID %in%WM_OH~"OH",
     ID %in%WM_OK~"OK",
     ID %in%WM_ON~"ON",
     ID %in%WM_OR~"OR",
     
     ID %in%WM_PA~"PA",
     ID %in%WM_PR~"PR",
     ID %in%WM_QUE~"QUE",
     ID %in%WM_RS~"RS",
     ID %in%WM_SC~"SC",
     ID %in%WM_SIN~"SIN",
     ID %in%WM_SK~"SK",
     ID %in%WM_WI~"WI",
     ID %in%WM_WA~"WA"
     
   )
   ) %>% select(-c(new)) %>% 
   mutate(Country = as.factor(Country),
          Year_Field = as.factor(Year_Field),
          State= as.factor(State) ) %>% 
   arrange(ID)
 
 
 
 
 ###For comparisons between hosts I had to take out the ones coming from baseline otherwise would not be fair
 
 boscalid.complete.host <- boscalid.complete %>% 
   filter(!source=="Baseline")
 
 ###For isolates from that have reduced sensitivity, assess sensitivity for isolates in the same field to account for heterogeneicity
 
 
 
 boscalid.complete.fields.cristian <- boscalid.complete.cris%>% 
   filter(  Field== "Baseline"|
            State== "CO"|
            State== "NE"|
              State== "ND"| 
              State== "WA"|
              State== "MI"
              ) %>% 
   arrange(source) %>% 
   filter( !source == "Fungicide Field Trials") %>% 
   
    mutate(State= as.character(State), 
      State =  ifelse(Field == "Baseline", "Baseline", State),
      State = as.factor(State)) %>% 
   filter(Year_Field == "2007_production_1"|
            Year_Field == "2007_production_2"|
            Year_Field == "2007_production_3"|
            Year_Field == "2007_production_4"|
            Year_Field == "2007_production_5"|
            Year_Field == "2007_production_6"|
            Year_Field == "2007_production_7"|
            Year_Field == "2007_production_8"|
            Year_Field == "2007_production_9"|
            Year_Field == "2008_production_10"|
            Year_Field == "2008_production_11"|
            Year_Field == "2008_production_12"|
            Year_Field == "2008_production_13"|
            Year_Field == "2009_production_17"|
            Year_Field == "2009_production_19"|
            Year_Field == "2009_production_20"|
            Year_Field == "2010_production_21"|
            Year_Field == "2007_production_22"|
            Year_Field == "2007_production_23"|
            Year_Field == "2007_production_24"|
            Year_Field == "2007_production_25"|
            Year_Field == "2007_production_26"|
            Year_Field == "2007_production_27"|
            Year_Field == "2007_production_28"|
            Year_Field == "2007_production_29"|
            Year_Field == "2007_production_30"|
             State == "Baseline")
 
 ###
 # boscalid.complete.2 <- boscalid.complete%>% 
 #   filter( !Field == "potato", !Field == "other_fields")
 # 
 # boscalid.complete.3 <- boscalid.complete%>%
 #   filter(
 #     !Country =="Brazil" & !Country == "Mexico"  
 #   )
 
 
 plot.bos.3.cristian <- boscalid.complete.fields.cristian %>%
   ggplot(aes(x = State, y = EC50DC)) +
   geom_jitter(
     width = .1,
     height = 0,
     shape = 21,
     color = "black",
     fill = "chartreuse2",
     size = 2,
     alpha = 3 / 4
   ) + stat_summary(
     fun.y = mean,
     geom = "point",
     shape = 95,
     size = 15,
     color = "black"
   ) +
   labs(x = "State", y = expression(bold(EC[bold("50") ~ (bold("D"))]) ~
                                      (bold(ppm ~ bold(
                                        "a.i."
                                      ))))) + theme(
                                        panel.border = element_rect(
                                          colour = "black",
                                          fill = NA,
                                          size = 1
                                        ),
                                        axis.title = element_text(size = 10, face = "bold", hjust = 0.5),
                                        axis.text = element_text(
                                          face = "bold",
                                          size = 10,
                                          family = "Arial"
                                        ),
                                        panel.background = element_rect(fill = "white", colour = "grey50")
                                      ) +  expand_limits(y = c(0.06, 0.18)) + scale_y_continuous(breaks = c(0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18))
 
 plot.bos.3.cristian
 
 
 boscalid.complete.fields.cristian$State <-
   factor(boscalid.complete.fields.cristian$State,
          c("Baseline", "WA", "ND", "CO", "MI", "NE"))
 
 
 #by field but because the name actually is "State"
 
 ##NORMALITY TEST-Shapiro_test
 
 shapiro.test.boscal.cristian <- boscalid.complete.fields.cristian %>%
   do(tidy(shapiro.test(.$EC50DC)))
 shapiro.test.boscal.cristian[[2]][[1]]
 # No normals
 object.1.cristian <- kruskal.test(EC50DC ~ State, data = boscalid.complete.fields.cristian)
 object.1.cristian[[3]][[1]]
 
 #There is difference p-value = 0.003064295
 
 object.2.cristian <-  DunnTest(EC50DC ~ State, data = boscalid.complete.fields.cristian, method = "bonferroni")
 
 # mean.rank.diff   pval    
 # WA-Baseline     -3.2077922 1.0000    
 # ND-Baseline    -14.7142857 1.0000    
 # CO-Baseline    -15.1269841 1.0000    
 # MI-Baseline    -16.7252747 0.3825    
 # NE-Baseline    -29.4047619 0.0019 ** 
 #   ND-WA          -11.5064935 1.0000    
 # CO-WA          -11.9191919 1.0000    
 # MI-WA          -13.5174825 1.0000    
 # NE-WA          -26.1969697 0.0465 *  
 #   CO-ND           -0.4126984 1.0000    
 # MI-ND           -2.0109890 1.0000    
 # NE-ND          -14.6904762 1.0000    
 # MI-CO           -1.5982906 1.0000    
 # NE-CO          -14.2777778 1.0000    
 # NE-MI          -12.6794872 1.0000
 # 
 

 
 ### picoxystrobin
 
 picoxystrobin.data.cristian <-
   read.csv("data/picoxystrobin_cristian_01.15.20.csv")
 glimpse(picoxystrobin.data.cristian)
 
 
 
 picoxystrobin.data.cristian.2 <-picoxystrobin.data.cristian %>%
   rename(ecuatorial= Radial_growth1..mm., polar= Radial_growth2..mm.) %>% 
   mutate(polar= replace(polar, polar == 0, 0.6)) %>%  #replacing 0 cm growth for the size of plug that is 0.6 cm or 6 mm
   mutate(ecuatorial= replace(ecuatorial, ecuatorial == 0, 0.6)) %>% #replacing 0 cm growth for the size of plug that is 0.6 or 6 mm
   group_by(Isolate_ID, Treatment, Rep) %>%
   mutate(growth = ((
     ecuatorial + polar
   ) / 2) / 10) %>% # getting average and transforming in cm
   select(-c(ecuatorial, polar, Notes)) %>%
   ungroup() %>%
   mutate(experimental_replicate = 1) %>%
   rename(repeats = Rep, ID = Isolate_ID) %>%
   select(ID, repeats,  experimental_replicate, growth, Treatment) %>%
   group_by(Treatment) %>% 
   mutate(grouped_id = row_number()) %>%
   spread(Treatment, growth) %>%
   select(-grouped_id) %>% 
   ungroup() %>%
   rename(control= Control, response= '0.01 ppm')
 

 #  Using the model to estimate the EC50DC
 picoxystrobin.cristian.filtered. <-picoxystrobin.data.cristian.2 %>%
   group_by(ID) %>%
   summarise(
     mean_response = mean(response, na.rm = TRUE),
     mean_control = mean(control, na.rm = TRUE)
   ) %>%
   mutate(RG = (mean_response / mean_control) * 100) %>%
   
   mutate(Estimate.50DC = exp # exponential is the opposite of log
          #using the model the intercept and the coefficient of the model
          (finalRG0.01[[1]][[1]] + finalRG0.01[[1]][[2]] * RG)) %>%
   ungroup() %>%
   rename(EC50DC= Estimate.50DC) %>% 
   mutate(ID = as.numeric(ID)) %>%
   select(c(ID, EC50DC)) %>% 
   mutate(new = ifelse(
     ID %in% picoxystrobin.filtered.survey.yes,
     "Yes","No")) 
 
 #Selecting the EC50D from the serial dilution
 picoxystrobin.complete.cristian <- final.picoxystrobin.DC %>% 
   ####WATCH OUT, Based on EC50DC
   rename(EC50DC = Estimate.50DC) %>% 
   select(ID, EC50DC) %>% 
   bind_rows(picoxystrobin.filtered. ) %>% 
   #mutate(new= as.character(new)) %>% 
   bind_rows(picoxystrobin.cristian.filtered. ) %>% 
   mutate(new= as.factor(new)) %>% 
   ungroup()
 
 # keeping juts the highest value of the repeated ID value
 picoxystrobin.complete.cristian<-picoxystrobin.complete.cristian%>%  
   mutate( ID= as.numeric(ID)) %>% 
   group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() 
 #Filterring specific fields by field scale
 picoxystrobin.complete.cris <-picoxystrobin.complete.cristian%>% 
   mutate(source = ifelse(
     ID %in% baseline_isolates_2,
     "Baseline",
     ifelse(
       ID %in% Farmer_fields,
       "Producer Fields","Fungicide Field Trials"
     )), Host = ifelse(
       ID %in% drybean,
       "Drybean", ifelse(
         ID %in% soybean,
         "Soybean","different_host")),
     Field = ifelse(
       ID %in% MI_East_Lansing,
       "MI:East Lansing",
       ifelse(
         ID %in% IA_Nashua,
         "IA:Nashua",
         ifelse(
           ID %in% WI_Hancock,
           "WI:Hancock",
           ifelse(
             ID %in% MI_Montcalm_2017,
             "MI:Montcalm 2",
             ifelse(
               ID %in% MI_Montcalm_2015,
               "MI: Montcalm 1",
               ifelse(
                 ID %in% WI_Cerny,
                 "WI:Cerny",
               ifelse(ID %in% potato, 
                      "potato", ifelse(ID %in% baseline_isolates_2, 
                                       "Baseline","other_fields"
                      )                 
                      )
               )
             )
           )
         )
       )
     )
   )%>%
   mutate(
     source = as.factor(source),
     Host = as.factor(Host),
     Field = as.factor(Field)
   ) %>%    
   #taking out "different host" becuse does not belong to drybean nor soybean
   filter(!ID == 8, !ID == 129)%>%
   mutate(Country = ifelse(
     ID %in% mexican.isolates,
     "Mexico",
     ifelse(ID %in% brazilian.isolates,
            "Brazil", "USA")
   )) 
 
 picoxystrobin.complete.cris <- picoxystrobin.complete.cris %>% 
   mutate(Year_Field= case_when(
     ID %in%WM_1963~"1963",
     ID %in%WM_1965~"1965",
     ID %in%WM_1969~"1969",
     ID %in%WM_1973~"1973",
     ID %in%WM_1974~"1974",
     ID %in%WM_1975~"1975",
     ID %in%WM_1976~"1976",
     ID %in%WM_1977~"1977",
     ID %in%WM_1978~"1978",
     ID %in%WM_1979~"1979",
     ID %in%WM_1980~"1980",
     ID %in%WM_1981~"1981",
     ID %in%WM_1982~"1982",
     ID %in%WM_1983~"1983",
     ID %in%WM_1984~"1984",
     ID %in%WM_1985~"1985",
     ID %in%WM_1987~"1987",
     ID %in%WM_1989~"1989",
     ID %in%WM_1990~"1990",
     ID %in%WM_1991~"1991",
     ID %in%WM_1992~"1992",
     ID %in%WM_1993~"1993",
     ID %in%WM_1994~"1994",
     ID %in%WM_1995~"1995",
     ID %in%WM_1996~"1996",
     ID %in%WM_1997~"1997",
     ID %in%WM_1997_WM_production_31~"1997_production_31",
     ID %in%WM_1997_WM_production_32~"1997_production_32",
     ID %in%WM_1988~"1988",
     ID %in%WM_1998~"1998",
     ID %in%WM_1999~"1999",
     ID %in%WM_2000~"2000",
     ID %in%WM_2002~"2002",
     ID %in%WM_2003~"2003",
     ID %in%WM_1997_WM_production_33~"1997_production_33",
     ID %in%WM_2003_WM_nursery_1~"2003_nursery_1",
     ID %in%WM_2003_WM_nursery_2~"2003_nursery_2",
     ID %in%WM_2003_WM_nursery_3~"2003_nursery_3",
     ID %in%WM_1996_WM_nursery_4~"1996_nursery_4",
     ID %in%WM_2003_WM_nursery_5~"2003_nursery_5",
     ID %in%WM_2004~"2004",
     ID %in%WM_2005~"2005",
     ID %in%WM_2006~"2006",
     ID %in%WM_2004_WM_nursery_6~"2004_nursery_6",
     ID %in%WM_2004_WM_nursery_7~"2004_nursery_7",
     ID %in%WM_2004_WM_nursery_8~"2004_nursery_8",
     ID %in%WM_2004_WM_nursery_9~"2004_nursery_9",
     ID %in%WM_2004_WM_nursery_10~"2004_nursery_10",
     ID %in%WM_2004_WM_nursery_11~"2004_nursery_11",
     ID %in%WM_2004_WM_nursery_12~"2004_nursery_12",
     ID %in%WM_2005_WM_nursery_13~"2005_nursery_13",
     ID %in%WM_2005_WM_nursery_14~"2005_nursery_14",
     ID %in%WM_2005_WM_nursery_15~"2005_nursery_15",
     ID %in%WM_2005_WM_nursery_16~"2005_nursery_16",
     ID %in%WM_2005_WM_nursery_17~"2005_nursery_17",
     ID %in%WM_2005_WM_nursery_18~"2005_nursery_18",
     ID %in%WM_2006_WM_nursery_19~"2006_nursery_19",
     ID %in%WM_2007~"2007",
     ID %in%WM_2007_WM_production_1~"2007_production_1",
     ID %in%WM_2007_WM_production_2~"2007_production_2",
     ID %in%WM_2007_WM_production_3~"2007_production_3",
     ID %in%WM_2007_WM_production_4~"2007_production_4",
     ID %in%WM_2007_WM_production_5~"2007_production_5",
     ID %in%WM_2007_WM_production_6~"2007_production_6",
     ID %in%WM_2007_WM_production_7~"2007_production_7",
     ID %in%WM_2007_WM_production_8~"2007_production_8",
     ID %in%WM_2007_WM_production_9~"2007_production_9",
     ID %in%WM_2008~"2008",
     ID %in%WM_2008_WM_production_10~"2008_production_10",
     ID %in%WM_2008_WM_nursery_20~"2008_nursery_20",
     ID %in%WM_2008_WM_production_11~"2008_production_11",
     ID %in%WM_2008_WM_production_12~"2008_production_12",
     ID %in%WM_2008_WM_production_13~"2008_production_13",
     ID %in%WM_2009~"2009",
     ID %in%WM_2008_WM_nursery_21~"2008_nursery_21",
     ID %in%WM_2009_WM_nursery_22~"2009_nursery_22",
     ID %in%WM_2009_WM_nursery_23~"2009_nursery_23",
     ID %in%WM_2009_WM_production_14~"2009_production_14",
     ID %in%WM_2009_WM_production_15~"2009_production_15",
     ID %in%WM_2009_WM_production_16~"2009_production_16",
     ID %in%WM_2009_WM_production_17~"2009_production_17",
     ID %in%WM_2009_WM_production_18~"2009_production_18",
     ID %in%WM_2009_WM_production_19~"2009_production_19",
     ID %in%WM_2009_WM_production_20~"2009_production_20",
     ID %in%WM_2010_WM_production_21~"2010_production_21",
     ID %in%WM_2010_WM_production_22~"2010_production_22",
     ID %in%WM_2010_WM_production_23~"2010_production_23",
     ID %in%WM_2010_WM_production_24~"2010_production_24",
     ID %in%WM_2010_WM_production_25~"2010_production_25",
     ID %in%WM_2010_WM_production_26~"2010_production_26",
     ID %in%WM_2010_WM_production_27~"2010_production_27",
     ID %in%WM_2010_WM_production_28~"2010_production_28",
     ID %in%WM_2010_WM_production_29~"2010_production_29",
     ID %in%WM_2010_WM_production_30~"2010_production_30",
     ID %in%WM_2011~"2011",
     ID %in%WM_2013~"2013",
     ID %in%WM_2014~"2014",
     ID %in%WM_2016~"2016",
     ID %in%WM_2016_po~"2016_po",
     ID %in%WM_2016_flo~"2016_flo",
     ID %in%WM_2016_lan~"2016_lan",
     ID %in%WM_2017_ho~"2017_ho",
     ID %in%WM_2017_an~"2017_an",
     ID %in%WM_2017_dod~"2017_dod",
     ID %in%WM_2017_va~"2017_va",
     ID %in%WM_2017_na~"2017_na",
     ID %in%WM_2017_sa~"2017_sa",
     ID %in%WM_2017_ha~"2017_ha",
     ID %in%WM_2017_cu~"2017_cu",
     ID %in%WM_2017_mon~"2017_mon",
     ID %in%WM_2017_ant~"2017_ant",
     ID %in%WM_2017_holt~"2017_holt",
     ID %in%WM_2017_arr~"2017_arr",
     ID %in%WM_2017_ini~"2017_ini",
     ID %in%WM_2017_ley~"2017_ley",
     ID %in%WM_2017_pre~"2017_pre",
     ID %in%WM_2017_ru~"2017_ru",
     ID %in%WM_2017_ca~"2017_ca",
     ID %in%WM_2017_st~"2017_st",
     ID %in%WM_2015_6~"2015_6",
     ID %in%WM_2015_8~"2015_8",
     ID %in%WM_2015_9~"2015_9",
     ID %in%WM_2015_10~"2015_10",
     ID %in%WM_2015_11~"2015_11",
     ID %in%WM_2015_18~"2015_18",
     ID %in%WM_2015_19~"2015_19",
     ID %in%WM_2015_20~"2015_20",
     ID %in%WM_2015_21~"2015_21",
     ID %in%WM_2015_22~"2015_22",
     ID %in%WM_2015_23~"2015_23",
     ID %in%WM_2015_24~"2015_24",
     ID %in%WM_2015_25~"2015_25",
     ID %in%WM_2015_34~"2015_34",
     ID %in%WM_2015_35~"2015_35",
     ID %in%WM_2015_45~"2015_45",
     ID %in%WM_2015_46~"2015_46",
     ID %in%WM_2015_47~"2015_47",
     ID %in%WM_2015_48~"2015_48",
     ID %in%WM_2015_60~"2015_60",
     ID %in%WM_2015_62~"2015_62",
     ID %in%WM_2015_66~"2015_66",
     ID %in%WM_2015_67~"2015_67",
     ID %in%WM_2015_68~"2015_68",
     ID %in%WM_2015_69~"2015_69",
     ID %in%WM_2015_75~"2015_75",
     ID %in%WM_2015_78~"2015_78",
     ID %in%WM_2015_97~"2015_97",
     ID %in%WM_2015_98~"2015_98",
     ID %in%WM_2015_99~"2015_99",
     ID %in%WM_2015_100~"2015_100",
     ID %in%WM_2015_101~"2015_101",
     ID %in%WM_2015_102~"2015_102",
     ID %in%WM_2015_103~"2015_103",
     ID %in%WM_2015_104~"2015_104",
     ID %in%WM_2015_Al~"2015_Al",
     ID %in%WM_2015_Hi~"2015_Hi",
     ID %in%WM_2015_In~"2015_In",
     ID %in%WM_2015_mont~"2015_mont",
     ID %in%WM_2015_Sand~"2015_Sand",
     ID %in%WM_2015_Sani~"2015_Sani",
     ID %in%WM_2015_St~"2015_St",
     ID %in%WM_2015_wau~"2015_wau",
     ID %in%WM_2015_la~"2015_la",
     ID %in%WM_2015_co~"2015_co",
     ID %in%WM_2015_cern~"2015_cern",
     ID %in%WM_2015_cer~"2015_cer",
     ID %in%WM_2010~"2010",
     ID %in%WM_2017~"2017",
     TRUE ~ "no specific field"
   ),
   State= case_when(
     ID %in%WM_AB~"AB",
     ID %in%WM_AZ~"AZ",
     ID %in%WM_BA~"BA",
     ID %in%WM_CA~"CA",
     ID %in%WM_CO~"CO",
     ID %in%WM_DE~"DE",
     ID %in%WM_DF~"DF",
     ID %in%WM_FL~"FL",
     ID %in%WM_GA~"GA",
     ID %in%WM_GO~"GO",
     
     ID %in%WM_IA~"IA",
     ID %in%WM_ID~"ID",
     ID %in%WM_KS~"KS",
     ID %in%WM_KY~"KY",
     ID %in%WM_LA~"LA",
     ID %in%WM_MD~"MD",
     ID %in%WM_MG~"MG",
     ID %in%WM_MI~"MI",
     ID %in%WM_MN~"MN",
     ID %in%WM_MO~"MO",
     
     ID %in%WM_MS~"MS",
     ID %in%WM_MT~"MT",
     ID %in%WM_NC~"NC",
     ID %in%WM_ND~"ND",
     ID %in%WM_NE~"NE",
     ID %in%WM_NY~"NY",
     ID %in%WM_OH~"OH",
     ID %in%WM_OK~"OK",
     ID %in%WM_ON~"ON",
     ID %in%WM_OR~"OR",
     
     ID %in%WM_PA~"PA",
     ID %in%WM_PR~"PR",
     ID %in%WM_QUE~"QUE",
     ID %in%WM_RS~"RS",
     ID %in%WM_SC~"SC",
     ID %in%WM_SIN~"SIN",
     ID %in%WM_SK~"SK",
     ID %in%WM_WI~"WI",
     ID %in%WM_WA~"WA"
     
   )
   ) %>% select(-c(new)) %>% 
   mutate(Country = as.factor(Country),
          Year_Field = as.factor(Year_Field),
          State= as.factor(State) ) %>% 
   arrange(ID)
 
 ##
 
 ###For comparisons between hosts I had to take out the ones coming from baseline otherwise would not be fair
 
 picoxystrobin.complete.host <-picoxystrobin.complete %>% 
   filter(!source=="Baseline")
 
 ###For isolates from that have reduced sensitivity, assess sensitivity for isolates in the same field to account for heterogeneicity
 ##WITH MEXICO & BRAZIL together
 # picoxystrobin.complete.fields.cristian <-picoxystrobin.complete.cris%>%
 #   filter(  Field== "Baseline"|
 #              State== "CO"|
 #              State== "NE"|
 #              State== "ND"|
 #              State== "WA"|
 #              State== "MI"|
 #              Country== "Mexico"|
 #              Country== "Brazil"
 #   ) %>%
 #   arrange(source) %>%
 #    #
 #   mutate(State= as.character(State),
 #          State =  ifelse(Field == "Baseline", "Baseline", State),
 #          State = as.factor(State)) %>%
 #   filter(Year_Field == "2007_production_1"|
 #            Year_Field == "2007_production_2"|
 #            Year_Field == "2007_production_3"|
 #            Year_Field == "2007_production_4"|
 #            Year_Field == "2007_production_5"|
 #            Year_Field == "2007_production_6"|
 #            Year_Field == "2007_production_7"|
 #            Year_Field == "2007_production_8"|
 #            Year_Field == "2007_production_9"|
 #            Year_Field == "2008_production_10"|
 #            Year_Field == "2008_production_11"|
 #            Year_Field == "2008_production_12"|
 #            Year_Field == "2008_production_13"|
 #            Year_Field == "2009_production_17"|
 #            Year_Field == "2009_production_19"|
 #            Year_Field == "2009_production_20"|
 #            Year_Field == "2010_production_21"|
 #            Year_Field == "2007_production_22"|
 #            Year_Field == "2007_production_23"|
 #            Year_Field == "2007_production_24"|
 #            Year_Field == "2007_production_25"|
 #            Year_Field == "2007_production_26"|
 #            Year_Field == "2007_production_27"|
 #            Year_Field == "2007_production_28"|
 #            Year_Field == "2007_production_29"|
 #            Year_Field == "2007_production_30"|
 #            State == "Baseline"|
 #            Country== "Mexico"|
 #            Country== "Brazil") %>% 
 #           mutate(State = fct_recode(State, "Mexico & Brazil"= "SIN" , "Mexico & Brazil"= "GO", "Mexico & Brazil" = "MG"))
 #          
 # 
 picoxystrobin.complete.fields.cristian <-picoxystrobin.complete.cris%>% 
   filter(  Field== "Baseline"|
              State== "CO"|
              State== "NE"|
              State== "ND"| 
              State== "WA"|
              State== "MI"
              
   ) %>% 
   arrange(source) %>% 
   filter( !source == "Fungicide Field Trials") %>% 
   
   mutate(State= as.character(State), 
          State =  ifelse(Field == "Baseline", "Baseline", State),
          State = as.factor(State)) %>% 
   filter(Year_Field == "2007_production_1"|
            Year_Field == "2007_production_2"|
            Year_Field == "2007_production_3"|
            Year_Field == "2007_production_4"|
            Year_Field == "2007_production_5"|
            Year_Field == "2007_production_6"|
            Year_Field == "2007_production_7"|
            Year_Field == "2007_production_8"|
            Year_Field == "2007_production_9"|
            Year_Field == "2008_production_10"|
            Year_Field == "2008_production_11"|
            Year_Field == "2008_production_12"|
            Year_Field == "2008_production_13"|
            Year_Field == "2009_production_17"|
            Year_Field == "2009_production_19"|
            Year_Field == "2009_production_20"|
            Year_Field == "2010_production_21"|
            Year_Field == "2007_production_22"|
            Year_Field == "2007_production_23"|
            Year_Field == "2007_production_24"|
            Year_Field == "2007_production_25"|
            Year_Field == "2007_production_26"|
            Year_Field == "2007_production_27"|
            Year_Field == "2007_production_28"|
            Year_Field == "2007_production_29"|
            Year_Field == "2007_production_30"|
            State == "Baseline")
 
 plot.pico.3.cristian <- picoxystrobin.complete.fields.cristian %>%
   ggplot(aes(x = State, y = EC50DC)) +
   geom_jitter(
     width = .1,
     height = 0,
     shape = 21,
     color = "black",
     fill = "chocolate1",
     size = 2,
     alpha = 3 / 4
   ) + stat_summary(
     fun.y = mean,
     geom = "point",
     shape = 95,
     size = 15,
     color = "black"
   ) +
   labs(x = "State", y = expression(bold(EC[bold("50") ~ (bold("D"))]) ~
                                      (bold(ppm ~ bold(
                                        "a.i."
                                      ))))) + theme(
                                        panel.border = element_rect(
                                          colour = "black",
                                          fill = NA,
                                          size = 1
                                        ),
                                        axis.title = element_text(size = 10, face = "bold", hjust = 0.5),
                                        axis.text = element_text(
                                          face = "bold",
                                          size = 10,
                                          family = "Arial"
                                        ),
                                        panel.background = element_rect(fill = "white", colour = "grey50")
                                      ) + scale_y_continuous(breaks = c(0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.020,0.022,0.024))
 
 
 
 
 plot.pico.3.cristian
 
 
 picoxystrobin.complete.fields.cristian$State <-
   factor(picoxystrobin.complete.fields.cristian$State,
          c( "WA", "ND", "MI", "NE", "CO", "Baseline"))


 shapiro.test.pico.fields.cristian <-picoxystrobin.complete.fields.cristian %>%
   do(tidy(shapiro.test(.$EC50DC)))
 shapiro.test.pico.fields.cristian 
 # NO NORMALILTY  0.000142
 # is by Field but it is written State
 # No normals
 object.3.cristian <- kruskal.test(EC50DC ~ State, data = picoxystrobin.complete.fields.cristian)
 object.3.cristian[[3]][[1]]
 
 #There is difference p-value = 0.044
 
 object.4.cristian <-  DunnTest(EC50DC ~ State, data = picoxystrobin.complete.fields.cristian, method = "bonferroni")
 object.4.cristian
 # 
 # mean.rank.diff   pval    
 # ND-WA           -0.6103896 1.0000    
 # MI-WA           -2.1048951 1.0000    
 # NE-WA          -13.6818182 1.0000    
 # CO-WA          -18.8484848 0.7677    
 # Baseline-WA    -19.3181818 0.2249    
 # MI-ND           -1.4945055 1.0000    
 # NE-ND          -13.0714286 1.0000    
 # CO-ND          -18.2380952 1.0000    
 # Baseline-ND    -18.7077922 0.6751    
 # NE-MI          -11.5769231 1.0000    
 # CO-MI          -16.7435897 1.0000    
 # Baseline-MI    -17.2132867 0.3321    
 # CO-NE           -5.1666667 1.0000    
 # Baseline-NE     -5.6363636 1.0000    
 # Baseline-CO     -0.4696970 1.0000    
 # ---
 
 ###tetraconazole
 
 tetraconazole.data.cristian <-
   read.csv("data/tetraconazole_cristian_01.15.20.csv")
 glimpse(tetraconazole.data.cristian)
 
 tetraconazole.data.cristian.2 <-tetraconazole.data.cristian %>%
   rename(ecuatorial= Radial_growth1..mm., polar= Radial_growth2..mm.) %>% 
   mutate(polar= replace(polar, polar == 0, 0.6)) %>%  #replacing 0 cm growth for the size of plug that is 0.6 cm or 6 mm
   mutate(ecuatorial= replace(ecuatorial, ecuatorial == 0, 0.6)) %>% #replacing 0 cm growth for the size of plug that is 0.6 or 6 mm
   group_by(Isolate_ID, Treatment, Rep) %>%
   mutate(growth = ((
     ecuatorial + polar
   ) / 2) / 10) %>% # getting average and transforming in cm
   select(-c(ecuatorial, polar, Notes)) %>%
   ungroup() %>%
   mutate(experimental_replicate = 1) %>%
   rename(repeats = Rep, ID = Isolate_ID) %>%
   select(ID, repeats,  experimental_replicate, growth, Treatment) %>%
   group_by(Treatment) %>% 
   mutate(grouped_id = row_number()) %>%
   spread(Treatment, growth) %>%
   select(-grouped_id) %>% 
   ungroup() %>%
   rename(control= Control, response= '2 ppm')
 
 #  Using the model to estimate the EC50DC
 tetraconazole.cristian.filtered. <-tetraconazole.data.cristian.2 %>%
   group_by(ID) %>%
   summarise(
     mean_response = mean(response, na.rm = TRUE),
     mean_control = mean(control, na.rm = TRUE)
   ) %>%
   mutate(RG = (mean_response / mean_control) * 100) %>%
   
   mutate(Estimate.50DC = exp (finalRG2[[1]][[1]] + finalRG2[[1]][[2]]*RG)) %>%
   ungroup() %>%
   rename(EC50DC= Estimate.50DC) %>% 
   mutate(ID = as.numeric(ID)) %>%
   select(c(ID, EC50DC)) %>% 
   mutate(new = ifelse(
     ID %in% picoxystrobin.filtered.survey.yes,
     "Yes","No")) 
 
 
 #Selecting the EC50D from the serial dilution
 tetraconazole.complete.cristian <- final.tetraconazole.DC %>% 
   ####WATCH OUT, Based on EC50DC
   rename(EC50DC = Estimate.50DC) %>% 
   select(ID, EC50DC) %>% 
   bind_rows(tetraconazole.filtered. ) %>% 
   #mutate(new= as.character(new)) %>% 
   bind_rows(tetraconazole.cristian.filtered. ) %>% 
   mutate(new= as.factor(new)) %>% 
   ungroup()
 
 # keeping juts the highest value of the repeated ID value
 tetraconazole.complete.cristian<-tetraconazole.complete.cristian%>%  
   mutate( ID= as.numeric(ID)) %>% 
   group_by(ID) %>%
   top_n(1, EC50DC) %>% #keeping the highest value of each ID value
   ungroup() 
 
 
 
 
 #Filterring specific fields by field scale
 tetraconazole.complete.cris <-tetraconazole.complete.cristian%>% 
   mutate(source = ifelse(
     ID %in% baseline_isolates_2,
     "Baseline",
     ifelse(
       ID %in% Farmer_fields,
       "Producer Fields","Fungicide Field Trials"
     )), Host = ifelse(
       ID %in% drybean,
       "Drybean", ifelse(
         ID %in% soybean,
         "Soybean","different_host")),
     Field = ifelse(
       ID %in% NE_Holt,
       "NE_Holt",
       ifelse(
         ID %in% NE_Field_34,
         "NE:Dodge",
         ifelse(
           ID %in% NE_Field_67,
           "NE:Cedar 1",
           ifelse(
             ID %in% NE_Field_75,
             "NE:Antelope",
             ifelse(
               ID %in% NE_Field_97,
               "NE:Pierce",
               ifelse(
                 ID %in% NE_Field_98,
                 "NE:Cedar 3",
                 ifelse(
                   ID %in% NE_Field_101,
                   "NE:Cedar 2",
                   ifelse(ID %in% potato, 
                          "potato", ifelse(ID %in% baseline_isolates_2, 
                                           "Baseline","other_fields"
                                           
                          )))
               )
             )
           )
         )
       )
     )
   )%>%
   mutate(
     source = as.factor(source),
     Host = as.factor(Host),
     Field = as.factor(Field)
   ) %>%    
   #taking out "different host" becuse does not belong to drybean nor soybean
   filter(!ID == 8, !ID == 129)%>%
   mutate(Country = ifelse(
     ID %in% mexican.isolates,
     "Mexico",
     ifelse(ID %in% brazilian.isolates,
            "Brazil", "USA")
   )) 
 
 tetraconazole.complete.cris <- tetraconazole.complete.cris %>% 
   mutate(Year_Field= case_when(
     ID %in%WM_1963~"1963",
     ID %in%WM_1965~"1965",
     ID %in%WM_1969~"1969",
     ID %in%WM_1973~"1973",
     ID %in%WM_1974~"1974",
     ID %in%WM_1975~"1975",
     ID %in%WM_1976~"1976",
     ID %in%WM_1977~"1977",
     ID %in%WM_1978~"1978",
     ID %in%WM_1979~"1979",
     ID %in%WM_1980~"1980",
     ID %in%WM_1981~"1981",
     ID %in%WM_1982~"1982",
     ID %in%WM_1983~"1983",
     ID %in%WM_1984~"1984",
     ID %in%WM_1985~"1985",
     ID %in%WM_1987~"1987",
     ID %in%WM_1989~"1989",
     ID %in%WM_1990~"1990",
     ID %in%WM_1991~"1991",
     ID %in%WM_1992~"1992",
     ID %in%WM_1993~"1993",
     ID %in%WM_1994~"1994",
     ID %in%WM_1995~"1995",
     ID %in%WM_1996~"1996",
     ID %in%WM_1997~"1997",
     ID %in%WM_1997_WM_production_31~"1997_production_31",
     ID %in%WM_1997_WM_production_32~"1997_production_32",
     ID %in%WM_1988~"1988",
     ID %in%WM_1998~"1998",
     ID %in%WM_1999~"1999",
     ID %in%WM_2000~"2000",
     ID %in%WM_2002~"2002",
     ID %in%WM_2003~"2003",
     ID %in%WM_1997_WM_production_33~"1997_production_33",
     ID %in%WM_2003_WM_nursery_1~"2003_nursery_1",
     ID %in%WM_2003_WM_nursery_2~"2003_nursery_2",
     ID %in%WM_2003_WM_nursery_3~"2003_nursery_3",
     ID %in%WM_1996_WM_nursery_4~"1996_nursery_4",
     ID %in%WM_2003_WM_nursery_5~"2003_nursery_5",
     ID %in%WM_2004~"2004",
     ID %in%WM_2005~"2005",
     ID %in%WM_2006~"2006",
     ID %in%WM_2004_WM_nursery_6~"2004_nursery_6",
     ID %in%WM_2004_WM_nursery_7~"2004_nursery_7",
     ID %in%WM_2004_WM_nursery_8~"2004_nursery_8",
     ID %in%WM_2004_WM_nursery_9~"2004_nursery_9",
     ID %in%WM_2004_WM_nursery_10~"2004_nursery_10",
     ID %in%WM_2004_WM_nursery_11~"2004_nursery_11",
     ID %in%WM_2004_WM_nursery_12~"2004_nursery_12",
     ID %in%WM_2005_WM_nursery_13~"2005_nursery_13",
     ID %in%WM_2005_WM_nursery_14~"2005_nursery_14",
     ID %in%WM_2005_WM_nursery_15~"2005_nursery_15",
     ID %in%WM_2005_WM_nursery_16~"2005_nursery_16",
     ID %in%WM_2005_WM_nursery_17~"2005_nursery_17",
     ID %in%WM_2005_WM_nursery_18~"2005_nursery_18",
     ID %in%WM_2006_WM_nursery_19~"2006_nursery_19",
     ID %in%WM_2007~"2007",
     ID %in%WM_2007_WM_production_1~"2007_production_1",
     ID %in%WM_2007_WM_production_2~"2007_production_2",
     ID %in%WM_2007_WM_production_3~"2007_production_3",
     ID %in%WM_2007_WM_production_4~"2007_production_4",
     ID %in%WM_2007_WM_production_5~"2007_production_5",
     ID %in%WM_2007_WM_production_6~"2007_production_6",
     ID %in%WM_2007_WM_production_7~"2007_production_7",
     ID %in%WM_2007_WM_production_8~"2007_production_8",
     ID %in%WM_2007_WM_production_9~"2007_production_9",
     ID %in%WM_2008~"2008",
     ID %in%WM_2008_WM_production_10~"2008_production_10",
     ID %in%WM_2008_WM_nursery_20~"2008_nursery_20",
     ID %in%WM_2008_WM_production_11~"2008_production_11",
     ID %in%WM_2008_WM_production_12~"2008_production_12",
     ID %in%WM_2008_WM_production_13~"2008_production_13",
     ID %in%WM_2009~"2009",
     ID %in%WM_2008_WM_nursery_21~"2008_nursery_21",
     ID %in%WM_2009_WM_nursery_22~"2009_nursery_22",
     ID %in%WM_2009_WM_nursery_23~"2009_nursery_23",
     ID %in%WM_2009_WM_production_14~"2009_production_14",
     ID %in%WM_2009_WM_production_15~"2009_production_15",
     ID %in%WM_2009_WM_production_16~"2009_production_16",
     ID %in%WM_2009_WM_production_17~"2009_production_17",
     ID %in%WM_2009_WM_production_18~"2009_production_18",
     ID %in%WM_2009_WM_production_19~"2009_production_19",
     ID %in%WM_2009_WM_production_20~"2009_production_20",
     ID %in%WM_2010_WM_production_21~"2010_production_21",
     ID %in%WM_2010_WM_production_22~"2010_production_22",
     ID %in%WM_2010_WM_production_23~"2010_production_23",
     ID %in%WM_2010_WM_production_24~"2010_production_24",
     ID %in%WM_2010_WM_production_25~"2010_production_25",
     ID %in%WM_2010_WM_production_26~"2010_production_26",
     ID %in%WM_2010_WM_production_27~"2010_production_27",
     ID %in%WM_2010_WM_production_28~"2010_production_28",
     ID %in%WM_2010_WM_production_29~"2010_production_29",
     ID %in%WM_2010_WM_production_30~"2010_production_30",
     ID %in%WM_2011~"2011",
     ID %in%WM_2013~"2013",
     ID %in%WM_2014~"2014",
     ID %in%WM_2016~"2016",
     ID %in%WM_2016_po~"2016_po",
     ID %in%WM_2016_flo~"2016_flo",
     ID %in%WM_2016_lan~"2016_lan",
     ID %in%WM_2017_ho~"2017_ho",
     ID %in%WM_2017_an~"2017_an",
     ID %in%WM_2017_dod~"2017_dod",
     ID %in%WM_2017_va~"2017_va",
     ID %in%WM_2017_na~"2017_na",
     ID %in%WM_2017_sa~"2017_sa",
     ID %in%WM_2017_ha~"2017_ha",
     ID %in%WM_2017_cu~"2017_cu",
     ID %in%WM_2017_mon~"2017_mon",
     ID %in%WM_2017_ant~"2017_ant",
     ID %in%WM_2017_holt~"2017_holt",
     ID %in%WM_2017_arr~"2017_arr",
     ID %in%WM_2017_ini~"2017_ini",
     ID %in%WM_2017_ley~"2017_ley",
     ID %in%WM_2017_pre~"2017_pre",
     ID %in%WM_2017_ru~"2017_ru",
     ID %in%WM_2017_ca~"2017_ca",
     ID %in%WM_2017_st~"2017_st",
     ID %in%WM_2015_6~"2015_6",
     ID %in%WM_2015_8~"2015_8",
     ID %in%WM_2015_9~"2015_9",
     ID %in%WM_2015_10~"2015_10",
     ID %in%WM_2015_11~"2015_11",
     ID %in%WM_2015_18~"2015_18",
     ID %in%WM_2015_19~"2015_19",
     ID %in%WM_2015_20~"2015_20",
     ID %in%WM_2015_21~"2015_21",
     ID %in%WM_2015_22~"2015_22",
     ID %in%WM_2015_23~"2015_23",
     ID %in%WM_2015_24~"2015_24",
     ID %in%WM_2015_25~"2015_25",
     ID %in%WM_2015_34~"2015_34",
     ID %in%WM_2015_35~"2015_35",
     ID %in%WM_2015_45~"2015_45",
     ID %in%WM_2015_46~"2015_46",
     ID %in%WM_2015_47~"2015_47",
     ID %in%WM_2015_48~"2015_48",
     ID %in%WM_2015_60~"2015_60",
     ID %in%WM_2015_62~"2015_62",
     ID %in%WM_2015_66~"2015_66",
     ID %in%WM_2015_67~"2015_67",
     ID %in%WM_2015_68~"2015_68",
     ID %in%WM_2015_69~"2015_69",
     ID %in%WM_2015_75~"2015_75",
     ID %in%WM_2015_78~"2015_78",
     ID %in%WM_2015_97~"2015_97",
     ID %in%WM_2015_98~"2015_98",
     ID %in%WM_2015_99~"2015_99",
     ID %in%WM_2015_100~"2015_100",
     ID %in%WM_2015_101~"2015_101",
     ID %in%WM_2015_102~"2015_102",
     ID %in%WM_2015_103~"2015_103",
     ID %in%WM_2015_104~"2015_104",
     ID %in%WM_2015_Al~"2015_Al",
     ID %in%WM_2015_Hi~"2015_Hi",
     ID %in%WM_2015_In~"2015_In",
     ID %in%WM_2015_mont~"2015_mont",
     ID %in%WM_2015_Sand~"2015_Sand",
     ID %in%WM_2015_Sani~"2015_Sani",
     ID %in%WM_2015_St~"2015_St",
     ID %in%WM_2015_wau~"2015_wau",
     ID %in%WM_2015_la~"2015_la",
     ID %in%WM_2015_co~"2015_co",
     ID %in%WM_2015_cern~"2015_cern",
     ID %in%WM_2015_cer~"2015_cer",
     ID %in%WM_2010~"2010",
     ID %in%WM_2017~"2017",
     TRUE ~ "no specific field"
   ),
   State= case_when(
     ID %in%WM_AB~"AB",
     ID %in%WM_AZ~"AZ",
     ID %in%WM_BA~"BA",
     ID %in%WM_CA~"CA",
     ID %in%WM_CO~"CO",
     ID %in%WM_DE~"DE",
     ID %in%WM_DF~"DF",
     ID %in%WM_FL~"FL",
     ID %in%WM_GA~"GA",
     ID %in%WM_GO~"GO",
     
     ID %in%WM_IA~"IA",
     ID %in%WM_ID~"ID",
     ID %in%WM_KS~"KS",
     ID %in%WM_KY~"KY",
     ID %in%WM_LA~"LA",
     ID %in%WM_MD~"MD",
     ID %in%WM_MG~"MG",
     ID %in%WM_MI~"MI",
     ID %in%WM_MN~"MN",
     ID %in%WM_MO~"MO",
     
     ID %in%WM_MS~"MS",
     ID %in%WM_MT~"MT",
     ID %in%WM_NC~"NC",
     ID %in%WM_ND~"ND",
     ID %in%WM_NE~"NE",
     ID %in%WM_NY~"NY",
     ID %in%WM_OH~"OH",
     ID %in%WM_OK~"OK",
     ID %in%WM_ON~"ON",
     ID %in%WM_OR~"OR",
     
     ID %in%WM_PA~"PA",
     ID %in%WM_PR~"PR",
     ID %in%WM_QUE~"QUE",
     ID %in%WM_RS~"RS",
     ID %in%WM_SC~"SC",
     ID %in%WM_SIN~"SIN",
     ID %in%WM_SK~"SK",
     ID %in%WM_WI~"WI",
     ID %in%WM_WA~"WA"
     
   )
   ) %>% select(-c(new)) %>% 
   mutate(Country = as.factor(Country),
          Year_Field = as.factor(Year_Field),
          State= as.factor(State) ) %>% 
   arrange(ID)
 
 ###For comparisons between hosts I had to take out the ones coming from baseline otherwise would not be fair
 
 tetraconazole.complete.host <-tetraconazole.complete %>% 
   filter(!source=="Baseline")
 
 ###For isolates from that have reduced sensitivity, assess sensitivity for isolates in the same field to account for heterogeneicity
 ##WITH MEXICO & BRAZIL together
 # tetraconazole.complete.fields.cristian <-tetraconazole.complete.cris%>%
 #   filter(  Field== "Baseline"|
 #              State== "CO"|
 #              State== "NE"|
 #              State== "ND"|
 #              State== "WA"|
 #              State== "MI"|
 #              Country== "Mexico"|
 #              Country== "Brazil"
 #   ) %>%
 #   arrange(source) %>%
 #    #
 #   mutate(State= as.character(State),
 #          State =  ifelse(Field == "Baseline", "Baseline", State),
 #          State = as.factor(State)) %>%
 #   filter(Year_Field == "2007_production_1"|
 #            Year_Field == "2007_production_2"|
 #            Year_Field == "2007_production_3"|
 #            Year_Field == "2007_production_4"|
 #            Year_Field == "2007_production_5"|
 #            Year_Field == "2007_production_6"|
 #            Year_Field == "2007_production_7"|
 #            Year_Field == "2007_production_8"|
 #            Year_Field == "2007_production_9"|
 #            Year_Field == "2008_production_10"|
 #            Year_Field == "2008_production_11"|
 #            Year_Field == "2008_production_12"|
 #            Year_Field == "2008_production_13"|
 #            Year_Field == "2009_production_17"|
 #            Year_Field == "2009_production_19"|
 #            Year_Field == "2009_production_20"|
 #            Year_Field == "2010_production_21"|
 #            Year_Field == "2007_production_22"|
 #            Year_Field == "2007_production_23"|
 #            Year_Field == "2007_production_24"|
 #            Year_Field == "2007_production_25"|
 #            Year_Field == "2007_production_26"|
 #            Year_Field == "2007_production_27"|
 #            Year_Field == "2007_production_28"|
 #            Year_Field == "2007_production_29"|
 #            Year_Field == "2007_production_30"|
 #            State == "Baseline"|
 #            Country== "Mexico"|
 #            Country== "Brazil") %>%
 #           mutate(State = fct_recode(State, "Mexico & Brazil"= "SIN" , "Mexico & Brazil"= "GO", "Mexico & Brazil" = "MG"))
 # 

 tetraconazole.complete.fields.cristian <-tetraconazole.complete.cris%>% 
   filter(  Field== "Baseline"|
              State== "CO"|
              State== "NE"|
              State== "ND"| 
              State== "WA"|
              State== "MI"
   ) %>% 
   arrange(source) %>% 
   filter( !source == "Fungicide Field Trials") %>% 
   
   mutate(State= as.character(State), 
          State =  ifelse(Field == "Baseline", "Baseline", State),
          State = as.factor(State)) %>% 
   filter(Year_Field == "2007_production_1"|
            Year_Field == "2007_production_2"|
            Year_Field == "2007_production_3"|
            Year_Field == "2007_production_4"|
            Year_Field == "2007_production_5"|
            Year_Field == "2007_production_6"|
            Year_Field == "2007_production_7"|
            Year_Field == "2007_production_8"|
            Year_Field == "2007_production_9"|
            Year_Field == "2008_production_10"|
            Year_Field == "2008_production_11"|
            Year_Field == "2008_production_12"|
            Year_Field == "2008_production_13"|
            Year_Field == "2009_production_17"|
            Year_Field == "2009_production_19"|
            Year_Field == "2009_production_20"|
            Year_Field == "2010_production_21"|
            Year_Field == "2007_production_22"|
            Year_Field == "2007_production_23"|
            Year_Field == "2007_production_24"|
            Year_Field == "2007_production_25"|
            Year_Field == "2007_production_26"|
            Year_Field == "2007_production_27"|
            Year_Field == "2007_production_28"|
            Year_Field == "2007_production_29"|
            Year_Field == "2007_production_30"|
            State == "Baseline") %>% 
   arrange(EC50DC)
 
 plot.tetra.3.cristian <- tetraconazole.complete.fields.cristian %>%
   ggplot(aes(x = State, y = EC50DC)) +
   geom_jitter(
     width = .1,
     height = 0,
     shape = 21,
     color = "black",
     fill = "deepskyblue",
     size = 2,
     alpha = 3 / 4
   ) + stat_summary(
     fun.y = mean,
     geom = "point",
     shape = 95,
     size = 15,
     color = "black"
   ) +
   labs(x = "State", y = expression(bold(EC[bold("50") ~ (bold("D"))]) ~
                                      (bold(ppm ~ bold(
                                        "a.i."
                                      ))))) + theme(
                                        panel.border = element_rect(
                                          colour = "black",
                                          fill = NA,
                                          size = 1
                                        ),
                                        axis.title = element_text(size = 10, face = "bold", hjust = 0.5),
                                        axis.text = element_text(
                                          face = "bold",
                                          size = 10,
                                          family = "Arial"
                                        ),
                                        panel.background = element_rect(fill = "white", colour = "grey50")
                                      ) + scale_y_continuous(breaks = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0))
 
 
 
 
 plot.tetra.3.cristian
 
 
 tetraconazole.complete.fields.cristian$State <-
   factor(tetraconazole.complete.fields.cristian$State,
 c("WA", "MI", "Baseline", "ND", "CO", "NE"))

 #by field but because the name actually is "State"
 
 ##NORMALITY TEST-Shapiro_test
 
 shapiro.test.tetra.cristian <- tetraconazole.complete.fields.cristian %>%
   do(tidy(shapiro.test(.$EC50DC)))
 shapiro.test.tetra.cristian[[2]][[1]]
 #1.159902e-05
 # No normals
 object.5.cristian <- kruskal.test(EC50DC ~ State, data = tetraconazole.complete.fields.cristian)
 object.5.cristian[[3]][[1]]
 
 #There is difference p-value = 0.007581091
 
 object.6.cristian <-  DunnTest(EC50DC ~ State, data = tetraconazole.complete.fields.cristian, method = "bonferroni")
 object.6.cristian
 
 
 # mean.rank.diff   pval    
 # MI-WA             6.258741 1.0000    
 # Baseline-WA      -1.865801 1.0000    
 # ND-WA            -9.818182 1.0000    
 # CO-WA           -10.818182 1.0000    
 # NE-WA           -23.318182 0.0908 .  
 # Baseline-MI      -8.124542 1.0000    
 # ND-MI           -16.076923 1.0000    
 # CO-MI           -17.076923 1.0000    
 # NE-MI           -29.576923 0.0042 ** 
 #   ND-Baseline      -7.952381 1.0000    
 # CO-Baseline      -8.952381 1.0000    
 # NE-Baseline     -21.452381 0.0537 .  
 # CO-ND            -1.000000 1.0000    
 # NE-ND           -13.500000 1.0000    
 # NE-CO           -12.500000 1.0000
 # 