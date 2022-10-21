libraries = c("dplyr", "tidyverse", "ggpubr", "rworldmap")
for(x in libraries) {library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE)}

theme_set(theme_bw())

## MPX incidence data
read.csv("../data/incidence/owid-monkeypox-data.csv") -> df_inci
df_inci %<>% filter(!(location %in% c("World"))) %>%  
dplyr::select(location, iso_code, date, new_cases, total_cases)

df_inci %>% arrange(date) %>% group_by(iso_code) %>% filter(rank(date)==1) %>% 
dplyr::select(iso_code, date) %>% rename(date_import = date) -> temp
merge(df_inci, temp, by=c("iso_code"), all.x=TRUE) -> df_inci

write.csv(df_inci, "../data/incidence/df_inci.csv")

## MSM population
read.csv("../data/MSM_pop/df_MSM_UNAIDS.csv") -> df_MSM1
read.csv("../data/MSM_pop/df_full_list_MSMsize_Fumi.csv") -> df_MSM2

## revising the pop2022 column from df_MSM2
read.csv("../data/MSM_pop/df_pop_raw.csv") %>% dplyr::select(cca2, pop2022) %>% rename(iso_2=cca2) -> df_pop1
read.csv("../data/MSM_pop/df_region.csv") %>% dplyr::select(alpha.2,alpha.3, name) %>%
rename(iso_2=alpha.2, iso_code=alpha.3, location=name) -> df_pop2
merge(df_pop1, df_pop2, by=c("iso_2")) -> df_pop
merge(df_MSM2 %>% dplyr::select(-pop2022), 
      df_pop %>% dplyr::select(iso_code, pop2022), by=c("iso_code"), all.x=TRUE) -> df_MSM2

## giving a priority to UNAIDS dashboard
merge(df_MSM1, df_MSM2, by=c("location"),all=TRUE) %>% 
dplyr::select(location, iso_code, estimate, MSM_size, region, sub_region, pop2022) %>%
mutate(prop=estimate/pop2022) %>%
mutate(imputed=case_when(!is.na(estimate) ~ estimate, is.na(estimate) & !is.na(MSM_size) ~ MSM_size), 
       prop=imputed/pop2022) -> temp

## imputation of missing values
temp %>% group_by(sub_region) %>% drop_na(imputed) %>% summarise(prop_avg = mean(prop)) -> prop_avg
merge(temp, prop_avg, by=c("sub_region"), all.x=TRUE) %>% drop_na(pop2022) %>%
mutate(prop=case_when(!is.na(prop)~prop, is.na(prop)~prop_avg),
       imputed=case_when(!is.na(imputed)~imputed, is.na(imputed)~prop*pop2022)) -> df_MSM_imputed


## modifying Samoa with the UNAIDS report
df_MSM_imputed %>% mutate(imputed=case_when(location==c("Samoa")~MSM_size, TRUE~imputed),
                          prop=case_when(location==c("Samoa")~MSM_size/pop2022, TRUE~prop)) -> df_MSM_imputed

df_MSM_imputed %>% head()
write.csv(df_MSM_imputed, "../data/MSM_pop/df_MSM_imputed.csv")

## merging incidence and MSM population data
merge(df_inci, df_MSM_imputed %>% dplyr::select(iso_code, imputed, pop2022, region, sub_region), 
      by=c("iso_code"), all.x=TRUE) -> df_all_inci
write.csv(df_all_inci, "../data/MSM_pop/df_all_inci.csv")

## figure for the proportion of MSM among total population
options(repr.plot.width=10,repr.plot.height=8)
df_MSM_imputed %>% dplyr::select(iso_code, prop) -> fig_MSM
Map = joinCountryData2Map(fig_MSM, joinCode = "ISO3", nameJoinColumn ="iso_code", mapResolution = "coarse") 
mapCountryData(Map, nameColumnToPlot="prop")

## adding data with zero case from day 0
day_0 <- as.Date("2022-04-20"); day_0_real <- as.Date("2022-05-01")

as.Date(df_all_inci$date) -> df_all_inci$date
as.data.frame(seq(day_0, max(df_all_inci$date),1)) -> temp_cal; colnames(temp_cal) <- c("date")
unique(df_all_inci$location) -> MPX_country_list

df_inci_list <- list()
for(i in 1:length(MPX_country_list)){
    df_all_inci %>% filter(location==MPX_country_list[i]) -> temp_inci
    merge(temp_cal, temp_inci, by=c("date"), all.x=TRUE) -> temp_inci_all

    temp_inci_all$location[is.na(temp_inci_all$location)] <- unique(temp_inci$location)
    temp_inci_all$iso_code[is.na(temp_inci_all$iso_code)] <- unique(temp_inci$iso_code)
    temp_inci_all$region[is.na(temp_inci_all$region)] <- unique(temp_inci$region)
    temp_inci_all$sub_region[is.na(temp_inci_all$sub_region)] <- unique(temp_inci$sub_region)
    temp_inci_all$new_cases[is.na(temp_inci_all$new_cases)] <- 0
    temp_inci_all$total_cases[is.na(temp_inci_all$total_cases)] <- 0
    temp_inci_all$date_import[is.na(temp_inci_all$date_import)] <- unique(temp_inci$date_import)
    temp_inci_all$pop2022[is.na(temp_inci_all$pop2022)] <- unique(temp_inci$pop2022)
    temp_inci_all$imputed[is.na(temp_inci_all$imputed)] <- unique(temp_inci$imputed)
    
    temp_inci_all %>% dplyr::select(iso_code, location, date, new_cases, total_cases, date_import, 
                                    imputed, pop2022, region, sub_region) -> temp_inci_all

    temp_inci_all -> df_inci_list[[i]]    
}

do.call("rbind", df_inci_list) %>% as.data.frame() %>% arrange(location, date) -> df_all_inci
as.Date(df_all_inci$date_import) -> df_all_inci$date_import

## applying a 14-day moving average for the newly reported and cumulative case counts
df_all_inci %>% group_by(location) %>% arrange(date) %>%
mutate(MA_new_cases = zoo::rollmean(new_cases, k=14, align="right", fill=NA),
       MA_total_cases = zoo::rollmean(total_cases, k=14, align="right", fill=NA)) -> temp
temp$MA_new_cases[is.na(temp$MA_new_cases)]<-0; temp$MA_total_cases[is.na(temp$MA_total_cases)]<-0; 

## calculating the cumlative incidence proportion (i.e., cumulative incidecne/MSM population)
temp %>% rename(MSM_pop=imputed) %>% mutate(MSM_pop=MSM_pop*1000) %>% 
group_by(location) %>% mutate(cum_icni_prop=MA_total_cases/(MSM_pop)) %>%
dplyr::select(iso_code, location, region, sub_region, date, new_cases, total_cases, MA_new_cases, MA_total_cases,
              date_import, MSM_pop, pop2022, cum_icni_prop) %>%
filter(date >= day_0_real) -> df_inci_MA

write.csv(df_inci_MA, "../data/df_inci_final.csv")
df_inci_MA %>% head()

unique(df_inci_MA$location)

## checking the smoothened incidence with the 14-day moving average
options(repr.plot.width=7,repr.plot.height=5)
df_inci_MA %>% filter(location==c("United States")) %>%
    ggplot() + 
    geom_bar(aes(x=date, y=MA_new_cases), fill="#FAAB18" ,stat='identity', width=0.7) +
    ggtitle("United States") +
    labs(x="\n Calendar date", y="Incidence \n") +
    theme(text = element_text(size=15, family="sans",color="black"),
          axis.text = element_text(size=12, family="sans",color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(size=15, family="sans",color="black"),
          axis.text.x=element_text(angle=60, hjust=1)) +
scale_x_date(date_labels = "%b %Y", date_breaks = "2 weeks", expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0)) +
coord_cartesian(ylim = c(0,600)) -> Fig1

df_inci_MA %>% filter(location==c("United Kingdom")) %>%
    ggplot() + 
    geom_bar(aes(x=date, y=MA_new_cases), fill="#FAAB18" ,stat='identity', width=0.7) +
    ggtitle("United Kingdom") +
    labs(x="\n Calendar date", y="Incidence \n") +
    theme(text = element_text(size=15, family="sans",color="black"),
          axis.text = element_text(size=12, family="sans",color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(size=15, family="sans",color="black"),
          axis.text.x=element_text(angle=60, hjust=1)) +
scale_x_date(date_labels = "%b %Y", date_breaks = "2 weeks", expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0)) +
coord_cartesian(ylim = c(0,100)) -> Fig2

df_inci_MA %>% filter(location==c("Belgium")) %>%
    ggplot() + 
    geom_bar(aes(x=date, y=MA_new_cases), fill="#FAAB18" ,stat='identity', width=0.7) +
    ggtitle("Belgium") +
    labs(x="\n Calendar date", y="Incidence \n") +
    theme(text = element_text(size=15, family="sans",color="black"),
          axis.text = element_text(size=12, family="sans",color="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(size=15, family="sans",color="black"),
          axis.text.x=element_text(angle=60, hjust=1)) +
scale_x_date(date_labels = "%b %Y", date_breaks = "2 weeks", expand = c(0, 0)) +
scale_y_continuous(expand = c(0, 0)) +
coord_cartesian(ylim = c(0,30)) -> Fig3

options(repr.plot.width=21,repr.plot.height=5)
ggarrange(Fig1, Fig2, Fig3, ncol = 3, nrow = 1, labels = c("A", "B", "C"), font.label = list(size = 20))
