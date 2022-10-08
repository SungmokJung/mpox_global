library(dplyr)
library(tidyverse)

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
write.csv(df_MSM_imputed, "../data/MSM_pop/df_MSM_imputed.csv")

merge(df_inci, df_MSM_imputed %>% dplyr::select(iso_code, imputed, pop2022, prop_avg), 
      by=c("iso_code"), all.x=TRUE) -> df_all_inci
write.csv(df_all_inci, "../data/MSM_pop/df_all_inci.csv")