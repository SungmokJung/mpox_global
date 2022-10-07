library(dplyr)

read.csv("../../data/incidence/owid-monkeypox-data.csv") -> df_inci
df_inci %<>% filter(!(location %in% c("World"))) %>%  
  dplyr::select(location, iso_code, date, new_cases, total_cases)

df_inci %>% arrange(date) %>% group_by(iso_code) %>% filter(rank(date)==1) %>% 
  dplyr::select(iso_code, date) %>% rename(date_import = date) -> temp
merge(df_inci, temp, by=c("iso_code"), all.x=TRUE) -> df_inci

df_inci %>% dplyr::select(iso_code, location) %>% distinct() -> country_list

write.csv(country_list, "../../data/incidence/country_list.csv")
write.csv(df_inci, "../../data/incidence/df_inci.csv")
