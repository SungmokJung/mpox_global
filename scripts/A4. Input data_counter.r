libraries = c("dplyr", "tidyverse", "ggpubr")
for(x in libraries) {library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE)}

theme_set(theme_bw())

read.csv("../data/df_inci_final_WHO_backproj.csv") -> df_inci
read.csv("../data/SAR_cip_Reff_excess.csv") -> df_Reff
read.csv("../data/flight/flight_matrix.csv") -> flight_matrix

## fixed parameters
w <- 14
SAR <- 0.1
time_0 <- as.Date("2022-04-17")
time_end <- as.Date("2022-10-01")

## Reff_i & G_i (fixing the country name issue)
df_inci %>% mutate(Reff_i=SAR, G_i=Reff_i*MA_new_cases) %>%
mutate(location=case_when(location==c("Democratic Republic of the Congo")~c("Congo, Democratic Republic of the"),
                          location==c("Curaçao")~c("Curacao"),
                          location==c("Czechia")~c("Czech Republic"),
                          location==c("Türkiye")~c("Turkey"),
                          location==c("Iran")~c("Iran, Islamic Republic of"),
                          location==c("Republic of Korea")~c("Korea, Republic of"),
                          location==c("United States")~c("United States of America"),
                          location==c("Venezuela (Bolivarian Republic of)")~c("Venezuela, Bolivarian Republic of"),
                          location==c("Republic of Moldova")~c("Moldova, Republic of"),
                          location==c("Russia")~c("Russian Federation"),
                          location==c("Bolivia")~c("Bolivia, Plurinational State of"),
                          location==c("Hong Kong")~c("Hong Kong, China"),
                          location==c("Republic of Congo") ~ c("Congo"),
                          location==c("Bosnia And Herzegovina") ~ c("Bosnia and Herzegovina"),
                          location==c("The United Kingdom") ~c("United Kingdom"),
                          location==c("Russia") ~ c("Russian Federation"),
                          TRUE~location)) %>%
filter(!(location %in% c("Gibraltar", "Guadeloupe", "Greenland", "Saint Martin"))) -> input

as.Date(input$date) -> input$date
input %>% mutate(time=as.numeric(date-time_0+1), censoring=0) %>% arrange(date) %>% dplyr::select(-X) -> input

## calculating Reff_i & G_i in countries without MPX importation
path <- "../data/flight/all_region/"; list.files(path = path, pattern = "*xlsx") -> file_list
substr(file_list,1,nchar(file_list)-5) -> flight_list
unique(input$location) -> country_list
sort(flight_list) -> flight_list_sort; sort(country_list) -> country_list_sort

setdiff(flight_matrix$destination, country_list_sort) -> country_no_list_sort
sort(country_no_list_sort) -> country_no_list_sort

as.data.frame(country_no_list_sort) %>% rename(location=country_no_list_sort) -> country_no_import

read.csv("../data/MSM_pop/df_MSM_imputed.csv") %>%
mutate(location=case_when(location==c("Côte d\'Ivoire")~c("Cote d'Ivoire"),
                          location==c("Macao")~c("Macao, China"),
                          location==c("Micronesia (Federated States of)")~c("Micronesia, Federated States of"),
                          location==c("Saint Vincent and the Grenadines")~c("Saint Vincent and The Grenadines"),
                          location==c("Virgin Islands (U.S.)")~c("United States Virgin Islands"),
                          location==c("Réunion")~c("Reunion"),
                          TRUE~location)) -> df_MSM_imputed

merge(country_no_import, df_MSM_imputed %>% dplyr::select(location, iso_code, imputed, pop2022, region, sub_region), 
      by=c("location"), all.x=TRUE) -> country_no_import_pop

df_inci_no_list <- list()

as.data.frame(seq(min(input$date), max(input$date),1)) -> temp_cal
colnames(temp_cal) <- c("date")
unique(country_no_import_pop$location) -> no_list

for(i in 1:length(no_list)){
    country_no_import_pop %>% filter(location==no_list[i]) %>% mutate(date=min(input$date))-> temp_inci
    merge(temp_cal, temp_inci, by=c("date"), all.x=TRUE) %>% 
    mutate(new_cases=NA, total_cases=NA, MA_new_cases=NA, MA_total_cases=NA, 
           date_import=NA, Reff_i=NA, G_i=NA, cum_icni_prop=NA) -> temp_inci_all

    temp_inci_all$location[is.na(temp_inci_all$location)] <- unique(temp_inci$location)
    temp_inci_all$iso_code[is.na(temp_inci_all$iso_code)] <- unique(temp_inci$iso_code)
    temp_inci_all$region[is.na(temp_inci_all$region)] <- unique(temp_inci$region)
    temp_inci_all$sub_region[is.na(temp_inci_all$sub_region)] <- unique(temp_inci$sub_region)
    temp_inci_all$new_cases[is.na(temp_inci_all$new_cases)] <- 0
    temp_inci_all$total_cases[is.na(temp_inci_all$total_cases)] <- 0
    temp_inci_all$MA_new_cases[is.na(temp_inci_all$MA_new_cases)] <- 0
    temp_inci_all$MA_total_cases[is.na(temp_inci_all$MA_total_cases)] <- 0
    temp_inci_all$Reff_i[is.na(temp_inci_all$Reff_i)] <- 0
    temp_inci_all$G_i[is.na(temp_inci_all$G_i)] <- 0
    temp_inci_all$cum_icni_prop[is.na(temp_inci_all$cum_icni_prop)] <- 0
    temp_inci_all$pop2022[is.na(temp_inci_all$pop2022)] <- unique(temp_inci$pop2022)
    temp_inci_all$imputed[is.na(temp_inci_all$imputed)] <- unique(temp_inci$imputed)
    
    temp_inci_all %>% dplyr::select(iso_code, location, date, new_cases, total_cases, 
                                    MA_new_cases, MA_total_cases, date_import, 
                                    imputed, pop2022, region, sub_region, cum_icni_prop, Reff_i, G_i) %>%
    rename(MSM_pop=imputed) -> temp_inci_all

    temp_inci_all -> df_inci_no_list[[i]]    
}

do.call("rbind", df_inci_no_list) %>% as.data.frame() %>% arrange(location, date) %>%
mutate(time=as.numeric(date-time_0+1), censoring=1) -> input_no

rbind(input, input_no) -> input_all
input_all %>% filter(date <= time_end) -> input_all
input_all %>% filter(is.na(pop2022)) %>% dplyr::select(location) %>% unique()

## countries without travel volume data
path <- "../data/flight/all_region/"; list.files(path = path, pattern = "*xlsx") -> file_list
substr(file_list,1,nchar(file_list)-5) -> flight_list
unique(input_all$location) -> country_list
country_list <- country_list[!country_list %in% c("Sudan", "Ghana", "Liberia", "Congo", "Nigeria" ,
                                                  "Congo, Democratic Republic of the",
                                                  "Venezuela, Bolivarian Republic of", "South Sudan", 
                                                  "Taiwan")] ## Taiwan can be added later
sort(flight_list) -> flight_list_sort; sort(country_list) -> country_list_sort
setdiff(country_list_sort, flight_matrix$destination)

## calculating F_i
F_i_country <- list(); F_i_time_list <- list()

for(i in 1:length(country_list_sort)){
    flight_matrix %>% filter(destination==country_list_sort[i]) %>% t() -> temp
    temp[4:nrow(temp),] %>% as.matrix() -> temp

    flight_list %>% as.matrix() -> flight_list_matrix
    cbind(temp, flight_list_matrix) %>% as.data.frame() %>% rename(value=V1, location=V2) -> V_ij_temp
    rownames(V_ij_temp) <- NULL

    input_all %>% dplyr::select(location, pop2022) %>% distinct() -> N_j_temp

    merge(V_ij_temp, N_j_temp, by=c("location"), all.x=TRUE) %>% na.omit() -> temp

    for(g in 1:max(input_all$time)) {
        input_all %>% filter(time==g) %>% dplyr::select(location, G_i) -> G_i_time
        merge(temp, G_i_time, by=c("location"), all.x=TRUE) -> temp_G_i_time
        temp_G_i_time$value <- as.numeric(temp_G_i_time$value)
        temp_G_i_time %>% mutate(F_i_t = w/365*value/pop2022*G_i, time=g) -> temp_G_i_time 
        sum(temp_G_i_time$F_i_t) -> F_i_time_list[[g]]
    }
    
    do.call("rbind", F_i_time_list) %>% as.data.frame() %>% 
    mutate(time=1:max(input_all$time), location=country_list_sort[i]) %>% rename(F_i=V1) -> F_i_country[[i]]
}

do.call("rbind", F_i_country) %>% as.data.frame() -> F_i_country_all
merge(input_all, F_i_country_all, by=c("location", "time"), all.x=TRUE) -> input_final_all

input_final_all$date_import <- as.Date(input_final_all$date_import)

## removing countires without travel volume data
input_final <- input_final_all[!is.na(input_final_all$F_i),]
write.csv(input_final, "../data/input_WHO_backproj_counter.csv")

input_final %>% group_by(location) %>% summarise(F_i_all=sum(F_i)) %>% ungroup() -> temp
temp%>% filter(F_i_all==0)






