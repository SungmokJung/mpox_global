libraries = c("dplyr", "tidyverse", "data.table", "readxl")
for(x in libraries) {library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE)}

theme_set(theme_bw())

path <- "../data/flight/all_region/"

options(warn=-1)
list.files(path = path, pattern = "*xlsx") -> file_list
read_xlsx_files <- function(x){
    df <- suppressMessages(read_xlsx(path = paste0(path, x), skip = 5))
    df[-1, c(1,2,27)] -> df_2019
    colnames(df_2019) <- c("destination", "series", "v2019")
    return(df_2019)}
lapply(file_list, read_xlsx_files) -> df_list

for(i in 1:length(df_list)){
    colnames(df_list[[i]]) <- c("destination", "series", substr(file_list[i],1,nchar(file_list[i])-5))
}

options(warn=0)

merge.all <- function(x, y) {merge(x, y, all=TRUE, by=c("destination", "series"))}
output <- Reduce(merge.all, df_list)

output %>% filter(!is.na(destination)) -> output

output %>% filter(destination %in% c("Czechia", "Czech Republic")) -> temp
unique(temp$series) -> series_list
Czech_list <- list()
for(i in 1:length(series_list)){
    temp %>% filter(series==series_list[i]) -> temp; temp[is.na(temp)] <- 0
    temp %>% dplyr::select(-c("destination", "series")) -> temp_col
    as.data.frame(colSums(temp_col))%>% t() %>% as.data.frame() %>%
    mutate(destination=c("Czech Republic"), series=series_list[i]) -> temp_all; rownames(temp_all) <- NULL
    temp_all -> Czech_list[[i]]
}

do.call(rbind, Czech_list) -> Czech_output

output %>% filter(destination %in% c("Türkiye", "Turkey")) -> temp
unique(temp$series) -> series_list
Turkey_list <- list()
for(i in 1:length(series_list)){
    temp %>% filter(series==series_list[i]) -> temp; temp[is.na(temp)] <- 0
    temp %>% dplyr::select(-c("destination", "series")) -> temp_col
    as.data.frame(colSums(temp_col))%>% t() %>% as.data.frame() %>%
    mutate(destination=c("Turkey"), series=series_list[i]) -> temp_all; rownames(temp_all) <- NULL
    temp_all -> Turkey_list[[i]]
}

do.call(rbind, Turkey_list) -> Turkey_output

output %>% filter(!(destination %in% c("END", "Türkiye", "Turkey", "Czechia", "Czech Republic"))) %>%
mutate(destination=case_when(destination==c("Curaçao")~c("Curacao"), 
                             destination==c("Côte d\'Ivoire")~c("Cote d'Ivoire"), 
                             destination==c("Taiwan Province of China")~c("Taiwan, Province of China"), 
                             destination==c("Saint Vincent and the Grenadines")~c("Saint Vincent and The Grenadines"),
                             destination==c("Palestine")~c("Palestine, State of"),
                             TRUE~destination)) -> output

rbind(output, Turkey_output, Czech_output) -> output

## coverting NAs to zero in the flight data
output[is.na(output)] <- 0
output %>% mutate(total_vol = rowSums(output[,3:ncol(output)])) %>% 
filter(total_vol > 0) %>% dplyr::select(-total_vol) -> output

## selecting the biggest value throughout the series
colMax <- function(data) sapply(data, max, na.rm = TRUE)

new_list <- list()
unique(output$destination) -> country_list


for (g in 1:length(country_list)){
    output %>% filter(destination==country_list[g]) -> temp
    as.data.frame(t(colMax(temp)))-> new_list[[g]]
}

            
do.call("rbind", new_list) -> flight_matrix
write.csv(flight_matrix, "../data/flight/flight_matrix.csv")     

flight_matrix %>% filter(destination==("Palestine, State of")) -> temp
temp[, 3:215] <- sapply(temp[, 3:215], as.numeric)
temp %>% dplyr::select_if(~ !is.numeric(.) || sum(.) != 0)

# ## selecting the series of data with prioritization
# new_list <- list()
# unique(output$destination) -> country_list

# for (g in 1:length(country_list)){
#     output %>% filter(destination==country_list[g]) -> temp
#     if(nrow(temp)>=2){
#         if(any(temp$series==c("VFR"))){temp %>% filter(series==c("VFR")) -> new_list[[g]]}
#         else if(any(temp$series==c("VFN")) && !(temp$series %in% c("VFR")))
#         {temp %>% filter(series==c("VFN")) -> new_list[[g]]}
#         else if(any(temp$series==c("TFR")) && !(temp$series %in% c("VFR", "VFN")))
#         {temp %>% filter(series==c("TFR")) -> new_list[[g]]}
#         else if(any(temp$series==c("TFN")) && !(temp$series %in% c("VFR", "VFR", "TFR")))
#         {temp %>% filter(series==c("TFN")) -> new_list[[g]]}
#         else if(any(temp$series==c("TCER")) && !(temp$series %in% c("VFR", "VFR", "TFR", "TFN")))
#         {temp %>% filter(series==c("TCER")) -> new_list[[g]]}
#         else if(any(temp$series==c("TCEN")) && !(temp$series %in% c("VFR", "VFR", "TFR", "TFN", "TCER")))
#         {temp %>% filter(series==c("TCEN")) -> new_list[[g]]}
#         else if(any(temp$series==c("THSR")) && !(temp$series %in% c("VFR", "VFR", "TFR", "TFN", "TCER", "TCEN")))
#         {temp %>% filter(series==c("THSR")) -> new_list[[g]]}
#         else if(any(temp$series==c("THSN")) && !(temp$series %in% c("VFR", "VFR", "TFR", "TFN", "TCER", "TCEN", "THSR")))
#         {temp %>% filter(series==c("THSN")) -> new_list[[g]]}
#     }
#     else{temp -> new_list[[g]]}
# }
            
# do.call("rbind", new_list) -> flight_matrix
# # flight_matrix %>% 
# write.csv(flight_matrix, "../data/flight/flight_matrix.csv")     

setdiff(flight_matrix$destination, colnames(flight_matrix[,3:ncol(flight_matrix)]))
setdiff(colnames(flight_matrix[,3:ncol(flight_matrix)]), flight_matrix$destination)


