libraries = c("dplyr", "tidyverse", "magrittr")
for(x in libraries) {library(x,character.only=TRUE,warn.conflicts=FALSE,quietly=TRUE)}
theme_set(theme_bw())

#### setting
start_date <- as.Date("2022-04-17") ## symptom onset date of the initial case in the UK

#### data with depletion effect
read.csv("../data/input_WHO_backproj.csv") -> df_input
df_input$date <- as.Date(df_input$date)
df_input$date_import <- as.Date(df_input$date_import)
censor_import <- max(df_input$date) ## for countries without any importation

## excluding the endemic countries (along with the UK)
df_input %<>% rename(censor=censoring, country=location) %>%
mutate(date_imp=case_when(censor==1~censor_import, censor==0~date_import)) %>%
filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                         "Congo","Congo, Democratic Republic of the","Ghana","Palestine, State of",
                         "United Kingdom")))

df_input %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> df_input_all

## regional specific data
data_list <- list()
df_input %>% filter(region == c("Europe")) %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list[[1]]

df_input %>% filter(region == c("Africa")) %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list[[2]]

df_input %>% filter(region == c("Americas")) %>%  dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list[[3]]

df_input %>% filter(region == c("Asia")) %>% filter(!sub_region %in% c("Central Asia", "Western Asia")) %>% 
dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list[[4]]

df_input %>% filter(sub_region %in% c("Central Asia", "Western Asia")) %>% 
dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list[[5]]

df_input %>% filter(region == c("Oceania")) %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list[[6]]

#### data without depletion effect
read.csv("../data/input_WHO_backproj_counter.csv") -> df_input_counter
df_input_counter$date <- as.Date(df_input_counter$date)
df_input_counter$date_import <- as.Date(df_input_counter$date_import)

## excluding the endemic countries (along with the UK)
df_input_counter %<>% rename(censor=censoring, country=location) %>%
mutate(date_imp=case_when(censor==1~censor_import, censor==0~date_import)) %>% 
filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                         "Congo","Congo, Democratic Republic of the","Ghana","Palestine, State of",
                         "United Kingdom")))

df_input_counter %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> df_input_counter_all

## regional specific data
data_list_counter <- list()
df_input_counter %>% filter(region == c("Europe")) %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list_counter[[1]]

df_input_counter %>% filter(region == c("Africa")) %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list_counter[[2]]

df_input_counter %>% filter(region == c("Americas")) %>%dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list_counter[[3]]

df_input_counter %>% filter(region == c("Asia")) %>% filter(!sub_region %in% c("Central Asia", "Western Asia")) %>% 
dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list_counter[[4]]

df_input_counter %>% filter(sub_region %in% c("Central Asia", "Western Asia")) %>% 
dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list_counter[[5]]

df_input_counter %>% filter(region == c("Oceania")) %>% dplyr::select(date, country, F_i, date_imp, censor) %>% 
filter(date >= start_date) %>% arrange(date, country) -> data_list_counter[[6]]

#### likelihood without random effect
LogL_i <- function(data, country_i){
    data_i <- data %>% filter(country==country_i)
    date_start_i <- data_i[1,1]
    date_import_i <- data_i[length(data_i[,1]), 4]
    surv_days_i <- as.numeric(date_import_i - date_start_i)+1
    F_i_vec <- data_i$F_i
    cens_i <- data_i[1,5] 
    
    function(alpha){
        return(
            (1-cens_i) * (log(alpha * F_i_vec[surv_days_i]) + (-sum(alpha * F_i_vec[1:surv_days_i]))) + 
            cens_i * (-sum(alpha * F_i_vec[1:surv_days_i])) 
        )
    }
}

LogL_full <- function(data, country_list){
  function(alpha){
      return(sum(sapply(country_list, FUN = function(x){LogL_i(data=data, country_i = x)(alpha=alpha)})))
  }
}

#### MLE for the model with depletion effect
## with global scaling factor
options(warn=-1)
optim(fn=LogL_full(data=df_input_all, country_list=(unique(df_input_all$country))), 
      par=c(0.001), method="Brent", lower=(0), upper=(1000), control = list(fnscale = -1, maxit=1000000)) -> est_all

## with regional specific scaling factors (separate MLE by region)
par_list <- list(); value_list <- list()
for(i in 1:length(data_list)){
    optim(fn=LogL_full(data=data_list[[i]], country_list=(unique(data_list[[i]]$country))), 
          par=c(0.001), method="Brent", lower=(0), upper=(1000), control = list(fnscale = -1, maxit=1000000)) -> est
    est$par -> par_list[[i]]; est$value -> value_list[[i]]
    
}

#### MLE for the model without depletion effect
## with global scaling factor
optim(fn=LogL_full(data=df_input_counter_all, country_list=(unique(df_input_counter_all$country))), 
      par=c(0.001), method="Brent", lower=(0), upper=(1000), 
      control = list(fnscale = -1, maxit=1000000)) -> est_all_counter

## with regional specific scaling factors (separate MLE by region)
par_list_counter <- list(); value_list_counter <- list()
for(i in 1:length(data_list_counter)){
    optim(fn=LogL_full(data=data_list_counter[[i]], country_list=(unique(data_list_counter[[i]]$country))), 
          par=c(0.001), method="Brent", lower=(0), upper=(1000), 
          control = list(fnscale = -1, maxit=1000000)) -> est
    est$par -> par_list_counter[[i]]; est$value -> value_list_counter[[i]]
    
}
options(warn=0)

saveRDS(par_list, "par_list.RDS")
saveRDS(par_list_counter, "par_list_counter.RDS")

readRDS("par_list.RDS") -> par_list
readRDS("par_list_counter.RDS") -> par_list_counter

#### likelihood with regional specific scaling factor
logL_region <- function(params){
    llk_all <- rep(0, length(data_list))
    
    for(k in 1:length(data_list)){
        llk <- rep(0, length(unique(data_list[[k]]$country)))
        
        for(g in 1:length(unique(data_list[[k]]$country))){
            data_i <- data_list[[k]] %>% filter(country==unique(data_list[[k]]$country)[g])
            date_start_i <- data_i[1,1]
            date_import_i <- data_i[length(data_i[,1]), 4]
            surv_days_i <- as.numeric(date_import_i - date_start_i)+1
            F_i_vec <- data_i$F_i
            cens_i <- data_i[1,5] 
            
            llk[g] <- (1-cens_i) * 
                      (params[k] + log(F_i_vec[surv_days_i]) + (-sum(exp(params[k]) * F_i_vec[1:surv_days_i]))) + 
                      cens_i * (-sum(exp(params[k]) * F_i_vec[1:surv_days_i])) 
        }
        llk_all[k] <- sum(llk)
    }
    return(sum(llk_all))
}

options(warn=-1)
do.call(rbind, par_list) %>% as.vector() -> initial_list
optim(fn=logL_region, par=c(log(initial_list)), 
      method="BFGS", control = list(fnscale = -1, maxit=1000000)) -> est_region
options(warn=0)

#### likelihood with regional specific scaling factor without depletion effect
logL_region_counter <- function(params){
    llk_all <- rep(0, length(data_list_counter))
    
    for(k in 1:length(data_list_counter)){
        llk <- rep(0, length(unique(data_list_counter[[k]]$country)))
        
        for(g in 1:length(unique(data_list_counter[[k]]$country))){
            data_i <- data_list_counter[[k]] %>% filter(country==unique(data_list_counter[[k]]$country)[g])
            date_start_i <- data_i[1,1]
            date_import_i <- data_i[length(data_i[,1]), 4]
            surv_days_i <- as.numeric(date_import_i - date_start_i)+1
            F_i_vec <- data_i$F_i
            cens_i <- data_i[1,5] 
            
            llk[g] <- (1-cens_i) * 
                      (params[k] + log(F_i_vec[surv_days_i]) + (-sum(exp(params[k]) * F_i_vec[1:surv_days_i]))) + 
                      cens_i * (-sum(exp(params[k]) * F_i_vec[1:surv_days_i])) 
        }
        llk_all[k] <- sum(llk)
    }
    return(sum(llk_all))
}

options(warn=-1)
do.call(rbind, par_list_counter) %>% as.vector() -> initial_list
optim(fn=logL_region_counter, par=c(log(initial_list)), 
      method="BFGS", control = list(fnscale = -1, maxit=1000000)) -> est_region_counter
options(warn=0)

#### likelihood with random effect and global scaling factor
LogL_i_random <- function(data, country_i){
    data_i <- data %>% filter(country==country_i)
    date_start_i <- data_i[1,1]
    date_import_i <- data_i[length(data_i[,1]), 4]
    surv_days_i <- as.numeric(date_import_i - date_start_i)+1
    F_i_vec <- data_i$F_i
    cens_i <- data_i[1,5] 
    
    function(alpha, theta){
        return(
            (1-cens_i) * 
            (alpha + log(F_i_vec[surv_days_i]) +  
            -(theta+1)/theta * log(1+theta*sum(exp(alpha)*F_i_vec[1:surv_days_i]))) +
            cens_i * -1/theta * log(1+theta*sum(exp(alpha)*F_i_vec[1:surv_days_i]))
        )
    }
}

LogL_full_random <- function(data, country_list){
  function(params){
    return(sum(sapply(country_list, 
                      FUN = function(x){LogL_i_random(data=data, country_i = x)(params[1], params[2])})))
  }
}

#### MLE for the model with depletion effect
## with global scaling factor
options(warn=-1)
optim(fn=LogL_full_random(data=df_input_all, country_list=(unique(df_input_all$country))), 
      par=c(log(0.001), 1), method="BFGS", control = list(fnscale = -1, maxit=1000000)) -> est_all_random


#### MLE for the model without depletion effect
## with global scaling factor
optim(fn=LogL_full_random(data=df_input_counter_all, country_list=(unique(df_input_counter_all$country))), 
      par=c(log(0.001), 1), method="BFGS", control = list(fnscale = -1, maxit=1000000)) -> est_all_counter_random

options(warn=0)

saveRDS(est_all, "original_global.rds")
saveRDS(est_all_counter, "original_counter.rds")
saveRDS(est_region, "original_region.rds")
saveRDS(est_region_counter, "original_region.rds_counter")

readRDS("original_global.rds") -> est_all
readRDS("original_counter.rds") -> est_all_counter
readRDS("original_region.rds") -> est_region
readRDS("original_region.rds_counter") -> est_region_counter

#### comparing AIC values
2*1-2*est_all$value
2*1-2*est_all_counter$value
2*6-2*est_region$value
2*6-2*est_region_counter$value

#### comparing BIC values
log(nrow(df_input_all))*1-2*est_all$value
log(nrow(df_input_counter_all))*1-2*est_all_counter$value
log(nrow(df_input_all))*6-2*est_region$value
log(nrow(df_input_counter_all))*6-2*est_region_counter$value

#### Fitting result
Survf_region <- list()
for(k in 1:length(par_list)){
    data_list[[k]] %>% filter(censor==0) -> df_case
    sort(unique(df_case$country)) -> country_list

    Survf_i <- list()
    for(i in 1:length(country_list)){

        data_i <- df_case %>% filter(country==country_list[i])
        date_start <- min(data_i$date)
        date_end <- max(data_i$date)
        date_import_i <- data_i[length(data_i[,1]), 4]
        surv_days_i <- as.numeric(date_import_i - date_start)+1

        surv <- rep(0,(as.numeric(date_end-date_start)+1))
        for(t in 1:(as.numeric(date_end-date_start)+1)){
            surv[t] <- exp(-sum(exp(est_region$par)[k]*data_i$F_i[1:t]))
        }
        surv -> Survf_i[[i]]
    }

    do.call("rbind", Survf_i) -> Survf_all
    cbind(as.data.frame(country_list), Survf_all) %>% rename(country=country_list) -> Survf_all

    df_case %>% filter(date==min(df_case$date)) %>% mutate(date_imp_num=as.numeric(date_imp-date)+1) %>% 
    dplyr::select(country, date_imp_num) -> temp

    merge(temp, Survf_all, by=c("country"), all.y=TRUE) %>% 
    gather(3:170, key='day', value='value') -> Survf_region[[k]]
}

do.call("rbind", Survf_region) -> Survf_region_fig
Survf_region_fig$day <- as.numeric(Survf_region_fig$day)

#### list of countries touching 0
Survf_region_fig %>% filter(day==max(Survf_region_fig$day)) %>% filter(value <= 0.025) %>% dplyr::select(country) -> temp
country_list0 = temp[['country']]

## calculating quantiles
quan_list0 <- list()
for(k in 1:length(country_list0)){
    Survf_region_fig %>% filter(country==country_list0[k]) %>% filter(abs(value-0.975)==min(abs(value-0.975))) %>%
    mutate(quan=c("q975")) -> q_975
    Survf_region_fig %>% filter(country==country_list0[k]) %>% filter(abs(value-0.75)==min(abs(value-0.75))) %>%
    mutate(quan=c("q75")) -> q_75
    Survf_region_fig %>% filter(country==country_list0[k]) %>% filter(abs(value-0.5)==min(abs(value-0.5))) %>%
    mutate(quan=c("q5")) -> q_5
    Survf_region_fig %>% filter(country==country_list0[k]) %>% filter(abs(value-0.25)==min(abs(value-0.25))) %>%
    mutate(quan=c("q25")) -> q_25
    Survf_region_fig %>% filter(country==country_list0[k]) %>% filter(abs(value-0.025)==min(abs(value-0.025))) %>%
    mutate(quan=c("q025")) -> q_025
    rbind(q_975, q_75, q_5, q_25, q_025) -> quan_list0[[k]]
}
do.call(rbind, quan_list0) -> result_quan_list0

result_quan_list0 %>% mutate(date=day+start_date-1) %>% dplyr::select(-c(value, day)) %>% 
spread(key=quan,value=date) -> quan_list0_fig

#### list of countries not toching 0
Survf_region_fig %>% filter(day==max(Survf_region_fig$day)) %>% filter(value > 0.025) %>% dplyr::select(country) -> temp
country_list_non = temp[['country']]

## calculating quantiles
quan_list <- list()
for(g in 1:length(country_list_non)){
    Survf_region_fig %>% filter(country==country_list_non[g]) %>% filter(abs(value-0.975)==min(abs(value-0.975))) %>%
    mutate(quan=c("q975")) -> q_975
    Survf_region_fig %>% filter(country==country_list_non[g]) %>% filter(abs(value-0.75)==min(abs(value-0.75))) %>%
    mutate(quan=c("q75")) -> q_75
    Survf_region_fig %>% filter(country==country_list_non[g]) %>% filter(abs(value-0.5)==min(abs(value-0.5))) %>%
    mutate(quan=c("q5")) -> q_5
    Survf_region_fig %>% filter(country==country_list_non[g]) %>% filter(abs(value-0.25)==min(abs(value-0.25))) %>%
    mutate(quan=c("q25")) -> q_25
    Survf_region_fig %>% filter(country==country_list_non[g]) %>% filter(abs(value-0.025)==min(abs(value-0.025))) %>%
    mutate(quan=c("q025")) -> q_025

    rbind(q_975, q_75, q_5, q_25, q_025) -> temp_quan

    Survf_region_fig %>% filter(country==country_list_non[g]) -> temp

    if(min(temp$value) >= 0.25 && min(temp$value) < 0.5){
        temp_quan %>% mutate(value=case_when(quan==c("q025")~min(temp$value), 
                                             quan==c("q25")~min(temp$value), 
                                             TRUE~value)) -> quan_list[[g]]
    } else if(min(temp$value) >= 0.5 && min(temp$value) < 0.75){
        temp_quan %>% mutate(value=case_when(quan==c("q025")~min(temp$value), 
                                             quan==c("q25")~min(temp$value),
                                             quan==c("q5")~min(temp$value), 
                                             TRUE~value)) -> quan_list[[g]]
    } else if(min(temp$value) >= 0.75 && min(temp$value) < 0.975){
        temp_quan %>% mutate(value=case_when(quan==c("q025")~min(temp$value), 
                                             quan==c("q25")~min(temp$value),
                                             quan==c("q5")~min(temp$value),
                                             quan==c("q75")~min(temp$value), 
                                             TRUE~value)) -> quan_list[[g]]
    } else if(min(temp$value) >= 0.975){
        temp_quan %>% mutate(value=case_when(quan==c("q025")~min(temp$value), 
                                             quan==c("q25")~min(temp$value),
                                             quan==c("q5")~min(temp$value),
                                             quan==c("q75")~min(temp$value), 
                                             quan==c("q975")~min(temp$value), 
                                             TRUE~value)) -> quan_list[[g]]
    } else{
        temp_quan %>% mutate(value=case_when(quan==c("q025")~min(temp$value), TRUE~value)) -> quan_list[[g]]}
}
        
do.call(rbind, quan_list) -> result_quan_list

result_quan_list %>% mutate(date=day+start_date-1) %>% dplyr::select(-c(value, day)) %>% 
spread(key=quan,value=date) -> quan_list_fig

quan_list_all %>% head()

library(lubridate)
bimonthly <- function(x) {
  x_range <- range(x, na.rm = TRUE)
  
  date_range <- c(
    floor_date(x_range[1], "month"),
    ceiling_date(x_range[2], "month")
  )
  monthly <- seq(date_range[1], date_range[2], by = "1 month")
  
  sort(c(monthly, monthly + days(14)))
}

#### figure
rbind(quan_list0_fig, quan_list_fig) -> quan_list_all
Survf_region_fig %>% mutate(date_imp=date_imp_num+start_date-1) %>% 
dplyr::select(country, date_imp) %>% distinct() %>% arrange(date_imp) -> quan_imp_list0

df_input %>% dplyr::select(country, region, sub_region) %>% distinct() -> temp
merge(quan_imp_list0, temp, by=c("country")) %>%
mutate(group=case_when(sub_region %in% c("Central Asia", "Western Asia") ~ c("Middle East"), 
                       TRUE~region)) -> quan_imp_list0

merge(quan_list_all, temp, by=c("country")) %>%
mutate(group=case_when(sub_region %in% c("Central Asia", "Western Asia") ~ c("Middle East"), 
                       TRUE~region)) -> quan_list_all

options(repr.plot.width=10,repr.plot.height=20)
quan_list_all %>% arrange(date_imp_num) %>% 
filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                        "Congo","Congo, Democratic Republic of the","Ghana","United Kingdom"))) %>%
ggplot() +
geom_boxplot(aes(x = reorder(country,date_imp_num), 
                 ymin = q025, lower = q25, middle = q5, upper = q75, ymax = q975, fill=group), 
             stat = "identity", alpha=0.1) +
scale_fill_manual("Region", values = c("#588300", "#1380A1", "#990000", "#FAAB18", "#E75B64FF", "#4C413FFF")) +
geom_errorbar(aes(x = reorder(country,date_imp_num), ymin = q025, ymax = q975, color=group)) +
geom_point(data=quan_imp_list0 %>%
           filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                                   "Congo","Congo, Democratic Republic of the","Ghana","United Kingdom"))), 
           aes(x=reorder(country,date_imp), y=date_imp, color=group), size=3) +
scale_color_manual("Region", values = c("#588300", "#1380A1", "#990000", "#FAAB18", "#E75B64FF", "#4C413FFF")) +
labs(x="Country \n", y="\n Date") +
ggtitle("Hazard of Mpox over time by country") + 
  theme(text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=15, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=18, family="sans",color="black")) +
scale_y_date(date_labels = "%b", breaks="1 month") +
coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

ggsave("../figures/fitting_regional.png", width = 10, height = 20)


#### coverage within the 50% and 95% range
quan_list_all %>% 
filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                        "Congo","Congo, Democratic Republic of the","Ghana","United Kingdom"))) %>%
mutate(date_imp=date_imp_num+start_date-1) %>% 
mutate(count_50=case_when(date_imp <= q25 & date_imp >= q75 ~ 1, TRUE ~ 0),
       count_95=case_when(date_imp <= q025 & date_imp >= q975 ~ 1, TRUE ~ 0)) -> counts

sum(counts$count_50)/nrow(counts); sum(counts$count_50)
sum(counts$count_95)/nrow(counts); sum(counts$count_95)

df_input %>% filter(date >= start_date) %>% arrange(date, country) -> df_input_COI
df_input_COI %>% mutate(h_i=case_when(region == c("Europe") ~ F_i*exp(est_region$par)[1],
                                      region == c("Africa") ~ F_i*exp(est_region$par)[2],
                                      region == c("Americas") ~ F_i*exp(est_region$par)[3],
                                      region == c("Asia") ~ F_i*exp(est_region$par)[4],
                                      sub_region %in% c("Central Asia", "Western Asia") ~ F_i*exp(est_region$par)[5],
                                      region == c("Oceania") ~ F_i*exp(est_region$par)[6])) -> df_input_final
write.csv(df_input_final, "result.csv")

#### time-varying FOI
df_input %>% filter(date >= start_date) %>% arrange(date, country) -> df_input_COI
df_input_COI %>% mutate(h_i=case_when(region == c("Europe") ~ F_i*exp(est_region$par)[1],
                                      region == c("Africa") ~ F_i*exp(est_region$par)[2],
                                      region == c("Americas") ~ F_i*exp(est_region$par)[3],
                                      region == c("Asia") ~ F_i*exp(est_region$par)[4],
                                      sub_region %in% c("Central Asia", "Western Asia") ~ F_i*exp(est_region$par)[5],
                                      region == c("Oceania") ~ F_i*exp(est_region$par)[6]),
                        group=case_when(sub_region %in% c("Central Asia", "Western Asia") ~ c("Middle East"),
                                        TRUE ~ region)) -> df_input_final
write.csv(df_input_final, "result.csv")

raw_global_FoI_data <- read.csv("result.csv")
raw_global_FoI_data$date <- as.Date(raw_global_FoI_data$date)
raw_global_FoI_data$date_import <- as.Date(raw_global_FoI_data$date_import)

global_FoI_data <- raw_global_FoI_data %>% filter(censor==0)

vis_data <- global_FoI_data %>%
filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                        "Congo","Congo, Democratic Republic of the","Ghana","Palestine, State of", 
                        "United Kingdom")))

#### Time-varying FoI figure
options(repr.plot.width=10,repr.plot.height=20)
vis_data %>% arrange(date_import) %>%
ggplot() +
  geom_raster(aes(x=reorder(country, date_import), y = date, fill=h_i)) +
#   scale_fill_gradient(low="white", high="red", name="Force of infection") +
#   scale_fill_gradient2(low="white", mid="red", high="black", space = "Lab", name="Force of infection") +
#   viridis::scale_fill_viridis(option="rocket", direction=-1, name="Force of infection") +
  scale_fill_gradientn(name="Force of infection", colours = c("white","red","purple","darkorchid4","black")) +
  geom_point(aes(x=reorder(country, date_import), y=date_import, color=group), size=3)+
  theme(text = element_text(size=15, family="sans",color="black"),
        axis.text = element_text(size=15, family="sans",color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text(size=18, family="sans",color="black")) +
scale_color_manual("Region", values = c("#588300", "#1380A1", "#990000", "#FAAB18", "#E75B64FF", "#4C413FFF")) +
labs(x="Country \n", y="\n Date") +
ggtitle("Force of infection of Mpox over time by country") + 
scale_y_date(date_labels = "%b", breaks="1 month", expand=c(0,0)) +
coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

ggsave("../figures/time_varying_FoI.png", width = 10, height = 20)

#### checking the linearity of log plot
raw_global_FoI_data %>% filter(censor==0) %>% dplyr::select(country, group, date, date_import, h_i) %>%
arrange(country, date) %>% group_by(country) %>% mutate(cum_h_i=cumsum(h_i)) %>% 
filter(!(country %in% c("Cameroon","Liberia","Central African Republic","Nigeria",
                        "Congo","Congo, Democratic Republic of the","Ghana","Palestine, State of",
                        "United Kingdom"))) %>% ungroup() -> cum_global_FoI_data

## distribution of the CFOI at the date of importation
cum_global_FoI_data %>% mutate(diff=date_import-date) %>% filter(diff==0) -> temp
temp %>% group_by(cum_h_i) %>% summarise(n=n()) -> temp2
temp %>% group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp2$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure

options(repr.plot.width=8,repr.plot.height=6)
temp_figure %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
scale_y_continuous(limits=c(-4.6, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) + 
ggtitle("Model with the depletion effect") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig
fig

ggsave("../figures/global_logplot.png", width = 8, height = 6)

# library(ggrepel)

# options(repr.plot.width=8,repr.plot.height=6)
# temp_figure %>% filter(group==c("Europe")) %>%
# ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
# geom_point() +
# geom_line() +
# geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="sans",color="black"),
#       axis.text = element_text(size=18, family="sans",color="black"),
#      plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(-4.6, 0)) + 
# scale_x_continuous(limits=c(0, 5.6)) +
# ggtitle("Model with the depletion effect (Europe)") + 
# labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_europe


# temp_figure %>% filter(group==c("Africa")) %>%
# ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
# geom_point() +
# geom_line() +
# geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="sans",color="black"),
#       axis.text = element_text(size=18, family="sans",color="black"),
#      plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(-4.6, 0)) + 
# scale_x_continuous(limits=c(0, 5.6)) +
# ggtitle("Model with the depletion effect (Africas)") + 
# labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_africa


# temp_figure %>% filter(group==c("Americas")) %>%
# ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
# geom_point() +
# geom_line() +
# geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="sans",color="black"),
#       axis.text = element_text(size=18, family="sans",color="black"),
#      plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(-4.6, 0)) + 
# scale_x_continuous(limits=c(0, 5.6)) +
# ggtitle("Model with the depletion effect (Americas)") + 
# labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_americas


# temp_figure %>% filter(group==c("Asia")) %>%
# ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
# geom_point() +
# geom_line() +
# geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="sans",color="black"),
#       axis.text = element_text(size=18, family="sans",color="black"),
#      plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(-4.6, 0)) + 
# scale_x_continuous(limits=c(0, 5.6)) +
# ggtitle("Model with the depletion effect (Asia)") + 
# labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_asia


# temp_figure %>% filter(group==c("Middle East")) %>%
# ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
# geom_point() +
# geom_line() +
# geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="sans",color="black"),
#       axis.text = element_text(size=18, family="sans",color="black"),
#      plot.title = element_text(size=22, family="sans",color="black")) +  
# # scale_y_continuous(limits=c(floor_dec(min(temp_figure5$cum_prop_n2), 2), 0)) + 
# scale_y_continuous(limits=c(-4.6, 0)) + 
# scale_x_continuous(limits=c(0, 5.6)) +
# ggtitle("Model with the depletion effect (Middle East)") + 
# labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_middleeast


# temp_figure %>% filter(group==c("Oceania")) %>%
# ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
# geom_point() +
# geom_line() +
# geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
# theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#       axis.title = element_text(size=18, family="sans",color="black"),
#       axis.text = element_text(size=18, family="sans",color="black"),
#      plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(-4.6, 0)) + 
# scale_x_continuous(limits=c(0, 5.6)) +
# ggtitle("Model with the depletion effect (Oceania)") + 
# labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_oceania


# library(ggpubr)
# options(repr.plot.width=24,repr.plot.height=12)
# ggarrange(fig_europe, fig_africa, fig_americas, fig_asia, fig_middleeast, fig_oceania, ncol = 3, nrow = 2,
#           labels = c("A","B","C","D","E","F"), font.label = list(size = 30))

# ggsave("../figures/regional_logplot.png", width = 24, height = 12)

floor_dec <- function(x, sigDig=1) {
  mostSig = ceiling(log10(abs(x)))
  floor(x*10^(sigDig-mostSig))*10^-(sigDig-mostSig)
}

library(ggrepel)
cum_global_FoI_data %>% filter(group==c("Europe")) %>% mutate(diff=date_import-date) %>% filter(diff==0) %>% 
group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp_round$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure1

options(repr.plot.width=8,repr.plot.height=6)
temp_figure1 %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(floor_dec(min(temp_figure1$cum_prop_n2), 2), 0)) + 
scale_y_continuous(limits=c(-0.6, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) +
ggtitle("Model with the depletion effect (Europe)") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_europe


cum_global_FoI_data %>% filter(group==c("Africa")) %>% mutate(diff=date_import-date) %>% filter(diff==0) %>% 
group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp_round$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure2

options(repr.plot.width=8,repr.plot.height=6)
temp_figure2 %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(floor_dec(min(temp_figure2$cum_prop_n2), 2), 0)) + 
scale_y_continuous(limits=c(-0.6, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) +
ggtitle("Model with the depletion effect (Africas)") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_africa


cum_global_FoI_data %>% filter(group==c("Americas")) %>% mutate(diff=date_import-date) %>% filter(diff==0) %>% 
group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp_round$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure3

options(repr.plot.width=8,repr.plot.height=6)
temp_figure3 %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(floor_dec(min(temp_figure3$cum_prop_n2), 2), 0)) + 
scale_y_continuous(limits=c(-0.75, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) +
ggtitle("Model with the depletion effect (Americas)") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_americas


cum_global_FoI_data %>% filter(group==c("Asia")) %>% mutate(diff=date_import-date) %>% filter(diff==0) %>% 
group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp_round$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure4

options(repr.plot.width=8,repr.plot.height=6)
temp_figure4 %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(floor_dec(min(temp_figure4$cum_prop_n2), 2), 0)) + 
scale_y_continuous(limits=c(-0.6, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) +
ggtitle("Model with the depletion effect (Asia)") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_asia


cum_global_FoI_data %>% filter(group==c("Middle East")) %>% mutate(diff=date_import-date) %>% filter(diff==0) %>% 
group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp_round$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure5

options(repr.plot.width=8,repr.plot.height=6)
temp_figure5 %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(floor_dec(min(temp_figure5$cum_prop_n2), 2), 0)) + 
scale_y_continuous(limits=c(-0.6, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) +
ggtitle("Model with the depletion effect (Middle East)") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_middleeast


cum_global_FoI_data %>% filter(group==c("Oceania")) %>% mutate(diff=date_import-date) %>% filter(diff==0) %>% 
group_by(cum_h_i) %>% summarise(n=n(), across()) %>% mutate(prop_n=n/sum(temp_round$n)) %>% 
mutate(cum_prop_n=cumsum(prop_n), cum_prop_n2=log(1-cum_prop_n), log_x=log(cum_h_i)) -> temp_figure6

options(repr.plot.width=8,repr.plot.height=6)
temp_figure6 %>%
ggplot(aes(x=cum_h_i, y=cum_prop_n2, label=country)) +
geom_point() +
geom_line() +
geom_label_repel(box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.title = element_text(size=18, family="sans",color="black"),
      axis.text = element_text(size=18, family="sans",color="black"),
     plot.title = element_text(size=22, family="sans",color="black")) +  
# scale_y_continuous(limits=c(floor_dec(min(temp_figure6$cum_prop_n2), 2), 0)) + 
scale_y_continuous(limits=c(-0.6, 0)) + 
scale_x_continuous(limits=c(0, 5.6)) +
ggtitle("Model with the depletion effect (Oceania)") + 
labs(x="CFOI on the date of importation", y="Log scale of 1-CDF") -> fig_oceania


library(ggpubr)
options(repr.plot.width=24,repr.plot.height=12)
ggarrange(fig_europe, fig_africa, fig_americas, fig_asia, fig_middleeast, fig_oceania, ncol = 3, nrow = 2,
          labels = c("A","B","C","D","E","F"), font.label = list(size = 30))

ggsave("../figures/regional_logplot2.png", width = 24, height = 12)


