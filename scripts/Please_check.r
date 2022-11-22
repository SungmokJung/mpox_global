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

## region-specific data
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

## region-specific data
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

## with region-specific scaling factors (separate MLE by region)
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

## with region-specific scaling factors (separate MLE by region)
par_list_counter <- list(); value_list_counter <- list()
for(i in 1:length(data_list)){
    optim(fn=LogL_full(data=data_list_counter[[i]], country_list=(unique(data_list_counter[[i]]$country))), 
          par=c(0.001), method="Brent", lower=(0), upper=(1000), 
          control = list(fnscale = -1, maxit=1000000)) -> est
    est$par -> par_list_counter[[i]]; est$value -> value_list_counter[[i]]
    
}
options(warn=0)

par_list
par_list_counter

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
            (1-cens_i) * (-(theta+1)/theta * log(1+theta*sum(alpha*F_i_vec[1:surv_days_i])) +
                         -1/theta * log(1+theta*sum(alpha*F_i_vec[1:surv_days_i]))) +
            cens_i * -1/theta * log(1+theta*sum(alpha*F_i_vec[1:surv_days_i]))
        )
    }
}

LogL_full_random <- function(data, country_list){
  function(params){
    return(sum(sapply(country_list, 
                      FUN = function(x){LogL_i_random(data=data, country_i = x)(params[1], params[2])})))
  }
}

# #### MLE for the model the depletion effect
# ## with global scaling factor
# options(warn=-1)
# optim(fn=LogL_full_random(data=df_input_all, country_list=(unique(df_input_all$country))), 
#       par=c(0.001, 1), method="L-BFGS-B", lower=c(0.00005,1e-6), control = list(fnscale = -1, maxit=1000000)) -> est_all_random
# options(warn=0)

# est_all_random

#### MLE for the model the depletion effect
## with global scaling factor
options(warn=-1)
optim(fn=LogL_full_random(data=df_input_all, country_list=(unique(df_input_all$country))), 
      par=c(0.001, 1), method="BFGS", control = list(fnscale = -1, maxit=1000000)) -> est_all_random
options(warn=0)

est_all_random

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
            (1-cens_i) * (-(exp(-theta)+1)/exp(-theta) * log(1+exp(-theta)*sum(exp(-alpha)*F_i_vec[1:surv_days_i])) +
                         -1/exp(-theta) * log(1+exp(-theta)*sum(exp(-alpha)*F_i_vec[1:surv_days_i]))) +
            cens_i * -1/exp(-theta) * log(1+exp(-theta)*sum(exp(-alpha)*F_i_vec[1:surv_days_i]))
        )
    }
}

LogL_full_random <- function(data, country_list){
  function(params){
    return(sum(sapply(country_list, 
                      FUN = function(x){LogL_i_random(data=data, country_i = x)(params[1], params[2])})))
  }
}

options(warn=-1)
optim(fn=LogL_full_random(data=df_input_all, country_list=(unique(df_input_all$country))), 
      par=c(-log(1), -log(1)), method="BFGS", control = list(fnscale = -1, maxit=1000000)) -> est_all_random
options(warn=0)

exp(-est_all_random$par[1])
exp(-est_all_random$par[2])

est_all_random

#### likelihood with random effect and global scaling factor
data=df_input_all
country_i=unique(df_input_all$country)[1]
    data_i <- data %>% filter(country==country_i)
    date_start_i <- data_i[1,1]
    date_import_i <- data_i[length(data_i[,1]), 4]
    surv_days_i <- as.numeric(date_import_i - date_start_i)+1
    F_i_vec <- data_i$F_i
    cens_i <- data_i[1,5] 
    
alpha=exp(-est_all_random$par[1])
theta=exp(-est_all_random$par[2])
        return(
            (1-cens_i) * (-(theta+1)/theta * log(1+theta*sum(alpha*F_i_vec[1:surv_days_i])) +
                         -1/theta * log(1+theta*sum(alpha*F_i_vec[1:surv_days_i]))) +
            cens_i * -1/theta * log(1+theta*sum(alpha*F_i_vec[1:surv_days_i]))
        )

-1/theta
log(1+theta*sum(alpha*F_i_vec[1:surv_days_i]))


