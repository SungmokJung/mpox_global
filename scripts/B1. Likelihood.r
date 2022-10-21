library(tidyverse); library(dplyr)

(censor_import <- max(df_input$date))

read.csv("../data/input.csv") -> df_input
df_input$date <- as.Date(df_input$date)
df_input$date_import <- as.Date(df_input$date_import)

df_input %>% rename(censor=censoring, country=location) %>%
mutate(date_imp=case_when(censor==1~censor_import, censor==0~date_import)) %>%
dplyr::select(date, country, F_i, date_imp, censor) -> df_input

df_input %>% head()

###2. Likelihood -----
## likelihood for country i 
LogL_i <- function(data, country_i){
  #data[1]=calender date; data[2]=country; data[3]=F_i; data[4]=date of importation; data[5]=dummy variable for censoring (0 observed, 1 censored)
  data_i <- data %>% filter(country==country_i)
  date_start_i <- data_i[1,1]
  date_import_i <- data_i[length(data_i[,1]), 4]
  surv_days_i <- as.numeric(date_import_i - date_start_i) + 1
  F_i_vec <- data_i$F_i
  cens_i <- data_i[1,5] #censored=1, observed=0
  #log-likelihood for country i 
  function(alpha){
    return(
      (1-cens_i) * (alpha * F_i_vec[surv_days_i]) + exp(-sum(alpha * F_i_vec[1:surv_days_i])) + #if observed, cens_i=0, the contribution to logL is log(h(t) * S(t))
        cens_i * exp(-sum(alpha * F_i_vec[1:surv_days_i])) #if censored, cens_i=1, the contribution to logL is log(S(t))
      )
  }
}

## full likelihood
LogL_full <- function(data, country_list){
  function(alpha){
    return(sum(sapply(country_list, FUN = function(x){LogL_i(data=data, country_i = x)(alpha=alpha)})))
  }
}

## MLE
options(warn=-1)
optim(fn=LogL_full(data=df_input, country_list=(unique(df_input$country))), 
      par=c(0.01), method="Brent", lower=(0), upper=(10000),
      control = list(fnscale = -1, maxit=1000000))
options(warn=0)


