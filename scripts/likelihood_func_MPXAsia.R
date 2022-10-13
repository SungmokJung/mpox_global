############################################################
#Aim: to construct a likelihood function for survival analysis
#Final edit: 13 October 2022
#Editor: Fumi Miura
############################################################
###Procedure
#0. Package 
#1. Test data
#2. Likelihood
############################################################

###0. Package -----
library(tidyverse)

###1. Test data -----
#hazard is set as constant: alpha_true
#survival prob: round(sapply(0:19, function(x){exp(-alpha_true*x)}),3)
alpha_true <- log(2)
no_of_country <- 20
no_of_days <- 20
#construct test data
test_data <- data.frame( #for the main analysis, we replace here by (1) read_csv("df_inci_final.csv") #"2022-05-01", "2022-10-03" and (2) select variables
  date = rep(seq(as.Date("2022-05-01"),as.Date("2022-05-01")+no_of_days-1,by=1), no_of_country),
  country = as.character(sapply(1:no_of_country,function(x){rep(x,no_of_days)})),
  F_i = rep(rep(alpha_true,no_of_days),no_of_country),
  date_imp = c(rep(rep(as.Date("2022-05-02"),no_of_days),10),
               rep(rep(as.Date("2022-05-03"),no_of_days),7),
               rep(rep(as.Date("2022-05-04"),no_of_days),2),
               rep(rep(as.Date("2022-05-20"),no_of_days),1)
  ), #note: I determined the numbers of countries that imported its first case by eye-judging
  censor = c(rep(rep(0, no_of_days),no_of_country-1), rep(1, no_of_days)) #censored=1, observed=0
)

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
      (1-cens_i) * (log(alpha * F_i_vec[surv_days_i]) + log(exp(-sum(alpha * F_i_vec[1:surv_days_i])))) + #if observed, cens_i=0, the contribution to logL is log(h(t) * S(t))
        cens_i * log(exp(-sum(alpha * F_i_vec[1:surv_days_i]))) #if censored, cens_i=1, the contribution to logL is log(S(t))
      )
  }
}

## full likelihood
LogL_full <- function(data, country_list){
  function(alpha){
    return(sum(sapply(country_list, FUN = function(x){LogL_i(data=data, country_i = x)(alpha=alpha)})))
  }
}
#check: LogL_full(data=test_data,country_list=(unique(test_data$country)))(alpha=log(1.1)) #-56.18391

## MLE
optim(fn=LogL_full(data=test_data, country_list=(unique(test_data$country))), 
      par=c(7), 
      control = list(fnscale = -1), 
      hessian = T) #estimated para = 0.397168
