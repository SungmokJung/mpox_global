###Survival analysis

### import data
alpha_true <- 0.2
data_global <- data.frame( #read_csv("df_inci_final.csv")
  date = seq(as.Date("2022-01-05"),as.Date("2022-01-19"),by=1),
  country = rep("Japan",15),
  F_i = exp(-alpha_true*(1:15))+rnorm(1,mean=0,sd=0.01)
)

### functions
F_i <- function(data, country_i, date_t){
  tmp_data <- data %>% filter(country==country_i) %>% filter(date==date_t)
  return(tmp_data$F_i)
}
#F_i(data=data_global,country="Japan",date_t="2022-01-07") #0.5488116

hazard_func <- function(alpha, data, country_i, date_t){
  return(alpha * F_i(data, country_i, date_t))
}
#hazard_func(alpha=0.1, data=data_global,country="Japan",date_t="2022-01-07") #0.05488116

survival_func <- function(alpha, data, country_i, date_t, date_0){ #date_0 = as.Date("2022-05-01")
  date_0 <- as.Date(date_0)
  date_t <- as.Date(date_t)
  int_date <- seq(date_0, date_t, by=1)
  #integration
  tmp_int <- sapply(int_date, FUN = function(x){hazard_func(alpha, data, country_i, date_t=x)})
  return(sum(tmp_int))
}
#survival_func(alpha=0.1, data=data_global,country_i="Japan",date_t="2022-01-10",date_0="2022-01-10") #0.1161125

LogL_i <- function(data, country_i, date_t, date_0){
  date_0 <- as.Date(date_0)
  date_t <- as.Date(date_t)
  int_date <- seq(date_0, date_t, by=1)
  function(alpha){
    tmp_int <- sapply(int_date, FUN = function(x){
      ###there must be mistake around here###
      log(hazard_func(alpha, data, country_i, date_t=x)) + log(survival_func(alpha, data, country_i, date_t=x, date_0))
      })
   return(sum(tmp_int)) 
  }
}
LogL_i(data=data_global,country_i="Japan",date_t="2022-01-12",date_0="2022-01-10")(alpha=10) #-20.10547

optim(fn= LogL_i(data=data_global,country_i="Japan",date_t="2022-01-15",date_0="2022-01-10"), par=c(0.4), control = list(fnscale = -1), hessian = T)