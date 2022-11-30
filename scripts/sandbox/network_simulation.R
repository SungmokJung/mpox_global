### set weibull parameters
alpha = 0.1; kappa=0.77; 
theta=(14/365)*((alpha/kappa)^(1/alpha));

dweibull_lt = function (x,shape,scale,u){
  if(min(x)<u){
    return(0)
    #stop('The value is less than left truncated point')
  }
  (dweibull(x,shape,scale)/(1-pweibull(u,shape,scale)))*heaviside(x,u)
}

aa <- rep(0,1000*365)
width <- 1/365
for(t in 1:(1000*365)){
  x <- t*width
  aa[t] <- x*dweibull_lt(x,alpha,theta,14/365)
}
mu1 <- sum(aa)*width ### weibull mean

width <- 1/365
N <- 1000000
upperl = 1001
S = rep(0,(upperl-1)/width)
for(t in 1:((upperl-1)/width)){
  x <- width*t 
  S[t] = width*dweibull_lt(x,alpha,theta,14/365)
}
wei_sum = sum(S)


### simulation code
Simulation_fn = function(n){
  upperl =1001
  fx = rep(0,(upperl-1)/width)
  sx = rep(0,(upperl-1)/width)
  meanx = rep(0,(upperl-1)/width)
  for(t in 1:((upperl-1)/width)){
    x <- width*t 
    fx[t] = width*x*max(0,x-1)*dweibull_lt(x,alpha,theta,14/365)*exp(-(x*n)/(mu1*N))
    sx[t] = width*dweibull_lt(x,alpha,theta,14/365)*exp(-(x*n)/(mu1*N))
    meanx[t] = width*x*dweibull_lt(x,alpha,theta,14/365)*exp(-(x*n)/(mu1*N))
  }
  beta = mu1/sum(fx[1:((upperl-1)/width)])   ### SAR
  m = wei_sum-sum(sx[1:((upperl-1)/width)])   ### CIP
  Reff1 = 0.10*sum(fx[1:((upperl-1)/width)])/mu1   ### Reff under SAR of 0.1
  Reff2 = 0.20*sum(fx[1:((upperl-1)/width)])/mu1   ### Reff under SAR of 0.2
  excess= sum(fx[1:((upperl-1)/width)])/sum(meanx[1:((upperl-1)/width)])   ### mean excess degree of susceptibles
  df = c(beta,m,Reff1,Reff2,excess) 
  return(df)
}

df <- mapply(Simulation_fn,seq(0,100*800,100))
result_df <- df %>% t() %>% as.data.frame()
colnames(result_df) <- c("SAR","Infections","Reff_1","Reff_2","Excess")
result_df <- cbind(result_df,seq(0,100*800,100))
colnames(result_df) <- c("SAR","Infections","Reff_1","Reff_2","Excess","n")
write.csv(result_df, "data/SAR_cip_Reff_excess.csv")