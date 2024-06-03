
setwd("/Users/peiqinhe/Desktop/Analysis Skills/R/")
load("Risk Analysis.RData", verbose=TRUE)

# The following packages and functions will be used for this assignment:
library(data.table)
library(ggplot2)
library(MASS)
library(moments)
library(parallel)
library(quantmod)
library(Quandl)
library(rugarch)
library(xts)
library(nloptr)

CalcVaRES <- function(r,alpha) {
  VaR <- quantile(r,1-alpha)
  ES  <- mean(r[r<VaR])
  VaR_ES <- c(VaR,ES)
  names(VaR_ES) <- c("VaR","ES")
  return(VaR_ES)
}
DATE <- function(yyyy,mm,dd) {
  dte  <- as.Date(sprintf("%i-%i-%i",yyyy,mm,dd),format="%Y-%m-%d")
  return(dte)
}
as.Date2 <- function(x) {
  tryfmt <- c("%Y-%m-%d","%m/%d/%Y","%Y/%m/%d","%b %d,%Y")
  return(as.Date(x,tryFormats=tryfmt))
}

# Option functions
# F is the price of the underlying asset
# X is the strike price of the option
# t is the time to maturity (in years)
# r is the tbill rate (in decimal form)
# sigma is the volatility of the underlying asset
# BFC is the price of the call option
# BFP is the price of the put option
# IVC is the implied volatility of the call option
# IVP is the implied volatility of the put option

BFC <- function(F,X,t,r,sigma) {
  d1 <- log(F/X) + 0.5*sigma^2 *t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  N1 <- pnorm(d1)
  N2 <- pnorm(d2)
  C <- exp(-r*t) * (F * N1 - X * N2 )
  return(C)
}

BFP <- function(F,X,t,r,sigma) {
  d1 <- log(F/X) + 0.5*sigma^2 *t
  d1 <- d1/(sigma*sqrt(t))
  d2 <- d1 - sigma*sqrt(t)
  NM1 <- pnorm(-d1)
  NM2 <- pnorm(-d2)
  P <- exp(-r*t) * (X * NM2 - F * NM1 )
  return(P)
}

IVC <- function(F,X,t,r,Call) {
  eval_f_C <- function(sigma) {
    return ( (Call-BFC(F,X,t,r,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_C,
                opts=opts)
  return(es$solution)
}
IVP <- function(F,X,t,r,Put) {
  eval_f_P <- function(sigma) {
    return ( (Put-BFP(F,X,t,r,sigma))^2 ) 
  }
  opts <- list("algorithm"="NLOPT_LN_COBYLA",
               "xtol_rel"=1.0e-8)
  xs <- 0.10
  es <- nloptr( x0=xs,
                eval_f=eval_f_P,
                opts=opts)
  return(es$solution)
}
#----------------------------------------------------

# Part 1 – R code

# Question 1

# The most appropriate statistical model to characterize these returns should be the GARCH~t
# Model Justification are shown below:

# logret_X Analysis
# The skewness of logret_X is:
round(skewness(logret_X.xts),2)
# -3.68 suggests logret_X is left-skewed (longer left tail), this is evidence against normality.
# The kurtosis of logret_X is:
round(kurtosis(logret_X.xts),2)
# 57.57 suggests logret_X has fat tails, this is evidence against normality.
# The JB statistic of logret_X is:
jb_logretX <- jarque.test(as.vector(logret_X.xts))
jb_logretX
round(jb_logretX$statistic,2)
# 159057.8 and p-value is almost 0, 
# meaning there is virtually zero probability of observing this JB static if the data were normally distributed,
# this is a strong rejection of normality.
# Run the acf() function to calculate the autocorrelation coefficients and save the output in acf.out:
acf.out_logretX <- acf(logret_X.xts,main="ACF of daily log returns")
acf.out_logretX
# The first order autocorrelation coefficient of logret_X is
rho1_logretX <- acf.out_logretX$acf[2]
round(rho1_logretX,3)
# -0.062, there is no evidence of autocorrelation.
# Run the acf function on |z|:
set.seed(123789)
rvec_abslogretX <- sample(abs(logret_X.xts),length(abs(logret_X.xts)),replace=FALSE)
rvec.acf.out <- acf(rvec_abslogretX)
rvec.acf.out
# The first order autocorrelation coefficient is
rvec.rho1_abslogretX <- rvec.acf.out$acf[2]
round(rvec.rho1_abslogretX,3)
# 0.483

# IID ~ Normal Distribution
IIDn <- fitdistr(logret_X.xts,"normal")
IIDn_mu <- IIDn$estimate[1]
IIDn_sig <- IIDn$estimate[2]
# The estimated mean of the normal distribution is
round(IIDn_mu,6)
# 0.000156
# The estimated standard deviation of the normal distribution is
round(IIDn_sig,6)
# 0.00647
# The AIC (Akaike Information Criterion) of this model is:
AIC_IIDn <- 2*length(IIDn$estimate)-2*IIDn$loglik
round(AIC_IIDn,0)
# -9115

# IID ~ t distribution
IIDt <- fitdistr(logret_X.xts,"t")
IIDt_mu <- IIDt$estimate[1]
IIDt_sig <- IIDt$estimate[2]
IIDt_tdf <- IIDt$estimate[3]
# The estimated mean of the t distribution is
round(IIDt_mu,6)
# 0.000663
#The estimated standard deviation of the t distribution is
round(IIDt_sig,6)
# 0.002667
# The estimated degree of freedom of the t distribution is
round(IIDt_tdf,6)
# 1.975892
# The AIC (Akaike Information Criterion) of this model is:
AIC_IIDt <- 2*length(IIDt$estimate)-2*IIDt$loglik
round(AIC_IIDt,0)
# -10093

# ARMA(1,1)
arma11 <- arima(logret_X.xts,c(1,0,1))
arma11
# The AIC (Akaike Information Criterion) of this model is:
AIC_arma11 <- arma11$aic
round(AIC_arma11,0)
# -9124

# GARCH~n model
uspec_N <- ugarchspec(variance.model = list(model = "sGARCH",garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0),include.mean = TRUE),         
                      distribution.model = "norm")
garch_N <- ugarchfit(spec = uspec_N, data = logret_X.xts, solver="hybrid")
round(garch_N@fit$coef,6)
# Save the standardized residuals in vector z:
z_garch_N <- garch_N@fit$z
# The mean of z is:
round(mean(z_garch_N),4)
# -0.0891, the mean of z is close to zero.  
# The standard deviation of z is:
round(sd(z_garch_N),4)
# 1.0067, the standard deviation of z is close to 1.
# The skewness of z is:
round(skewness(z_garch_N),2)
# -1.88, the skewness of z is slightly negative.  
# The kurtosis of z is:
round(kurtosis(z_garch_N),2)
# 11.86, the kurtosis of z is a lot bigger than 3.
# The JB statistic is:
jb_garch_N <- jarque.test(as.vector(z_garch_N))
jb_garch_N
round(jb_garch_N$statistic,2) 
# 4861.83 with p-value almost 0, there is evidence of non-normality.
# Run the acf function on z:
acf.out_garch_N_z <- acf(z_garch_N,main="acf of standardized residuals")
acf.out_garch_N_z
# The first order autocorrelation coefficient of z is
rho1_garch_N_z <- acf.out_garch_N_z$acf[2]
round(rho1_garch_N_z,3)
# 0.011, there is no evidence of autocorrelation.
# Run the acf function on |z|:
acf.out_garch_N_absz <- acf(abs(z_garch_N),main="acf of |standardized residuals|")
acf.out_garch_N_absz
# The first order autocorrelation coefficient of |z| is
rho1_garch_N_absz <- acf.out_garch_N_absz$acf[2]
round(rho1_garch_N_absz,3)
# 0.044, there is no evidence of volatility clustering.
# The AIC of the estimated model:
AIC_garchN <- 2 * length(garch_N@fit$coef) - 2 * garch_N@fit$LLH
round(AIC_garchN,0)
# -10032
# These diagnostics indicate that the AR(1)-GARCH(1,1)~n model captures all of the volatility clustering in the data. 
# But it is only able to capture some, but not all, of the non-normality in our data.

# GARCH~t model
uspec_t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
garch_t <- ugarchfit(spec = uspec_t, data = logret_X.xts, solver="hybrid")
round(garch_t@fit$coef,6)
# Save the standardized residuals in vector z:
z_garch_t <- garch_t@fit$z
# The mean of z is:
round(mean(z_garch_t),4)
# -0.1295, the mean of z is close to zero.
# The standard deviation of z is:
round(sd(z_garch_t),4)
# 1.0224, the standard deviation of z is close to 1.
# The skewness of z is:
round(skewness(z_garch_t),2)
# -1.9
# The kurtosis of z is:
round(kurtosis(z_garch_t),2)
# 12.01
# The JB statistic is:
jb_garch_t <- jarque.test(as.vector(z_garch_t))
jb_garch_t
round(jb_garch_t$statistic,2) 
# 5014.29
# Run the acf function on z:
acf.out_garch_t_z <- acf(z_garch_t,main="acf of standardized residuals")
acf.out_garch_t_z
# The first order autocorrelation coefficient of z is
rho1_garch_t_z <- acf.out_garch_t_z$acf[2]
round(rho1_garch_t_z,3)
# 0.018, there is no evidence of autocorrelation.
# Run the acf function on |z|:
acf.out_garch_t_absz <- acf(abs(z_garch_t),main="acf of |standardized residuals|")
acf.out_garch_t_absz
# The first order autocorrelation coefficient of |z| is
rho1_garch_t_absz <- acf.out_garch_t_absz$acf[2]
round(rho1_garch_t_absz,3)
# 0.055, there is no evidence of volatility clustering.
# The AIC of the estimated model:
AIC_garcht <- 2 * length(garch_t@fit$coef) - 2 * garch_t@fit$LLH
round(AIC_garcht,0)
# -10408

# Compare the models:
model_names <- c("IID~n", "IID~t", "ARMA(1,1)", "AR(1)-GARCH(1,1)~n", "AR(1)-GARCH(1,1)~t")
AIC_values <- c(AIC_IIDn, AIC_IIDt, AIC_arma11, AIC_garchN, AIC_garcht)
AIC_values_whole <- round(AIC_values,0)
model_AIC <- data.frame(Model = model_names, AIC = AIC_values_whole)
model_AIC
# Based on the above table, GARCH~t has the lowest AIC, thus this model is the most appropriate one

# five (trading) days simulation using GARCH~t
detectCores(logical=FALSE)
ncore <- max(1,floor(detectCores(logical=FALSE)/2))
ncore
cl <- makeCluster(ncore)
set.seed(123789)        # for reproducibility
nper <- 5              # for 5-period ahead simulation
nsim <- 100000          # number of simulations
boot_garch <- ugarchboot(fitORspec=garch_t,
                         method="Partial",
                         sampling="raw",
                         n.ahead=nper,
                         cluster=cl,
                         n.bootpred=nsim)
simmat <- boot_garch@fseries
str(simmat)
sim <- apply(simmat,1,sum)
str(sim)
alpha <- 0.95
VaR_ES <- CalcVaRES(sim,alpha)
round(VaR_ES,6)
# The value-at-risk (VaR) is:
round(VaR_ES[1],4) 
# -0.0469
# The expected shortfall (ES) is
round(VaR_ES[2],4)
# -0.0851
stopCluster(cl)

#----------------------------------------------------

# Question 2
holdings_X[, t := as.numeric(Expiration-Date)/365]
call_data<- holdings_X[`Call/Put`==1]
put_data <- holdings_X[`Call/Put`==2]

# find implied volatility of Call option
IVC_values <- numeric(nrow(call_data))
for (i in 1:nrow(call_data)) {
  row <- call_data[i, ]
  IVC_values[i] <- IVC(F = row$FutPrice,
                       X = row$Strike,
                       t = row$t,
                       r = row$r,
                       Call = row$OptPrice)
}
call_data$IV <- IVC_values

# find implied volatility of Put option
IVP_values <- numeric(nrow(put_data))
for (i in 1:nrow(put_data)) {
  row <- put_data[i, ]
  IVP_values[i] <- IVP(F = row$FutPrice,
                       X = row$Strike,
                       t = row$t,
                       r = row$r,
                       Put = row$OptPrice)
}
put_data$IV <- IVP_values

# Create a scatter plot of implied volatility against strike price
ggplot(call_data, aes(x = Strike, y = IV)) +
  geom_point() +
  labs(x = "Strike Price", y = "Implied Volatility") +
  ggtitle("Implied Volatility vs. Strike Price")
# The call options do not exhibit volatility skew. The trend is not clear.

ggplot(put_data, aes(x = Strike, y = IV)) +
  geom_point() +
  labs(x = "Strike Price", y = "Implied Volatility") +
  ggtitle("Implied Volatility vs. Strike Price")
# The implied volatility tends to decrease as the strike prices increases for the put options (downward-sloping).

#----------------------------------------------------

# Question 3

# Perform the following scenario analysis. Suppose the implied volatility of all options increase by 50%, 100%, 
# 200%, or 300% on Nov 1, 2017, what happens to the value of the fund’s investments, assuming everything else 
# remains unchanged?
data <- rbind(call_data, put_data)
data[, market.value := Contracts * OptPrice * 250]
cash_position = 628226078
calculate_total_fund_value = function(multiplier, cash_position) {
  for(i in 1:nrow(data)) {
    row_data = data[i]
    # Calculate new option prices
    if(row_data$`Call/Put` == 1) {
      # Call option
      data[i, new_price := BFC(row_data$FutPrice, row_data$Strike, row_data$t, row_data$r, row_data$IV*multiplier)]
    } else {
      # Put option
      data[i, new_price := BFP(row_data$FutPrice, row_data$Strike, row_data$t, row_data$r, row_data$IV*multiplier)]
    }
  }
  # we have full contracts, which means each contract contains 250 options
  data[, new_value := 250 * Contracts * new_price ]
  # Exclude the redemption fee and the annual fund operation fees adapted to two-month period
  total_options_old = sum(data$market_value) * (1-0.01 - 0.0324*(2/12))
  total_options_new <- sum(data$new_value) * (1- 0.01 - 0.0324*(2/12))
  
  return(total_options_new - total_options_old + cash_position)
}

## Scenario 1, increase by 50%
scenario_1 = calculate_total_fund_value(1.5, cash_position)

## Scenario 2, increase by 100%
scenario_2 = calculate_total_fund_value(2, cash_position)

## Scenario 3, increase by 200%
scenario_3 = calculate_total_fund_value(3, cash_position)

## Scenario 4, increase by 300%
scenario_4 = calculate_total_fund_value(4, cash_position)


scenarios = c("Increase by 50%", "Increase by 100%", "Increase by 200%", "Increase by 300%")
values = c(scenario_1, scenario_2, scenario_3, scenario_4)

scenario_analysis <- data.frame(Scenarios = scenarios, Total_Fund_Value = values)
scenario_analysis

# How high will the implied volatilities have to rise to reduce the value of the fund (cash + market value of options) 
# by 50%?
totalfv = cash_position + sum(data$market_value)
target_value = totalfv / 2

iv_multiplier <- 1
target_iv_multiplier <- NA
results <- data.frame(multiplier = numeric(), fund_value = numeric())

while (iv_multiplier <= 4) {
  fund_value <- calculate_total_fund_value(iv_multiplier, cash_position)
  
  if (is.na(target_iv_multiplier) && fund_value <= target_value) {
    target_iv_multiplier <- iv_multiplier
  }
  results <- rbind(results, data.frame(multiplier = iv_multiplier, fund_value = fund_value))
  iv_multiplier <- iv_multiplier + 0.05
}

IV_increase <- paste0((target_iv_multiplier -1) *100,"%")
IV_increase

# What other factor(s), besides the implied volatilities of options, would you take into
# consideration for your scenario analysis?

# In addition to implied volatilities of options, several other factors should be considered for scenario analysis:
# Underlying Asset Price Movement: Changes in the S&P 500 index can significantly impact option prices. Analyzing different potential movements in the underlying asset's price is essential.
# Interest Rates (particularly for options with longer expiration periods): Higher interest rates can increase option prices, while lower rates may decrease them.
# Market Sentiment: Bullish or bearish market sentiments may affect the demand for options and subsequently impact their prices.
# Market Liquidity: Liquidity in the options market can affect bid-ask spreads and execution prices. Lower liquidity can result in wider spreads and higher trading costs.
# Time to Expiration: Options with more time to expiration generally have higher premiums compared to options with less time.
# Dividends (especially for American-style options): Ex-dividend dates and dividend amounts should be considered in scenario analysis.

#----------------------------------------------------

# Part 2 – R code for Question 2

# Question 4
str(logret_G.xts)
tail(logret_G.xts)
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH",garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,0),include.mean = TRUE),    
                      distribution.model = "std")
garch_t <- ugarchfit(spec = uspec.t, data = logret_G.xts, solver='hybrid')
round(garch_t@fit$coef,6)

start_date <- as.Date("2020-12-31")
end_date <- as.Date("2021-06-30")


n2020 <- sum(index(logret_G.xts) <= as.Date("2020-12-31"))

roll.garch <- ugarchroll(spec=uspec.t,
                         data=logret_G.xts, 
                         n.ahead=1,
                         n.start=n2020,  
                         refit.every=1,
                         refit.window="recursive",
                         calculate.VaR=TRUE,
                         VaR.alpha=0.05)
dfvar <- roll.garch@forecast$VaR
names(dfvar) <- c("VaR","actual")
str(dfvar)

position_value <- 1000000
dfvar$MarginRequirement <- position_value * (exp(dfvar$VaR) - 1)
dfvar

#----------------------------------------------------

#Question 5
n_start<-n2020
forecast_dates <- index(logret_G.xts)[(n_start + 1):length(logret_G.xts)]
forecast_dates <- forecast_dates[1:nrow(dfvar)]
dfvar$date <- as.Date(forecast_dates)
dfvar
library(ggplot2)
dfvar$date <- as.Date(dfvar$date)

#Checking for NA values
sum(is.na(dfvar$VaR))     
sum(is.na(dfvar$actual))   
sum(is.na(dfvar$date))     


#Plotting the graphs
cols <- c("VaR" = "red", "Actual" = "blue")

ggplot(dfvar, aes(x = date)) +
  geom_line(aes(y = VaR, color = "VaR"), linewidth = 1) +  
  geom_line(aes(y = actual, color = "Actual"), linewidth = 1) + 
  scale_color_manual(values = cols) + 
  labs(title = "1-day 95% VaR from Dec 31, 2020 to Jun 30, 2021",
       x = "",
       y = "") +
  theme(legend.position = "bottom")


#How many days did the actual return fall below the VaR at the 95% confidence level?
finalans<-sum(dfvar$actual<dfvar$VaR)
print(finalans)
sprintf("The actual return fell below the VaR at the 95%% confidence level on %d days for the given time period.", finalans)
#----------------------------------------------------
