rm(list = ls())

############  VECTOR AUTO-REGRESSION MODEL (VAR MODEL) ###################

# LIBRARIES 

library(MASS)
library(strucchange)
library(zoo)
library(sandwich)
library(lmtest)
library(urca)
library(vars)
library(mFilter)
library(tseries)
library(forecast)
library(tidyverse)
library(ggplot2)
library(TSstudio)

# Loading the data set: VAR_geometric_data.xlsx


#Visualize data 


d <- ggplot(data = VAR_geometric_data)
d + geom_smooth(method = lm, mapping = aes(x = oil_yoy, y = inflation_yoy)) + geom_point(mapping = aes(x = oil_yoy, y = inflation_yoy))
d + geom_smooth(method = lm, mapping = aes(x = gas_yoy, y = inflation_yoy)) + geom_point(mapping = aes(x = gas_yoy, y = inflation_yoy))
d + geom_smooth(method = lm, mapping = aes(x = local_price_yoy, y = inflation_yoy)) + geom_point(mapping = aes(x = local_price_yoy, y = inflation_yoy))

graph <- cbind(cpi_geo_y,loc_geo_y,gas_geo_y,oil_geo_y)
colnames(graph) <- c("Inflacion", "Combustibles", "Gasolina", "Petr?leo")
plot.ts(graph)

############  DECLARE VARIABLES to TIME SERIES FOR THE VAR  ###################  

loc_geo_y <- ts(VAR_geometric_data$local_price_yoy, start = c(2001,2), frequency = 12)
gas_geo_y <- ts(VAR_geometric_data$gas_yoy, start = c(2001,2), frequency = 12)
oil_geo_y <- ts(VAR_geometric_data$oil_yoy, start = c(2001,2), frequency = 12)
cpi_geo_y <- ts(VAR_geometric_data$inflation_yoy, start = c(2001,2), frequency = 12)

#plot the series 

ts_plot(cpi_geo_y)
ts_plot(loc_geo_y)
ts_plot(gas_geo_y)
ts_plot(oil_geo_y)

autoplot(cbind(loc_geo_y, gas_geo_y, oil_geo_y, cpi_geo_y))

############  INDIVIDUAL ORDINARY LEAST SQUARES ################################

OLS1 <- lm(cpi_geo_y ~ loc_geo_y)
OLS2 <- lm(cpi_geo_y ~ gas_geo_y)
OLS3 <- lm(cpi_geo_y ~ oil_geo_y)
OLS4 <- lm(oil_geo_y ~ gas_geo_y)


plot(OLS1)
par(mfrow=c(2,2))

abline(OLS1)



resume_lm <- summary.lm(OLS1)
print(resume_lm, digits = max(3, getOption("digists") - 3))


resume_lm$coefficients

anova(OLS1)
summary(OLS1)



############  DETERMINE THE PERSISTANCE OF THE MODEL  ###########################


#Autocorrelation Function (ACF) and Partial Autocorrelation Function (PACF)
#to figure out the order of an AR model, you need to look at the PACF
#to figure out the order of an MA model, you need to look at the ACF

par(mar = c(1, 1, 1, 1))

acf(loc_geo_y, main = "ACF para Precios Locales de Combustibles")
acf(gas_geo_y, main = "ACF para Precios Externos de Gasolina")
acf(oil_geo_y, main = "ACF para Precios Externos de Petroleo")
acf(cpi_geo_y, main = "ACF para inflacion")

pacf(loc_geo_y, main = "PACF para Precios Locales de Combustibles")
pacf(gas_geo_y, main = "PACF para Precios Externos de Gasolina")
pacf(oil_geo_y, main = "PACF para Precios Externos de Petroleo")
pacf(cpi_geo_y, main = "PACF para inflacion")

#Augmented Dickey Fuller for STATIONARITY

adf_1 <- adf.test(loc_geo_y)
adf_2 <- adf.test(gas_geo_y)
adf_3 <- adf.test(oil_geo_y)
adf_4 <- adf.test(cpi_geo_y)

dickey <- c(adf_1$p.value, adf_2$p.value, adf_3$p.value, adf_4$p.value)
dickey

#Finding Optimal Lags and Group  time series together

arreglo1 <- cbind(cpi_geo_y, oil_geo_y) 
colnames(arreglo1) <- cbind("Inflacion","Petroleo") 

arreglo2 <- cbind(cpi_geo_y, oil_geo_y, gas_geo_y) 
colnames(arreglo2) <- cbind("Inflacion", "Petroleo", "Gasolina") 

arreglo3 <- cbind(cpi_geo_y, oil_geo_y, loc_geo_y) 
colnames(arreglo3) <- cbind("Inflacion", "Petroleo", "Combustibles") 

arreglo4 <- cbind(cpi_geo_y, oil_geo_y, gas_geo_y, loc_geo_y) 
colnames(arreglo4) <- cbind("Inflacion", "Petroleo", "Gasolina", "Combustibles") 

arreglo5 <- cbind(cpi_geo_y, oil_geo_y,  loc_geo_y, gas_geo_y) 
colnames(arreglo5) <- cbind("Inflacion", "Petroleo", "Combustibles", "Gasolina") 

#Results Akaike Hannan-Quin Swartz

lagselect1 <- VARselect(arreglo1, lag.max = 10, type = "const")
lagselect2 <- VARselect(arreglo2, lag.max = 10, type = "const")
lagselect3 <- VARselect(arreglo3, lag.max = 10, type = "const")
lagselect4 <- VARselect(arreglo4, lag.max = 10, type = "const")
lagselect5 <- VARselect(arreglo5, lag.max = 10, type = "const")

lagselection <- c(lagselect1$selection[1], lagselect2$selection[1],
                  lagselect3$selection[1], lagselect4$selection[1],
                  lagselect5$selection[1])
lagselection

############  REDUCED MODEL ESTIMATION #########################################

#Reduced form of the model,a priori knowledge of the variables included
# 'cuz structural can?t be estimated directly due effects of the contemporaneous variables

reduced1 <- VAR(arreglo1, p = 7, type = "const", season = NULL, exogen = NULL)
reduced2 <- VAR(arreglo2, p = 6, type = "const", season = NULL, exogen = NULL)
reduced3 <- VAR(arreglo3, p = 2, type = "const", season = NULL, exogen = NULL)
reduced4 <- VAR(arreglo4, p = 2, type = "const", season = NULL, exogen = NULL)
reduced5 <- VAR(arreglo5, p = 2, type = "const", season = NULL, exogen = NULL)


summary(reduced1)
summary(reduced2)
summary(reduced3)
summary(reduced4)
summary(reduced5)

############  TESTING THE reduced MODEL ROBUSTESS ##############################

#Serial correlation: p-value has to be GREATER than 5% for NO Serial Correlation

serial1 <- serial.test(reduced1, lags.bg = 24, type = "BG")
serial2 <- serial.test(reduced2, lags.bg = 24, type = "BG")
serial3 <- serial.test(reduced3, lags.bg = 24, type = "BG")
serial4 <- serial.test(reduced4, lags.bg = 24, type = "BG")
serial5 <- serial.test(reduced5, lags.bg = 24, type = "BG")

serial1

#Heteroscedasticity: p-value has to be GREATER than 5% for No heteroscedasticity

arch1 <- arch.test(reduced1, lags.multi = 12, multivariate.only = TRUE)
arch2 <- arch.test(reduced2, lags.multi = 12, multivariate.only = TRUE)
arch3 <- arch.test(reduced3, lags.multi = 12, multivariate.only = TRUE)
arch4 <- arch.test(reduced4, lags.multi = 12, multivariate.only = TRUE)
arch5 <- arch.test(reduced5, lags.multi = 12, multivariate.only = TRUE)

arch1

#Normal Distribution of the Residual: p-value has to be GREATER than 5% for Normality

norm1 <- normality.test(reduced1, multivariate.only = TRUE)
norm2 <- normality.test(reduced2, multivariate.only = TRUE)
norm3 <- normality.test(reduced3, multivariate.only = TRUE)
norm4 <- normality.test(reduced4, multivariate.only = TRUE)
norm5 <- normality.test(reduced5, multivariate.only = TRUE)

#Structural Breaks in the Residuals

stability1 <- stability(reduced1, type = "OLS-CUSUM")
stability2 <- stability(reduced2, type = "OLS-CUSUM")
stability3 <- stability(reduced3, type = "OLS-CUSUM")
stability4 <- stability(reduced4, type = "OLS-CUSUM")
stability5 <- stability(reduced5, type = "OLS-CUSUM")


par(mar = c(1, 1, 1, 1))
#plot(1:30)
plot(stability4)


############  GRANGER CASUALITY - IMPULSE RESPONSE #############################

#GRANGER CAUSALITY searching for causality more than a simple correlation in both directions
#only to change the reduced model for every variable

granger <- causality(reduced3, cause = "Inflacion")
granger <- causality(reduced3, cause = "Combustibles")
granger <- causality(reduced3, cause = "Gasolina")
granger <- causality(reduced3, cause = "Petroleo")

granger

#IMPULSE - RESPONSE FUNCTION
#how the variable responds to a shock of the other variables


OILirf1 <- irf(reduced1, impulse = "Petroleo", response = "Inflacion", n.ahead = 36, boot = TRUE)
OILirf2 <- irf(reduced1, impulse = "Combustibles", response = "Inflacion", n.ahead = 36, boot = TRUE)
OILirf3 <- irf(reduced1, impulse = "Gasolina", response = "Inflacion", n.ahead = 36, boot = TRUE)


OILirf <- irf(reduced1, impulse = "Petroleum", response = "Inflacion", n.ahead = 36, boot = TRUE)


plot(OILirf, ylab = "Inflacion", main = "Shock from OIL")
plot(COMBirf, ylab = "Inflacion", main = "Shock from Petroleum")
#Variance Decomposition

vadc_1 <- fevd(reduced1, n.ahead = 10)
plot(vadc_1)

############  FORECAST AND VARIANCE DECOMPOSITION ##############################

model_f <- predict(reduced3, n.ahead = 12, ci = 0.95)

fanchart(model_f, xlab = "periodos", names = "Inflacion")
fanchart(model_f, names = "Combustibles")
fanchart(model_f, names = "Petroleo")

############  STRUCTURAL VAR W/ RESTRICTIONS #######################

#Restrictions for Structural Var

amat2 <- diag(2)
amat2[1, 2] <- NA

amat3 <- diag(3)
amat3[1, 2] <- NA
amat3[1, 3] <- NA
amat3[2, 3] <- NA

#amat3[2, 1] <- NA
#amat3[3, 1] <- NA
#amat3[3, 2] <- NA

amat4 <- diag(4)
amat4[1, 2] <- NA
amat4[1, 3] <- NA
amat4[1, 4] <- NA
amat4[2, 3] <- NA
amat4[2, 4] <- NA
amat4[3, 4] <- NA

#amat4[2, 1] <- NA
#amat4[3, 1] <- NA
#amat4[3, 2] <- NA
#amat4[4, 1] <- NA
#amat4[4, 2] <- NA
#amat4[4, 3] <- NA

bmat2 <- diag(2)
diag(bmat2) <- NA

bmat3 <- diag(3)
diag(bmat3) <- NA
bmat3

bmat4 <- diag(4)
diag(bmat4) <- NA



structural1 <- SVAR(reduced1, Amat = NULL, Bmat = bmat2, 
                    hessian = TRUE, estmethod = "scoring")
structural2 <- SVAR(reduced2, Amat = amat3, Bmat = bmat3, 
                    hessian = TRUE, estmethod = c("scoring", "direct"))

structural4 <- SVAR(reduced4, Amat = amat4, Bmat = NULL, 
                    hessian = TRUE, estmethod = c("scoring", "direct"))
structural5 <- SVAR(reduced5, Amat = amat4, Bmat = bmat4, 
                    hessian = TRUE, estmethod = c("scoring", "direct"))

############  CONCLUSIONS ###########

#BEST MODEL w/ no autocorrelation between variables asumption 
structural3 <- SVAR(reduced3, Amat = amat3, Bmat = bmat3, 
                    hessian = TRUE, estmethod = c("scoring", "direct"))

#Forecast

modelsvar <- predict(structural1, n.ahead = 24, ci = 0.95)

sCOMBirf <- irf(structural3, impulse = "Petroleo", response = "Inflacion", n.ahead = 36, boot = TRUE)
plot(sCOMBirf, ylab = "Inflacion", main = "Shock from Petroleum")
