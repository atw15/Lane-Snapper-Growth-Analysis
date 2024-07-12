# Load necessary libraries
library("tidyverse")
library("EnvStats") # Library for truncated normal distribution
library("rjags") # Interface to JAGS for Bayesian analysis
library("R2jags") # Another interface to JAGS
library("ggplot2")
library("dplyr")
library("magrittr")
library("ggridges")
theme_set(theme_bw()) # Set theme for ggplot2

# Read in data from CSV file
LaneAll = read.csv("Age Master Final_5_12.csv")
summary(LaneAll) # Summarize the dataset

# Convert Fork Length to numeric and remove rows with NA values in FFL.mm and Final.Age
LaneAll$FFL.mm = as.numeric(LaneAll$FFL.mm)
LaneAll = LaneAll[!is.na(LaneAll$FFL.mm) & !is.na(LaneAll$Final.Age),]
summary(LaneAll) # Summarize cleaned data
table(LaneAll$FFL.mm) # Display frequency table for FFL.mm

# Plot mean age by gear type and combined
ggplot(LaneAll) +
  stat_summary(aes(x = trunc(Final.Age), y = FFL.mm, color = Fishery)) +
  stat_summary(aes(x = trunc(Final.Age), y = FFL.mm), alpha = 0.5) +
  xlab("Mean Fractional Age (years)") + 
  ylab("Fork Length (mm)") +
  ggtitle("Lane Snapper: Mean Fork Length", subtitle = "974 samples")

# Plot histograms of lengths above and below size limit by Fishery
ggplot(filter(LaneAll, Final.Age < 9), aes(x = FFL.mm, color = Fishery, fill = Fishery)) +
  geom_histogram(alpha = 0.6, binwidth = 10) +
  geom_vline(xintercept = 203.2) +
  xlab("Fork Length (mm)") + 
  ylab("Count") +
  ggtitle("Lane Snapper: Lengths Above and Below Size Limit by Fishery", subtitle = "974 samples\nFractional Age < 9 years")

# Check for fish below size limit in the FD data
LaneAll = mutate(LaneAll, BelowLimit = ifelse(FFL.mm < 203.2 & Fishery != "FI", 1, 0))
table(LaneAll$BelowLimit) # Display table of BelowLimit

# Filter out fish below the size limit
LaneAll = filter(LaneAll, BelowLimit == 0)

# Define function to return the sum of the negative log likelihood for truncated normal
vonBertTrunc = function(params) {
  Linf = params[1]
  K = params[2]
  t0 = params[3]
  sigma = params[4]
  age = LaneAll$Final.Age
  L = LaneAll$FFL.mm
  Lpred = Linf * (1 - exp(-K * (age - t0)))
  SizeLimit = ifelse(LaneAll$Fishery == "FI", -Inf, 203.2)
  LL = log(dnormTrunc(L, mean = Lpred, sd = sigma, min = SizeLimit, max = Inf))
  -sum(LL)
}

# Define function to return the sum of the negative log likelihood for regular normal
vonBertNormal = function(params) {
  Linf = params[1]
  K = params[2]
  t0 = params[3]
  sigma = params[4]
  age = LaneAll$Final.Age
  L = LaneAll$FFL.mm
  Lpred = Linf * (1 - exp(-K * (age - t0)))
  LL2 = dnorm(L, mean = Lpred, sd = sigma, log = TRUE)
  -sum(LL2)
}

# Fit both models and compare
ModFit1 = nlminb(start = c(Linf = 400, K = 0.5, t0 = -0.5, sigma = 100), objective = vonBertTrunc)
ModFit1
ModFit2 = nlminb(start = c(Linf = 400, K = 0.5, t0 = -0.5, sigma = 10), objective = vonBertNormal)
ModFit2
# The parameters are only a little different
ModFit1$par # Truncated Normal parameters
ModFit2$par # Regular Normal parameters

# Define von Bertalanffy growth function
vonBertFunc = function(age, Linf, K, t0) Linf * (1 - exp(-K * (age - t0)))

# Make predictions from both fits
Predicted = data.frame(Age = seq(0, 12, 0.1)) %>% 
  mutate(Truncated = vonBertFunc(Age, ModFit1$par[1], ModFit1$par[2], ModFit1$par[3]),
         Normal = vonBertFunc(Age, ModFit2$par[1], ModFit2$par[2], ModFit2$par[3]))

# Plot the data with both fits
ggplot(Predicted) + 
  geom_line(aes(x = Age, y = Truncated), color = "red", lwd = 1) +
  geom_line(aes(x = Age, y = Normal), color = "blue", lwd = 1) +
  geom_point(data = LaneAll, aes(x = Final.Age, y = FFL.mm, color = Fishery)) +
  xlab("Fractional Age (years)") + 
  ylab("Fork Length (mm)") +
  ggtitle("Lane Snapper: Size at Age by Fishery", subtitle = "974 samples\nRed is Truncated Normal, Blue is Untruncated")

## Bayesian version to get credible intervals
# This also fits truncated, ordinary normal, and one that excludes FD data from ages that might be influenced by the size limit

write("model{
  # Priors for the parameters
  Linf ~ dunif(0, 1000)    # Uniform prior for Linf
  K ~ dunif(0, 3)          # Uniform prior for K
  t0 ~ dunif(-10, 0)       # Uniform prior for t0
  sigma ~ dunif(0, 200)    # Uniform prior for sigma
  tau <- 1 / (sigma * sigma)  

  # Likelihood for fishery-independent data
  for (i in 1:N.FI) {
    Lpred.FI[i] <- Linf * (1 - exp(-K * (age.FI[i] - t0)))  # Predicted length
    L.FI[i] ~ dnorm(Lpred.FI[i], tau)  # Likelihood with normal distribution
    resid.FI[i] <- L.FI[i] - Lpred.FI[i]  # Residuals
  }

  # Likelihood for fishery-dependent data with truncation
  for (i in 1:N.FD) {
    Lpred.FD[i] <- Linf * (1 - exp(-K * (age.FD[i] - t0)))  # Predicted length
    L.FD[i] ~ dnorm(Lpred.FD[i], tau) T(203.2,)  # Truncated normal distribution
    resid.FD[i] <- L.FD[i] - Lpred.FD[i]  # Residuals
  }
}", file = "vonbertTrunc.txt")

# Prepare data for JAGS
data1 <- list(
  age.FI = LaneAll$Final.Age[LaneAll$Fishery == "FI"],  # Ages for fishery-independent data
  L.FI = LaneAll$FFL.mm[LaneAll$Fishery == "FI"],  # Lengths for fishery-independent data
  N.FI = nrow(LaneAll[LaneAll$Fishery == "FI",]),  # Number of fishery-independent samples
  age.FD = LaneAll$Final.Age[LaneAll$Fishery != "FI"],  # Ages for fishery-dependent data
  L.FD = LaneAll$FFL.mm[LaneAll$Fishery != "FI"],  # Lengths for fishery-dependent data
  N.FD = nrow(LaneAll[LaneAll$Fishery != "FI",])  # Number of fishery-dependent samples
)

# Fit the model using JAGS
model1 <- jags(data = data1,
               parameters.to.save = c("Linf", "K", "t0", "sigma", "Lpred.FI", "Lpred.FD", "resid.FI", "resid.FD"),
               model.file = "vonbertTrunc.txt",
               n.iter = 21000, n.burnin = 1000, n.thin = 1)
round(model1$BUGSoutput$summary[c("Linf", "K", "t0", "sigma"),], 3)

# Extract residuals and predicted values
residsFImodel1 <- model1$BUGSoutput$summary[paste0("resid.FI[", 1:data1$N.FI, "]"), "mean"]
residsFDmodel1 <- model1$BUGSoutput$summary[paste0("resid.FD[", 1:data1$N.FD, "]"), "mean"]
LaneAll$residualmodel1[LaneAll$Fishery == "FI"] <- residsFImodel1
LaneAll$residualmodel1[LaneAll$Fishery != "FI"] <- residsFDmodel1
summary(LaneAll$residualmodel1)

LpredFImodel1 <- model1$BUGSoutput$summary[paste0("Lpred.FI[", 1:data1$N.FI, "]"), "mean"]
LpredFDmodel1 <- model1$BUGSoutput$summary[paste0("Lpred.FD[", 1:data1$N.FD, "]"), "mean"]
LaneAll$predictedmodel1[LaneAll$Fishery == "FI"] <- LpredFImodel1
LaneAll$predictedmodel1[LaneAll$Fishery != "FI"] <- LpredFDmodel1

# Plot residuals
plot(LaneKnownSex$predictedmodel1, LaneKnownSex$residualmodel1)
ggplot(LaneAll) +
  geom_point(aes(x = predictedmodel1, y = residualmodel1)) +
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Combined (sexed and unsexed) Model Truncated: Residuals") +
  xlab("Residuals") +
  ylab("Predicted")

# QQ plot of residuals
ggplot(LaneAll, aes(sample = residualmodel1)) +
  geom_qq() +
  geom_qq_line() +
  ggtitle("Combined (sexed and unsexed) Model Truncated QQNormal of Residuals")

# Ordinary normal model without truncation
write("model{
  # Priors for the parameters
  Linf ~ dunif(0, 1000)
  K ~ dunif(0, 3)
  t0 ~ dunif(-10, 0)
  sigma ~ dunif(0, 200)
  tau <- 1 / (sigma * sigma)

  # Likelihood for all data
  for (i in 1:N) {
    Lpred[i] <- Linf * (1 - exp(-K * (age[i] - t0)))  # Predicted length
    L[i] ~ dnorm(Lpred[i], tau)  # Likelihood with normal distribution
    resid[i] <- L[i] - Lpred[i]  # Residuals
  }
}", file = "vonbertNorm.txt")

# Prepare data for JAGS
data2 <- list(
  age = LaneAll$Final.Age,
  L = LaneAll$FFL.mm,
  N = nrow(LaneAll)
)

# Fit the model using JAGS
model2 <- jags(data = data2,
               parameters.to.save = c("Linf", "K", "t0", "sigma", "Lpred", "resid"),
               model.file = "vonbertNorm.txt",
               n.iter = 21000, n.burnin = 1000, n.thin = 1)
round(model2$BUGSoutput$summary[c("Linf", "K", "t0", "sigma"),], 3)

# Extract residuals and predicted values
residsmodel2 <- model2$BUGSoutput$summary[paste0("resid[", 1:data2$N, "]"), "mean"]
LaneAll$residualmodel2 <- residsmodel2

Lpredmodel2 <- model2$BUGSoutput$summary[paste0("Lpred[", 1:data2$N, "]"), "mean"]
LaneAll$predictedmodel2 <- Lpredmodel2

# Plot residuals
plot(LaneAll$predictedmodel2, LaneAll$residualmodel2)

ggplot(LaneAll) +
  geom_point(aes(x = predictedmodel2, y = residualmodel2)) +
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Combined (sexed and unsexed) Model Normal: Residuals") +
  xlab("Residuals") +
  ylab("Predicted")

# QQ plot of residuals
ggplot(LaneAll, aes(sample = residualmodel2)) +
  geom_qq() +
  geom_qq_line() +
  ggtitle("Combined (sexed and unsexed) Model Normal QQNormal of Residuals")

# Get predictions from both model results for plots
predvalsTrunc <- bind_rows(
  data.frame(model1$BUGSoutput$summary[paste0("Lpred.FI[", 1:data1$N.FI, "]"), c("mean", "2.5%", "97.5%")]),
  data.frame(model1$BUGSoutput$summary[paste0("Lpred.FD[", 1:data1$N.FD, "]"), c("mean", "2.5%", "97.5%")])
)
predvalsTrunc$Age <- c(data1$age.FI, data1$age.FD)

predvalsNormal <- data.frame(model2$BUGSoutput$summary[paste0("Lpred[", 1:nrow(LaneAll), "]"), c("mean", "2.5%", "97.5%")])
predvalsNormal$Age <- LaneAll$Final.Age

# Plot predictions from both models
ggplot() + 
  geom_line(data = predvalsTrunc, aes(x = Age, y = mean), color = "red", lwd = 1) +
  geom_ribbon(data = predvalsTrunc, aes(x = Age, y = mean, ymin = X2.5., ymax = X97.5.), fill = "red", alpha = 0.25) +
  geom_line(data = predvalsNormal, aes(x = Age, y = mean), color = "blue", lwd = 1) +
  geom_ribbon(data = predvalsNormal, aes(x = Age, y = mean, ymin = X2.5., ymax = X97.5.), fill = "blue", alpha = 0.25)

# Define function to return sum of negative log likelihood for truncated normal
vonBertTrunc<-function(params)  {
  Linf<-params[1]
  K<-params[2]
  t0<-params[3]
  sigma<-params[4]
  age<-LaneAll$Final.Age
  L<-LaneAll$FFL.mm
  Lpred<-Linf*(1-exp(-K*(age-t0))) # Von Bertalanffy growth function
  SizeLimit<-ifelse(LaneAll$Fishery=="FI",-Inf,203.2) # Size limit for fishery-dependent data
  LL<-log(dnormTrunc(L, mean = Lpred, sd = sigma, min = SizeLimit, max = Inf)) # Truncated normal log-likelihood
  -sum(LL)
}

# Define function to return sum of negative log likelihood for regular normal
vonBertNormal<-function(params)  {
  Linf<-params[1]
  K<-params[2]
  t0<-params[3]
  sigma<-params[4]
  age<-LaneAll$Final.Age
  L<-LaneAll$FFL.mm
  Lpred<-Linf*(1-exp(-K*(age-t0))) # Von Bertalanffy growth function
  LL2<-dnorm(L, mean = Lpred, sd = sigma,  log = TRUE) # Normal log-likelihood
  -sum(LL2)
}

# Check for fish below size limit in the FD data and filter them out
LaneKnownSex<-LaneAll[LaneAll$Macro.sex !="U",] # Remove unknown sex samples
LaneKnownSex$sexnum<-ifelse(LaneKnownSex$Macro.sex=="F",1,2) # Assign numeric values to sex
table(LaneKnownSex$sexnum,LaneKnownSex$Macro.sex) # Check sex distribution
table(LaneKnownSex$Macro.sex) # Check sex distribution
sex=LaneKnownSex$Macro.sex
# 1 is female, 2 is male
LaneKnownSex<-mutate(LaneKnownSex,BelowLimit=ifelse(FFL.mm<203.2 & Fishery!="FI",1,0)) # Mark samples below size limit
table(LaneAll$BelowLimit) # Check count of samples below limit

LaneKnownSex<-filter(LaneKnownSex,BelowLimit==0) # Filter out samples below size limit

# Prepare data list for JAGS model with sex-specific K and t0 values, single Linf, and fishery FD vs FI
data2K<-list(
  age.FI=LaneKnownSex$Final.Age[LaneKnownSex$Fishery=="FI"],
  L.FI=LaneKnownSex$FFL.mm[LaneKnownSex$Fishery=="FI"],
  N.FI=nrow(LaneKnownSex[LaneKnownSex$Fishery=="FI",]),
  age.FD=LaneKnownSex$Final.Age[LaneKnownSex$Fishery!="FI"],
  L.FD=LaneKnownSex$FFL.mm[LaneKnownSex$Fishery!="FI"],
  N.FD=nrow(LaneKnownSex[LaneKnownSex$Fishery!="FI",]),
  sex.FD=LaneKnownSex$sexnum[LaneKnownSex$Fishery!="FI"],
  sex.FI=LaneKnownSex$sexnum[LaneKnownSex$Fishery=="FI"]
)

# Write JAGS model to file with sex-dependent truncation
write("model{
  # priors
  Linf~dunif(0,1000)
  K[1]~dunif(0,3)
  K[2]~dunif(0,3)
  t0[1]~dunif(-10,0)
  t0[2]~dunif(-10,0)
  sigma~dunif(0,200)
  tau<-1/(sigma*sigma)
  for(i in 1:N.FI) {
   Lpred.FI[i]<-Linf*(1-exp(-K[sex.FI[i]]*(age.FI[i]-t0[sex.FI[i]])))
   L.FI[i]~dnorm(Lpred.FI[i],tau)
   resid.FI[i]<-L.FI[i]-Lpred.FI[i]
  }
  for(i in 1:N.FD) {
   Lpred.FD[i]<-Linf*(1-exp(-K[sex.FD[i]]*(age.FD[i]-t0[sex.FD[i]])))
   L.FD[i]~dnorm(Lpred.FD[i],tau)T(203.2,)
   resid.FD[i]<-L.FD[i]-Lpred.FD[i]
  }
}",file="vonbertTrunc2K.txt")

# Model with sex-dependent truncation using sexed samples only and male and female parameters
model2KTrunc<-jags(data=data2K,
                   parameters.to.save = c("Linf","K[1]","K[2]","t0[1]","t0[2]","sigma","Lpred.FI","Lpred.FD","resid.FI", "resid.FD"),
                   model.file = "vonbertTrunc2K.txt",
                   n.iter=21000,n.burnin = 1000,n.thin=1) 

# Display summary of model results
round(model2KTrunc$BUGSoutput$summary[c("Linf","K[1]","K[2]","t0[1]","t0[2]","sigma"),],3)

# Extract residuals and predicted values from model results
residsFImodel2KT<-model2KTrunc$BUGSoutput$summary[paste0("resid.FI[",1:data2K$N.FI,"]"),"mean"]
residsFDmodel2KT<-model2KTrunc$BUGSoutput$summary[paste0("resid.FD[",1:data2K$N.FD,"]"),"mean"]
LaneKnownSex$residualmodel2KT[LaneKnownSex$Fishery=="FI"]<-residsFImodel2KT
LaneKnownSex$residualmodel2KT[LaneKnownSex$Fishery!="FI"]<-residsFDmodel2KT
summary(LaneKnownSex$residualmodel2KT)

LpredFImodel2KT<-model2KTrunc$BUGSoutput$summary[paste0("Lpred.FI[",1:data2K$N.FI,"]"),"mean"]
LpredFDmodel2KT<-model2KTrunc$BUGSoutput$summary[paste0("Lpred.FD[",1:data2K$N.FD,"]"),"mean"]
LaneKnownSex$predictedmodel2KT[LaneKnownSex$Fishery=="FI"]<-LpredFImodel2KT
LaneKnownSex$predictedmodel2KT[LaneKnownSex$Fishery!="FI"]<-LpredFDmodel2KT
model2KTRP<-plot(LaneKnownSex$predictedmodel2KT,LaneKnownSex$residualmodel2KT)

# Plot residuals for sex-specific model with truncation
ggplot(LaneKnownSex)+geom_point(aes(x=predictedmodel2KT,y=residualmodel2KT))+geom_abline(intercept=0,slope=0)+
  ggtitle("Sex-Specific Model Truncated: Residuals")+
  xlab("Residuals")+
  ylab("Predicted")

# Plot QQ plot of residuals for sex-specific model with truncation
ggplot(LaneKnownSex, aes(sample=residualmodel2KT)) + 
  geom_qq() + 
  geom_qq_line() +
  ggtitle("Sex-Specific Model Truncated QQNormal of Residuals")

# Calculate DIC for model2KTrunc
model2KTrunc$BUGSoutput$DIC # Use DIC for model comparison; smaller DIC indicates a better model

# Model without sex comparison, using only sexed samples combined
model1KTrunc <- jags(data = data2K,
                     parameters.to.save = c("Linf", "K", "t0", "sigma", "Lpred.FI", "Lpred.FD", "resid.FI", "resid.FD"),
                     model.file = "vonbertTrunc.txt",
                     n.iter = 21000, n.burnin = 1000, n.thin = 1)

# Display summary of model1KTrunc results
round(model1KTrunc$BUGSoutput$summary[c("Linf", "K", "t0", "sigma"), ], 3)

# Extract residuals and predicted values from model1KTrunc results
residsFImodel1KT <- model1KTrunc$BUGSoutput$summary[paste0("resid.FI[", 1:data2K$N.FI, "]"), "mean"]
residsFDmodel1KT <- model1KTrunc$BUGSoutput$summary[paste0("resid.FD[", 1:data2K$N.FD, "]"), "mean"]
LaneKnownSex$residualmodel1KT[LaneKnownSex$Fishery == "FI"] <- residsFImodel1KT
LaneKnownSex$residualmodel1KT[LaneKnownSex$Fishery != "FI"] <- residsFDmodel1KT
summary(LaneKnownSex$residualmodel1KT)

LpredFImodel1KT <- model1KTrunc$BUGSoutput$summary[paste0("Lpred.FI[", 1:data2K$N.FI, "]"), "mean"]
LpredFDmodel1KT <- model1KTrunc$BUGSoutput$summary[paste0("Lpred.FD[", 1:data2K$N.FD, "]"), "mean"]
LaneKnownSex$predictedmodel1KT[LaneKnownSex$Fishery == "FI"] <- LpredFImodel1KT
LaneKnownSex$predictedmodel1KT[LaneKnownSex$Fishery != "FI"] <- LpredFDmodel1KT

model1KTRP <- plot(LaneKnownSex$predictedmodel1KT, LaneKnownSex$residualmodel1KT)

# Plot residuals for combined (sexed) model with truncation
ggplot(LaneKnownSex) + 
  geom_point(aes(x = predictedmodel1KT, y = residualmodel1KT)) + 
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Combined (sexed) Model Truncated: Residuals") +
  xlab("Residuals") +
  ylab("Predicted")

# Plot QQ plot of residuals for combined (sexed) model with truncation
ggplot(LaneKnownSex, aes(sample = residualmodel1KT)) + 
  geom_qq() + 
  geom_qq_line() +
  ggtitle("Combined (sexed) Model Truncated QQNormal of Residuals")

# Calculate DIC for model1KTrunc
model1KTrunc$BUGSoutput$DIC # Use DIC for model comparison; smaller DIC indicates a better model

# Model without sex comparison, normal likelihood, one set of growth parameters
write("model{
  # priors
  Linf~dunif(0,1000)
  K~dunif(0,3)
  t0~dunif(-10,0)
  sigma~dunif(0,200)
  tau<-1/(sigma*sigma)
  for(i in 1:N.FI) {
   Lpred.FI[i]<-Linf*(1-exp(-K*(age.FI[i]-t0)))
   L.FI[i]~dnorm(Lpred.FI[i],tau)
   resid.FI[i]<-L.FI[i]-Lpred.FI[i]
  }
  for(i in 1:N.FD) {
   Lpred.FD[i]<-Linf*(1-exp(-K*(age.FD[i]-t0)))
   L.FD[i]~dnorm(Lpred.FD[i],tau)
   resid.FD[i]<-L.FD[i]-Lpred.FD[i]
  }
}", file = "vonbert1kNorm.txt")

# Model without sex comparison, normal likelihood
model1KNorm <- jags(data = data2K,
                    parameters.to.save = c("Linf", "K", "t0", "sigma", "Lpred.FI", "Lpred.FD", "resid.FI", "resid.FD"),
                    model.file = "vonbert1kNorm.txt",
                    n.iter = 21000, n.burnin = 1000, n.thin = 1)

# Display summary of model1KNorm results
round(model1KNorm$BUGSoutput$summary[c("Linf", "K", "t0", "sigma"), ], 3)

# Extract residuals and predicted values from model1KNorm results
residsFImodel1KN <- model1KNorm$BUGSoutput$summary[paste0("resid.FI[", 1:data2K$N.FI, "]"), "mean"]
residsFDmodel1KN <- model1KNorm$BUGSoutput$summary[paste0("resid.FD[", 1:data2K$N.FD, "]"), "mean"]
LaneKnownSex$residualmodel1KN[LaneKnownSex$Fishery == "FI"] <- residsFImodel1KN
LaneKnownSex$residualmodel1KN[LaneKnownSex$Fishery != "FI"] <- residsFDmodel1KN
summary(LaneKnownSex$residualmodel1KN)

LpredFImodel1KN <- model1KNorm$BUGSoutput$summary[paste0("Lpred.FI[", 1:data2K$N.FI, "]"), "mean"]
LpredFDmodel1KN <- model1KNorm$BUGSoutput$summary[paste0("Lpred.FD[", 1:data2K$N.FD, "]"), "mean"]
LaneKnownSex$predictedmodel1KN[LaneKnownSex$Fishery == "FI"] <- LpredFImodel1KN
LaneKnownSex$predictedmodel1KN[LaneKnownSex$Fishery != "FI"] <- LpredFDmodel1KN

model1KNRP <- plot(LaneKnownSex$predictedmodel1KN, LaneKnownSex$residualmodel1KN)

# Plot residuals for combined (sexed) model with normal likelihood
ggplot(LaneKnownSex) + 
  geom_point(aes(x = predictedmodel1KN, y = residualmodel1KN)) + 
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Combined (sexed) Model Normal: Residuals") +
  xlab("Residuals") +
  ylab("Predicted")

# Plot QQ plot of residuals for combined (sexed) model with normal likelihood
ggplot(LaneKnownSex, aes(sample = residualmodel1KN)) + 
  geom_qq() + 
  geom_qq_line() +
  ggtitle("Combined (sexed) Model Normal QQNormal of Residuals")

# Calculate DIC for model1KNorm
model1KNorm$BUGSoutput$DIC # Use DIC for model comparison; smaller DIC indicates a better model

# Write JAGS model with sex-specific parameters
write("model{
  # priors
  Linf~dunif(0,1000) #females
  K[1]~dunif(0,3)
  K[2]~dunif(0,3)
  t0[1]~dunif(-10,0)
  t0[2]~dunif(-10,0)
  sigma~dunif(0,200)
  tau<-1/(sigma*sigma)
  for(i in 1:N.FI) {
   Lpred.FI[i]<-Linf*(1-exp(-K[sex.FI[i]]*(age.FI[i]-t0[sex.FI[i]])))
   L.FI[i]~dnorm(Lpred.FI[i],tau)
   resid.FI[i]<-L.FI[i]-Lpred.FI[i]
  }
  for(i in 1:N.FD) {
   Lpred.FD[i]<-Linf*(1-exp(-K[sex.FD[i]]*(age.FD[i]-t0[sex.FD[i]])))
   L.FD[i]~dnorm(Lpred.FD[i],tau)
   resid.FD[i]<-L.FD[i]-Lpred.FD[i]
  }
}", file = "vonbertnormal2K.txt")

table(LaneAll$Fishery) # Display the count of samples by Fishery type

# Run the JAGS model for only sexed samples with two sets of growth parameters, not truncated
model2KNorm <- jags(data = data2K,
                    parameters.to.save = c("Linf", "K", "t0", "sigma", "Lpred.FI", "Lpred.FD", "resid.FI", "resid.FD"),
                    model.file = "vonbertnormal2K.txt",
                    n.iter = 21000, n.burnin = 1000, n.thin = 1) 

# Display summary of model2KNorm results
round(model2KNorm$BUGSoutput$summary[c("Linf", "K[1]", "K[2]", "t0[1]", "t0[2]", "sigma"), ], 3)

# Extract residuals and predicted values from model2KNorm results
residsFImodel2KN <- model2KNorm$BUGSoutput$summary[paste0("resid.FI[", 1:data2K$N.FI, "]"), "mean"]
residsFDmodel2KN <- model2KNorm$BUGSoutput$summary[paste0("resid.FD[", 1:data2K$N.FD, "]"), "mean"]
LaneKnownSex$residualmodel2KN[LaneKnownSex$Fishery == "FI"] <- residsFImodel2KN
LaneKnownSex$residualmodel2KN[LaneKnownSex$Fishery != "FI"] <- residsFDmodel2KN
summary(LaneKnownSex$residualmodel2KN)

LpredFImodel2KN <- model2KNorm$BUGSoutput$summary[paste0("Lpred.FI[", 1:data2K$N.FI, "]"), "mean"]
LpredFDmodel2KN <- model2KNorm$BUGSoutput$summary[paste0("Lpred.FD[", 1:data2K$N.FD, "]"), "mean"]
LaneKnownSex$predictedmodel2KN[LaneKnownSex$Fishery == "FI"] <- LpredFImodel2KN
LaneKnownSex$predictedmodel2KN[LaneKnownSex$Fishery != "FI"] <- LpredFDmodel2KN

model2KNRP <- plot(LaneKnownSex$predictedmodel2KN, LaneKnownSex$residualmodel2KN)

# Plot residuals for sex-specific model with normal likelihood
ggplot(LaneKnownSex) + 
  geom_point(aes(x = predictedmodel2KN, y = residualmodel2KN)) + 
  geom_abline(intercept = 0, slope = 0) +
  ggtitle("Sex-Specific Model Normal: Residuals") +
  xlab("Residuals") +
  ylab("Predicted")

# Plot QQ plot of residuals for sex-specific model with normal likelihood
ggplot(LaneKnownSex, aes(sample = residualmodel2KN)) + 
  geom_qq() + 
  geom_qq_line() +
  ggtitle("Sex-Specific Model Normal QQNormal of Residuals")

# Calculate DIC for each model for comparison
model2$BUGSoutput$DIC # All samples combined, normal likelihood
model1$BUGSoutput$DIC # All samples combined, truncated normal likelihood

# The model with all samples combined and truncated has the lowest DIC
# 10063.28 vs 10133.7

# Model with only sexed samples without sex comparison, normal likelihood
model1KNorm$BUGSoutput$DIC # Use DIC for each, whichever has the smallest is the best
round(model1KNorm$BUGSoutput$summary[c("Linf", "K", "t0", "sigma"), ], 3)

# Model with only sexed samples with one set of growth parameters, truncated normal likelihood
model1KTrunc$BUGSoutput$DIC # Use DIC for each, whichever has the smallest is the best
model4 <- round(model1KTrunc$BUGSoutput$summary[c("Linf", "K", "t0", "sigma"), ], 3)

# Model with sexed samples, two sets of growth parameters, and normal likelihood
model2KNorm$BUGSoutput$DIC # Use DIC for each, whichever has the smallest is the best
model5 <- round(model2KNorm$BUGSoutput$summary[c("Linf", "K[1]", "K[2]", "t0[1]", "t0[2]", "sigma"), ], 3)

# Model with only sexed samples, two sets of growth parameters, and truncated normal likelihood
model2KTrunc$BUGSoutput$DIC # Use DIC for each, whichever has the smallest is the best
round(model2KTrunc$BUGSoutput$summary[c("Linf", "K[1]", "K[2]", "t0[1]", "t0[2]", "sigma"), ], 3)
