# Clear all existing objects from the R environment to ensure a clean workspace.
rm(list=ls(all=TRUE))

# Load the dataset from a CSV file.
LaneAll = read.csv("Age Master Final_5_12.csv")

# Generate a summary of the 'Final.Age' variable in the dataset.
summary(LaneAll$Final.Age)

# Ensure necessary libraries are installed and loaded. Install 'devtools' if not already installed.
if (!require('devtools')) install.packages('devtools'); require('devtools')

# Install the 'FSAmisc' package from GitHub, forcing the newest version.
devtools::install_github('droglenc/FSAmisc', force = TRUE)

# Display citation information for R.
citation()

# Load additional required libraries.
library(magrittr)
library(FSA)
library(dplyr)
library(car)
library(nlstools)
library(patchwork)

# View a summary of the data and a frequency table of 'Macro.sex'.
summary(LaneAll)
table(LaneAll$Macro.sex)

# Convert specific columns to numeric data types and remove rows with NA values in these columns.
LaneAll$FFL.mm = as.numeric(LaneAll$FFL.mm)
LaneAll$Final.Age = as.numeric(LaneAll$Final.Age)
LaneAll = LaneAll[!is.na(LaneAll$FFL.mm) & !is.na(LaneAll$Final.Age),]

# Confirm the changes to data types.
is.numeric(LaneAll$Final.Age)
is.numeric(LaneAll$FFL.mm)
summary(LaneAll)

# DATA STATISTICS SECTION

# Subset data to only include rows where 'Macro.sex' is known (not 'U').
LaneKnownSex <- LaneAll[LaneAll$Macro.sex != "U",]
LaneKnownSex$sexnum <- ifelse(LaneKnownSex$Macro.sex == "F", 1, 2) # Assign numeric values to sexes: 1 for female, 2 for male.
table(LaneKnownSex$sexnum, LaneKnownSex$Macro.sex)
table(LaneKnownSex$Macro.sex)
summary(LaneKnownSex$FFL.mm)

# Get initial parameter estimates for nonlinear modeling using the von Bertalanffy growth function.
svTyp = vbStarts(FFL.mm ~ Final.Age, data = LaneAll, plot = TRUE)
svTyp = vbStarts(FFL.mm ~ Final.Age, data = LaneKnownSex, plot = TRUE)

# Display help information for 'vbStarts'.
?vbStarts

# Display the maximum value of 'FFL.mm'.
max(LaneAll$FFL.mm)

# Define the von Bertalanffy Growth Function (VBGF) to predict mean length at a given age.
vbTyp = function(Final.Age, Linf, K, t0) Linf * (1 - exp(-K * (Final.Age - t0)))

# Examples of predicting lengths using the VBGF.
vbTyp(2, Linf = 495, K = 0.67, t0 = 0)
vbTyp(3, Linf = 495, K = 0.67, t0 = 0)
vbTyp

# Fit the VBGF model to data, estimate parameters, and calculate confidence intervals.
svTyp
fitTyp = nls(FFL.mm ~ vbTyp(Final.Age, Linf, K, t0), data = LaneAll, start = svTyp)
coef(fitTyp)
confint(fitTyp)

# Perform bootstrapping for robustness check of the nonlinear model.
bootTyp = nlsBoot(fitTyp)
headtail(bootTyp$coefboot, n = 2)
confint(bootTyp, plot = TRUE)

# Display a detailed summary of the model fit, including correlations.
summary(fitTyp, correlation = TRUE)

## Predictions 12.3.4 Derek Ogle

# Predict lengths at specific ages using the fitted model.
nd = data.frame(Final.Age = c(1, 3, 5, 9, 12, 14))
predict(fitTyp, nd) # Predictions for the specified ages using the nonlinear model.

# Predict lengths at specific ages using the VBGF function and coefficients from the model.
vbTyp(c(1, 3, 5, 9, 12, 14), coef(fitTyp))
vbTyp(3, bootTyp$coefboot[1,]) # Prediction at age 3 using the first set of bootstrapped coefficients.

# Apply the VBGF function to each set of bootstrapped coefficients at age 3.
p3Typ = apply(bootTyp$coefboot, MARGIN = 1, FUN = vbTyp, t = 3)
p3Typ[1:6] # Display the first six predictions.

# Calculate the 2.5th and 97.5th percentiles of the bootstrapped predictions.
quantile(p3Typ, c(0.025, 0.975))

## Visualizing Model Fit 12.3.5

# Retrieve the coefficients from the fitted model.
coef(fitTyp)

# Generate a sequence of ages for predictions.
x = seq(0, 15, length.out = 199)
# Predicted lengths using the VBGF and the model coefficients.
pTyp = vbTyp(x, Linf = coef(fitTyp))

# Determine the range for the x and y axes.
xlmts = range(c(x, LaneAll$Final.Age))
ylmts = range(c(pTyp, LaneAll$FFL.mm))

# Plot the observed data and the predicted lengths.
plot(FFL.mm ~ Final.Age, data = LaneAll, xlab = "Age", ylab = "Full Fork Length (mm)",
     xlim = xlmts, ylim = ylmts, pch = 19, col = rgb(0, 0, 0, 1/3))
lines(pTyp, lwd = 2)

# Initialize vectors to store the lower and upper confidence intervals.
LCI = UCI = numeric(length(x))

# Calculate the confidence intervals for each age.
for (i in 1:length(x)) {
  tmp = apply(bootTyp$coefboot, MARGIN = 1, FUN = vbTyp, t = x[i])
  LCI[i] = quantile(tmp, 0.025)
  UCI[i] = quantile(tmp, 0.975)
}

# Update the y-axis limits to include the confidence intervals.
ylmts = range(c(pTyp, LCI, UCI, LaneAll$FFL.mm))

# Plot the observed data, predicted lengths, and confidence intervals.
plot(FFL.mm ~ Final.Age, data = LaneAll, xlab = "Age", ylab = "Full Fork Length (mm)",
     xlim = xlmts, ylim = ylmts, pch = 19, col = rgb(0, 0, 0, 1/3))
lines(pTyp ~ x, lwd = 2)
lines(UCI ~ x, lwd = 2, lty = "dashed")
lines(LCI ~ x, lwd = 2, lty = "dashed")

# Plot the residuals of the nonlinear model fit.
FSAmisc::residPlot(fitTyp)

# Add a log-transformed Full Fork Length variable to the dataset.
LaneAll %<>% mutate(logFFL = log(FFL.mm))

# Fit a nonlinear model to the log-transformed data.
fitTypM = nls(logFFL ~ log(vbTyp(Final.Age, Linf, K, t0)),
              data = LaneAll, start = svTyp)

# Plot the residuals of both the original and log-transformed model fits.
FSAmisc::residPlot(fitTyp)
FSAmisc::residPlot(fitTypM)

# Display the summary of the log-transformed model fit.
fitTypM

## Model Fitting 12.4.2 Derek Ogle

# Subset data to exclude unknown sex and create a numeric variable for sex.
LaneKnownSex <- LaneAll[LaneAll$Macro.sex != "U",]
LaneKnownSex$sexnum <- ifelse(LaneKnownSex$Macro.sex == "F", 1, 2) # 1 for female, 2 for male.
table(LaneKnownSex$sexnum, LaneKnownSex$Macro.sex)
sex = LaneKnownSex$Macro.sex
table(LaneKnownSex$sexnum, LaneKnownSex$Macro.sex)
table(sex)
FFL = LaneKnownSex$FFL.mm
Age = LaneKnownSex$Final.Age

# Summarize 'Catch.Year', 'State', and 'GulfSide' variables.
summary(LaneAll$Catch.Year)
table(LaneAll$State)
table(LaneAll$GulfSide)

# Define various von Bertalanffy growth function (VBGF) equations with different parameters.
vbLKt = FFL ~ Linf[sexnum] * (1 - exp(-K[sexnum] * (Age - t0[sexnum])))
vbLK = FFL ~ Linf[sexnum] * (1 - exp(-K[sexnum] * (Age - t0)))
vbLt = FFL ~ Linf[sexnum] * (1 - exp(-K * (Age - t0[sexnum])))
vbKt = FFL ~ Linf * (1 - exp(-K[sexnum] * (Age - t0[sexnum])))
vbL = FFL ~ Linf[sexnum] * (1 - exp(-K * (Age - t0)))
vbK = FFL ~ Linf * (1 - exp(-K[sexnum] * (Age - t0)))
vbt = FFL ~ Linf * (1 - exp(-K * (Age - t0[sexnum])))
vb0 = FFL ~ Linf * (1 - exp(-K * (Age - t0)))

# Get initial parameter estimates for nonlinear models.
vbStarts(FFL.mm ~ Final.Age, data = LaneKnownSex, vbStartsDP = TRUE)
sv0 = vbStarts(FFL ~ Age, data = LaneKnownSex)
sv0

# Replicate starting values for different groups (change '2' to the number of groups).
svLKt = Map(rep, sv0, c(2, 2, 2))
svLKt

# Fit the nonlinear model with the specified starting values.
fitLKt = nls(vbLKt, data = LaneKnownSex, start = svLKt)

# Plot the residuals of the fitted model.
FSAmisc::residPlot(fitLKt, col = rgb(0, 0, 0, 1/3))
fit0 = nls(vb0, data = LaneKnownSex, start = sv0)

# Perform a likelihood ratio test to compare models.
lrt(fit0, com = fitLKt, com.name = "All pars differ", sim.names = "No pars differ")
extraSS(fit0, com = fitLKt, com.name = "All pars diff", sim.names = "No pars diff")

# Interpretation of likelihood ratio and extra sum of squares tests.
# Likelihood Ratio (X^2 = 18.95, p = 0.00028) indicates significant difference (p < 0.05).
# Extra sum of squares test (F = 6.35, p = 0.0003037) also indicates significant difference.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
svLK = Map(rep, sv0, c(2, 2, 1))
svLK

# Starting values for other nested models.
svLt = Map(rep, sv0, c(2, 1, 2))
svLt
svKt = Map(rep, sv0, c(1, 2, 2))
svKt

# Fit the nested models.
fitLK = nls(vbLK, data = LaneKnownSex, start = svLK)
fitLt = nls(vbLt, data = LaneKnownSex, start = svLt)
fitKt = nls(vbKt, data = LaneKnownSex, start = svKt)

# Compare the nested models using likelihood ratio tests.
lrt(fitLK, fitLt, fitKt, com = fitLKt, com.name = "All pars diff", 
    sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: Strong evidence that all models fit equally well. The K,t0 model has the greatest log-likelihood value.

# Further starting values for additional comparisons.
svL = Map(rep, sv0, c(2, 1, 1))
svt = Map(rep, sv0, c(1, 1, 2))
svK = Map(rep, sv0, c(1, 2, 1))

# Fit additional models.
fitL = nls(vbL, data = LaneKnownSex, start = svL)
fitt = nls(vbt, data = LaneKnownSex, start = svt)
fitK = nls(vbK, data = LaneKnownSex, start = svK)

# Compare additional models.
lrt(fitt, fitK, com = fitKt, com.name = "K, t0 dif", sim.names = c("t0 dif", "K dif"))
# Interpretation: Results suggest that K (p = 0.06) fits as well as (K,t0), but not t0 (p = 0.0008).

lrt(fit0, com = fitK, com.name = "K dif", sim.names = "No pars dif")
# Interpretation: Comparison confirms a significant difference between groups for K (p < 0.5).

# Model Selection with AIC or BIC

# Fit model with starting values.
svK = Map(rep, sv0, c(1, 2, 1))
fitK = nls(vbK, data = LaneKnownSex, start = svK)

# Compare models using AIC and BIC.
cbind(AIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0),
      BIC(fitLKt, fitLK, fitLt, fitKt, fitL, fitK, fitt, fit0))

# Interpretation: The fitKt model has the lowest AIC (6424.84) and fitK has the lowest BIC (6448.50), indicating the best model according to these criteria.

# Plot best fit model
predict(fitKt) # Predict values using the fitKt model.
coef(fitKt) # Extract coefficients from the fitKt model.

# Plot observed data.
plot(LaneKnownSex$Final.Age, LaneKnownSex$FFL.mm)

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Plot the VBGF for females (red line) and males (green line).
lines(seq(0, 14, 0.1), vonBert(coef(fitKt)[1], coef(fitKt)[4], coef(fitKt)[2], seq(0, 14, 0.1)), col = 2) # Female
lines(seq(0, 14, 0.1), vonBert(coef(fitKt)[1], coef(fitKt)[5], coef(fitKt)[3], seq(0, 14, 0.1)), col = 3) # Male

# Display the distribution of 'Fishery' in the dataset.
table(LaneKnownSex$Fishery)

# Create a plot with customized title and labels.
plot(LaneKnownSex$Final.Age, LaneKnownSex$FFL.mm, main = "Fork Length at Fractional Age for Lane Snapper", font.main = 1,
     xlab = "Fractional Age (years)", ylab = "Fork Length (mm)", col = "white")

# Plot observed data points for females and males.
points(FFL.mm ~ Final.Age, data = lsf, pch = 1, col = rgb(0, 0, 0, 1/2), cex = 0.8) # Female points
points(FFL.mm ~ Final.Age, data = lsm, pch = 8, col = rgb(0, 0, 0, 1/2), cex = 0.8) # Male points

# Plot VBGF lines for females and males.
lines(seq(0, 14, 0.1), vonBert(coef(fitKt)[1], coef(fitKt)[4], coef(fitKt)[2], seq(0, 14, 0.1)), col = "black", lwd = 3, lty = 2) # Female
lines(seq(0, 14, 0.1), vonBert(coef(fitKt)[1], coef(fitKt)[5], coef(fitKt)[3], seq(0, 14, 0.1)), col = "black", lwd = 2, lty = 1) # Male

# Add a legend to the plot.
legend("bottomright", c("Female", "Male", "Female Predicted", "Male Predicted"),
       pch = c(1, 8, NA, NA), lty = c(NA, NA, 2, 1), bty = "n", cex = 0.8, lwd = c(1, 1, 2, 2))

# Summarizing the Model Fit 12.4.3
# The following code creates a function for the typical VBGF, creates separate subsets for females and males,
# fits the model, and extracts parameter estimates and bootstrap confidence intervals (CI) for each sex.

vbTyp = vbFuns("typical")

## Females
lsf = dplyr::filter(LaneKnownSex, Macro.sex == "F") # Subset data for females.
fitf = nls(FFL.mm ~ vbTyp(Final.Age, Linf, K, t0), data = lsf, start = sv0) # Fit the model.
bcf = nlsBoot(fitf) # Bootstrap the model.
cbind(coef(fitf), confint(bcf)) # Combine coefficients and confidence intervals.

summary(lsf$FFL.mm) # Summary of female fork lengths.
summary(lsm$FFL.mm) # Summary of male fork lengths.

## Males
lsm = dplyr::filter(LaneKnownSex, Macro.sex == "M") # Subset data for males.
fitm = nls(FFL.mm ~ vbTyp(Final.Age, Linf, K, t0), data = lsm, start = sv0) # Fit the model.
bcm = nlsBoot(fitm) # Bootstrap the model.
cbind(coef(fitm), confint(bcm)) # Combine coefficients and confidence intervals.

# Constructing a visual of separate fits

# Generate predictions for females.
xf = seq(min(lsf$Final.Age), max(lsf$Final.Age), length.out = 199)
pf = vbTyp(xf, Linf = coef(fitf))

# Generate predictions for males.
xm = seq(min(lsm$Final.Age), max(lsm$Final.Age), length.out = 199)
pm = vbTyp(xm, Linf = coef(fitm))

# Determine the range for x and y axes.
xlmts = range(c(xf, xm))
ylmts = range(c(lsf$FFL.mm, lsm$FFL.mm))

# Plot the observed data for females and males.
plot(FFL.mm ~ Final.Age, data = LaneKnownSex, xlab = "Age", ylab = "FFL (mm)",
     xlim = xlmts, ylim = ylmts, col = "white")
points(FFL.mm ~ Final.Age, data = lsf, pch = 1, col = rgb(0, 0, 0, 1/2), cex = 0.8)
points(FFL.mm ~ Final.Age, data = lsm, pch = 8, col = rgb(0, 0, 0, 1/2), cex = 0.8)

# Add lines for the predicted growth for females and males.
lines(pf ~ xf, lwd = 2, lty = "solid")
lines(pm ~ xm, lwd = 2, lty = "dashed")

# Add a legend to the plot.
legend("bottomright", c("Female", "Male"), pch = c(1, 8), lwd = 2, lty = c("solid", "dashed"), bty = "n", cex = 0.8)

# Display summaries of the fitted models.
summary(fitf)
summary(fitm)

# Logarithmic transformations and models

# Define various logarithmic VBGF equations with different parameters.
lvbLKt = log(FFL) ~ log(Linf[sexnum] * (1 - exp(-K[sexnum] * (Age - t0[sexnum]))))
lvbLK = log(FFL) ~ log(Linf[sexnum] * (1 - exp(-K[sexnum] * (Age - t0))))
lvbLt = log(FFL) ~ log(Linf[sexnum] * (1 - exp(-K * (Age - t0[sexnum]))))
lvbKt = log(FFL) ~ log(Linf * (1 - exp(-K[sexnum] * (Age - t0[sexnum]))))
lvbL = log(FFL) ~ log(Linf[sexnum] * (1 - exp(-K * (Age - t0))))
lvbK = log(FFL) ~ log(Linf * (1 - exp(-K[sexnum] * (Age - t0))))
lvbt = log(FFL) ~ log(Linf * (1 - exp(-K * (Age - t0[sexnum]))))
lvb0 = log(FFL) ~ log(Linf * (1 - exp(-K * (Age - t0))))

# Add log-transformed variables to the dataset.
LaneKnownSex %<>% mutate(logFFL = log(FFL.mm)) # Add log FFL variable.
LaneKnownSex %<>% mutate(logFinalAge = log(Final.Age)) # Add log Final Age variable.

# Fit a nonlinear model to the log-transformed data.
fitTypM = nls(logFFL ~ log(vbTyp(Final.Age, Linf, K, t0)), data = LaneKnownSex, start = svTyp)
FSAmisc::residPlot(fitTypM) # Plot the residuals of the model fit.

# Fit the logarithmic VBGF model with different starting values.
lfitLKt = nls(lvbLKt, data = LaneKnownSex, start = svLKt)
FSAmisc::residPlot(lfitLKt, col = rgb(0, 0, 0, 1/3)) # Plot the residuals.

# Fit the simplest model with logarithmic transformations.
lfit0 = nls(lvb0, data = LaneKnownSex, start = sv0)
lfit0
lfitLKt
fitLKt
fit0

# Perform likelihood ratio test to compare models.
lrt(lfit0, com = lfitLKt, com.name = "All pars differ", sim.names = "No pars differ")

# Perform extra sum of squares test to compare models.
extraSS(lfit0, com = lfitLKt, com.name = "All pars diff", sim.names = "No pars diff")
# Likelihood Ratio (X^2 = 21.419, p < 0.05) indicates significant difference (p < 0.05).
# Extra sum of squares test (F = 7.195, p < 0.5) also indicates significant difference.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
svLK = Map(rep, sv0, c(2, 2, 1))
svLK

# Starting values for other nested models.
svLt = Map(rep, sv0, c(2, 1, 2))
svLt
svKt = Map(rep, sv0, c(1, 2, 2))
svKt
svLK

# Fit the nested models with logarithmic transformations.
lfitLK = nls(lvbLK, data = LaneKnownSex, start = svLK)
lfitLt = nls(lvbLt, data = LaneKnownSex, start = svLt)
lfitKt = nls(lvbKt, data = LaneKnownSex, start = svKt)

# Compare each of the nested models using likelihood ratio test.
lrt(lfitLK, lfitLt, lfitKt, com = lfitLKt,
    com.name = "All pars diff", 
    sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: There is strong evidence that all models fit the data well, but 
# the K,t0 model has the highest log-likelihood ratio and will be used for further model reduction.

# Define starting values for the reduced models.
svt = Map(rep, sv0, c(1, 1, 2))
svt
svK = Map(rep, sv0, c(1, 2, 1))
svK
svL = Map(rep, sv0, c(2, 1, 1))
svL

# Fit the reduced models.
lfitt = nls(lvbt, data = LaneKnownSex, start = svt)
lfitK = nls(lvbK, data = LaneKnownSex, start = svK)
lfitL = nls(lvbL, data = LaneKnownSex, start = svL)

# Compare the reduced models to the full model.
lrt(lfitK, lfitt, com = lfitKt, com.name = "K,t0 dif", sim.names = c("K dif", "t0 dif"))
# Interpretation: Neither K (p = 0.0149) nor t0 (p < 0.05) fit as well as the full model.

lrt(lfit0, com = lfitL, com.name = "Linf dif", sim.names = c("No pars dif"))
# Interpretation: The comparison of Linf to the full model confirms that a difference between groups exists for Linf (p < 0.05).

# Model Selection with AIC or BIC
model_comparison <- cbind(
  AIC(lfitLKt, lfitLK, lfitLt, lfitKt, lfitL, lfitK, lfitt, lfit0),
  BIC(lfitLKt, lfitLK, lfitLt, lfitKt, lfitL, lfitK, lfitt, lfit0)
)
model_comparison
# Interpretation: The lfitKt model has the lowest AIC (-364.5202), and the lfitL model has the lowest BIC (-341.3835), indicating the best models using these criteria.

# Plot the best fit model
# Predict values using the lfitKt model.
predict(lfitKt)
summary(LaneKnownSex$Final.Age)

# Plot observed data with the best fit model.
plot(LaneKnownSex$Final.Age, LaneKnownSex$FFL.mm, main = "Fork Length at Fractional Age for Lane Snapper", font.main = 1,
     xlab = "Fractional Age (years)", ylab = "Fork Length (mm)", col = "white")

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Plot observed data points for females and males.
points(FFL.mm ~ Final.Age, data = lsf, pch = 1, col = rgb(0, 0, 0, 1/2), cex = 0.8) # Female points 
points(FFL.mm ~ Final.Age, data = lsm, pch = 8, col = rgb(0, 0, 0, 1/2), cex = 0.8) # Male points

# Add VBGF lines for females and males.
lines(seq(0, 14, 0.1), vonBert(coef(lfitKt)[1], coef(lfitKt)[4], coef(lfitKt)[2], seq(0, 14, 0.1)), col = "black", lwd = 3, lty = 2) # Female
lines(seq(0, 14, 0.1), vonBert(coef(lfitKt)[1], coef(lfitKt)[5], coef(lfitKt)[3], seq(0, 14, 0.1)), col = "black", lwd = 2, lty = 1) # Male

# Add a legend to the plot.
legend("bottomright", c("Female", "Male", "Female Predicted", "Male Predicted"),
       pch = c(1, 8, NA, NA), lty = c(NA, NA, 2, 1), bty = "n", cex = 0.8, lwd = c(1, 1, 2, 2))

# Display the coefficients of the best fit model.
coef(lfitKt)

# 1 is female, 2 is male
predict(lfitKt) # Predict values using the lfitKt model.

# Plot observed data.
plot(LaneKnownSex$Final.Age, LaneKnownSex$FFL.mm)

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Add VBGF lines for females (red) and males (blue).
lines(seq(0, 14, 0.1), vonBert(coef(lfitKt)[1], coef(lfitKt)[4], coef(lfitKt)[2], seq(0, 14, 0.1)), col = 2) # Female
lines(seq(0, 14, 0.1), vonBert(coef(lfitKt)[1], coef(lfitKt)[5], coef(lfitKt)[3], seq(0, 14, 0.1)), col = "blue") # Male
coef(lfitL)

# Summarizing the Model Fit 12.4.3 
# The code below creates a function for the typical VBGF, creates separate subsets of
# female and male Lane Snapper, and follows methods from 12.3.3 to fit the model and 
# extract parameter estimates and bootstrap CI for each sex.

vbTyp = vbFuns("typical")

## Females
lsf = dplyr::filter(LaneKnownSex, Macro.sex == "F") # Subset data for females.
fitf = nls(logFFL ~ log(vbTyp(Final.Age, Linf, K, t0)), data = lsf, start = sv0) # Fit the model.
bcf = nlsBoot(fitf) # Bootstrap the model.
cbind(coef(fitf), confint(bcf)) # Combine coefficients and confidence intervals.

## Males
lsm = dplyr::filter(LaneKnownSex, Macro.sex == "M") # Subset data for males.
fitm = nls(logFFL ~ vbTyp(Final.Age, Linf, K, t0), data = lsm, start = sv0) # Fit the model.
bcm = nlsBoot(fitm) # Bootstrap the model.
cbind(coef(fitm), confint(bcm)) # Combine coefficients and confidence intervals.

# Constructing a visual of separate fits

# Generate predictions for females.
xf = seq(min(lsf$Final.Age), max(lsf$Final.Age), length.out = 199)
pf = vbTyp(xf, Linf = coef(fitf))

# Generate predictions for males.
xm = seq(min(lsm$Final.Age), max(lsm$Final.Age), length.out = 199)
pm = vbTyp(xm, Linf = coef(fitm))

# Determine the range for x and y axes.
xlmts = range(c(xf, xm))
ylmts = range(c(lsf$FFL.mm, lsm$FFL.mm))

# Plot the observed data for females and males.
plot(FFL.mm ~ Final.Age, data = LaneKnownSex, xlab = "Age", ylab = "FFL (mm)",
     xlim = xlmts, ylim = ylmts, col = "white")
points(FFL.mm ~ Final.Age, data = lsf, pch = 1, col = rgb(0, 0, 0, 1/2), cex = 0.8)
points(FFL.mm ~ Final.Age, data = lsm, pch = 8, col = rgb(0, 0, 0, 1/2), cex = 0.8)

# Add lines for the predicted growth for females and males.
lines(pf ~ xf, lwd = 2, lty = "solid")
lines(pm ~ xm, lwd = 2, lty = "dashed")

# Add a legend to the plot.
legend("bottomright", c("Female", "Male"), pch = c(1, 8), lwd = 2, lty = c("solid", "dashed"), bty = "n", cex = 0.8)

# Display summaries of the fitted models.
summary(fitf)
summary(fitm)

# Examining Differences Between Eastern vs Western GOM

# Get a list of starting values.
svTyp = vbStarts(FFL.mm ~ Final.Age, data = LaneAll) # Generate starting values.
svTyp
summary(LaneAll$FFL.mm)

vbStarts(FFL.mm ~ Final.Age, data = LaneAll, vbStartsDP = TRUE)

max(LaneAll$FFL.mm)

# Define the VBGF function.
vbTyp = function(Final.Age, Linf, K, t0) Linf * (1 - exp(-K * (Final.Age - t0)))
vbTyp(2, Linf = 319, K = 0.63, t0 = -0.08)

# Create the VBGF function using vbFuns().
vbTyp = vbFuns()

vbTyp(3, Linf = 495, K = 0.67, t0 = -0.88) # Example prediction with specific parameters.

# Fit the VBGF model to data and estimate parameters and confidence intervals.
fitTyp = nls(FFL.mm ~ vbTyp(Final.Age, Linf, K, t0), data = LaneAll, start = svTyp)
fitTyp
coef(fitTyp)
confint(fitTyp)

# Bootstrapping the nonlinear model.
bootTyp = nlsBoot(fitTyp)
headtail(bootTyp$coefboot, n = 2)
confint(bootTyp, plot = TRUE)

# Display more results for the model fit.
summary(fitTyp, correlation = TRUE)

## Predictions 12.3.4

# Predicted lengths at specific ages using the fitted model.
nd = data.frame(Final.Age = c(1, 3, 5, 9, 12, 14))
predict(fitTyp, nd) # Predict lengths for the specified ages using the fitTyp model.

# Predict lengths using the VBGF function and coefficients from the model.
vbTyp(c(1, 3, 5, 9, 12, 14), coef(fitTyp))
vbTyp(3, bootTyp$coefboot[1,]) # Prediction at age 3 using the first set of bootstrapped coefficients.

# Apply the VBGF function to each set of bootstrapped coefficients at age 3.
p3Typ = apply(bootTyp$coefboot, MARGIN = 1, FUN = vbTyp, t = 3)
p3Typ[1:6] # Display the first six predictions.

# Calculate the 2.5th and 97.5th percentiles of the bootstrapped predictions.
quantile(p3Typ, c(0.025, 0.975))

## Visualizing Model Fit 12.3.5

# Generate a sequence of ages for predictions.
x = seq(0, 15, length.out = 199)

# Predicted lengths using the VBGF and the model coefficients.
pTyp = vbTyp(x, Linf = coef(fitTyp))

# Determine the range for the x and y axes.
xlmts = range(c(x, LaneAll$Final.Age))
ylmts = range(c(pTyp, LaneAll$FFL.mm))

# Plot the observed data and the predicted lengths.
plot(FFL.mm ~ Final.Age, data = LaneAll, xlab = "Age", ylab = "Full Fork Length (mm)",
     xlim = xlmts, ylim = ylmts, pch = 19, col = rgb(0, 0, 0, 1/3))
lines(pTyp, lwd = 2)

# Initialize vectors to store the lower and upper confidence intervals.
LCI = UCI = numeric(length(x))

# Calculate the confidence intervals for each age.
for (i in 1:length(x)) {
  tmp = apply(bootTyp$coefboot, MARGIN = 1, FUN = vbTyp, t = x[i])
  LCI[i] = quantile(tmp, 0.025)
  UCI[i] = quantile(tmp, 0.975)
}

# Update the y-axis limits to include the confidence intervals.
ylmts = range(c(pTyp, LCI, UCI, LaneAll$FFL.mm))

# Plot the observed data, predicted lengths, and confidence intervals.
plot(FFL.mm ~ Final.Age, data = LaneAll, xlab = "Age", ylab = "Full Fork Length (mm)",
     xlim = xlmts, ylim = ylmts, pch = 19, col = rgb(0, 0, 0, 1/3))
lines(pTyp ~ x, lwd = 2)
lines(UCI ~ x, lwd = 2, lty = "dashed")
lines(LCI ~ x, lwd = 2, lty = "dashed")

# Add a log-transformed Full Fork Length variable to the dataset.
LaneAll %<>% mutate(logFFL = log(FFL.mm))

# Fit a nonlinear model to the log-transformed data.
fitTypM = nls(logFFL ~ log(vbTyp(Final.Age, Linf, K, t0)), data = LaneAll, start = svTyp)

# Plot the residuals of both the original and log-transformed model fits.
FSAmisc::residPlot(fitTypM)
FSAmisc::residPlot(fitTyp)

## Model Fitting 12.4.2

# Subset data to exclude unknown sex and create a numeric variable for sex.
LaneKnownSex <- LaneAll[LaneAll$Macro.sex != "U",]
LaneKnownSex$sexnum <- ifelse(LaneKnownSex$Macro.sex == "F", 1, 2)
table(LaneKnownSex$sexnum, LaneKnownSex$Macro.sex)
sex = LaneKnownSex$Macro.sex

# Subset data to exclude unknown GulfSide and Atlantic.
LKGS = LaneAll[!(LaneAll$GulfSide == "U" | LaneAll$GulfSide == "Atlantic"),]
LKGS$GSnum <- ifelse(LKGS$GulfSide == "East", 1, 2)
table(LKGS$GulfSide)
table(LKGS$State)

# Extract age and length data.
Age = as.numeric(LKGS$Final.Age)
FFL = as.numeric(LKGS$FFL.mm)

# Define various von Bertalanffy growth function (VBGF) equations with different parameters.
GSvbLKt = FFL ~ Linf[GSnum] * (1 - exp(-K[GSnum] * (Age - t0[GSnum])))
GSvbLK = FFL ~ Linf[GSnum] * (1 - exp(-K[GSnum] * (Age - t0)))
GSvbLt = FFL ~ Linf[GSnum] * (1 - exp(-K * (Age - t0[GSnum])))
GSvbKt = FFL ~ Linf * (1 - exp(-K[GSnum] * (Age - t0[GSnum])))
GSvbL = FFL ~ Linf[GSnum] * (1 - exp(-K * (Age - t0)))
GSvbK = FFL ~ Linf * (1 - exp(-K[GSnum] * (Age - t0)))
GSvbt = FFL ~ Linf * (1 - exp(-K * (Age - t0[GSnum])))
GSvb0 = FFL ~ Linf * (1 - exp(-K * (Age - t0)))

# Get initial parameter estimates for nonlinear models.
GSsv0 = vbStarts(FFL ~ Age, data = LKGS)
GSsv0

# Replicate starting values for different groups.
GSsvLKt = Map(rep, GSsv0, c(2, 2, 2))
GSsvLKt

# Fit the nonlinear model with the specified starting values.
GSfitLKt = nls(GSvbLKt, data = LKGS, start = GSsvLKt)

# Plot the residuals of the fitted model.
FSAmisc::residPlot(GSfitLKt, col = rgb(0, 0, 0, 1/3))

# Fit the simplest model.
GSfit0 = nls(GSvb0, data = LKGS, start = GSsv0)

# Perform a likelihood ratio test to compare models.
lrt(GSfit0, com = GSfitLKt, com.name = "All pars differ", sim.names = "No pars differ")
extraSS(GSfit0, com = GSfitLKt, com.name = "All pars diff", sim.names = "No pars diff")
# Interpretation: Likelihood Ratio (X^2 = 4.62, p = 0.2015) and Extra sum of squares test (F = 1.53, p = 0.2059) indicate no significant difference.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
GSsvLK = Map(rep, GSsv0, c(2, 2, 1))
GSsvLK

# Starting values for other nested models.
GSsvLt = Map(rep, GSsv0, c(2, 1, 2))
GSsvLt
GSsvKt = Map(rep, GSsv0, c(1, 2, 2))
GSsvKt

# Fit the nested models.
GSfitLK = nls(GSvbLK, data = LKGS, start = GSsvLK)
GSfitLt = nls(GSvbLt, data = LKGS, start = GSsvLt)
GSfitKt = nls(GSvbKt, data = LKGS, start = GSsvKt)

# Compare each of the nested models using likelihood ratio test.
lrt(GSfitLK, GSfitLt, GSfitKt, com = GSfitLKt, com.name = "All pars diff", sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: Evidence that (Linf, K) (p = 0.1321) and (Linf, t0) (p = 0.062) fit the data as well as (Linf, K, t0). The (Linf, K) model is used as the basis due to the greatest log-likelihood value (-2456.6192).

# Further starting values for additional comparisons.
GSsvL = Map(rep, GSsv0, c(2, 1, 1))
GSsvL
GSsvK = Map(rep, GSsv0, c(1, 2, 1))
GSsvK

# Fit additional models.
GSfitL = nls(GSvbL, data = LKGS, start = GSsvL)
summary(GSfitL)
GSfitK = nls(GSvbK, data = LKGS, start = GSsvK)

# Compare additional models.
lrt(GSfitL, GSfitK, com = GSfitLK, com.name = "Linf,K diff", sim.names = c("Linf diff", "K diff"))
# Interpretation: Results suggest that Linf (p = 0.129) and K (p = 0.099) fit as well as (Linf, K).

lrt(GSfit0, com = GSfitL, com.name = "Linf diff", sim.names = "No pars diff")
# Interpretation: Comparison confirms that a difference between groups is not evident for Linf (p = 0.49).

lrt(GSfit0, com = GSfitK, com.name = "K diff", sim.names = "No pars diff")
# Interpretation: Comparison confirms that a difference between groups is not evident for K (p = 0.82).

# Model Selection with AIC or BIC

# Replicate starting values for the final model.
GSsvt = Map(rep, GSsv0, c(1, 1, 2))
GSfitt = nls(GSvbt, data = LKGS, start = GSsvt)

# Compare models using AIC and BIC.
model_comparison <- cbind(
  AIC(GSfitLKt, GSfitLK, GSfitLt, GSfitKt, GSfitL, GSfitK, GSfitt, GSfit0),
  BIC(GSfitLKt, GSfitLK, GSfitLt, GSfitKt, GSfitL, GSfitK, GSfitt, GSfit0)
)
model_comparison
# Interpretation: The Omega model has the lowest AIC (4924.004) and BIC (4940.624), indicating the best model using either criterion.

## LOG Model Fitting 12.4.2 for Gulf Comparisons

# Subset data to exclude unknown sex and create a numeric variable for sex.
LaneKnownSex <- LaneAll[LaneAll$Macro.sex != "U",]
LaneKnownSex$sexnum <- ifelse(LaneKnownSex$Macro.sex == "F", 1, 2)
table(LaneKnownSex$sexnum, LaneKnownSex$Macro.sex)
sex = LaneKnownSex$Macro.sex

# Subset data to exclude unknown GulfSide and Atlantic.
LKSaGS = LaneKnownSex[!(LaneKnownSex$GulfSide == "U" | LaneKnownSex$GulfSide == "Atlantic"),]
LKSaGS$GSnum <- ifelse(LKSaGS$GulfSide == "East", 1, 2)
table(LKSaGS$GSnum)

# Extract age and length data.
Age = as.numeric(LKSaGS$Final.Age)
FFL = as.numeric(LKSaGS$FFL.mm)

# Define various logarithmic von Bertalanffy growth function (VBGF) equations with different parameters.
lGSvbLKt = log(FFL) ~ log(Linf[GSnum] * (1 - exp(-K[GSnum] * (Age - t0[GSnum]))))
lGSvbLK = log(FFL) ~ log(Linf[GSnum] * (1 - exp(-K[GSnum] * (Age - t0))))
lGSvbLt = log(FFL) ~ log(Linf[GSnum] * (1 - exp(-K * (Age - t0[GSnum]))))
lGSvbKt = log(FFL) ~ log(Linf * (1 - exp(-K[GSnum] * (Age - t0[GSnum]))))
lGSvbL = log(FFL) ~ log(Linf[GSnum] * (1 - exp(-K * (Age - t0))))
lGSvbK = log(FFL) ~ log(Linf * (1 - exp(-K[GSnum] * (Age - t0))))
lGSvbt = log(FFL) ~ log(Linf * (1 - exp(-K * (Age - t0[GSnum]))))
lGSvb0 = log(FFL) ~ log(Linf * (1 - exp(-K * (Age - t0))))

# Plot the residuals of the original log-transformed model.
FSAmisc::residPlot(fitTypM)

# Get initial parameter estimates for nonlinear models.
GSsv0 = vbStarts(FFL ~ Age, data = LKSaGS)
GSsv0

# Replicate starting values for different groups.
lGSsvLKt = Map(rep, GSsv0, c(2, 2, 2))
lGSsvLKt

# Fit the nonlinear model with the specified starting values.
lGSfitLKt = nls(lGSvbLKt, data = LKSaGS, start = lGSsvLKt)

# Plot the residuals of the fitted model.
FSAmisc::residPlot(lGSfitLKt, col = rgb(0, 0, 0, 1/3))

# Fit the simplest model.
lGSfit0 = nls(lGSvb0, data = LKSaGS, start = GSsv0)

# Perform a likelihood ratio test to compare models.
lrt(lGSfit0, com = lGSfitLKt, com.name = "All pars differ", sim.names = "No pars differ")
extraSS(lGSfit0, com = lGSfitLKt, com.name = "All pars diff", sim.names = "No pars diff")
# Interpretation: Likelihood Ratio (X^2 = 6.26, p = 0.09944) and Extra sum of squares test (F = 2.0758, p = 0.1025) indicate no significant difference.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
GsvLK = Map(rep, GSsv0, c(2, 2, 1))
GsvLK

# Starting values for other nested models.
GsvLt = Map(rep, GSsv0, c(2, 1, 2))
GsvLt
GsvKt = Map(rep, GSsv0, c(1, 2, 2))
GsvKt

# Fit the nested models.
lGSfitLK = nls(lGSvbLK, data = LKSaGS, start = GsvLK)
lGSfitLt = nls(lGSvbLt, data = LKSaGS, start = GsvLt)
lGSfitKt = nls(lGSvbKt, data = LKSaGS, start = GsvKt)

# Compare each of the nested models using likelihood ratio test.
lrt(lGSfitLK, lGSfitLt, lGSfitKt, com = lGSfitLKt, com.name = "All pars diff", sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: Evidence that (Linf, K) (p = 0.03971), (Linf, t0) (p = 0.02114), and (K, t0) (p = 0.03044) are all statistically different from (Linf, K, t0) and not considered better. The (Linf, t0) model has the greatest log-likelihood value.

# Further starting values for additional comparisons.
GSsvL = Map(rep, GSsv0, c(2, 1, 1))
GSsvL
GSsvt = Map(rep, GSsv0, c(1, 1, 2))
GSsvt
GSsvK = Map(rep, GSsv0, c(1, 2, 1))
GSsvK

# Fit additional models.
lGSfitK = nls(lGSvbK, data = LKSaGS, start = GSsvK)
lGSfitt = nls(lGSvbt, data = LKSaGS, start = GSsvt)
lGSfitL = nls(lGSvbL, data = LKSaGS, start = GSsvL)

# Compare additional models.
lrt(lGSfitL, lGSfitt, com = lGSfitLt, com.name = "Linf,t dif", sim.names = c("Linf dif", "K dif"))
# Interpretation: Results suggest that Linf (p = 0.46) and K (p = 0.79) fit as well as (Linf, t).

lrt(lGSfit0, com = lGSfitL, com.name = "Linf dif", sim.names = "No pars dif")
# Interpretation: Comparison confirms that a difference between groups is not evident for Linf (p = 0.52).

lrt(lGSfit0, com = lGSfitK, com.name = "K dif", sim.names = "No pars dif")
# Interpretation: Comparison confirms that a difference between groups is not evident for K (p = 0.28).

# Model Selection with AIC or BIC

# Replicate starting values for the final model.
GSsvt = Map(rep, GSsv0, c(1, 1, 2))
lGSfitt = nls(lGSvbt, data = LKSaGS, start = GSsvt)

# Compare models using AIC and BIC.
model_comparison <- cbind(
  AIC(lGSfitLKt, lGSfitLK, lGSfitLt, lGSfitKt, lGSfitL, lGSfitK, lGSfitt, lGSfit0),
  BIC(lGSfitLKt, lGSfitLK, lGSfitLt, lGSfitKt, lGSfitL, lGSfitK, lGSfitt, lGSfit0)
)
model_comparison
# Interpretation: The lGSfitLKt model has the lowest AIC (-495.8272) and the lGSfit0 has the lowest BIC (-478.7853), indicating the best models using either criterion.

# Display the coefficients of the best fit models.
coef(lGSfitLKt)
coef(lGSfit0) # No parameters different.

## Regular Examining Differences Among Capture Years

# Display frequency tables for Catch Year and State.
table(LaneAll$Catch.Year)
table(LaneAll$State)

# Create a copy of the dataset and convert Catch Year to numeric.
LaneKnownCY <- LaneAll
LaneKnownCY$YRnum = as.numeric(factor(LaneKnownCY$Catch.Year)) # Converts years to numeric and creates a new column.

# Combine 2020 and 2022 samples with 2021 catch year.
LaneKnownCY[LaneKnownCY$Catch.Year == 2022,]$Catch.Year <- 2021
LaneKnownCY[LaneKnownCY$Catch.Year == 2020,]$Catch.Year <- 2021
LaneKnownCY$YRnum = as.numeric(factor(LaneKnownCY$Catch.Year)) # Recreate the numeric column after changes.

# Display frequency tables for the modified dataset.
table(LaneKnownCY$YRnum)
table(LaneKnownCY$Catch.Year)
table(LaneKnownCY$State)

# Get a list of starting values for the von Bertalanffy growth function (VBGF) model.
YRsvTyp = vbStarts(FFL.mm ~ Final.Age, data = LaneKnownCY)
YRsvTyp = list(Linf = max(LaneKnownCY$FFL.mm, na.rm = TRUE), K = 0.67, t0 = -0.88)
YRsvTyp

## Model Fitting 12.4.2

# Extract age and length data.
Age = as.numeric(LaneKnownCY$Final.Age)
FFL = as.numeric(LaneKnownCY$FFL.mm)

# Define various VBGF equations with different parameters for each capture year.
YRvbLKt = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - t0[YRnum])))
YRvbLK = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - t0)))
YRvbLt = FFL ~ Linf[YRnum] * (1 - exp(-K * (Age - t0[YRnum])))
YRvbKt = FFL ~ Linf * (1 - exp(-K[YRnum] * (Age - t0[YRnum])))
YRvbL = FFL ~ Linf[YRnum] * (1 - exp(-K * (Age - t0)))
YRvbK = FFL ~ Linf * (1 - exp(-K[YRnum] * (Age - t0)))
YRvbt = FFL ~ Linf * (1 - exp(-K * (Age - t0[YRnum])))
YRvb0 = FFL ~ Linf * (1 - exp(-K * (Age - t0)))

# Get initial parameter estimates for nonlinear models.
YRsv0 = vbStarts(FFL ~ Age, data = LaneKnownCY)
YRsv0

# Replicate starting values for different groups.
YRsvLKt = Map(rep, YRsv0, c(4, 4, 4))
YRsvLKt

# Fit the nonlinear model with the specified starting values.
YRfitLKt = nls(YRvbLKt, data = LaneKnownCY, start = YRsvLKt)

# Plot the residuals of the fitted model.
FSAmisc::residPlot(YRfitLKt, col = rgb(0, 0, 0, 1/3))

# Fit the simplest model.
YRfit0 = nls(YRvb0, data = LaneKnownCY, start = YRsv0)

# Perform a likelihood ratio test to compare models.
lrt(YRfit0, com = YRfitLKt, com.name = "All pars differ", sim.names = "No pars differ")
extraSS(YRfit0, com = YRfitLKt, com.name = "All pars diff", sim.names = "No pars diff")
# Interpretation: Likelihood Ratio (X^2 = 159.16, p < 2.2e-16) and Extra sum of squares test (F = 19.81, p < 2.2e-16) indicate significant differences.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
YRsvLK = Map(rep, YRsv0, c(4, 4, 1))
YRsvLK

# Starting values for other nested models.
YRsvLt = Map(rep, YRsv0, c(4, 1, 4))
YRsvLt
YRsvKt = Map(rep, YRsv0, c(1, 4, 4))
YRsvKt

# Fit the nested models.
YRfitLK = nls(YRvbLK, data = LaneKnownCY, start = YRsvLK)
YRfitLt = nls(YRvbLt, data = LaneKnownCY, start = YRsvLt)
YRfitKt = nls(YRvbKt, data = LaneKnownCY, start = YRsvKt)

# Compare each of the nested models using likelihood ratio test.
lrt(YRfitLK, YRfitLt, YRfitKt, com = YRfitLKt, com.name = "All pars diff", sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: Evidence that (K,t0) (p = 0.06) and (Linf,K) (p = 0.06) fit the data as well as (Linf, K, t0), but not (Linf,t0) (p = 0.003). The (K,t0) model is used for further reduction.

# Further starting values for additional comparisons.
YRsvK = Map(rep, YRsv0, c(1, 4, 1))
YRsvK
YRsvt = Map(rep, YRsv0, c(1, 1, 4))
YRsvt

# Fit additional models.
YRfitK = nls(YRvbK, data = LaneKnownCY, start = YRsvK)
summary(YRfitK)
YRfitt = nls(YRvbt, data = LaneKnownCY, start = YRsvt)

# Compare additional models.
lrt(YRfitL, YRfitK, com = YRfitLK, com.name = "L,K diff", sim.names = c("L diff", "K diff"))
# Interpretation: Results suggest that neither K nor t0 (p = 2.844e-06) fit as well as (K,t0).

lrt(YRfit0, com = YRfitK, com.name = "K diff", sim.names = "No pars diff")
# Interpretation: Comparison confirms that a difference between groups exists and that difference is only evident for K (p < 2.2e-16).

# Model Selection with AIC or BIC

# Replicate starting values for the final model.
YRsvL = Map(rep, YRsv0, c(4, 1, 1))
YRsvL
YRfitL = nls(YRvbL, data = LaneKnownCY, start = YRsvL)

# Compare models using AIC and BIC.
model_comparison <- cbind(
  AIC(YRfitLKt, YRfitLK, YRfitLt, YRfitKt, YRfitL, YRfitK, YRfitt, YRfit0),
  BIC(YRfitLKt, YRfitLK, YRfitLt, YRfitKt, YRfitL, YRfitK, YRfitt, YRfit0)
)
model_comparison
# Interpretation: The YRfitKt model has the lowest AIC (6083.285) and YRfitK has the lowest BIC (6115.396).

# Plot for best fit based on AIC
predict(YRfitKt) # Predict values using the YRfitKt model.
coef(YRfitKt) # Extract coefficients from the YRfitKt model.

# Plot observed data.
plot(LKSaYR$Final.Age, LKSaYR$FFL.mm)

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Add VBGF lines for each year group using coefficients from the YRfitKt model.
# Group 1
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[6], coef(YRfitKt)[2], seq(0, 14, 0.1)), col = 2) 
# Group 2
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[7], coef(YRfitKt)[3], seq(0, 14, 0.1)), col = 3) 
# Group 3
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[8], coef(YRfitKt)[4], seq(0, 14, 0.1)), col = 4) 
# Group 4
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[9], coef(YRfitKt)[5], seq(0, 14, 0.1)), col = 5) 

# Display frequency table for year groups.
table(LKSaYR$YRnum)

# Plot for best fit BIC
predict(YRfitK) # Predict values using the YRfitK model.
coef(YRfitK) # Extract coefficients from the YRfitK model.

# Plot observed data.
plot(LKSaYR$Final.Age, LKSaYR$FFL.mm)

# Add VBGF lines for each year group using coefficients from the YRfitK model.
# Group 1
lines(seq(0, 14, 0.1), vonBert(coef(YRfitK)[1], coef(YRfitK)[6], coef(YRfitK)[2], seq(0, 14, 0.1)), col = 2) 
# Group 2
lines(seq(0, 14, 0.1), vonBert(coef(YRfitK)[1], coef(YRfitK)[6], coef(YRfitK)[3], seq(0, 14, 0.1)), col = 3) 
# Group 3
lines(seq(0, 14, 0.1), vonBert(coef(YRfitK)[1], coef(YRfitK)[6], coef(YRfitK)[4], seq(0, 14, 0.1)), col = 4) 
# Group 4
lines(seq(0, 14, 0.1), vonBert(coef(YRfitK)[1], coef(YRfitK)[6], coef(YRfitK)[5], seq(0, 14, 0.1)), col = 5) 

## Regular Examining Differences Among Capture Years with FIXED t0

# Display frequency tables for Catch Year.
table(LaneAll$Catch.Year)

# Create a copy of the dataset and convert Catch Year to numeric.
LaneKnownCY <- LaneAll
LaneKnownCY$YRnum = as.numeric(factor(LaneKnownCY$Catch.Year)) # Converts years to numeric and creates a new column.

# Combine 2020 and 2022 samples with 2021 catch year.
LaneKnownCY[LaneKnownCY$Catch.Year == 2022,]$Catch.Year <- 2021
LaneKnownCY[LaneKnownCY$Catch.Year == 2020,]$Catch.Year <- 2021
LaneKnownCY$YRnum = as.numeric(factor(LaneKnownCY$Catch.Year)) # Recreate the numeric column after changes.

# Subset for only adult hook and line samples.
LaneKnownCYAdults = LaneKnownCY[LaneKnownCY$Final.Age > 1,]  # Subset for samples older than 1 year.
LaneKnownCYAdultsHL = LaneKnownCYAdults[LaneKnownCYAdults$Gear.Code == "HL",]

# Display frequency tables for the modified dataset.
table(LaneKnownCYAdults$YRnum)
table(LaneKnownCYAdults$Catch.Year)
table(LaneKnownCYAdults$State)

# Get a list of starting values for the VBGF model with fixed t0.
YRsvTyp = vbStarts(FFL.mm ~ Final.Age, data = LaneKnownCYAdults, fixed = list(t0 = 0))
YRsvTyp = list(Linf = max(LaneKnownCYAdults$FFL.mm, na.rm = TRUE), K = 0.67, t0 = 0)
YRsvTyp

## Model Fitting 12.4.2

# Extract age and length data.
Age = as.numeric(LaneKnownCYAdults$Final.Age)
FFL = as.numeric(LaneKnownCYAdults$FFL.mm)

# Define various VBGF equations with different parameters for each capture year with fixed t0.
YRvbLKt0 = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbLK0 = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbLt0 = FFL ~ Linf[YRnum] * (1 - exp(-K * (Age - 0)))
YRvbKt0 = FFL ~ Linf * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbL0 = FFL ~ Linf[YRnum] * (1 - exp(-K * (Age - 0)))
YRvbK0 = FFL ~ Linf * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbt0 = FFL ~ Linf * (1 - exp(-K * (Age - 0)))
YRvb0 = FFL ~ Linf * (1 - exp(-K * (Age - 0)))
YRvbt0 = FFL ~ Linf * (1 - exp(-K * (Age - 0)))

# Get initial parameter estimates for nonlinear models with fixed t0.
YRsv0 = vbStarts(FFL ~ Age, data = LaneKnownCYAdults, fixed = list(t0 = 0))
YRsv0

# Remove t0 from the starting values list.
YRsvt0 = YRsv0
YRsvt0$t0 = NULL

# Replicate starting values for different groups.
YRsvLKt0 = Map(rep, YRsv0, c(4, 4, 4))
YRsvLKt0$t0 = NULL

# Fit the nonlinear model with the specified starting values and fixed t0.
YRfitLKt0 = nls(YRvbLKt0, data = LaneKnownCYAdults, start = YRsvLKt0)
YRfitLKt0
FSAmisc::residPlot(YRfitLKt0, col = rgb(0, 0, 0, 1/3))

# Fit the nonlinear model with the specified starting values and fixed t0.
YRfitt0 = nls(YRvbt0, data = LaneKnownCYAdults, start = YRsvt0)
YRfitt0

# Perform a likelihood ratio test to compare models.
lrt(YRfitt0, com = YRfitLKt0, com.name = "All pars differ", sim.names = "No pars differ")
extraSS(YRfitt0, com = YRfitLKt0, com.name = "All pars diff", sim.names = "No pars diff")
# Interpretation: Likelihood Ratio (X^2 = 159.16, p < 2.2e-16) and Extra sum of squares test (F = 19.81, p < 2.2e-16) indicate significant differences.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
YRsvLK = Map(rep, YRsv0, c(4, 4, 1))
YRsvLK

# Remove t0 from the starting values list for specific models.
YRsvLK0$t0 = NULL
YRsvLK0

# Starting values for other nested models.
YRsvLt = Map(rep, YRsv0, c(4, 1, 4))
YRsvLt
YRsvLt0$t0 = NULL
YRsvLt0

YRsvKt = Map(rep, YRsv0, c(1, 4, 4))
YRsvKt
YRsvKt$t0 = NULL
YRsvLKt0

# Fit the nested models.
YRfitLK0 = nls(YRvbLK0, data = LaneKnownCYAdults, start = YRsvLK0)
YRfitLt0 = nls(YRvbLt0, data = LaneKnownCYAdults, start = YRsvLt0)
YRfitKt0 = nls(YRvbKt0, data = LaneKnownCYAdults, start = YRsvKt0)
YRfitKt0

# Compare each of the nested models using likelihood ratio test.
lrt(YRfitLK, YRfitLt, YRfitKt, com = YRfitLKt0, com.name = "All pars diff", sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: Evidence that (K,t0) (p = 0.06) and (Linf,K) (p = 0.06) fit the data as well as (Linf, K, t0), but not (Linf,t0) (p = 0.003). The (K,t0) model is used for further reduction.

# Further starting values for additional comparisons.
YRsvK = Map(rep, YRsv0, c(1, 4, 1))
YRsvK
YRsvt = Map(rep, YRsv0, c(1, 1, 4))
YRsvt

# Fit additional models.
YRfitK = nls(YRvbK, data = LaneKnownCYAdults, start = YRsvK)
summary(YRfitK)
YRfitt = nls(YRvbt, data = LaneKnownCYAdults, start = YRsvt)

# Compare additional models.
lrt(YRfitL, YRfitK, com = YRfitLK, com.name = "L,K diff", sim.names = c("L diff", "K diff"))
# Interpretation: Results suggest that neither K nor t0 (p = 2.844e-06) fit as well as (K,t0).

lrt(YRfit0, com = YRfitK, com.name = "K diff", sim.names = "No pars diff")
# Interpretation: Comparison confirms that a difference between groups exists and that difference is only evident for K (p < 2.2e-16).

# Model Selection with AIC or BIC

# Replicate starting values for the final model.
YRsvL = Map(rep, YRsv0, c(4, 1, 1))
YRsvL
YRfitL = nls(YRvbL, data = LaneKnownCYAdults, start = YRsvL)

# Compare models using AIC and BIC.
model_comparison <- cbind(
  AIC(YRfitLKt, YRfitLK, YRfitLt, YRfitKt, YRfitL, YRfitK, YRfitt, YRfit0),
  BIC(YRfitLKt, YRfitLK, YRfitLt, YRfitKt, YRfitL, YRfitK, YRfitt, YRfit0)
)
model_comparison
# Interpretation: The YRfitKt model has the lowest AIC (6083.285) and YRfitK has the lowest BIC (6115.396).

# Plot for best fit based on AIC
predict(YRfitKt) # Predict values using the YRfitKt model.
coef(YRfitKt) # Extract coefficients from the YRfitKt model.

# Plot observed data.
plot(LKSaYR$Final.Age, LKSaYR$FFL.mm)

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Add VBGF lines for each year group using coefficients from the YRfitKt model.
# Group 1
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[6], coef(YRfitKt)[2], seq(0, 14, 0.1)), col = 2) 
# Group 2
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[7], coef(YRfitKt)[3], seq(0, 14, 0.1)), col = 3) 
# Group 3
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[8], coef(YRfitKt)[4], seq(0, 14, 0.1)), col = 4) 
# Group 4
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[9], coef(YRfitKt)[5], seq(0, 14, 0.1)), col = 5) 

# Display frequency table for year groups.
table(LKSaYR$YRnum)

# Plot for best fit BIC
predict(YRfitK) # Predict values using the YRfitK model.
coef(YRfitK) # Extract coefficients from the YRfitK model.

# Plot observed data.
plot(LKSaYR$Final.Age, LKSaYR$FFL.mm)

# Add VBGF lines for each year group using coefficients from the YRfitK model.
# Group 1
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[6], coef(YRfitKt)[2], seq(0, 14, 0.1)), col = 2) 
# Group 2
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[6], coef(YRfitKt)[3], seq(0, 14, 0.1)), col = 3) 
# Group 3
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[6], coef(YRfitKt)[4], seq(0, 14, 0.1)), col = 4) 
# Group 4
lines(seq(0, 14, 0.1), vonBert(coef(YRfitKt)[1], coef(YRfitKt)[6], coef(YRfitKt)[5], seq(0, 14, 0.1)), col = 5)

## LOG-Examining Differences Among Capture Years
LKSaYR <- LaneAll[LaneAll$Macro.sex != "U",] # Subset for samples with known sex.

LKSaYR$YRnum = as.numeric(factor(LaneKnownSex$Catch.Year)) # Convert years to numeric and create a new column.
LKSaYR[LKSaYR$Catch.Year == 2022,]$Catch.Year <- 2021 # Combine 2022 samples with 2021 catch year.
LKSaYR[LKSaYR$Catch.Year == 2020,]$Catch.Year <- 2021 # Combine 2020 samples with 2021 catch year.
LKSaYR$YRnum = as.numeric(factor(LKSaYR$Catch.Year)) # Recreate the numeric column after changes.

# Display frequency tables for the modified dataset.
table(LKSaYR$Catch.Year)
table(LKSaYR$YRnum)

# Get a list of starting values for the VBGF model.
YRsvTyp = vbStarts(FFL.mm ~ Final.Age, data = LKSaYR)

## Model Fitting 12.4.2

# Extract age and length data.
Age = as.numeric(LKSaYR$Final.Age)
FFL = as.numeric(LKSaYR$FFL.mm)

# Define various logarithmic von Bertalanffy growth function (VBGF) equations with different parameters for each capture year.
lYRvbLKt = log(FFL) ~ log(Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - t0[YRnum]))))
lYRvbLK = log(FFL) ~ log(Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - t0))))
lYRvbLt = log(FFL) ~ log(Linf[YRnum] * (1 - exp(-K * (Age - t0[YRnum]))))
lYRvbKt = log(FFL) ~ log(Linf * (1 - exp(-K[YRnum] * (Age - t0[YRnum]))))
lYRvbL = log(FFL) ~ log(Linf[YRnum] * (1 - exp(-K * (Age - t0))))
lYRvbK = log(FFL) ~ log(Linf * (1 - exp(-K[YRnum] * (Age - t0))))
lYRvbt = log(FFL) ~ log(Linf * (1 - exp(-K * (Age - t0[YRnum]))))
lYRvb0 = log(FFL) ~ log(Linf * (1 - exp(-K * (Age - t0))))

# Get initial parameter estimates for nonlinear models.
YRsv0 = vbStarts(FFL ~ Age, data = LKSaYR)
YRsv0

# Replicate starting values for different groups.
YRsvLKt = Map(rep, YRsv0, c(4, 4, 4))
YRsvLKt

# Fit the nonlinear model with the specified starting values.
lYRfitLKt = nls(lYRvbLKt, data = LKSaYR, start = YRsvLKt)

# Plot the residuals of the fitted model.
FSAmisc::residPlot(lYRfitLKt, col = rgb(0, 0, 0, 1/3))

# Fit the simplest model.
lYRfit0 = nls(lYRvb0, data = LKSaYR, start = YRsv0)

# Perform a likelihood ratio test to compare models.
lrt(lYRfit0, com = lYRfitLKt, com.name = "All pars differ", sim.names = "No pars differ")
extraSS(lYRfit0, com = lYRfitLKt, com.name = "All pars diff", sim.names = "No pars diff")
# Interpretation: Likelihood Ratio (X^2 = 243.65, p < 2.2e-16) and Extra sum of squares test (F = 14.34, p < 2.2e-16) indicate significant differences.

# Finding the Best Subset Model

# Replicate starting values for other models to compare.
lYRsvLK = Map(rep, YRsv0, c(4, 4, 1))
lYRsvLK

# Starting values for other nested models.
lYRsvLt = Map(rep, YRsv0, c(4, 1, 4))
lYRsvLt
lYRsvKt = Map(rep, YRsv0, c(1, 4, 4))
lYRsvKt

# Fit the nested models.
lYRfitLK = nls(lYRvbLK, data = LKSaYR, start = lYRsvLK)
lYRfitLt = nls(lYRvbLt, data = LKSaYR, start = lYRsvLt)
lYRfitKt = nls(lYRvbKt, data = LKSaYR, start = lYRsvKt)

# Compare each of the nested models using likelihood ratio test.
lrt(lYRfitLK, lYRfitLt, lYRfitKt, com = lYRfitLKt, com.name = "All pars diff", sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Interpretation: Evidence that only (K, t0) (p = 0.2) fit the data as well as (Linf, K, t0). The model (K, t0) is used as the basis because it has the greatest log-likelihood value.

# Further starting values for additional comparisons.
lYRsvK = Map(rep, YRsv0, c(1, 4, 1))
lYRsvK
lYRsvt = Map(rep, YRsv0, c(1, 1, 4))
lYRsvt

# Fit additional models.
lYRfitK = nls(lYRvbK, data = LKSaYR, start = lYRsvK)
summary(lYRfitK)
lYRfitt = nls(lYRvbt, data = LKSaYR, start = lYRsvt)

# Compare additional models.
lrt(lYRfitK, lYRfitt, com = lYRfitKt, com.name = "K,t0 diff", sim.names = c("K diff", "t0 diff"))
# Interpretation: Results suggest that neither K (p = 0.002) nor t0 (p = 0.005) fit as well as (K,t0).

# Model Selection with AIC or BIC

# Replicate starting values for the final model.
lYRsvL = Map(rep, YRsv0, c(4, 1, 1))
lYRsvL
lYRfitL = nls(lYRvbL, data = LKSaYR, start = lYRsvL)

# Compare models using AIC and BIC.
model_comparison <- cbind(
  AIC(lYRfitLKt, lYRfitLK, lYRfitLt, lYRfitKt, lYRfitL, lYRfitK, lYRfitt, lYRfit0),
  BIC(lYRfitLKt, lYRfitLK, lYRfitLt, lYRfitKt, lYRfitL, lYRfitK, lYRfitt, lYRfit0)
)
model_comparison
# Interpretation: The lYRfitKt model has the lowest AIC (-562.81) and lYRfitK has the lowest BIC (-583.5675).

# Plot for best fit based on AIC
predict(lYRfitKt) # Predict values using the lYRfitKt model.
coef(lYRfitKt) # Extract coefficients from the lYRfitKt model.

# Plot observed data.
plot(LKSaYR$Final.Age, LKSaYR$FFL.mm)

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Add VBGF lines for each year group using coefficients from the lYRfitKt model.
# Group 1
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitKt)[1], coef(lYRfitKt)[6], coef(lYRfitKt)[2], seq(0, 14, 0.1)), col = 2) 
# Group 2
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitKt)[1], coef(lYRfitKt)[7], coef(lYRfitKt)[3], seq(0, 14, 0.1)), col = 3) 
# Group 3
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitKt)[1], coef(lYRfitKt)[8], coef(lYRfitKt)[4], seq(0, 14, 0.1)), col = 4) 
# Group 4
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitKt)[1], coef(lYRfitKt)[9], coef(lYRfitKt)[5], seq(0, 14, 0.1)), col = 5) 

# Display frequency table for year groups.
table(LKSaYR$YRnum)

# Plot for best fit based on BIC
predict(lYRfitK) # Predict values using the lYRfitK model.
coef(lYRfitK) # Extract coefficients from the lYRfitK model.

# Plot observed data.
plot(LKSaYR$Final.Age, LKSaYR$FFL.mm)

# Define the von Bertalanffy growth function (VBGF).
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))

# Add VBGF lines for each year group using coefficients from the lYRfitK model.
# Group 1
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitK)[1], coef(lYRfitK)[6], coef(lYRfitK)[2], seq(0, 14, 0.1)), col = 2) 
# Group 2
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitK)[1], coef(lYRfitK)[6], coef(lYRfitK)[3], seq(0, 14, 0.1)), col = 3) 
# Group 3
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitK)[1], coef(lYRfitK)[6], coef(lYRfitK)[4], seq(0, 14, 0.1)), col = 4) 
# Group 4
lines(seq(0, 14, 0.1), vonBert(coef(lYRfitK)[1], coef(lYRfitK)[6], coef(lYRfitK)[5], seq(0, 14, 0.1)), col = 5)

# SUMMARY STATS

# Subset by sex for only fish with final ages and lengths.
LaneAll = LaneAll[!is.na(LaneAll$FFL.mm) & !is.na(LaneAll$Final.Age),]

femaleLS = LaneAll[LaneAll$Macro.sex == "F",] # Subset for females.
maleLS = LaneAll[LaneAll$Macro.sex == "M",]  # Subset for males.
unkLS = LaneAll[LaneAll$Macro.sex == "U",]   # Subset for unknown sex.

# Display frequency table for sex.
table(LaneAll$Macro.sex)
length(unique(LaneAll$Macro.sex))

# Display max age by sex.
max(maleLS$Final.Age)
max(femaleLS$Final.Age)
max(unkLS$Final.Age)

# Mean heights and standard deviations by sex.
mheight.mean = mean(maleLS$FFL.mm)
mheight.mean
mheightsd = sd(maleLS$FFL.mm)
mheightsd

fheight.mean = mean(femaleLS$FFL.mm)
fheight.mean
fheightsd = sd(femaleLS$FFL.mm)
fheightsd

# Plot females only: length vs. age.
plot(FFL.mm ~ Final.Age, data = femaleLS, xlim = c(0, 15), ylim = c(100, 500), pch = 21, bg = "pink2", main = "Female: Length at Age", xlab = "Age (years)", ylab = "Fork Length (mm)")

# Plot males only: length vs. age.
plot(FFL.mm ~ Final.Age, data = maleLS, xlim = c(0, 15), ylim = c(100, 500), pch = 21, bg = "slateblue4", main = "Male: Length at Age", xlab = "Age (years)", ylab = "Fork Length (mm)")

# Subset by fishery type.
table(LaneAll$Fishery)
length(unique(LaneAll$Fishery))

RecLS = LaneAll[LaneAll$Fishery == "REC",] # Subset for recreational.
FILS = LaneAll[LaneAll$Fishery == "FI",]   # Subset for fishery independent.
CM = LaneAll[LaneAll$Fishery == "CM",]     # Subset for commercial.

# Display frequency tables for gear types by fishery type.
table(RecLS$Gear.Name)
table(FILS$Gear.Name)
table(CM$Gear.Name)
table(LaneAll$Gear.Name)

# Plot only recreational samples: length vs. age.
plot(RecLS$FFL.mm ~ RecLS$Final.Age, pch = 1, main = "Lane Snapper Fishery: Length At Age", xlab = 'Fractional Age (years)', ylab = 'Fork Length (mm)', xlim = c(0, 20), ylim = c(100, 500))

# Summary statistics for commercial samples.
summary(CM$FFL.mm)
mean(CM$FFL.mm)
sd(CM$FFL.mm)
mean(CM$Final.Age)
sd(CM$Final.Age)
min(CM$FFL.mm)
max(CM$FFL.mm)
min(CM$Final.Age)
max(CM$Final.Age)

# Summary statistics for all samples.
summary(LaneAll$FFL.mm)

# Summary statistics for recreational samples.
summary(RecLS$FFL.mm)
mean(RecLS$FFL.mm)
sd(RecLS$FFL.mm)
mean(RecLS$Final.Age)
sd(RecLS$Final.Age)
min(RecLS$FFL.mm)
max(RecLS$FFL.mm)
min(RecLS$Final.Age)
max(RecLS$Final.Age)

# Summary statistics for fishery independent samples.
summary(FILS$FFL.mm)
mean(FILS$FFL.mm)
sd(FILS$FFL.mm)
mean(FILS$Final.Age)
sd(FILS$Final.Age)
min(FILS$FFL.mm)
max(FILS$FFL.mm)
min(FILS$Final.Age)
max(FILS$Final.Age)

# Add other points to the plot.
points(FILS$FFL.mm ~ FILS$Final.Age, pch = 3) # Adds fishery independent observed on the same graph.
points(CM$FFL.mm ~ CM$Final.Age, pch = 9)    # Adds commercial observed on the same graph.

# Add legend to the plot.
legend("bottomright", legend = c("Recreational", "Fishery Independent", "Commercial"), lty = 0, pch = c(1, 3, 9))

# Add horizontal lines to the plot.
abline(h = x)
abline(h = 203.2, col = "blue")

# Summary statistics for commercial samples.
mean(CM$FFL.mm)
sd(CM$FFL.mm)
min(CM$FFL.mm)
max(CM$FFL.mm)
mean(CM$Final.Age)
sd(CM$Final.Age)
min(CM$Final.Age)
max(CM$Final.Age)

# Summary statistics for recreational samples.
mean(RecLS$FFL.mm)
sd(RecLS$FFL.mm)
min(RecLS$FFL.mm)
max(RecLS$FFL.mm)
mean(RecLS$Final.Age)
sd(RecLS$Final.Age)
min(RecLS$Final.Age)
max(RecLS$Final.Age)

# Summary statistics for fishery independent samples.
mean(FILS$FFL.mm)
sd(FILS$FFL.mm)
min(FILS$FFL.mm)
max(FILS$FFL.mm)
mean(FILS$Final.Age)
sd(FILS$Final.Age)
min(FILS$Final.Age)
max(FILS$Final.Age)

# Create points for state landed.
as.data.frame(table(LaneAll$State))
sum(65 + 536 + 3 + 15 + 19 + 215 + 1)

# Subset by state.
FL = LaneAll[LaneAll$State == "FL",]  # Subset for Florida.
AL = LaneAll[LaneAll$State == "AL",]  # Subset for Alabama.
TX = LaneAll[LaneAll$State == "TX",]  # Subset for Texas.
UNK = LaneAll[LaneAll$State == "UNK",]  # Subset for Unknown.
LA = LaneAll[LaneAll$State == "LA",]  # Subset for Louisiana.
MS = LaneAll[LaneAll$State == "MS",]  # Subset for Mississippi.
NC = LaneAll[LaneAll$State == "NC",]  # Subset for North Carolina.

# Display frequency table for state.
table(LaneAll$State)

# Plot and points for state landed.

# Plot only Florida samples.
plot(FL$FFL.mm ~ FL$Final.Age, pch = 21, bg = 'light green', 
     main = "State Landed: Fork Length at Fractional Age for Lane Snapper",
     xlab = 'Fractional Age (years)', 
     ylab = 'Fork Length (mm)', xlim = c(0, 15), ylim = c(50, 500))

# Add points for other states.
points(AL$FFL.mm ~ AL$Final.Age, pch = 21, bg = "orange")  # Adds Alabama observed on same graph.
points(TX$FFL.mm ~ TX$Final.Age, pch = 21, bg = "purple2")  # Adds Texas observed on same graph.
points(UNK$FFL.mm ~ UNK$Final.Age, pch = 21, bg = "red")  # Adds Unknown observed on same graph.
points(LA$FFL.mm ~ LA$Final.Age, pch = 21, bg = "yellow")  # Adds Louisiana observed on same graph.
points(MS$FFL.mm ~ MS$Final.Age, pch = 21, bg = "dark green")  # Adds Mississippi observed on same graph.
points(NC$FFL.mm ~ NC$Final.Age, pch = 21, bg = "blue")  # Adds North Carolina observed on same graph.

# Add legend to the plot.
legend("bottomright", legend = c("AL", "FL", "LA", "MS", "NC", "TX", "Unk"), col = c("orange", "light green", "yellow", "dark green", "blue", "purple2", "red"), lty = 0, pch = 20, 21)

# Create points for edge type.
E2 = LaneAll[LaneAll$Edge.type == "2",]  # Subset for edge type 2.
E4 = LaneAll[LaneAll$Edge.type == "4",]  # Subset for edge type 4.
E6 = LaneAll[LaneAll$Edge.type == "6",]  # Subset for edge type 6.

# Plot by edge type.
plot(FFL.mm ~ Final.Age, data = LaneAll, pch = 21, bg = 'gray', main = "Edge Type: Length At Age", xlab = 'Age (years)', 
     ylab = 'Fork Length (mm)', xlim = c(0, 15), ylim = c(100, 500))

# Add points for different edge types.
points(E2$FFL.mm ~ E2$Final.Age, pch = 21, bg = "orange")  # Adds edge type 2 observed on same graph.
points(E4$FFL.mm ~ E4$Final.Age, pch = 21, bg = "purple2")  # Adds edge type 4 observed on same graph.
points(E6$FFL.mm ~ E6$Final.Age, pch = 21, bg = "red")  # Adds edge type 6 observed on same graph.

# Display frequency table for edge type.
table(LaneAll$Edge.type)

# Plot males and females.

plot(FFL.mm ~ Final.Age, data = LaneKnownSexF, pch = 21, bg = 'pink', main = "Edge Type: Length At Age", xlab = 'Age (years)', 
     ylab = 'Fork Length (mm)', xlim = c(0, 18), ylim = c(100, 500))

# Add points for males.
points(LaneKnownSexM$FFL.mm ~ LaneKnownSexM$Final.Age, pch = 21, bg = "blue")  # Adds males observed on same graph.

# Subset data by catch year.
as.data.frame(table(LaneAll$Catch.Year))

Y12 = LaneAll[LaneAll$Catch.Year == "2012",]  # Subset for 2012 samples.
Y13 = LaneAll[LaneAll$Catch.Year == "2013",]  # Subset for 2013 samples.
Y14 = LaneAll[LaneAll$Catch.Year == "2014",]  # Subset for 2014 samples.
Y15 = LaneAll[LaneAll$Catch.Year == "2015",]  # Subset for 2015 samples.
Y16 = LaneAll[LaneAll$Catch.Year == "2016",]  # Subset for 2016 samples.
Y17 = LaneAll[LaneAll$Catch.Year == "2017",]  # Subset for 2017 samples.
Y21 = LaneAll[LaneAll$Catch.Year == "2021",]  # Subset for 2021 samples.
Y22 = LaneAll[LaneAll$Catch.Year == "2022",]  # Subset for 2022 samples.

# Plot only 2015 samples.
plot(Y15$FFL.mm ~ Y15$Final.Age, pch = 21, bg = 'dark green', main = "Year-Landed: Length At Age", xlab = 'Age (years)',
     ylab = 'Fork Length (mm)', xlim = c(0, 18), ylim = c(100, 500))

# Add points for other years.
points(Y16$FFL.mm ~ Y16$Final.Age, pch = 21, bg = "orange")  # Adds 2016 observed on same graph.
points(Y17$FFL.mm ~ Y17$Final.Age, pch = 21, bg = "purple")  # Adds 2017 observed on same graph.
points(Y21$FFL.mm ~ Y21$Final.Age, pch = 21, bg = "yellow")  # Adds 2021 observed on same graph.

# Add legend to the plot.
legend("bottomright", legend = c("2015", "2016", "2017", "2021"), col = c("dark green", "orange", "purple", "yellow"), lty = 0, pch = 20, 21)

# Subset for SEAMAP.
SP = LaneAll[LaneAll$Source.Code == "SEAMAP",]  # Subset for SEAMAP.

# Add SEAMAP points to the plot.
points(SP$FFL.mm ~ SP$Final.Age, pch = 21, bg = "blue")  # Adds SEAMAP observed on same graph.

##############VonB#########

## Both sexes combined ##
age = as.numeric(LaneAll$Final.Age)  # Extracts age data (years old).
cfl = as.numeric(LaneAll$FFL.mm)  # Extracts length data.

options(na.action = na.exclude)  # Tells R to ignore missing variables.

# Create a nonlinear model using least squares (lsvb: least square von B).
alllsvb = nls(FFL.mm ~ Linf * (1 - exp(-k * (Final.Age - tnot))), data = LaneAll, 
              start = list(Linf = 500, k = 0.4, tnot = 0.0), trace = TRUE)  # Trace shows parameter search.

summary(alllsvb)  # Extract parameter estimates.
attributes(summary(alllsvb))  # List attribute names.
testsd = summary(alllsvb)$par  # Select parameter estimates with p-values and SD.
testsd
nlsparams = summary(alllsvb)$par[,1]  # Get just the estimates for NLS.
lsvbSE = summary(alllsvb)$par[,2]  # SE to estimate 95% confidence intervals (+/- 2 SD are 95%).
lsvbSE
nlsparams

# Plot combined fit model.
vonBert <- function(Linf, t0, K, age) Linf * (1 - exp(-K * (age - t0)))  # Calling the VB function.
plot(LaneAll$Final.Age, LaneAll$FFL.mm,
     xlab = "Fractional Age (years)", ylab = "Fork Length (mm)")
lines(seq(0, 14, 0.1), vonBert(coef(alllsvb)[1], coef(alllsvb)[3],
                               coef(alllsvb)[2], seq(0, 14, 0.1)), col = "purple", lwd = 4)

# Subset for specific ages.
old = LaneAll[LaneAll$Final.Age == "18.45722",]  # Subset for only 18yo.
points(old$FFL.mm ~ old$Final.Age, pch = 21, bg = "orange")
old2 = LaneAll[LaneAll$Final.Age == "15.64956",]  # Subset for only 15yo.
points(old2$FFL.mm ~ old2$Final.Age, pch = 21, bg = "red")

# Add legend to the plot.
legend("bottomright", inset = .05, c("Females, Males and Unknown Sex"),
       bty = "n", cex = 0.8, fill = c("purple"))

# Subset for older age groups.
olderthan14 = LaneAll[LaneAll$Final.Age > 14,]  # Subset for older than 14yo.
olderthan10 = LaneAll[LaneAll$Final.Age > 10,]  # Subset for older than 10yo.
table(olderthan10$Fishery)

# Display coefficients of the fitted model.
coef(alllsvb)

# Confidence Intervals for sexes combined parameters.
sample.n = length(LaneAll$FFL.mm)
sample.n
alpha = 0.05
degrees.freedom = sample.n - 1
degrees.freedom
t.score = qt(p = alpha / 2, df = degrees.freedom, lower.tail = FALSE)
print(t.score)

# Calculate margins of error and confidence intervals.
Lse = 3.28 
Lmargin.error = t.score * Lse
Kse = 0.0249023
Kmargin.error = t.score * Kse
tse = 0.1082970 
tmargin.error = t.score * tse

Llower.bound <- 359.5398911 - Lmargin.error
Lupper.bound <- 359.5398911 + Lmargin.error
print(c(Llower.bound, Lupper.bound))

Klower.bound <- 0.4300123 - Kmargin.error
Kupper.bound <- 0.4300123 + Kmargin.error
print(c(Klower.bound, Kupper.bound))

tlower.bound <- -0.8015067 - tmargin.error
tupper.bound <- -0.8015067 + tmargin.error
print(c(tlower.bound, tupper.bound))

Lmargin.error
Kmargin.error
tmargin.error

# Display frequency table for state.
table(LaneAll$State)

# NORMAL Confidence Intervals for FEMALE and MALE parameters.
coef(fitKt)
summary(fitKt)  # Extract parameter estimates.
attributes(summary(fitKt))  # List attribute names.
ftestsd = summary(fitKt)$par  # Select parameter estimates with p-values and SD.
ftestsd
fnlsparams = summary(fitKt)$par[,1]  # Get just the estimates for NLS.
flsvbSE = summary(fitKt)$par[,2]  # SE to estimate 95% confidence intervals (+/- 2 SD are 95%).
flsvbSE
fnlsparams

# Calculate confidence intervals for the female group.
fsample.n = length(lsf$FFL.mm)
fsample.n
falpha = 0.05
fdegrees.freedom = fsample.n - 1
fdegrees.freedom
ft.score = qt(p = falpha / 2, df = fdegrees.freedom, lower.tail = FALSE)
print(ft.score)

# Linfsd
# Calculate margins of error and confidence intervals for the female group.
fLse = 5.50682947
fLmargin.error = ft.score * fLse
fKse = 0.3243378
fKmargin.error = ft.score * fKse
ftse = -1.0822993
ftmargin.error = ft.score * ftse

# Display margin of error for Linf.
fLmargin.error

# Calculate and display confidence intervals for Linf, K, and t0.
fLlower.bound <- 373.1817015 - fLmargin.error
fLupper.bound <- 373.1817015 + fLmargin.error
print(c(fLlower.bound, fLupper.bound))

fKlower.bound <- 0.3243378 - fKmargin.error
fKupper.bound <- 0.3243378 + fKmargin.error
print(c(fKlower.bound, fKupper.bound))

ftlower.bound <- -1.0822993 - ftmargin.error
ftupper.bound <- -1.0822993 + ftmargin.error
print(c(ftlower.bound, ftupper.bound))

fLmargin.error
fKmargin.error
ftmargin.error

# Male Group 2
# Calculate sample size, degrees of freedom, and t-score for males.
msample.n = 306
msample.n
malpha = 0.05
mdegrees.freedom = sample.n - 1
mdegrees.freedom
mt.score = qt(p = malpha / 2, df = mdegrees.freedom, lower.tail = FALSE)
print(mt.score)

# Display SE and NLS parameters.
flsvbSE
fnlsparams

# Calculate margins of error and confidence intervals for the male group.
mLse = 6.51339051
mLmargin.error = mt.score * mLse
mKse = 0.03623843
mKmargin.error = mt.score * mKse
mtse = 0.10207888
mtmargin.error = mt.score * mtse

# Display margin of error for Linf.
mLmargin.error

# Calculate and display confidence intervals for Linf, K, and t0.
mLlower.bound <- 366.5629548 - mLmargin.error
mLupper.bound <- 366.5629548 + mLmargin.error
print(c(mLlower.bound, mLupper.bound))

mKlower.bound <- 0.4465838 - mKmargin.error
mKupper.bound <- 0.4465838 + mKmargin.error
print(c(mKlower.bound, mKupper.bound))

mtlower.bound <- -0.4703319 - mtmargin.error
mtupper.bound <- -0.4703319 + mtmargin.error
print(c(mtlower.bound, mtupper.bound))

mLmargin.error
mKmargin.error
mtmargin.error

# Fishery Gear subsets
FITR = LaneAll[LaneAll$Gear.Name == "Trawl",]  # Subset for only Fishery Independent Trawl.

# Display summary statistics for gear names.
summary(LaneAll$Gear.Name)

# LOG Confidence Intervals for FEMALE and MALE parameters
coef(lfitKt)
summary(lfitKt)  # Extract parameter estimates.
attributes(summary(lfitKt))  # List attribute names.
ftestsd = summary(lfitKt)$par  # Select parameter estimates with p-values and SD.
ftestsd
fnlsparams = summary(lfitKt)$par[,1]  # Get just the estimates for NLS.
flsvbSE = summary(lfitKt)$par[,2]  # SE to estimate 95% confidence intervals (+/- 2 SD are 95%).
flsvbSE
fnlsparams

# Calculate sample size, degrees of freedom, and t-score for the female group.
fsample.n = length(lsf$FFL.mm)
fsample.n
falpha = 0.05
fdegrees.freedom = sample.n - 1
fdegrees.freedom
ft.score = qt(p = falpha / 2, df = fdegrees.freedom, lower.tail = FALSE)
print(ft.score)

# Linfsd
# Calculate margins of error and confidence intervals for the female group.
fLse = 5.79974751
fLmargin.error = ft.score * fLse
fKse = 0.02415698
fKmargin.error = ft.score * fKse
ftse = 0.19473366
ftmargin.error = ft.score * ftse

# Display margin of error for Linf.
fLmargin.error

# Calculate and display confidence intervals for Linf, K, and t0.
fLlower.bound <- 375.4189617 - fLmargin.error
fLupper.bound <- 375.4189617 + fLmargin.error
print(c(fLlower.bound, fLupper.bound))

fKlower.bound <- 0.3024100 - fKmargin.error
fKupper.bound <- 0.3024100 + fKmargin.error
print(c(fKlower.bound, fKupper.bound))

ftlower.bound <- -1.0547568 - ftmargin.error
ftupper.bound <- -1.0547568 + ftmargin.error
print(c(ftlower.bound, ftupper.bound))

fLmargin.error
fKmargin.error
ftmargin.error
flsvbSE
fnlsparams

# Male Group 2
# Calculate sample size, degrees of freedom, and t-score for males.
msample.n = 305
msample.n
malpha = 0.05
mdegrees.freedom = sample.n - 1
mdegrees.freedom
mt.score = qt(p = alpha / 2, df = mdegrees.freedom, lower.tail = FALSE)
print(mt.score)

# Display SE and NLS parameters.
flsvbSE  # 1 is female, 2 is male.

# Calculate margins of error and confidence intervals for the male group.
mLse = 5.79974751
mLmargin.error = mt.score * mLse
mKse = 0.03715673
mKmargin.error = mt.score * mKse
mtse = 0.16544387
mtmargin.error = mt.score * mtse

# Display margin of error for Linf.
mLmargin.error

# Calculate and display confidence intervals for Linf, K, and t0.
mLlower.bound <- 375.4189617 - mLmargin.error
mLupper.bound <- 375.4189617 + mLmargin.error
print(c(mLlower.bound, mLupper.bound))

mKlower.bound <- 0.4101349 - mKmargin.error
mKupper.bound <- 0.4101349 + mKmargin.error
print(c(mKlower.bound, mKupper.bound))

mtlower.bound <- -0.6217888 - mtmargin.error
mtupper.bound <- -0.6217888 + mtmargin.error
print(c(mtlower.bound, mtupper.bound))

mLmargin.error
mKmargin.error
mtmargin.error

# Display coefficients of the fitted model.
coef(fitKt)
fnlsparams

# Fishery Gear subsets
FITR = LaneAll[LaneAll$Gear.Name == "Trawl",]  # Subset for only Fishery Independent Trawl.

# Display summary statistics for gear names.
summary(LaneAll$Gear.Name)

# Display frequency tables for catch year and Gulf side.
table(LaneKnownSex$Catch.Year)
table(LKSaGS$Catch.Year)
table(LKSaYR$Catch.Year)
table(LaneKnownSex$GulfSide)
table(LKSaGS$GulfSide)

# Display frequency table for sex.
table(LaneKnownSex$Macro.sex)

# Create and display sex table.
Sex.Table <- table(LaneKnownSex$Macro.sex)
Sex.Table

# Perform chi-square test for sex ratio.
chisq.test(x = c(305, 306), p = c(.5, .5))
binom.test(305, 611)  # Not significant, ratio is 1:1.
# X^2 = 0.0016367, degrees of freedom = 1, p-value = 0.9677. Since p-value is above 0.05, we accept the null hypothesis of a 50/50 sex ratio.

305 / (305 + 306)
0.4991

## Regular Examining Differences Among Capture Years with FIXED t0 ONLY Hand Line

# Display frequency table for catch year.
table(LaneAll$Catch.Year)

# Create numeric year column and adjust years.
LaneKnownCY <- LaneAll
LaneKnownCY$YRnum = as.numeric(factor(LaneKnownCY$Catch.Year))  # Makes all years a number and creates a column with it.
LaneKnownCY[LaneKnownCY == 2022] <- 2021  # Add 2022 samples to 2021 catch year.
LaneKnownCY[LaneKnownCY == 2020] <- 2021  # Add 2020 samples to 2021 catch year.
LaneKnownCY$YRnum = as.numeric(factor(LaneKnownCY$Catch.Year))  # Makes all years a number and creates a column with it.

# Subset for only adult hook and line samples.
LaneKnownCYAdults = LaneKnownCY[LaneKnownCY$Final.Age > 1,]  # Subset for older than 1 year old.
LaneKnownCYAdultsHL = LaneKnownCYAdults[LaneKnownCYAdults$Gear.Code == "HL",]  # Subset for hand line samples.

# Remove 2021 samples due to insufficient sample size.
LaneKnownCYAdultsHL = LaneKnownCYAdultsHL[LaneKnownCYAdultsHL$Catch.Year < 2021,]

# Display frequency tables for the subset.
table(LaneKnownCYAdultsHL$YRnum)
table(LaneKnownCYAdultsHL$Catch.Year)
table(LaneKnownCYAdultsHL$State)

# Get a list of starting values.
YRsvTyp = vbStarts(FFL.mm ~ Final.Age, data = LaneKnownCYAdultsHL, fixed = list(t0 = 0))  # Formula for generating starting values and assigning to an object for later use.

# Display starting values.
vbStarts(FFL.mm ~ Final.Age, data = LaneKnownCYAdultsHL, vbStartsDP = TRUE, fixed = list(t0 = 0))
YRsvTyp = list(Linf = 319, K = 0.67, t0 = 0)
YRsvTyp = list(Linf = max(LaneKnownCYAdultsHL$FFL.mm, na.rm = TRUE), K = 0.67, t0 = 0)
YRsvTyp

## Model Fitting 12.4.2

# Convert age and length data to numeric.
Age = as.numeric(LaneKnownCYAdultsHL$Final.Age)  # Pulls the age data (years old).
FFL = as.numeric(LaneKnownCYAdultsHL$FFL)  # Pulls the length data.

# Define model formulas.
YRvbLKt = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - t0[YRnum])))
YRvbLKt0 = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbLK0 = FFL ~ Linf[YRnum] * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbLt0 = FFL ~ Linf[YRnum] * (1 - exp(-K * (Age - 0[YRnum])))
YRvbKt0 = FFL ~ Linf * (1 - exp(-K[YRnum] * (Age - 0[YRnum])))
YRvbL0 = FFL ~ Linf[YRnum] * (1 - exp(-K * (Age - 0)))
YRvbK0 = FFL ~ Linf * (1 - exp(-K[YRnum] * (Age - 0)))
YRvbt0 = FFL ~ Linf * (1 - exp(-K * (Age - 0[YRnum])))
YRvb0 = FFL ~ Linf * (1 - exp(-K * (Age - 0)))
YRvbt0 = FFL ~ Linf * (1 - exp(-K * (Age - 0)))

# Get starting values.
YRsv0 = vbStarts(FFL ~ Age, data = LaneKnownCYAdultsHL, fixed = list(t0 = 0))
YRsv0

YRsvt0 = YRsv0
YRsvt0$t0 = NULL

YRsvLKt = Map(rep, YRsv0, c(3, 3, 3))
YRsvLKt

YRsvLKt0 = YRsvLKt
YRsvLKt0$t0 = NULL
YRsvLKt0

# Fit the nonlinear models.
YRfitLKt0 = nls(YRvbLKt0, data = LaneKnownCYAdultsHL, start = YRsvLKt0)
YRfitLKt0

# Residual plot.
FSAmisc::residPlot(YRfitLKt, col = rgb(0, 0, 0, 1/3))

# Fit another nonlinear model.
YRfitt0 = nls(YRvbt0, data = LaneKnownCYAdultsHL, start = YRsvt0)
YRfitt0

# Likelihood ratio test.
lrt(YRfitt0, com = YRfitLKt0, com.name = "All pars differ",
    sim.names = "No pars differ")
extraSS(YRfitt0, com = YRfitLKt0, com.name = "All pars diff",
        sim.names = "No pars diff")
# Likelihood Ratio (X^2 = 159.16, p < 2.2e-16): Yes, significantly different.
# Extra sum of squares test (F = 19.81, p < 2.2e-16): Yes, significantly different.

# Finding the Best Subset Model.
YRsvLK = Map(rep, YRsv0, c(3, 3, 1))
YRsvLK

YRsvLK0 = YRsvLK
YRsvLK0$t0 = NULL
YRsvLK0

# Starting values for other models.
YRsvLt = Map(rep, YRsv0, c(3, 1, 3))
YRsvLt

YRsvLt0 = YRsvLt
YRsvLt0$t0 = NULL
YRsvLt0

YRsvKt = Map(rep, YRsv0, c(1, 3, 3))
YRsvKt

YRsvKt0 = YRsvKt
YRsvKt0$t0 = NULL
YRsvKt0

# Fit the nested models.
YRfitLK0 = nls(YRvbLK0, data = LaneKnownCYAdultsHL, start = YRsvLK0)
YRfitLt0 = nls(YRvbLt0, data = LaneKnownCYAdultsHL, start = YRsvLt0)
YRfitKt0 = nls(YRvbKt0, data = LaneKnownCYAdultsHL, start = YRsvKt0)

# Likelihood ratio test for nested models.
lrt(YRfitLK0, YRfitLt0, YRfitKt0, com = YRfitLKt0,
    com.name = "All pars diff", 
    sim.names = c("Linf,K diff", "Linf,t0 diff", "K,t0 diff"))
# Evidence that (K,t0) (0.06) and (Linf,K) (p = .06) fit the data as well as (Linf, K, t0),
# but not Linf,t0 (0.003).
# The model (Kinf, t0) is used for further reduction.

# Fit another model.
YRfitK = nls(YRvbK, data = LaneKnownCYAdultsHL, start = YRsvK)
summary(YRfitK)

YRfitt = nls(YRvbt, data = LaneKnownCYAdultsHL, start = YRsvt)

# Likelihood ratio test.
lrt(YRfitL, YRfitK, com = YRfitLK, com.name = "L,K dif", sim.names = c("L dif", "K dif"))
# These results suggest that neither K nor t0 (p = 2.844e-06) fit as well as (K,t0).

lrt(YRfit0, com = YRfitK, com.name = "K dif", sim.names = "No pars dif")
# The comparison of K to omega confirms that a difference between groups exists,
# and that difference is only evident for K (p < 2.2e-16).

# Model Selection with AIC or BIC.
YRsvL = Map(rep, YRsv0, c(4, 1, 1))
YRsvL
YRfitL = nls(YRvbL, data = LaneKnownCYAdultsHL, start = YRsvL)

# Compare models using AIC and BIC.
cbind(AIC(YRfitLKt, YRfitLK, YRfitLt, YRfitKt, YRfitL, YRfitK, YRfitt, YRfit0),
      BIC(YRfitLKt, YRfitLK, YRfitLt, YRfitKt, YRfitL, YRfitK, YRfitt, YRfit0))

# Subset for only recreational fishery.
RecLS = LaneAll[LaneAll$Fishery == "REC",]

# Summarize data by fishery.
my_sum <- LaneAll %>%
  group_by(Fishery) %>%
  summarise(
    n = n(),
    mean = mean(FFL.mm),
    sd = sd(FFL.mm)
  ) %>%
  mutate(se = sd / sqrt(n)) %>%
  mutate(ic = se * qt((1 - 0.05) / 2 + .5, n - 1))
my_sum
