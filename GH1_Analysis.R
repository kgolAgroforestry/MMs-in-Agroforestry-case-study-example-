# A companion R script for the publication:
# A system-scale approach to perennial agroecosystem research - 
# sampling strategies and best practice statistical methods

# Authors: KG, H-PP, E-MM, WB, AG-S, JO, LB, AG, and SJ
# Journal: Methods in Ecology and Evolution
# Date of submission: 28 October 2022

# This code is designed to act as an introduction to 
# multilevel models. It is not suitable for analysis 
# of other data without substantial modifications.
# It is recommended to run it until Step 9 first (starting with custom-made
# functions in lines 604), before implementing any changes.  

############################################################################x

# Step 1: Selection and classification of variables--------------------------
# (a) Identify the research question
# (b) Identify the variables:
#        - Dependent variable
#        - Explanatory variable(s)
# (c) Consider relevant spatial/temporal scales for the dependent variable
# (d) Identify the random effect (the grouping variable for 
#     multiple observations that are correlated)

# Load (or install) the required packages 
library(tidyverse) # used to summarize the raw data and for plotting
library(nlme) # used for mixed effects modeling
library(lattice) # used for plotting
library(ggmap) # used to show sampling locations in a graphical form
# ggmap works in wgs84 projection (convert coords if needed)
library(maptools) # used for mapping of the sampling locations
library(sp) # used for mapping of the sampling locations
library(rgdal) # used for mapping of the sampling locations
library(gstat) # used for variogram modelling 
library(outliers) # identification of outliers (tests)
library(EnvStats) # identification of outliers (tests, more detailed)
library(stargazer) # tabularizing the results for reporting
library(partR2) # used for estimation of R^2 values, not recommended here; 
# issues with MerMod objects

# Step 2: Data investigation and linear regression fit (with fixed eff.)-----

# Import and check the data
yield <- read.csv("~/Downloads/yield.csv")

head(yield)
names(yield)

# Ensure the variables of interest are coded as factors
sampleID<-as.factor(yield$sample)
trans<-as.factor(yield$transect)
row<-as.factor(yield$row)
direction<-as.factor(yield$direction)
dist<-as.factor(yield$distance)
steepness<-as.factor(yield$slope)
steepness.deg<-as.numeric(yield$slope.deg)
wheat<-as.numeric(yield$grain.t.ha)
land.prod<-as.factor(yield$land.prod.basic)
NDVI<-as.numeric(yield$ndvi)
SM<-as.numeric(yield$sm)
x<-as.numeric(yield$xcoord)
y<-as.numeric(yield$ycoord)

# Coerce into a data frame or models won't run
my.yield=data.frame(wheat,dist,direction,row, trans, sampleID, land.prod, SM,
                    NDVI, steepness.deg, x, y) 
my.yield # use to run models lm and lme

my.yield %>% group_by(dist) %>%         # summary statistics for distance
             summarize(min = min(wheat), 
             max = max(wheat),
             mean = mean(wheat),
             median = median(wheat),
             sd = sd(wheat),
             n = length(wheat))

my.yield %>% group_by(direction) %>%     # summary statistics for direction
             summarize(min = min(wheat), 
             max = max(wheat),
             mean = mean(wheat),
             median = median(wheat),
             sd = sd(wheat),
             n = length(wheat))

# Understand your data, inspect with plots, boxplots, check for outliers
boxplot(wheat~row, data = my.yield) # the variation might be a good indicator 
# that this variable should be fitted as fixed or random effect
boxplot(wheat~dist, data = my.yield)
boxplot(wheat~direction, data = my.yield)
plot(wheat~steepness.deg, data = my.yield)
boxplot(wheat~steepness, data = my.yield)
boxplot(wheat~land.prod, data = my.yield)
plot(wheat~SM, data = my.yield)
plot(wheat~NDVI, data = my.yield)

xyplot(wheat ~ dist | row, data = my.yield)
xyplot(wheat ~ row, groups=direction, data = my.yield)

# Looking for outliers (you can use PMCMRplus & ggstatsplot for 
# excellent graphics but check if the 'insight' package is updated or 
# ggstatsplot might not run as of May 2022)

boxplot(my.yield$wheat) # no visible outliers

boxplot(my.yield$wheat ~ my.yield$row) # no visible outliers

boxplot(my.yield$wheat ~ my.yield$direction) # no visible outliers

# suspicious value in row 4

outliers<-boxplot(my.yield$wheat ~ my.yield$row, plot=FALSE)$out
outliers # no outliers 

outlier.test <- tapply(my.yield$wheat, my.yield$row, rosnerTest, k=2)

outlier.test # no statistically sig. outliers

# Develop an lm model
M.lm<-lm(wheat~dist*direction, data=my.yield)
M.lm
summary(M.lm)$coef
drop1(M.lm, test='F') # check if the interaction is significant

# Extract the residuals
E.GLM <- resid(M.lm)

# Inspect your residuals 
# Are there any patterns, if so, violation of assumptions

par(mfrow = c(2, 2))
plot(M.lm)
par(mfrow = c(1, 1))

# There are patterns in the residuals (observed values - fitted values)
# as we move from left to right, the prediction errors increase

# Are the residuals normally distributed? 
# If not, transform to perform (e.g., log)

qqnorm(E.GLM)
qqline(E.GLM, col = "dark green", lwd = 2)
shapiro.test(E.GLM) # normal distribution

# Is there equality of variance?

plot(dist, E.GLM, xlab = "Distance", ylab = "Residuals")
plot(direction, E.GLM, xlab = "Direction", ylab = "Residuals")

# We can accept the null hypothesis that the variances are equal
# for both predictors

# Looking at residual independence in some more detail
# Spatial dependence in transect sampling is a feature of the design

E<-rstandard(M.lm)
spatial.errors<-data.frame(E, 
                           yield$xcoord,
                           yield$ycoord)

spatial.errors # use the headings from the data frame spatial.errors
# for coordinates() function

coordinates(spatial.errors)<-c("yield.xcoord","yield.ycoord")

B1<-bubble(spatial.errors,"E",col=c("dark green","orange"),
       main="Residuals",xlab="X-coordinates",
       ylab="Y-coordinates") # clustering present

B1 # the bubble plot code from Zuur et al. (2009),
# apart from showing the spatial distribution of residuals,
# it is very well-suited for displaying results from transects
# i.e., at what end of the transects are the errors higher?

# Plot variogram 

vario1=variogram(E~1,spatial.errors)

plot(vario1) # variability as a function of distance,
# plots are spatially autocorrelated  

# Plot variogram (multi directional)

vario2=variogram(E~1,
                 spatial.errors,
                 alpha=c(0,45,90,135))

plot(vario2) # spatial dependence in the same direction (isotropy)
# one of the conditions to be met 

# In summary:
# normality - check 
# equality of variance - check 
# independence of residual error terms - check 
# be ware of high leverage points (see Cook's distance plot below)

# Outliers - revisited

cooksdist <- cooks.distance(M.lm)
nrow(my.yield)

plot(cooksdist, pch="*", cex=1, 
     main="High leverage points")
abline(h = 4/144, col="dark green") # following the 4/total sample size rule 
# based on Bollen and Jackman (1990)

text(x=1:length(cooksdist)+1, y=cooksdist, 
     labels=ifelse(cooksdist>4/144, names(cooksdist),""), col="dark green")

high.leverage<-as.numeric(names(cooksdist)[(cooksdist > (4/144))])

wheat.out<-my.yield[-high.leverage,] # taking out 8 values 
count(wheat.out)

wheat.out

wheat.out %>% group_by(row) %>% count(dist)

# The high leverage values are found at a distance of 4.5m; 
# instead of removing (multiple) outliers, considering how would the 
# removal of these values change the parameter estimates is essential for 
# reaching final conclusions

# Fit linear regression with GLS
# Fitting with GLS is needed to compare lm with lme; the output is the same

f<-formula(wheat~dist*direction, data=my.yield) 

M.gls<-gls(f,method = "ML") # refit lm with gls to allow for
# comparison against mixed effects models

# Step 3: Fitting a Marginal Model------------------------------------------
# Investigate a Marginal Model
# (GLS that accounts for residual dependencies)
# benefit: no issues with REML vs. ML bias
# benefit: allows for the variances in the covariance structure to be negative
# which is also possible in lme but needs to be specified (addressed below)

# Visual representation of the actual crop yield distribution

yield.on.map <- ggplot(data = my.yield,
                mapping = aes(x = x, y = y, color = wheat)) +
                geom_point(size = 3) +
                scale_color_gradientn(colors = c("red", "yellow", "dark green"))
yield.on.map

MM <-gls(f, method="ML", 
            correlation = corExp(form=~x+y),
            data=my.yield)

anova(M.gls, MM)

MM.Gau <- gls(f, method="ML", 
              correlation = corGaus(form=~x+y),
              data=my.yield)
MM.Lin <- gls(f, method="ML", 
              correlation = corLin(form=~x+y),
              data=my.yield)
MM.Rat <- gls(f, method="ML", 
              correlation = corRatio(form=~x+y),
              data=my.yield)
MM.Sph <- gls(f, method="ML", 
              correlation = corSpher(form=~x+y),
              data=my.yield)

AIC(MM,MM.Gau,MM.Lin,MM.Rat,MM.Sph)
plot(Variogram(MM, form=~x+y)) 
plot(Variogram(MM.Rat, form=~x+y)) 

# include nugget

MM.nugget <- update(MM, correlation = corExp(form=~x+y, nugget=T)) 

anova(MM, MM.nugget) # not significantly different, MM is superior

summary(MM)

anova(MM) # anova is sequential, mind the order in which the explanatory 
# variables are arranged

# Investigating patterns in the residuals 

plot(Variogram(MM, form=~x+y,resType = "n")) # no patterns

# Model diagnostics (run the functions in the line 604 first)
# Refit the model with REML prior to running model diagnostics

diagnostics.simple(MM)  # the spread of residuals increases with distance
bubble.plotting(MM)     # the bubble plot shows less clustering than in the 
# M.lm 

# Allowing for different variances per 'direction'

MM.2 <-gls(f, method="ML", 
              correlation = corExp(form=~x+y),
              weights =     varIdent(form=~1|direction),
              data =        my.yield)

anova(MM,MM.2) # not significantly better

# Allowing for different spread per 'distance' and 
# 'direction' and per 'distance'

MM.2.1 <-gls(f, method="ML", 
                correlation = corExp(form=~x+y),
                weights =     varIdent(form = ~ 1 | dist),
                data =        my.yield)

MM.3 <-gls(f, method="ML", 
              correlation = corExp(form=~x+y),
              weights =     varComb(varIdent(form = ~ 1 | dist), 
                                 varIdent(form = ~ 1 | direction)),
              data =        my.yield)

anova(MM, MM.2.1, MM.3)

bubble.plotting(MM.3) # lower errors
diagnostics.simple(MM.3) # residuals continue to display heterogeneity,
# residual spread might not be related to the measured explanatory variables
# fitting an alternative model i.e., with random effects is recommended

# Keep in mind:
# GLS weaknesses: effect size for GLS is usually weaker 
# especially when transects (nesting) is concerned 

# Step 4: Fitting LME: Selecting an error structure-----------------------------

# Random intercept-only model 
# [fitted with 'ML' to allow for using log likelihood ratios i.e., anova]

M0.lme=lme(wheat~1,random = ~1|row, 
                   method = "ML", data=my.yield)

summary(M0.lme) # Different fixed structure to other models; 
# don't compare with GLS

# Random intercept-only model with nested random effects
# Can be used e.g., to quantify the amount of variance explained by 
# different random effects which can lead to model simplification

M1.lme=lme(wheat~1,random = ~1|row/trans, 
                   method = "ML", 
                   data=my.yield)
summary(M1.lme)

anova(M0.lme,M1.lme) # M1 performs better than M0; p-value < 0.5, 
# the complex model has lower AIC

# Calculating variance components
SD  <- c(0.1880902, 0.3959564) # for M1.lme fitted with REML, not ML; 
# taking SD of the intercept for row and transect

variance <- SD ^ 2

100 * variance / sum(variance) 

# Transect accounts for 82% of variance in the response variable; 
# row accounts for further 18%

# Individual plots are not used as a random effect in this example
# to demonstrate the variance partitioning because there is one
# yield measurement per plot (random effects should have > 5 levels)

# Step 5: Fitting LME: A random intercept model------------------------------
# Fixed effects are included

f

M3.lme=lme(f,random = ~1|row/trans, 
             method = "ML", data=my.yield)

anova(M.gls, M3.lme) # the models are significantly different;
# lme performs better

summary(M3.lme) # t-values higher than 2 and -2 for all predictors
anova(M3.lme)

# Step 6: Fitting LME: Random slope and intercept model-----------------------

# We will consider terrain steepness as equivalent to 'direction'
# and fit a random slope model

M.4.lme=lme(wheat ~ dist * steepness.deg, 
            random = ~ 1 + steepness.deg|row, 
            method = "ML", data=my.yield)

M.5.lme=lme(wheat ~ dist * steepness.deg, 
            random = ~ 1 + steepness.deg|row/trans,
            method = "ML",
            data=my.yield)
 
# Overparameterization in model 5; hence no convergence 
# a potential sample size issue (random slope models increase the 
# sample size requirement), 
# it's possible to increase the number of model iterations with e.g.
# lmeControl(niterEM = 5000, msMaxIter = 5000).

# This model would have to go through rigorous testing, the same to what has 
# been for the other GLS and LME models, and will not be explored further
# because data exploration does not support
plot(wheat ~ steepness.deg, col=factor(row)) 
# including terrain steepness as a significant predictor i.e. there is
# no evidence for linear or nonlinear relationships between variables.

# Step 7: LME error structure parameterization ---------------------------------

# Fit the covariance structure and compare against the parameterized GLS model
# lme has default assumptions that
# 1. the covariance structure for random effects is exchangeable and
# 2. the correlation structure for residuals is independent
# These assumptions can lead to incorrect results

# Model with a specified (corExp) covariance structure
M3.Exp.lme=lme(f,random = ~1|row/trans, 
                 corr=corExp(form=~x+y|row/trans),
                 method = "ML", data=my.yield)

summary(M3.Exp.lme)

anova(M.gls, M3.Exp.lme) # M3.Exp.lme performs better
anova(M3.lme, M3.Exp.lme) # M3.Exp.lme performs better
anova(MM.2.1,M3.Exp.lme) # marginal model performs better
anova(MM.3,M3.Exp.lme) # marginal model performs better (lowest AIC)

# include other correlation structures 

M3.Gaus.lme=lme(f,random = ~1|row/trans, 
                  corr=corGaus(form=~x+y|row/trans),
                  method = "ML", data=my.yield)
M3.Lin.lme=lme(f,random = ~1|row/trans, 
                 corr=corLin(form=~x+y|row/trans),
                 method = "ML", data=my.yield)
M3.Ratio.lme=lme(f,random = ~1|row/trans, 
                   corr=corRatio(form=~x+y|row/trans),
                   method = "ML", data=my.yield) # refit with REML
M3.Spher.lme=lme(f,random = ~1|row/trans, 
                   corr=corSpher(form=~x+y|row/trans),
                   method = "ML", data=my.yield)

AIC(M3.Gaus.lme,M3.Lin.lme,M3.Ratio.lme,M3.Spher.lme)

anova(M3.lme, M3.Ratio.lme)

plot(Variogram(M3.Ratio.lme, form=~x+y)) 

# This model for the residual (in this package)
# restricts spatial covariance to points within the same transect,  
# whereas points in different transects are assumed uncorrelated

anova(MM.3, M3.Ratio.lme) # significantly different, MM.3 has lower AIC

# Allowing for different variances per explanatory variable

M3.1.lme=lme(f,   random = ~1|row/trans, 
                  weights = varIdent(form =~ 1|direction),
                  method = "ML", 
                  data = my.yield)

anova(M3.lme, M3.1.lme) # M3.1.lme is marginally better

M3.2.lme=lme(f,   random = ~1|row/trans, 
                  weights = varIdent(form = ~ 1 | dist),
                  method = "ML", 
                  data = my.yield)

anova(M3.1.lme, M3.2.lme) #  M3.2.lme is significantly better 

M4.1.Ratio.lme=lme(f,   random = ~1|row/trans, 
                 corr =    corRatio(form=~x+y|row/trans),
                 weights = varIdent(form = ~ 1 | direction),
                 method = "ML", 
                 data = my.yield)

M4.Ratio.lme=lme(f,   random = ~1|row/trans, 
                 corr =    corRatio(form=~x+y|row/trans),
                 weights = varIdent(form = ~ 1 | dist),
                 method = "ML", 
                 data = my.yield)

anova(M4.1.Ratio.lme, M4.Ratio.lme) 

anova(M3.Ratio.lme, M4.Ratio.lme) # M4.Ratio.lme is significantly better

anova(M3.2.lme, M4.Ratio.lme) # M4.Ratio.lme is significantly better

# Step 8: Selection of the best performing model and final reporting----------

# Do we have homogeneity of variance and independence of residuals?

# Extract (standardized) residuals
# Inspect the residuals 
# Are there any patterns (heteroscedasticity) , if so, improvements needed
# Residuals per predictor should be evenly distributed & with mean around 0

# Diagnostics (refit with REML before running model diagnostics)

# Original model

M.GLS <- update(M.gls, method = "REML")

# 3 best performing models

LME.1 <- update(M3.2.lme, method = "REML")
LME.2 <- update(M4.Ratio.lme, method = "REML")
MM.4  <- update(MM.3, method = "REML")

diagnostics.simple(M.GLS) # known violation of assumptions
bubble.plotting(M.GLS) # clustered errors

diagnostics.simple(LME.1) # good residual spread
bubble.plotting(LME.1) # less clustered errors

diagnostics.simple(LME.2) # some patterns in the residuals but still good
# spread per predictor
bubble.plotting(LME.2) # less clustered errors

diagnostics.simple(MM.4) # residuals vs. fitted values look okay but 
# but the residual spread per distance is problematic
bubble.plotting(MM.4) # # less clustered errors

# Residual plots should be included in e.g. Appendices

# Tip for reporting

stargazer(LME.1,LME.2,MM.4, type='text',digits=2)

# The estimates are very similar. MM.3 has the lowest AIC value but the
# diagnostic plots indicate that it should be improved (e.g., by
# adding covariates). In our case, we reject MM.3 based on diagnostic plots
# in lieu of LME.2 which has an AIC value of 359 and fairly acceptable 
# diagnostic plots. The LME.2 is a model with a parameterized covariance
# structure.

# Further checks 
# Fitting actual vs predicted wheat yields for LME.1, LME.2 and MM.4

plot(LME.1, wheat~fitted(.))
plot(LME.1, wheat~fitted(.)|row) # looks not very promising but in this example,
# hypothesis testing and elimination of any violation to assumptions is of
# primary interest

plot(LME.2, wheat~fitted(.)) 
plot(LME.2, wheat~fitted(.)|row)

plot(MM.4, wheat~fitted(.)) # residual standard error of 0.90; 
# percentage error of 25.5%
plot(MM.4, wheat~fitted(.)|row)

# In order to improve the predictive power of the model and further
# reduce the bubble plot clustering, we could: 
#                       1. Add covariates through backward or forward selection
#                          (GH1 has variable topography and site conditions)
#                       2. Add a  missing interaction term 
#                       3. Investigate a non-linear relationship 
#                          (mind we only have 3 levels for distance and
#                           2 levels for direction)

# Step 9: Avoiding hasty conclusions with visualization  -------------

# Distance from tree strips is a significant predictor (removing 'dist' makes
# the model weaker). However, concluding that it is the effect of the trees as  
# opposed to the locations of the sampling plots which drives the patterns
# observed in the data can be a hasty conclusion. 

lattice::xyplot(wheat~dist | row, groups=row, 
                data=my.yield, type=c('p','r'), auto.key=F)

# Wheat yield is highest at 4.5m - this is odd, if the trees were negatively
# affecting the crop yield (only 1 year after planting!),
# the yield should be highest in the middle of the field 

lattice::xyplot(wheat~dist | direction, groups=row, 
                data=my.yield, type=c('p','r'), auto.key=T)

# There is some evidence for the yield increasing further away from tree strips
# but not for both directions. However, the high yield at 4.5 remains a 
# consistent pattern.

# In conclusion, we fitted some extremely complex models but based them
# on a research question that does not necessarily fit our observations.
# Land management (i.e., a change in the width of the machinery to 
# accommodate tree strips is likley to have resulted in more seeds being 
# delivered to sampling locations at 4.5m away from the tree strips)
# is a more likely driver of the observed trends than the impact of trees
# in this very young agroforestry system. This highlights the importance of 
# preliminary data collection and has implications for
# future sampling schemes.

# We hope this example has served as a good introduction to using
# a subset of multilevel models in agroforestry research. 
# This script highlights the main aspects of data exploration, analysis and
# diagnostics which might be helpful to researchers who are new to
# multilevel models. The data set was provided to encourage readers of this
# Opinion piece to modify the script (e.g., re-run the models with two levels
# of the distance variable), ask alternative research questions,
# include other variance-covariance structures,
# data distributions, non-linear relationships, etc. and to think about 
# potential pitfalls of analyzing agroforestry data. 
# We are looking forward to seeing how the data analysis methods in  
# agroforestry research evolve in the upcoming years!

# Speeding up the workflow

bubble.plotting<-function(mltlvl.model) {
  resid.mltlvl<-resid(mltlvl.model, type="normalized")
  spatial.errors<-data.frame(resid.mltlvl, 
                             yield$xcoord,
                             yield$ycoord)
  coordinates(spatial.errors)<-c("yield.xcoord","yield.ycoord")
  plot.errors<- bubble(spatial.errors,"resid.mltlvl",
                       col=c("dark green","orange"),
                       main="Residuals", xlab="X-coordinates",
                       ylab="Y-coordinates") 
  return(plot.errors)
}

diagnostics.simple<-function(mltlvl.model) {
  mltlvl.errors  <- resid(mltlvl.model, type="normalized")
  mltlvl.fitted  <- fitted(mltlvl.model)
  op<-par(mfrow=c(2,2), mar=c(5,5,3,2), cex.lab = 1.5, cex.axis = 1)
  MyYlab="Residuals"
  plot(x=mltlvl.fitted,y=mltlvl.errors,xlab="Fitted values", ylab=MyYlab)
  boxplot(mltlvl.errors ~ dist, data=my.yield, 
          xlab="Distance from tree rows", 
          ylab=MyYlab)
  boxplot(mltlvl.errors ~ direction, data=my.yield,
          xlab="Direction",
          ylab=MyYlab)
  qqnorm(mltlvl.errors)
  qqline(mltlvl.errors, col = "dark green", lwd = 2)
  shapiro.test(mltlvl.errors)
}

