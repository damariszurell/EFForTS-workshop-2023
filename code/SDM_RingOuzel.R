# Example SDM workflow
# Upscaling workshop CRC-990 EFForTS Göttingen 14-Feb-2023
# (c) Damaris Zurell, Univ. Potsdam
# More elaborate materials: https://damariszurell.github.io/SDM-Intro/ and https://damariszurell.github.io/EEC-MGC/ 

#--------------------------------------------------------------------------------
#
#         Load packages
# 
#--------------------------------------------------------------------------------

library(terra)   # Manipulate geographic data
library(corrplot)   # Graphical displays of correlation matrices
library(PresenceAbsence)    # package providing performance measures for SDM
library(AUC)  # package providing performance measures for SDM
library(randomForest)   # package for calibrating random forests
library(RColorBrewer)   # package for defining colour palettes
library(lattice)    # package for advanced visualisation
library(PresenceAbsence) # package for assessing SDM performance statistics
library(mgcv) # for fitting generalised additive models, used in select07() functions

# Source additional helper functions. Also available in the R package mecofun (meant just for teaching purposes): https://gitup.uni-potsdam.de/macroecology/mecofun
source('code/SDM_helper_functions.R')




#--------------------------------------------------------------------------------
#
#         SET WORKING DRECTORY
# 
#--------------------------------------------------------------------------------

# Set your working directory to the workshop folder, e.g 
setwd('EFForTS_workshop/SDM')


#--------------------------------------------------------------------------------
#
#         TYPICAL SDM DATA
# 
#--------------------------------------------------------------------------------

# DATA

# Data stem from Citizen Science Atlas UK (https://doi.org/10.1111/geb.12906). 
# Data set contains presence-absence data of the Ring Ouzel in UK for breeding period 2008-2011. Climate data stem from worldclim. 

# Read in presence-absence data
sp_dat <- read.table('data/ATLAS_RingOuzel.txt',header=T)

# Inspect data
summary(sp_dat)

# Map data within Britsh Isle
# Read in background data defining British Isle land mass:
bg <- terra::rast('data/UK_mask.tif')

# Plot GB land mass:
plot(bg,col='grey',axes=F,legend=F)

# Plot presences in red and absences in black:
plot(extend(terra::rast(sp_dat[,1:3], crs=crs(bg), type='xyz'), bg), col=c('black','red'), legend=F,add=T)


#--------------------------------------------------------------------------------
#
#         SIMPLE SDM FITTING (GLM)
# 
#--------------------------------------------------------------------------------


#---------------
# GLM FORMULAS
#---------------

# We first fit a GLM for the bio11 variable assuming a linear relationship:
m1 <- glm(Turdus_torquatus ~ bio11, family="binomial", data= sp_dat)

# We can get a summary of the model:
summary(m1) 

# Some options for more complex model specifications
# Fit a quadratic relationship with bio11:
summary( glm(Turdus_torquatus ~ bio11 + I(bio11^2), family="binomial", data= sp_dat))

# Or use the poly() function:
summary( glm(Turdus_torquatus ~ poly(bio11,2) , family="binomial", data= sp_dat) )

# Fit two variables with second-order polynomials:
summary( glm(Turdus_torquatus ~ poly(bio11,2) + poly(bio8,2), family="binomial", data= sp_dat) )


#---------------
# COLLINEARITY
#---------------

# We first estimate a correlation matrix from the predictors. 
# We use Spearman rank correlation coefficient, as we do not know 
# whether all variables are normally distributed.
cor_mat <- cor(sp_dat[,-c(1:3)], method='spearman')

# We can visualise this correlation matrix. For better visibility, 
# we plot the correlation coefficients as percentages.
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Use select07 method to identify all pairs of variables that have correlation |rho|>0.7 and remove the less important variable
# Function described in Dormann et al. (2013): http://dx.doi.org/10.1111/j.1600-0587.2012.07348.x

# Run select07() (function contained in SDM_helper_functions.R and mecofun package):
var_sel <- select07(X=sp_dat[,-c(1:3)], 
                    y=sp_dat$Turdus_torquatus, 
                    threshold=0.7)

# Check out the structure of the resulting object:
str(var_sel)

# We extract the names of the weakly correlated predictors ordered by the univariate variable importance in terms of AIC:
pred_sel <- var_sel$pred_sel

# How many presence points do we have? Rule of thumb: you need 10 presences per parameter in the model
sum(sp_dat$Turdus_torquatus)


#----------------
# MODEL SELECTION
#----------------

# Fit the full model with 4 parameters:
m_full <- glm( Turdus_torquatus ~ bio11 + I(bio11^2) + bio8 + I(bio8^2), 
               family='binomial', data=sp_dat)

# Inspect the model:
summary(m_full)

# Explained deviance (function contained in SDM_helper_functions.R and mecofun package):
expl_deviance(obs = sp_dat$Turdus_torquatus,
              pred = m_full$fitted)

# Simplify model by AIC-based stepwise variable selection (by default in both directions)
m_step <- step(m_full) 
summary(m_step)

# Explained deviance:
expl_deviance(obs = sp_dat$Turdus_torquatus,
              pred = m_step$fitted)


# The final model only selected the linear terms for bio11 and bio8. The explained deviance is slightly lower than for the quadratic model, but the linear model is more parsimonious.




#--------------------------------------------------------------------------------
#
#         SDM ASSESSMENT
# 
#--------------------------------------------------------------------------------

#----------------
# RESPONSE CURVES
#----------------

# Partial response plots plot the predicted response along one environmental gradient while keeping all other gradients constant at their mean. 
# You can easily construct partial response plots yourself when you know how to make predictions. For simplicity, I include a function partial_response() in the SDM_helper_functions.R (also available in mecofun package).

# Assign a vector with names of our variables:
my_preds <- c('bio11', 'bio8')

# We want two panels next to each other:
par(mfrow=c(1,2))

# Plot the partial responses
partial_response(m_step, predictors = sp_dat[,my_preds], ylab="Occurrence probability", ylim=c(0,1))

# Switch back to single panel:
par(mfrow=c(1,1))

#--------------------------------------
# CROSS-VALIDATION & PREDICTIVE ACCURACY
#--------------------------------------

# We assess predictive accuracy of the model using 5-fold cross-validation.
# We split our data into 5 folds, re-calibrate the model using only 4/5 of the original data and predict the model to the hold-out 1/5 of the data. Again, for simplicity I have implemented a function crossvalSDM() in the SDM_helper_functions.R (also available in teaching package "mecofun").

# Run 5-fold cross-validation. Output is a numeric vector of cross-predictions.
preds_cv <- crossvalSDM(m_step, traindat = sp_dat, colname_species = 'Turdus_torquatus', colname_pred = my_preds)

# These cross-predictions can be used now to calculate how well the model predicts to hold-out data

#--------------------------------------
# Threshold-dependent performance measures
#--------------------------------------

# We first have find an optimal threshold to binarise the predictions = convert the predicted probabilities into predicted presences and absences. Thresholding approaches described in Liu et al. (2005): https://doi.org/10.1111/j.0906-7590.2005.03957.x

# Prepare cross-validated predictions:
thresh_dat <- data.frame(
  ID = seq_len(nrow(sp_dat)), 
  obs = sp_dat$Turdus_torquatus,
  pred = preds_cv)

# Then, we find the optimal thresholds: 
(thresh_cv <- PresenceAbsence::optimal.thresholds(DATA= thresh_dat))

# Here, we use the threshold that maximies the sum of sensitivity and specificity (=maximises the true skill statistic TSS).
# We construct the contingency table:
(cmx_maxSSS <- PresenceAbsence::cmx(DATA= thresh_dat, threshold=thresh_cv[3,2]))

# From the contingency table, we can calculate different performance measures:
# Proportion of correctly classified observations
PresenceAbsence::pcc(cmx_maxSSS, st.dev=F)

# Sensitivity = true positive rate
PresenceAbsence::sensitivity(cmx_maxSSS, st.dev=F)

# Specificity = true negative rate
PresenceAbsence::specificity(cmx_maxSSS, st.dev=F)

# Kappa
PresenceAbsence::Kappa(cmx_maxSSS, st.dev=F)

# True skill statistic (function contained in SDM_helper_functions.R and mecofun package)
TSS(cmx_maxSSS) 


#--------------------------------------
# Threshold-independent performance measures
#--------------------------------------

# Here, we only want to calculate the AUC - the area under ROC (receiver operating characteristic curve). Different packages allow calculating AUC. We here use the "AUC" package as it allows plotting the ROC curve, which may facilitate understanding.

# Let's have a look a the ROC curve. It plots the sensitivity and 1-specificity for all possible threshold values. The area under the curve is the predictive accuracy measure AUC (AUC=1 perfect discrimination, AUC=0.5 random, AUC>0.7 fair predictions)
roc_cv <- AUC::roc(preds_cv, as.factor(sp_dat$Turdus_torquatus))
plot(roc_cv, col = "grey70", lwd = 2)

# Compute the AUC:
AUC::auc(roc_cv)

# Hint: for rank-based equivalent of AUC, look at the c-index in the rcorr.cens() function in the Hmisc package.

#--------------------------------------
# All performance measures
#--------------------------------------

# For convenience, I have included a function evalSDM() in the SDM_helper_functions.R (and in the teaching package "mecofun") to produce an output with all measures mentioned above
evalSDM(sp_dat$Turdus_torquatus, preds_cv)


#--------------------------------------------------------------------------------
#
#         SDM PREDICTIONS
# 
#--------------------------------------------------------------------------------


# We can make predictions using the function predict() and the argument newdata. This function is generic and works for almost all SDM algorithms in a similar way.
# The newdata argument expects a data.frame of the environmental data. Thus, if we want to make predictions to a specific landscape, we first have to make a data frame with the environmental data for all locations. Here, we need climate layers.

# We read climate layers from file - for current and future climatic conditions. (Note, for simplicity only a single future climate layer is used here. It stems from the GCM "NorESM1-M" using RCP4.5 for 2050.)
bio_curr <- terra::rast('data/UK_bio_curr.tif')
bio_fut <- terra::rast('data/UK_bio_fut.tif')


# It is now straight forward to make continuous predictions to the current and the future climate:
# Prepare data frames
bio_curr_df <- data.frame(crds(bio_curr),as.points(bio_curr))
bio_fut_df <- data.frame(crds(bio_fut),as.points(bio_fut))

# Make continuous predictions:
bio_curr_df$pred_glm <- predict(m_step, newdata= bio_curr_df, type="response")
bio_fut_df$pred_glm <- predict(m_step, newdata= bio_fut_df, type="response")

# Map the continuous climate suitability predictions:
par(mfrow=c(1,2))

# Make raster of predictions to current environment:
r_pred_curr <- terra::rast(bio_curr_df[,c('x','y','pred_glm')], type='xyz', crs=crs(bg))
plot(r_pred_curr, axes=F, main='Occ. prob. - today')

# Make raster stack of predictions to future environment:
r_pred_fut <- terra::rast(bio_fut_df[,c('x','y','pred_glm')], type='xyz', crs=crs(bg))
plot(r_pred_fut, axes=F, main='Occ. prob. - 2050')


# translate the continuous predictions into binary predictions and plot the resulting maps

# Make binary predictions:
bio_curr_df$bin_glm <- ifelse(bio_curr_df$pred_glm >= thresh_cv[3,2], 1, 0)
bio_fut_df$bin_glm <- ifelse(bio_fut_df$pred_glm >= thresh_cv[3,2], 1, 0)

# Make raster stack of predictions to current environment:
r_pred_curr <- terra::rast(bio_curr_df[,c('x','y','pred_glm','bin_glm')], type='xyz', crs=crs(bg))
plot(r_pred_curr, axes=F)

# Make raster stack of predictions to future environment:
r_pred_fut <- terra::rast(bio_fut_df[,c('x','y','pred_glm','bin_glm')], type='xyz', crs=crs(bg))
plot(r_pred_fut, axes=F)


#--------------------------------------------------------------------------------
#
#         OTHER SDM ALGORITHMS
# 
#--------------------------------------------------------------------------------

#----------------
# RANDOM FOREST
#----------------

# Fit RF with 1000 trees (a question/warning pops up whether we really want to do regression: YES, we want to. With binary data, you could also use classification by setting the argument y=as.factor(sp_dat$Turdus_torquatus) but from experience regression produces the better results.)
m_rf <- randomForest( x=sp_dat[,my_preds], y=sp_dat$Turdus_torquatus, 
                      ntree=1000, importance =T)

# Variable importance (estimated by a permutation procedure, which measures for each variable the drop in mean accuracy when this variables is permutated:
importance(m_rf,type=1)

# Look at single trees:
head(getTree(m_rf,1,T))

# Let's plot a 3D-response surface to get a better impression of the ruggedness of the predictions
# For the response surface, we first prepare the 3D-grid with environmental gradient and predictions
xyz <- expand.grid(
  seq(min(sp_dat[,my_preds[1]]),max(sp_dat[,my_preds[1]]),length=50),
  seq(min(sp_dat[,my_preds[2]]),max(sp_dat[,my_preds[2]]),length=50))
names(xyz) <- my_preds

# Define colour palette:
cls <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))(100)

# Make predictions to gradients:
xyz$z <- predict(m_rf, xyz, type='response')

# Plot response surface:
wireframe(z ~ bio11 + bio8, data = xyz, zlab = list("Occurrence prob.", rot=90), 
          drape = TRUE, col.regions = cls, scales = list(arrows = FALSE), 
          zlim = c(0, 1), main='Random Forest', xlab='bio11', ylab='bio8', 
          screen=list(z = -120, x = -70, y = 3))


# Plot partial response curves:
par(mfrow=c(1,2)) 
partial_response(m_rf, predictors = sp_dat[,my_preds], main='Random Forest')


# Make cross-validated predictions for RF:
preds_cv_rf <- crossvalSDM(m_rf, traindat = sp_dat, colname_species = 'Turdus_torquatus', colname_pred = my_preds)

# Performance measures of RF: 
(perf_rf <- evalSDM(sp_dat$Turdus_torquatus, preds_cv_rf))


# Map predictions:
# Make continuous predictions:
bio_curr_df$pred_rf <- predict(m_rf, newdata= bio_curr_df, type="response")
bio_fut_df$pred_rf <- predict(m_rf, newdata= bio_fut_df, type="response")

# Map the continuous climate suitability predictions:
par(mfrow=c(1,2))

# Make raster of predictions to current environment:
r_pred_curr_rf <- terra::rast(bio_curr_df[,c('x','y','pred_rf')], type='xyz', crs=crs(bg))
plot(r_pred_curr_rf, axes=F, main='RF: Occ. prob. - today')

# Make raster stack of predictions to future environment:
r_pred_fut_rf <- terra::rast(bio_fut_df[,c('x','y','pred_rf')], type='xyz', crs=crs(bg))
plot(r_pred_fut_rf, axes=F, main='RF: Occ. prob. - 2050')

# Make binary predictions:
bio_curr_df$bin_rf <- ifelse(bio_curr_df$pred_rf >= perf_rf$thresh, 1, 0)
bio_fut_df$bin_rf <- ifelse(bio_fut_df$pred_rf >= perf_rf$thresh, 1, 0)

# Make raster stack of predictions to current environment:
r_pred_curr_rf <- terra::rast(bio_curr_df[,c('x','y','pred_rf','bin_rf')], type='xyz', crs=crs(bg))
plot(r_pred_curr_rf, axes=F)

# Make raster stack of predictions to future environment:
r_pred_fut_rf <- terra::rast(bio_fut_df[,c('x','y','pred_rf','bin_rf')], type='xyz', crs=crs(bg))
plot(r_pred_fut_rf, axes=F)
