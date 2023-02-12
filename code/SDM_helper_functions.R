# helper functions for SDM workshop, (c) Damaris Zurell, damaris@zurell.de
# Upscaling workshop CRC-990 EFForTS Göttingen 14-Feb-2023
# Functions also available in R package mecofun (only for teaching purposes): https://gitup.uni-potsdam.de/macroecology/mecofun.git

# Select 07
# Select weakly correlated variables based on univariate importance based on Dormann et al. (2013). Univariate variable importance is based on AIC. Variable importance can also be pre-defined by hand.


select07 <- function (X, y, family = "binomial", univar = "glm2", threshold = 0.7, 
                      method = "spearman", sequence = NULL, weights = NULL) 
{
  require(mgcv)
  # selects variables based on removing correlations > 0.7, retaining those
  # variables more important with respect to y
  # Order of importance can be provided by the character vector 'sequence'
  
  # 1. step: cor-matrix
  # 2. step: importance vector
  # 3. step: identify correlated pairs
  # 4. step: in order of importance: remove collinear less important variable,
  #           recalculate correlation matrix a.s.f.
  
  var.imp <- function(variable, response, univar, family, weights) {
    m1 <- switch(univar, glm1 = glm(response ~ variable, 
                                    family = family, weights = weights), 
                         glm2 = glm(response ~ poly(variable, 2), family = family, weights = weights), 
                         gam = mgcv::gam(response ~ s(variable, k = 4), family = family, 
                                 weights = weights))
    AIC(m1)
  }
  cm <- cor(X, method = method)
  if (is.null(sequence)) {
    a <- try(var.imp(X[, 1], y, univar = univar, family = family, 
                     weights = weights))
    if (is.numeric(a) != 1) {
      stop("invalid univar method")
    }
    imp <- apply(X, 2, var.imp, response = y, family = family, 
                 univar = univar, weights = weights)
    sort.imp <- names(sort(imp))
  }
  else {
    sort.imp <- sequence
  }
  pairs <- which(abs(cm) >= threshold, arr.ind = T)
  index <- which(pairs[, 1] == pairs[, 2])
  pairs <- pairs[-index, ]
  exclude <- NULL
  for (i in 1:length(sort.imp)) {
    if ((sort.imp[i] %in% row.names(pairs)) & ((sort.imp[i] %in% 
                                                exclude) == F)) {
      cv <- cm[setdiff(row.names(cm), exclude), sort.imp[i]]
      cv <- cv[setdiff(names(cv), sort.imp[1:i])]
      exclude <- c(exclude, names(which((abs(cv) >= threshold))))
    }
  }
  pred_sel <- sort.imp[!(sort.imp %in% unique(exclude)), drop = F]
  return(list(AIC = sort(imp), cor_mat = cm, pred_sel = pred_sel))
}


#---------------------------

# Explained Deviance

expl_deviance <- function(obs, pred, family='binomial'){
  if (family == "binomial") {
    pred <- ifelse(pred < 1e-05, 1e-05, ifelse(pred > 0.9999, 0.9999, pred))
  }
  null_pred <- rep(mean(obs), length(obs))
  1 - (calc.deviance(obs, pred, family = family)/calc.deviance(obs, null_pred, family = family))
}


#---------------------

# calc.deviance function from dismo package
calc.deviance <- function (obs, pred, weights = rep(1, length(obs)), family = "binomial", 
          calc.mean = TRUE) 
{
  if (length(obs) != length(pred)) {
    stop("observations and predictions must be of equal length")
  }
  y_i <- obs
  u_i <- pred
  family = tolower(family)
  if (family == "binomial" | family == "bernoulli") {
    deviance.contribs <- (y_i * log(u_i)) + ((1 - y_i) * 
                                               log(1 - u_i))
    deviance <- -2 * sum(deviance.contribs * weights)
  }
  else if (family == "poisson") {
    deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) - 
      (y_i - u_i)
    deviance <- 2 * sum(deviance.contribs * weights)
  }
  else if (family == "laplace") {
    deviance <- sum(abs(y_i - u_i))
  }
  else if (family == "gaussian") {
    deviance <- sum((y_i - u_i) * (y_i - u_i))
  }
  else {
    stop("unknown family, should be one of: \"binomial\", \"bernoulli\", \"poisson\", \"laplace\", \"gaussian\"")
  }
  if (calc.mean) 
    deviance <- deviance/length(obs)
  return(deviance)
}

#---------------------

# True skill statistic

TSS = function(cmx){
  require(PresenceAbsence)
  PresenceAbsence::sensitivity(cmx, st.dev=F) + 
    PresenceAbsence::specificity(cmx, st.dev=F) - 1
}


#---------------------

# evaluation statistics
evalSDM <- function (observation, predictions, thresh.method = "MaxSens+Spec", 
                     req.sens = 0.85, req.spec = 0.85, FPC = 1, FNC = 1) 
  {
  thresh.dat <- data.frame(ID = seq_len(length(observation)), 
                           obs = observation, pred = predictions)
  thresh <- PresenceAbsence::optimal.thresholds(DATA = thresh.dat, 
                                                req.sens = req.sens, req.spec = req.spec, FPC = FPC, 
                                                FNC = FNC)
  cmx.opt <- PresenceAbsence::cmx(DATA = thresh.dat, threshold = thresh[thresh$Method == 
                                                                          thresh.method, 2])
  data.frame(AUC = PresenceAbsence::auc(thresh.dat, st.dev = F), 
             TSS = TSS(cmx.opt), Kappa = PresenceAbsence::Kappa(cmx.opt, st.dev = F), 
             Sens = PresenceAbsence::sensitivity(cmx.opt, st.dev = F), 
             Spec = PresenceAbsence::specificity(cmx.opt, st.dev = F), 
             PCC = PresenceAbsence::pcc(cmx.opt, st.dev = F), 
             D2 = expl_deviance(observation, predictions),thresh = thresh[thresh$Method == thresh.method, 2])
}

#-------------------------

# Make predictions
# Here, just implemented for GLM and RandomForest. More algorithms available in mecofun package (only for teaching purposes): https://gitup.uni-potsdam.de/macroecology/mecofun.git

predictSDM <- function(model, newdata) {
  require(randomForest)
  
  switch(class(model)[1],  
         glm = predict(model, newdata, type = "response"), 
         randomForest.formula = switch(model$type, 
                                       regression = predict(model, newdata, type = "response"),
                                       classification = predict(model, newdata, type = "prob")[,2]),
         randomForest = switch(model$type, 
                               regression = predict(model, newdata, type = "response"),
                               classification = predict(model, newdata, type = "prob")[, 2]))
}

#------------------------

# cross validation
# Here, just implemented for GLM and RandomForest. More algorithms available in mecofun package (only for teaching purposes): https://gitup.uni-potsdam.de/macroecology/mecofun.git


crossvalSDM <- function(model, kfold=5, traindat, colname_species, colname_pred,
                        env_r=NULL, colname_coord=NULL, weights=NULL)
{
  
  require(randomForest)
  
  weights.full <- weights
  
  if (length(kfold)==1) {
    # Make k-fold data partitions
    ks <- dismo::kfold(traindat, k = kfold, by = traindat[,colname_species])
  } else {
    ks <- kfold
    kfold <- length(unique(kfold))
  }
  
  cross_val_preds = numeric(length = nrow(traindat))
  
  for(i in seq_len(kfold)){
    cv_train <- traindat[ks!=i,]
    cv_test <- traindat[ks==i,]
    
    if (!is.null(weights)) {
      weights <- weights.full[ks!=i]
    }
    
    # We update the model for the new training data
    modtmp <- switch(class(model)[1],
                     glm = update(model, data=cv_train),
                     randomForest = update(model, data=cv_train),    
                     randomForest.formula = update(model, data=cv_train))
    
    # We make predictions for k-fold:
    cross_val_preds[ks==i] <- predictSDM(modtmp, cv_test[, colname_pred, drop=F])
  }
  cross_val_preds
  
}

#--------------------------

# Inflated response curves

inflated_response=function(object,predictors,select.columns=NULL,label=NULL,len=50,lhsample=100,lwd=1,
      ylab="Occurrence probabilities",method="stat3",disp="all",overlay.mean=T,
      col.curves='grey',col.novel='grey',col.mean='black',lwd.known=2,lwd.mean=2,ylim=c(0,1),...)
{
  
  require(lhs)
  
  if (is.null(select.columns)) select.columns=seq_len(ncol(predictors))
  
  for (i in select.columns)
  {
    summaries=data.frame(matrix(0,6,ncol(predictors)))
    for (iz in 1:ncol(predictors)) {
      summaries[,iz]=summary(predictors[,iz])
    }
    if (method=="stat3") {
      summaries.j=as.matrix(summaries[c(1,4,6),-i],ncol=(ncol(predictors)-1));comb=min(lhsample,3^(ncol(predictors)-1));nc=3
    } else
      if (method=="stat6") {
        summaries.j=as.matrix(summaries[,-i],ncol=(ncol(predictors)-1));comb=min(lhsample,6^(ncol(predictors)-1));nc=6
      } else
        if (method=="mean") {
          summaries.j=as.matrix(summaries[4,-i],ncol=(ncol(predictors)-1));comb=1;nc=1;overlay.mean=F
        }
    
    dummy.j=as.matrix(predictors[1:len,-i],ncol=(ncol(predictors)-1))
    
    if (comb<lhsample) {
      mat=vector("list",ncol(dummy.j))
      for (m in 1:ncol(dummy.j)) mat[[m]]=1:nc
      mat=expand.grid(mat)
    } else {
      mat=round(qunif(lhs::randomLHS(lhsample,ncol(dummy.j)),1,nrow(summaries.j)),0)
    }
    
    if (is.null(label)) {
      label=names(predictors)
    }
    
    for (r in 1:nrow(mat))
    {
      for (j in 1:ncol(dummy.j))
      {
        dummy.j[,j]=as.vector(rep(summaries.j[mat[r,j],j],len))
      }
      
      dummy=data.frame(seq(min(predictors[,i]),max(predictors[,i]),length=len),dummy.j)
      names(dummy)[-1]=names(predictors)[-i]
      names(dummy)[1]=names(predictors)[i]
      
      curves <- predictSDM(object, dummy)
      
      # display all lines in same type
      if (disp=='all')
      {
        if (r==1)
        {
          if (i==1) plot(dummy[,names(predictors)[i]],
                         curves,type="l",ylim=ylim,xlab=label[i],ylab=ylab,
                         lwd=lwd,col=col.curves,...)
          else plot(dummy[,names(predictors)[i]],
                    curves,type="l",ylim=ylim,xlab=label[i],ylab="",lwd=lwd,col=col.curves,...)
        }
        else lines(dummy[,names(predictors)[i]],
                   curves,lwd=lwd,col=col.curves,...)
      }
      
      # highlight extrapolation to novel environmental conditions
      if (disp=='eo.mask')
      {
        novel=eo.mask(predictors,dummy)
        curves.known=curves
        curves.known[novel==1]=NA
        curves.novel=curves
        curves.novel[novel==0]=NA
        
        if (r==1)
        {
          if (i==1) {plot(dummy[,names(predictors)[i]],
                          curves.known,type="l",ylim=ylim,xlab=label[i],ylab=ylab,
                          lwd=lwd.known,col=col.curves,...)
            lines(dummy[,names(predictors)[i]],
                  curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
          else {plot(dummy[,names(predictors)[i]],
                     curves.known,type="l",ylim=ylim,xlab=label[i],ylab="",lwd=lwd.known,col=col.curves,...)
            lines(dummy[,names(predictors)[i]],
                  curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
        }
        else {lines(dummy[,names(predictors)[i]],
                    curves.known,lwd=lwd.known,col=col.curves,...)
          lines(dummy[,names(predictors)[i]],
                curves.novel,lwd=lwd,col=col.novel,lty='dotted',...)}
      }
    }
    
    #-------------------------------------------------
    # now, this is for overlaying mean response curve
    if (overlay.mean==T)
    {
      dummy=predictors[1:len,]
      dummy[,i]=seq(min(predictors[,i]),max(predictors[,i]),length=len)
      for (j in 1:ncol(predictors))
      {
        if (j!=i) 
        {
          dummy[,j]=rep(mean(predictors[,j]),len)
        }
      }
      
      curves <- predictSDM(object, dummy)
      
      lines(dummy[,names(predictors)[i]],
            curves,lwd=lwd.mean,col=col.mean,...)
    }    
  }
}


#----------------------

# partial response plots

partial_response=function(object,predictors,select.columns=NULL, label=NULL, len=50,
                          ylab="Occurrence probability", col='black',...)
{
  inflated_response(object,predictors,select.columns,label,len,method='mean',col.curves=col, ylab=ylab, ...)
}
