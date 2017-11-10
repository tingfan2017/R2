


  # 1.train_partition
  train_partition <- function (data, proportion, outcome)
  {
    label = as.factor(get(outcome, data))
    classes = lapply(levels(label), function(x) which(label ==
                                                        x))
    train = lapply(classes, function(class) sample(class, proportion *
                                                     length(class), replace = F))
    train = unlist(train)
    test = (1:nrow(data))[-train]
    training = data[train, ]
    testing = data[test, ]
    return(list(train_set = training, test_set = testing))
  }

  # 2.result
  result <- function(model,
                     train_data,
                     test_data,
                     level_train,
                     level_test)
  {
    # no level column in train_data or test_data
    library(pROC)
    options(digits = 3)
    train_data <- as.data.frame(train_data)
    test_data  <- as.data.frame(test_data)
    level_train <- as.factor(level_train)
    level_test  <- as.factor(level_test)
    # - outcome
    pred_train <- predict(model,type="response")
    pred_test  <- predict(model,test_data,type="response")
    roc_train  <- roc(level_train,pred_train)
    roc_test   <- roc(level_test ,pred_test)
    # threshold/spe/auc
    c_train <- coords(roc_train,'best'); (auc_train  <- roc_train$auc)
    c_test  <- coords(roc_test ,'best'); (auc_test   <- roc_test$auc)
    threshold_train <- c_train[1]
    threshold_test  <- c_test [1]
    Specificity_train <- c_train[2]
    Specificity_test  <- c_test [2]
    Sensitivity_train <- c_train[3]
    Sensitivity_test  <- c_test [3]
    # acc
    z1  <- ifelse( pred_train < threshold_train,0,1)
    z2  <- ifelse( pred_test  < threshold_test ,0,1)
    zz1 <- table(level_train,z1,dnn =c("Actual","Predicted"))
    zz2 <- table(level_test ,z2,dnn =c("Actual","Predicted"))
    acc_train <- sum(diag(zz1))/sum(zz1)
    acc_test  <- sum(diag(zz2))/sum(zz2)

    df <- data.frame(Train_set = c(auc_train,
                                   acc_train,
                                   Specificity_train,
                                   Sensitivity_train,
                                   threshold_train) ,
                     Test_set  = c(auc_test,
                                   acc_test,
                                   Specificity_test,
                                   Sensitivity_test,
                                   threshold_test) ,
                     row.names = c("AUC",
                                   "ACC",
                                   "Specificity",
                                   "Sensitivity",
                                   "Threshold"))
    Confusion_Matrix_train <- zz1
    Confusion_Matrix_test  <- zz2
    out <- list(CM_Train = zz1,
                CM_Test  = zz2,
                More_information = df,
                feature_numbers = ncol(train_data))

    return(out)
  }

  glm_cv <- function()

  # 3.reduce_redundency
  reduce_redundency <- function (dat, threshold = 0.9, method = "spearman")
  {

    if (!("data.frame" %in% class(dat))) {
      stop("Input data must be class data.frame")
    }
    if (sum(is.na(dat)) != 0) {
      stop("Input data contain missing values")
    }
    feature_names = colnames(dat)
    types = sapply(dat, class)
    dataIndex = which(types == "numeric" | types == "integer")
    despIndex = which(types != "numeric" & types != "integer")
    data_numeric = dat[, dataIndex]
    if (length(despIndex) > 0) {
      desp = dat[, despIndex]
    }
    else {
      desp = NULL
    }
    cor = cor(data_numeric, method = "spearman")
    cor[upper.tri(cor)] = 0
    diag(cor) = 0
    dat.redd = data_numeric[, !apply(cor, 2, function(x) any(abs(x) >
                                                               threshold))]
    features_selected = colnames(dat.redd)
    if (length(despIndex) > 0) {
      dat.redd = data.frame(desp, dat.redd)
    }
    return(list(names = features_selected, dat.redd = dat.redd))
  }

  # 4.Radiomics score
  radscore <- function(model,train_set,test_set){
    {
      # model : logistic by AK features
      # train_set: AK features
      Data  <- data.frame(rbind(train_set,test_set))
      xn <- dim(Data)[1]      # row
      yn <- dim(Data)[2]      # column
      beta           <- coef(model,s="Coefficients")   # get beta = Coefficients
      betai_Matrix   <- as.matrix( beta[-1] )      # get betai
      beta0_Matrix   <- matrix(beta[1], xn, 1 )    # get beta0
      Rad_Matrix     <- as.matrix( Data )          # get X
      Radcore_Matrix <- Rad_Matrix %*% betai_Matrix + beta0_Matrix  # get Radscore
      radscore_all   <- as.numeric(Radcore_Matrix)
      Radscore_train <- radscore_all[1:dim(train_set)[1]]    # train rads
      Radscore_test  <- radscore_all[c( (dim(train_set)[1]+1) :xn)] # test rads

      return(list(Radscore_train = Radscore_train,
                  Radscore_test  = Radscore_test))

    }
  }

  # 5.DCA
  dca <- function(data, outcome, predictors,
                  xstart=0.01, xstop=0.99,
                  xby=0.01, ymin=-0.05,
                  probability=NULL,
                  harm=NULL,
                  graph=TRUE,
                  intervention=FALSE,
                  interventionper=100,
                  smooth=FALSE,
                  loess.span=0.10)
  {

    # LOADING REQUIRED LIBRARIES
    require(stats)

    # data MUST BE A DATA FRAME
    if (class(data)!="data.frame") {
      stop("Input data must be class data.frame")
    }

    #ONLY KEEPING COMPLETE CASES
    data=data[complete.cases(data[append(outcome,predictors)]),append(outcome,predictors)]

    # outcome MUST BE CODED AS 0 AND 1
    if (max(data[[outcome]])>1 | min(data[[outcome]])<0) {
      stop("outcome cannot be less than 0 or greater than 1")
    }
    # xstart IS BETWEEN 0 AND 1
    if (xstart<0 | xstart>1) {
      stop("xstart must lie between 0 and 1")
    }

    # xstop IS BETWEEN 0 AND 1
    if (xstop<0 | xstop>1) {
      stop("xstop must lie between 0 and 1")
    }

    # xby IS BETWEEN 0 AND 1
    if (xby<=0 | xby>=1) {
      stop("xby must lie between 0 and 1")
    }

    # xstart IS BEFORE xstop
    if (xstart>=xstop) {
      stop("xstop must be larger than xstart")
    }

    #STORING THE NUMBER OF PREDICTORS SPECIFIED
    pred.n=length(predictors)

    #IF probability SPECIFIED ENSURING THAT EACH PREDICTOR IS INDICATED AS A YES OR NO
    if (length(probability)>0 & pred.n!=length(probability)) {
      stop("Number of probabilities specified must be the same as the number of predictors being checked.")
    }

    #IF harm SPECIFIED ENSURING THAT EACH PREDICTOR HAS A SPECIFIED HARM
    if (length(harm)>0 & pred.n!=length(harm)) {
      stop("Number of harms specified must be the same as the number of predictors being checked.")
    }

    #INITIALIZING DEFAULT VALUES FOR PROBABILITES AND HARMS IF NOT SPECIFIED
    if (length(harm)==0) {
      harm=rep(0,pred.n)
    }
    if (length(probability)==0) {
      probability=rep(TRUE,pred.n)
    }


    #CHECKING THAT EACH probability ELEMENT IS EQUAL TO YES OR NO,
    #AND CHECKING THAT PROBABILITIES ARE BETWEEN 0 and 1
    #IF NOT A PROB THEN CONVERTING WITH A LOGISTIC REGRESSION
    for(m in 1:pred.n) {
      if (probability[m]!=TRUE & probability[m]!=FALSE) {
        stop("Each element of probability vector must be TRUE or FALSE")
      }
      if (probability[m]==TRUE & (max(data[predictors[m]])>1 | min(data[predictors[m]])<0)) {
        stop(paste(predictors[m],"must be between 0 and 1 OR sepcified as a non-probability in the probability option",sep=" "))
      }
      if(probability[m]==FALSE) {
        model=NULL
        pred=NULL
        model=glm(data.matrix(data[outcome]) ~ data.matrix(data[predictors[m]]), family=binomial("logit"))
        pred=data.frame(model$fitted.values)
        pred=data.frame(pred)
        names(pred)=predictors[m]
        data=cbind(data[names(data)!=predictors[m]],pred)
        print(paste(predictors[m],"converted to a probability with logistic regression. Due to linearity assumption, miscalibration may occur.",sep=" "))
      }
    }

    # THE PREDICTOR NAMES CANNOT BE EQUAL TO all OR none.
    if (length(predictors[predictors=="all" | predictors=="none"])) {
      stop("Prediction names cannot be equal to all or none.")
    }

    #########  CALCULATING NET BENEFIT   #########
    N=dim(data)[1]
    event.rate=colMeans(data[outcome])

    # CREATING DATAFRAME THAT IS ONE LINE PER THRESHOLD PER all AND none STRATEGY
    nb=data.frame(seq(from=xstart, to=xstop, by=xby))
    names(nb)="threshold"
    interv=nb

    nb["all"]=event.rate - (1-event.rate)*nb$threshold/(1-nb$threshold)
    nb["none"]=0

    # CYCLING THROUGH EACH PREDICTOR AND CALCULATING NET BENEFIT
    for(m in 1:pred.n){
      for(t in 1:length(nb$threshold)){
        # COUNTING TRUE POSITIVES AT EACH THRESHOLD
        tp=mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome])*sum(data[[predictors[m]]]>=nb$threshold[t])
        # COUNTING FALSE POSITIVES AT EACH THRESHOLD
        fp=(1-mean(data[data[[predictors[m]]]>=nb$threshold[t],outcome]))*sum(data[[predictors[m]]]>=nb$threshold[t])
        #setting TP and FP to 0 if no observations meet threshold prob.
        if (sum(data[[predictors[m]]]>=nb$threshold[t])==0) {
          tp=0
          fp=0
        }

        # CALCULATING NET BENEFIT
        nb[t,predictors[m]]=tp/N - fp/N*(nb$threshold[t]/(1-nb$threshold[t])) - harm[m]
      }
      interv[predictors[m]]=(nb[predictors[m]] - nb["all"])*interventionper/(interv$threshold/(1-interv$threshold))
    }

    # CYCLING THROUGH EACH PREDICTOR AND SMOOTH NET BENEFIT AND INTERVENTIONS AVOIDED
    for(m in 1:pred.n) {
      if (smooth==TRUE){
        lws=loess(data.matrix(nb[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(nb[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
        nb[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted

        lws=loess(data.matrix(interv[!is.na(nb[[predictors[m]]]),predictors[m]]) ~ data.matrix(interv[!is.na(nb[[predictors[m]]]),"threshold"]),span=loess.span)
        interv[!is.na(nb[[predictors[m]]]),paste(predictors[m],"_sm",sep="")]=lws$fitted
      }
    }

    # PLOTTING GRAPH IF REQUESTED
    if (graph==TRUE) {
      require(graphics)

      # PLOTTING INTERVENTIONS AVOIDED IF REQUESTED
      if(intervention==TRUE) {
        # initialize the legend label, color, and width using the standard specs of the none and all lines
        legendlabel <- NULL
        legendcolor <- NULL
        legendwidth <- NULL
        legendpattern <- NULL

        #getting maximum number of avoided interventions
        ymax=max(interv[predictors],na.rm = TRUE)

        #INITIALIZING EMPTY PLOT WITH LABELS
        plot(x=nb$threshold, y=nb$all, type="n" ,xlim=c(xstart, xstop), ylim=c(ymin, ymax),
             xlab="Threshold probability", ylab=paste("Net reduction in interventions per",interventionper,"patients"))

        #PLOTTING INTERVENTIONS AVOIDED FOR EACH PREDICTOR
        for(m in 1:pred.n) {
          if (smooth==TRUE){
            lines(interv$threshold,data.matrix(interv[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
          } else {
            lines(interv$threshold,data.matrix(interv[predictors[m]]),col=m,lty=2)
          }

          # adding each model to the legend
          legendlabel <- c(legendlabel, predictors[m])
          legendcolor <- c(legendcolor, m)
          legendwidth <- c(legendwidth, 1)
          legendpattern <- c(legendpattern, 2)
        }
      } else {
        # PLOTTING NET BENEFIT IF REQUESTED

        # initialize the legend label, color, and width using the standard specs of the none and all lines
        legendlabel <- c("None", "All")
        legendcolor <- c(17, 8)
        legendwidth <- c(2, 2)
        legendpattern <- c(1, 1)

        #getting maximum net benefit
        ymax=max(nb[names(nb)!="threshold"],na.rm = TRUE)

        # inializing new benfit plot with treat all option
        plot(x=nb$threshold, y=nb$all, type="l", col=8, lwd=2 ,
             xlim=c(xstart, xstop), ylim=c(ymin, ymax),
             xlab="Threshold probability", ylab="Net benefit",
             cex.lab = 1.2,cex.axis=1.2)
        axis(side = 1,at=seq(0,0.9,0.1),cex.axis=1.2)
        # adding treat none option
        lines(x=nb$threshold, y=nb$none,lwd=2,lty=1,col=9)
        #PLOTTING net benefit FOR EACH PREDICTOR
        for(m in 1:pred.n) {
          if (smooth==TRUE){
            lines(nb$threshold,data.matrix(nb[paste(predictors[m],"_sm",sep="")]),col=m,lty=2)
          } else {
            lines(nb$threshold,data.matrix(nb[predictors[m]]),lty=1,lwd=2,type="l",col=m+3)
          }
          # adding each model to the legend
          legendlabel <- c(legendlabel, predictors[m])
          legendcolor <- c(legendcolor, m+3)
          legendwidth <- c(legendwidth, 2)
          legendpattern <- c(legendpattern, 1)
        }
      }
      # then add the legend
      legend(0,0.4,
             legendlabel,
             cex = 1,
             col=legendcolor,
             lwd=legendwidth,
             lty=legendpattern)

    }

    #RETURNING RESULTS
    results=list()
    results$N=N
    results$predictors=data.frame(cbind(predictors,harm,probability))
    names(results$predictors)=c("predictor","harm.applied","probability")
    results$interventions.avoided.per=interventionper
    results$net.benefit=nb
    results$interventions.avoided=interv

    return(results)

  }
 # sample
  {

    # dca(data = dca_df,
    #     outcome = "Pathology",
    #     predictors = c("Radiomics nomogram",
    #                    "Clinicopathological nomogram"),
    #     xby=0.01, ymin=-0.1,
    #     xstart=0.01, xstop=0.7)


    # - DCA
    # dca_df <- data.frame(level      =  c(1,1,1,0,0,0,1,0,1,0),
    #                      WithoutRad =  c(0,0,0,0,0,0,1,1,1,1),
    #                      Rad        =  c(1,0,1,1,1,0,1,0,1,0)  )
    # dca(data = dca_df,
    #     outcome = "level",
    #     predictors = c("WithoutRad","Rad"),
    #     xby=0.01, ymin=-0.18)







    # dca_df <- data.frame(MVI  = as.numeric(response_test)-1,
    #                      Radiomics.logistic = pred_test,
    #                      Clinical.APE_NSTM_PHI = pred_test1,
    #                      Clinical.APE_NSTM_PHI_any_two_or_three = pred_test2)
    # names(dca_df) <- c("MVI",
    #                    "Radiomics model",
    #                    "Expertise APE NSTM PHI",
    #                    "Expertise any two of three")
    #
    #
    # p1 <- dca(data = dca_df,
    #           outcome = "MVI",
    #           predictors = c("Radiomics model",
    #                          "Expertise APE NSTM PHI",
    #                          "Expertise any two of three"),
    #           xby=0.01, ymin=-0.18)


  }













