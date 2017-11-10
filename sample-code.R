 .

 .

 install.packages("verification")
 install.packages("pROC")
 install.packages("rms")
 install.packages("glmnet")

 library(verification)
 library(pROC)    # ROC
 library(rms)     # nomogram
 library(glmnet)  # cv.glmnet


 .
 #--------------- Environment settings -----------------
 rm(list=ls())
 setwd("C:/Users/502711516/Desktop")
 getwd()
 #---------------- data import-------------
 df <- read.csv("sample1.csv",header = T)
 View(df)
 dim(df)   # 205  1+64

 #---------------- cleaning & distribution ----------------
 # deal NAN value
 {
   #Determine any NA value
   anyNA(df)
   length(which(is.na(df)==T))
   #find num of NA
   q <- c()
   for(i in 1:dim(df)[1]) { q[i] <- length(which(is.na(df[i,]==T))) }
   index <- which(q>0)
   NAN <- data.frame(index,q[index])
   names(NAN) <- c("Index","the num of NAN")
   NAN
   # write.csv(NAN,"NAN.csv")       # outcome
   # - Delete unrelated columns
   y0 <- 2
   yn <- dim(df)[2]
   # - Replace by median
   for(i in y0:yn){ df[,i][which( is.na(df[,i]) ) ] <- median(df[,i],na.rm = T) }
   # - check any NA again
   anyNA(df)  # no NA~~
   NAN
 }

 # data classification
 {
   set.seed( 1010 )
   t = train_partition(df, 0.7, names(df[1]) )
   train = t$train_set[,-1]
   test  = t$test_set [,-1]
   response_train <- as.factor(t$train_set[,1] )   # train level
   response_test  <- as.factor(t$test_set [,1] )   #  test level
   dim(train);dim(test);length(response_train);length(response_test)

   train_index <- row.names(t$train_set)
   test_index  <- row.names(t$test_set)
   data_classification <- list(train_num = length(train_index),
                               test_num  = length(test_index),
                               train_index = train_index,
                               test_index  = test_index)
   cs <- data.frame(sum = c(dim(df)[1],
                            length(response_train),
                            length(response_test)),
                    x0  = c(length(which(   df[,1]==0)),
                            length(which(response_train==0)),
                            length(which(response_test==0))),
                    x1  = c(length(which(   df[,1]==1)),
                            length(which(response_train==1)),
                            length(which(response_test==1))))
   row.names(cs) <- c("Data ","Train ","Test ")
   colnames(cs)  <- c("SUM","0","1")
   data_classification ; cs
 }

 # train data：train 136 * 328
 # test  data：test   59 * 328

 #-------------- dimensionality reduction ------------------

 #- revome ONE value & Factor
 {
   trn <- train
   trn_var <- diag( var(trn) )
   trn_rm <- trn[ , which(trn_var != 0) ]
   dim(trn_rm)
 } # trn_rm,v=328

 #- Anova/KW test
 {
   train_gg <- trn_rm
   # Normality test
   ps0  <- apply(train_gg[which(response_train == "0"),], 2,
                 function(x)  shapiro.test(   x  )$p.value)
   ps1  <- apply(train_gg[which(response_train == "1"),], 2,
                 function(x)  shapiro.test(   x  )$p.value)
   w0  <- which(ps0 >= 0.05 & ps1 >= 0.05);length(w0)# 两类同时满足正态性
   pb  <- apply(train_gg,2,function(x) bartlett.test(x,response_train)$p.value)
   w1  <- which(ps0 >= 0.05 & ps1 >= 0.05 & pb <= 0.05)
   l1  <- length(w1);l1   # 正态+方差不齐
   v1 <- names(w1)
   train_k  <- train_gg[!(names(train_gg) %in% names(w1))]

   wilcox.test(x = train_k[,1][which(response_train == "0")],
               y = train_k[,1][which(response_train == "1")])$p.value
   pw <- apply(train_k,2,function(x) wilcox.test(x[which(response_train == "0")],
                                                      x[which(response_train == "1")])$p.value)
   w2  <- which(pw <= 0.05)
   l2  <- length(w2);l2
   v2   <-  names(w2)
   vkw <-  names(sort(c(w1,w2)))
   train_kw   <- train_gg[vkw]
   dim(train_kw)
   # names(train_kw)
 } # train_kw,v=291

 #- GLM (binomial) test
 {
   glm_test <- train_kw
   yn <- dim(glm_test)[2]
   p_value <- c()
   for(n in 1:yn) { p_value[n] <- summary(glm(response_train ~ glm_test[,n],
                                              family= binomial))$coefficients[2,4]} # glm_test[,20]
   # p_value ;p_value[which(p_value <= 0.05)]
   trainname <- names(glm_test)[which(p_value <= 0.05)]
   length(trainname)
   train_glm <- glm_test[trainname]
   dim(train_glm)
   # names(train_glm)
 } # train_glm,v=254

 #- lasso test
 {
   set.seed( 33211 )
   cv.x <- as.matrix( train_glm )
   cv.y <- response_train
   lasso.model <- cv.glmnet(x = cv.x,
                            y = cv.y,
                            family="binomial",
                            type.measure = "deviance",
                            alpha=1)
   plot(lasso.model)
   coefPara <- coef(lasso.model,s="lambda.min")
   lasso_values <- as.data.frame( which(coefPara !=0,arr.ind = T) )
   lasso_values_spm <- rownames(lasso_values)[-1]
   train_glm_lasso <- train_kw[lasso_values_spm]
   dim(train_glm_lasso)
 } # train_glm_lasso,v=9

 #- SPM test:train_glm/train_glm_lasso
 {
   spm_test  <- train_glm
   train_spm <- data.frame(reduce_redundency(spm_test,
                                             threshold = 0.9)$dat.redd)
   dim(train_spm)
   # names(train_spm)
 } # train_spm,v = 36

 # - PCA: train_kw,train_glm,train_glm_lasso
 {
   train_p <- train_spm
   test_p  <- test[names(train_p)]         # get test_p
   dim(train_p);dim(test_p)
   mvi_pr  <- princomp(train_p,cor = T)    # use train_p to make pca model,get scores
   sum <- summary(mvi_pr,loadings = T)
   # screeplot(mvi_pr,type="lines")
   train_matrix <- as.matrix(scale(train_p))
   test_matrix  <- as.matrix(scale(test_p))                     # 37*23
   pca_coefficient_matrix <- as.matrix(sum$loadings)   # 23*12
   # function
   train_pca_matrix <- train_matrix %*% pca_coefficient_matrix
   test_pca_matrix  <- test_matrix  %*% pca_coefficient_matrix

   train_pca <- data.frame(train_pca_matrix)          # logistic train set
   test_pca   <- data.frame(test_pca_matrix)          # logistic test  set
   summary(mvi_pr,loadings = F)
 } # train_pca,test_pca,v=20,p:4



 #----------------- model building ---------------------------

 #- glm_spm model
 {
   train_g <- train_spm
   test_g  <- test[names(train_g)]
   glm_fit <- glm(response_train ~ .,
                  family= binomial,
                  data=train_g,
                  control=list(maxit=1000))
   sum_glm <- summary(glm_fit)  #p value
   r_spm <- result(glm_fit,
                   train_g,
                   test_g,
                   response_train,
                   response_test)
   r_spm
   # step
   glm_spm_step <- step(glm_fit)
   summary(glm_spm_step)
   r_spm_step <- result(glm_spm_step,
                        train_g,
                        test_g,
                        response_train,
                        response_test)
   r_spm_step
 }

 #- glm_lasso model
 {
   train_g_lasso <- train_glm_lasso
   test_g_lasso  <- test[names(train_g_lasso)]
   glm_fit_lasso <- glm(response_train ~ .,
                        family= binomial,
                        data=train_g_lasso,
                        control=list(maxit=1000))
   sum_glm <- summary(glm_fit_lasso)  #p value
   r_lasso <- result(glm_fit_lasso,
                     train_g_lasso,
                     test_g_lasso,
                     response_train,
                     response_test)
   r_lasso
   # step
   glm_lasso_step <- step(glm_fit_lasso)
   summary(glm_lasso_step)

   train_step <- train_g_lasso[names(glm_lasso_step$coefficients)[-1]]
   test_step  <- test_g_lasso[names(train_step)]
   r_lasso_step <- result(glm_lasso_step,
                          train_step,
                          test_step,
                          response_train,
                          response_test)
   r_lasso_step
 }

 #- glm_pca model
 {
   x = 12
   train_g_pca <- train_pca[1:x]
   test_g_pca  <- test_pca
   glm_fit <- glm(response_train ~ .,
                  family= binomial,
                  data=train_g_pca,
                  control=list(maxit=1000))
   sum_glm <- summary(glm_fit)  #p value
   r_pca <- result(glm_fit,
                   train_g_pca,
                   test_g_pca,
                response_train,
                response_test)
   r_pca
   # step
   glm_fit_step <- step(glm_fit)
   r_pca_step <- result(glm_fit_step,
                        train_g_pca,
                        test_g_pca,
                        response_train,
                        response_test)
   summary(glm_fit_step)
   r_pca_step
 }

 #------------------- Model diagnosis -------------------------
 # ROC
 {
   model <- glm_fit_lasso
   train_last <- train_glm_lasso
   test_last  <- test[names(train_last)]
   response_trn_pred <- predict(model,  data=train_last,type="response")
   response_tst_pred <- predict(model,newdata=test_last,type="response")
   par(mfrow=c(1,2))
   # train ROC
   plot.roc(response_train,
            response_trn_pred,
            main="The ROC of Train set",
            col="1",
            add=F,
            reuse.auc=T,
            axes=T,
            max.auc.polygon=T,
            legacy.axes=F,   # x:0 -> 1
            grid=c(0.2),
            #grid.col=c("green", "red"),
            # auc.polygon.col="skyblue",
            print.thres=T,
            print.auc=T,
            auc.polygon=T,
            identity.lty=1,
            identity.lwd=2,
            asp=1,
            mar=c(4, 4, 2, 2)+.1,
            mgp=c(2.5, 1, 0))
   # test ROC
   plot.roc(response_test,
            response_tst_pred,
            main="The ROC of Test set",
            col="1",
            add=F,
            reuse.auc=T,
            axes=T,
            max.auc.polygon=T,
            legacy.axes=F,   # x:0 -> 1
            grid=c(0.2),
            #grid.col=c("green", "red"),
            #auc.polygon.col="skyblue",
            print.thres=T,
            print.auc=T,
            auc.polygon=T,
            identity.lty=1,
            identity.lwd=2,
            asp=1,
            mar=c(4, 4, 2, 2)+.1,
            mgp=c(2.5, 1, 0))

   roc(response_train,response_trn_pred,ci=T)

 }
 # Radiomics Score
 {
   Data  <- data.frame(rbind(train_last,test_last))[ names(train_last) ]
   model_last <- glm_fit_lasso
   # View(Data)
   dim(Data)  # 126  10
   xn <- dim(Data)[1]      # row     126
   yn <- dim(Data)[2]      # column  10
   beta           <- coef(model_last,s="Coefficients")   # get beta = Coefficients
   betai_Matrix   <- as.matrix( beta[-1] )      # get betai
   beta0_Matrix   <- matrix(beta[1], xn, 1 )    # get beta0
   Rad_Matrix     <- as.matrix( Data )          # get X
   Radcore_Matrix <- Rad_Matrix %*% betai_Matrix + beta0_Matrix  # get Radscore
   radscore_all   <- as.numeric(Radcore_Matrix)
   Radscore_train  <- radscore_all[1:dim(train)[1]]    # train rads
   Radscore_test      <- radscore_all[c( (dim(train)[1]+1) :xn)] #test rads
   Radscore <- list(Radscore_train = Radscore_train,
                    Radscore_test  = Radscore_test)
   Radscore
 }

 # nomogram
 {
   dt <- data.frame(Radscore_train)
   ddist    <- datadist(dt)
   options(datadist="ddist")
   flrm     <- lrm(response_train ~ .,data=dt,x=T,y=T)

   # va <- validate(flrm,method="boot",B=150,dxy=T,pr=T)
   # cal  <- calibrate(flrm,method="boot",B=150)
   nomogplot <- nomogram(flrm,
                         fun=plogis,
                         fun.at=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1),
                         lp=F,
                         funlabel="Probability")
   par(mfrow=c(1,1))
   plot(nomogplot)
 }

 # Bar Chart
 {
   Data <- test_last
   Levels <-  response_test

   Rad.score <- Radscore_test
   Bar.Chart <- reorder(c(1:length(Levels)),Rad.score)
   p<-ggplot(Data, aes(x = Bar.Chart,
                       # xlab ="Bar Chart",
                       y = Rad.score,
                       fill = Levels,
                       #  color = c("1",'4'),
                       width=0.8) ) +
     geom_bar(stat="identity",color="white", alpha=0.5)
   p
 }

 # DCA
 {
   pred_test_spm   <- predict(glm_spm_step,test_g,type="response")
   pred_test_lasso <- predict(glm_lasso_step,test_g_lasso,type="response")
   pred_test_pca   <- predict(glm_fit_step,test_g_pca,type="response")
   dca_df <- data.frame(Group = as.numeric(response_test)-1 ,
                        SPM   = pred_test_spm,
                        LASSO = pred_test_lasso,
                        PCA   = pred_test_pca)

   dcaplot <- dca(data = dca_df,
                  outcome = "Group",
                  predictors = c("SPM",
                    "LASSO",
                    "PCA"),
                  xby=0.01, ymin=-0.42,
                  xstart=0.01, xstop=0.99)
   dcaplot
   net.benefit1 <- dcaplot$net.benefit[,4]
   net.benefit2 <- dcaplot$net.benefit[,5]
   net.benefit3 <- dcaplot$net.benefit[,6]

   kw <- kruskal.test(x = c(net.benefit1,net.benefit2,net.benefit3),
                      g = as.factor(c(rep(0,99),rep(1,99),rep(2,99))));kw

   s1 <- sum(net.benefit1)
   s2 <- sum(net.benefit2)
   s3 <- sum(net.benefit3)
   s1;s2;s3

 }
