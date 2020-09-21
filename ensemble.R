library(e1071)
library(randomForest)
library(RankAggreg)
library(MASS)
library(nnet)
library(class)
library(adabag)
library(rpart)
library(plsgenomics)
#library(penalized)
library(glmnet)
library(ROCR)

source("validation.R")
source("miscellaneous.R")

ensembleClassifier <- function(x, y, M=51, fit.individual=TRUE, varimp=FALSE, rfs=FALSE, nf=ceiling(sqrt(ncol(x))),
                               boost=NULL,
    algorithms=c("svm", "rf", "lda", "qda", "lr", "nnet", "pls_lda", "pls_qda", "pls_rf", "pls_lr",
                    "pca_lda", "pca_qda", "pca_rf", "pca_lr", "rpart", "adaboost", "plr", "pls"), 
    validation=c("accuracy", "sensitivity", "specificity"), progress=NULL,fold=NULL, step=NULL,
    ncomp=5, nunits=3, lambda1=5, kernel="radial",
    distance="Spearman", weighted=TRUE, verbose=TRUE, seed=NULL, ...){

    rownames(x) <- NULL # to suppress the warning message about duplicate rownames

    if(!is.null(seed))
        set.seed(seed)  
    
    # if(length(algorithms) < 2)
    #     stop("Ensemble classifier needs at least 2 classification algorithms")

    # error on p >> N for classical methods
    if(nrow(x) < ncol(x))
        if(length(intersect(algorithms, c("lda", "qda", "lr"))) > 0)
            stop("LDA, QDA and LR cannot be used if nrow(x) < ncol(x) (p > N)")

    n <- length(y)
    nalg <- length(algorithms)
    nvm <- length(validation)
    ly <- levels(y)

    if(length(ly) == 2 && all(ly != c("0","1")))
        stop("For binary classification, levels in y must be 0 and 1")

    fittedModels <- list() #to keep the fitted algorithms
    fittedBoost<-list()
	  varImportance <- NULL
	  if(length(ly)>2)
	    fam="multinomial"
	  else
	    fam="binomial"
#--------- for progress bar--------------
    for(k in 1:M){
      if(!is.null(progress))
      {
        ind<-(fold-1)*M+k
        if(step<0)
        {
          svalue(progress) <- svalue(progress)+abs(step)
        }
        else
        {
          if(ind%%step==0)
          {
            #ss<-svalue(progress)
            svalue(progress) <- ind%/%step
          }
        }
      }
#-----------------------------------      
        repeat{ # to make sure that at least two classes are present 
            s <- sample(n, replace=TRUE)
            if(length(table(y[s])) >= 2 & length(table(y[-s])) >= 2)
                break
        }
        
        if(rfs) # perform feature selection?
            fs <- sample(1:ncol(x), nf)  else
            fs <- 1:ncol(x)
            
        training <- x[s, fs]
        testing  <- x[-unique(s), fs]
        trainY <- y[s]  
        
        # construct PLS latent variables (store in plsX) if one of the PLS methods used
        # this method is taken from library CMA (bioconductor)
        if("pls_lda" %in% algorithms || "pls_qda" %in% algorithms || "pls_rf" %in% algorithms ||
                "pls_lr" %in% algorithms){
            if(ncomp >= ncol(x))
                stop("Decrease ncomp for the PLS models; must be smaller than ncol(x)")
            fpr <- pls.regression(training, transformy(trainY), ncomp=ncomp)
            plsX <- scale(training, scale=FALSE, center=fpr$meanX)%*%fpr$R
            plsTestX <- scale(testing, scale=FALSE, center=fpr$meanX)%*%fpr$R
        }
        
        # construct PCA latent variables (store in pcaX)
        if("pca_lda" %in% algorithms || "pca_qda" %in% algorithms || "pca_rf" %in% algorithms ||
                "pca_lr" %in% algorithms){
            if(ncomp >= ncol(x))
                stop("Decrease ncomp for the PLS models; must be smaller than ncol(x)")
            pcaX <- prcomp(training)$x[,1:ncomp]
            pcaTestX <- prcomp(testing)$x[,1:ncomp]
        }        
                

        # train all algorithms on the subset
        Res <- list()   
        for(j in 1:nalg){
            Res[[j]] <- switch(algorithms[j],
                "svm"      = svm(training, trainY, probability=TRUE, kernel=kernel, ...),
                "rf"       = randomForest(training, trainY, ...),
                "lda"      = lda(training, trainY, ...),
                "qda"      = qda(training, trainY, ...),
                "lr"       = glm(y~., family=binomial, data=data.frame(y=trainY, x=training), ...),
                "plr"      = {
                              out=createFolds(trainY, k = 10, list = F, returnTrain = F)
                              cv.lasso <- cv.glmnet(training, trainY, alpha = 1, family = fam, foldid=out)
                              # Fit the final model on the training data
                              glmnet(training, trainY, alpha = 1, family = fam,intercept = F,
                              lambda = cv.lasso$lambda.min)
                              },#penalized(trainY, penalized=training, trace=T, model="logistic", lambda1=lambda1, ...),
                "pls_lda"  = lda(plsX, trainY, ...),
                "pls_qda"  = qda(plsX, trainY, ...),
                "pls_rf"   = randomForest(plsX, trainY, ...),
                "pls_lr"   = glm(y~., family=binomial, data=data.frame(y=trainY, x=plsX), ...),
                "pca_lda"  = lda(pcaX, trainY, ...),
                "pca_qda"  = qda(pcaX, trainY, ...),
                "pca_rf"   = randomForest(pcaX, trainY, ...),
                "pca_lr"   = glm(y~., family=binomial, data=data.frame(y=trainY, x=pcaX), ...),              
                "nnet"     = nnet(training, class.ind(trainY), size=nunits, trace=FALSE, ...),
                "rpart"    = rpart(y~., data=data.frame(y=trainY, training), ...),
                "adaboost" = adaboost.M1(y~., data=data.frame(y=trainY, training), ...),
                "pls"      = pls.regression(training, transformy(trainY), ncomp=ncomp, ...)
            )
            attr(Res[[j]], "algorithm") <- algorithms[j]
            if("pls_lda" == algorithms[j] || "pls_qda" == algorithms[j] || "pls_rf" == algorithms[j] ||
                   "pls_lr" == algorithms[j]){
                attr(Res[[j]], "meanX") <- fpr$meanX
                attr(Res[[j]], "R") <- fpr$R
            }
            if("pca_lda" == algorithms[j] || "pca_qda" == algorithms[j] || "pca_rf" == algorithms[j] ||
                   "pca_lr" == algorithms[j])
                attr(Res[[j]], "ncomp") <- ncomp
            #cat(algorithms[j], "\n")
        }
        
        if(length(algorithms) >1)
        {
        # predict using fitted models on oob (class 1 probabilities for ROC curves)
        probabilities <- list()
        predicted <- list()
        for(j in 1:nalg){
            switch(algorithms[j],
                "svm"       = {predicted[[j]] <- predict(Res[[j]], testing, prob=TRUE)
                            probabilities[[j]] <- attr(predicted[[j]], "probabilities")[,
                                which(colnames(attr(predicted[[j]], "probabilities")) == "1")]},
                "rf"        = {predicted[[j]] <- predict(Res[[j]], testing, type="class")
                                temp <- predict(Res[[j]], testing, type="prob")
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "lda"       = {predicted[[j]] <- predict(Res[[j]], testing)$class
                                temp <- predict(Res[[j]], testing)$posterior
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "qda"       = {predicted[[j]] <- predict(Res[[j]], testing)$class
                                temp <- predict(Res[[j]], testing)$posterior
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "lr"        = {probabilities[[j]] <- predict(Res[[j]], data.frame(x=testing), type="response")
                                temp <- probabilities[[j]] > .5
                                temp[temp] <- 1
                                predicted[[j]] <- factor(temp)},
                "plr"       = {probabilities[[j]] <- predict(Res[[j]], testing,type="response")
                                temp <- probabilities[[j]] > .5
                                temp[temp] <- 1
                                predicted[[j]] <- factor(temp)},
                "pls_lda"   = {predicted[[j]] <- predict(Res[[j]], plsTestX)$class
                                temp <- predict(Res[[j]], plsTestX)$posterior
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "pls_qda"   = {predicted[[j]] <- predict(Res[[j]], plsTestX)$class
                                temp <- predict(Res[[j]], plsTestX)$posterior
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "pls_lr"    = {probabilities[[j]] <- predict(Res[[j]], data.frame(x=plsTestX), type="response")
                                temp <- probabilities[[j]] > .5
                                temp[temp] <- 1
                                predicted[[j]] <- factor(temp)},
                "pls_rf"    = {predicted[[j]] <- predict(Res[[j]], plsTestX, type="class")
                                temp <- predict(Res[[j]], plsTestX, type="prob")
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "pca_lda"   = {predicted[[j]] <- predict(Res[[j]], pcaTestX)$class
                                temp <- predict(Res[[j]], pcaTestX)$posterior
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "pca_qda"   = {predicted[[j]] <- predict(Res[[j]], pcaTestX)$class
                                temp <- predict(Res[[j]], pcaTestX)$posterior
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "pca_lr"    = {probabilities[[j]] <- predict(Res[[j]], data.frame(x=pcaTestX), type="response")
                                temp <- probabilities[[j]] > .5
                                temp[temp] <- 1
                                predicted[[j]] <- factor(temp)},
                "pca_rf"    = {predicted[[j]] <- predict(Res[[j]], pcaTestX, type="class")
                                temp <- predict(Res[[j]], pcaTestX, type="prob")
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},                         
                "nnet"      = {predicted[[j]] <- predictNNET(Res[[j]], testing) # custom fn defined above
                                temp <- predict(Res[[j]], testing)
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "rpart"     = {predicted[[j]] <- predict(Res[[j]], newdata=data.frame(testing), type="class")
                                temp <- predict(Res[[j]], newdata=data.frame(testing))
                                probabilities[[j]] <- temp[,which(colnames(temp) == "1")]},
                "adaboost"  = {predicted[[j]] <- predict.boost(Res[[j]], newdata=data.frame(testing))$class
                                probabilities[[j]] <- predict.boost(Res[[j]], newdata=data.frame(testing))$probs},
                "pls"       = {predicted[[j]] <- predict.pls(Res[[j]], newdata=testing)$class
                				probabilities[[j]] <- predict.pls(Res[[j]], newdata=testing)$probs}
            )
          #cat(algorithms[j], "\n")
          }

        # compute validation measures           
        scores <- matrix(0, nalg, nvm)
        rownames(scores) <- algorithms
        colnames(scores) <- validation

        truth <- y[-unique(s)]
        for(i in 1:nalg)
            for(j in 1:nvm)
                scores[i,j] <- switch(validation[j],
                    "accuracy"    = accuracy(truth, factor(predicted[[i]], levels=ly)),
                    "sensitivity" = sensitivity(truth, factor(predicted[[i]], levels=ly)),
                    "specificity" = specificity(truth, factor(predicted[[i]], levels=ly)),
                    "auc"         = AUC(truth, probabilities[[i]]))     
        

        # perform rank aggregation
        convScores <- convertScores(scores)
        }
        
        if(!is.null(boost))
        {
          fittedBoost[[k]] <- Res[[boost]]
        }
        if(nvm > 1 && nalg>1 && nalg <= 6)
          if(weighted)
            fittedModels[[k]] <- Res[[which(algorithms == BruteAggreg(convScores$ranks,
                                                                      nalg, convScores$weights, distance=distance)$top.list[1])]]
        else
          fittedModels[[k]] <- Res[[which(algorithms == BruteAggreg(convScores$ranks, nalg,
                                                                    distance=distance)$top.list[1])]]
        else if(nvm > 1 && nalg > 6)
          if(weighted)
            fittedModels[[k]] <- Res[[which(algorithms == RankAggreg(convScores$ranks,
                                                                     nalg, convScores$weights, distance=distance, verbose=FALSE)$top.list[1])]]
          else
            fittedModels[[k]] <- Res[[which(algorithms == RankAggreg(convScores$ranks, nalg,
                                                                   distance=distance, verbose=FALSE)$top.list[1])]]        
        else
          if(nalg<2)
            fittedModels[[k]] <- Res[[1]]
          else
            fittedModels[[k]] <- Res[[which.max(scores[,1])]]

		###############################################################################################
		# variable importance as in Random Forest
		if(varimp){
			predicted <- matrix(0, nrow(testing), ncol(testing)+1) # one extra for untouched data
			kalg <- attr(fittedModels[[k]], "algorithm")
			nts <- nrow(testing) #number of testing samples
			
			if(!kalg %in% c("pls_lda", "pls_qda", "pls_rf", "pls_lr","pca_lda", "pca_qda", "pca_rf", "pca_lr")){
				testingX <- testing
				testingY <- as.numeric(truth)
				
				for(i in 1:(ncol(testing)+1)){
					if(i != ncol(testing)+1){
						rsample <- sample(nts, nts)
						testing <- testingX
						testing[,i] <- testingX[rsample, i]
					}
					else
						testing <- testingX
					
					colnames(testing) <- colnames(testingX)
											
					switch(kalg,
							"svm"       = predicted[,i] <- predict(fittedModels[[k]], testing, prob=TRUE),
							"rf"        = predicted[,i] <- predict(fittedModels[[k]], testing, type="class"),
							"lda"       = predicted[,i] <- predict(fittedModels[[k]], testing)$class,
							"qda"       = predicted[,i] <- predict(fittedModels[[k]], testing)$class,
							"lr"        = {temp <- predict(fittedModels[[k]], data.frame(x=testing), type="response")
											temp <- temp > .5
											temp[temp] <- 1
											predicted[,i] <- factor(temp)},
							"plr"       = {temp <- predict(fittedModels[[k]], testing,type="response")
											temp <- temp > .5
											temp[temp] <- 1
											predicted[,i] <- factor(temp)},                     
							"nnet"      = predicted[,i] <- predictNNET(fittedModels[[k]], testing), # custom fn defined above
							"rpart"     = predicted[,i] <- predict(fittedModels[[k]], newdata=data.frame(testing), type="class"),
							"adaboost"  = predicted[,i] <- predict.boost(fittedModels[[k]], newdata=data.frame(testing))$class,
			                "pls"       = predicted[,i] <- predict.pls(fittedModels[[k]], newdata=testing)$class
					)
				}
				untouchedAccuracy <- sum(predicted[,ncol(predicted)] == testingY)/length(testingY)
				tempImp <- rep(0, ncol(testing))
				for(i in 1:ncol(testing))
					tempImp[i] <- untouchedAccuracy - sum(predicted[,i] == testingY)/length(testingY)
				varImportance <- rbind(varImportance, tempImp)
			}
		}
		##############################################################################################################
        #fittedModels[[k]]$selectedFeatures <- fs
        
        # some output
        if(verbose)
            cat("Iter ", k, "\n")
    } # loop 1:M

    # how many times each algorithms was the best?
    bestAlg <- unlist(sapply(fittedModels, FUN = function(x) attr(x, "algorithm")))
	rawImportance <- varImportance
	if(!is.null(varImportance)){
		varImportance <- matrix(c(apply(varImportance, 2, mean), apply(varImportance, 2, mean)/
				(apply(varImportance, 2, sd)/sqrt(nrow(varImportance)))), ncol=2)
		rownames(varImportance) <- colnames(x)
		colnames(varImportance) <- c("MeanDecreaseAcc", "StdMeanDecreaseAcc")
	}
    Res <- list()  
    
    ###########################################################################
    # train all classifiers individually on all training data
    ###########################################################################
    if(fit.individual){
        # construct PLS latent variables (store in plsX) if one of the PLS methods used
        # this method is taken from library CMA (bioconductor)
        if("pls_lda" %in% algorithms || "pls_qda" %in% algorithms || "pls_rf" %in% algorithms ||
                "pls_lr" %in% algorithms){
            if(ncomp >= ncol(x))
                stop("Decrease ncomp for the PLS models; must be smaller than ncol(x)")
            fpr <- pls.regression(x, transformy(y), ncomp=ncomp)
            plsX <- scale(x, scale=FALSE, center=fpr$meanX)%*%fpr$R
        }
        
        # construct PCA latent variables (store in pcaX)
        if("pca_lda" %in% algorithms || "pca_qda" %in% algorithms || "pca_rf" %in% algorithms ||
                "pca_lr" %in% algorithms){
            if(ncomp >= ncol(x))
                stop("Decrease ncomp for the PLS models; must be smaller than ncol(x)")
            pcaX <- prcomp(x)$x[,1:ncomp]
        }            
         
        training <- x #train using all data
        trainY <- y

        for(j in 1:nalg){
            Res[[j]] <- switch(algorithms[j],
                "svm"      = svm(training, trainY, probability=TRUE, kernel=kernel, ...),
                "rf"       = randomForest(training, trainY, ...),
                "lda"      = lda(training, trainY, ...),
                "qda"      = qda(training, trainY, ...),
                "lr"       = glm(y~., family=binomial, data=data.frame(y=trainY, x=training), ...),
                "plr"      = {
                            out=createFolds(trainY, k = 10, list = F, returnTrain = F)
                            cv.lasso <- cv.glmnet(training, trainY, alpha = 1, family = fam, foldid=out)
                            # Fit the final model on the training data
                            glmnet(training, trainY, alpha = 1, family = fam,intercept = F,
                            lambda = cv.lasso$lambda.min)},#penalized(trainY, penalized=training, trace=FALSE, model="logistic", lambda1=lambda1, ...),
                "pls_lda"  = lda(plsX, trainY, ...),
                "pls_qda"  = qda(plsX, trainY, ...),
                "pls_rf"   = randomForest(plsX, trainY, ...),
                "pls_lr"   = glm(y~., family=binomial, data=data.frame(y=trainY, x=plsX), ...),
                "pca_lda"  = lda(pcaX, trainY, ...),
                "pca_qda"  = qda(pcaX, trainY, ...),
                "pca_rf"   = randomForest(pcaX, trainY, ...),
                "pca_lr"   = glm(y~., family=binomial, data=data.frame(y=trainY, x=pcaX), ...),              
                "nnet"     = nnet(training, class.ind(trainY), size=nunits, trace=FALSE, ...),
                "rpart"    = rpart(y~., data=data.frame(y=trainY, training), ...),
                "adaboost" = adaboost.M1(y~., data=data.frame(y=trainY, training), ...),
                "pls"      = pls.regression(training, transformy(trainY), ncomp=ncomp, ...)                
            )
            attr(Res[[j]], "algorithm") <- algorithms[j]
            if("pls_lda" == algorithms[j] || "pls_qda" == algorithms[j] || "pls_rf" == algorithms[j] ||
                    "pls_lr" == algorithms[j]){
                attr(Res[[j]], "meanX") <- fpr$meanX
                attr(Res[[j]], "R") <- fpr$R
            }
            if("pca_lda" == algorithms[j] || "pca_qda" == algorithms[j] || "pca_rf" == algorithms[j] ||
                   "pca_lr" == algorithms[j])
                attr(Res[[j]], "ncomp") <- ncomp            
        }
    }       
    res <- list(models = fittedModels, indModels = Res, boost=fittedBoost,rawImportance=rawImportance, M = M, 
		bestAlg = bestAlg, levels=ly, importance=varImportance)#,convScores=convScores)
    class(res) <- "ensemble"
    res
}
            
predictEns <- function(EnsObject, newdata, y=NULL, flag=TRUE,plot=TRUE,graph=NULL,flagMult=FALSE){
    M <- EnsObject$M
    n <- nrow(newdata)
    predicted <- matrix(0, n, M)
    if (flag)
    {
      test_models=EnsObject$models
    }
    else
    {
      test_models=EnsObject$boost
    }
    for(i in 1:M){
        # constract components for PLS and PCA
        testing <- newdata #[,test_models[[i]]$selectedFeatures]
        if(attr(test_models[[i]], "algorithm") %in% c("pls_lda", "pls_qda", "pls_rf", "pls_lr")){
            R <- attr(test_models[[i]], "R")
            meanX <- attr(test_models[[i]], "meanX")
            plsTestX <- scale(testing, scale=FALSE, center=meanX)%*%R
        }
        if(attr(test_models[[i]], "algorithm") %in% c("pca_lda", "pca_qda", "pca_rf", "pca_lr")){
            pcaTestX <- prcomp(testing)$x[,1:as.numeric(attr(test_models[[i]], "ncomp"))]
        }        
        
        switch(attr(test_models[[i]], "algorithm"),
                "svm"       = predicted[,i] <- predict(test_models[[i]], testing, prob=TRUE),
                "rf"        = predicted[,i] <- predict(test_models[[i]], testing, type="class"),
                "lda"       = predicted[,i] <- predict(test_models[[i]], testing)$class,
                "qda"       = predicted[,i] <- predict(test_models[[i]], testing)$class,
                "lr"        = {temp <- predict(test_models[[i]], data.frame(x=testing), type="response")
                                temp <- temp > .5
                                temp[temp] <- 1
                                predicted[,i] <- factor(temp)},
                "plr"       = {temp <- predict(test_models[[i]], testing,type="response")
                                temp <- temp > .5
                                temp[temp] <- 1
                                predicted[,i] <- factor(temp)},
                "pls_lda"   = predicted[,i] <- predict(test_models[[i]], plsTestX)$class,
                "pls_qda"   = predicted[,i] <- predict(test_models[[i]], plsTestX)$class,
                "pls_lr"    = {temp <- predict(test_models[[i]], data.frame(x=plsTestX), type="response")
                                temp <- temp > .5
                                temp[temp] <- 1
                                predicted[,i] <- factor(temp)},
                "pls_rf"    = predicted[,i] <- predict(test_models[[i]], plsTestX, type="class"),
                "pca_lda"   = predicted[,i] <- predict(test_models[[i]], pcaTestX)$class,
                "pca_qda"   = predicted[,i] <- predict(test_models[[i]], pcaTestX)$class,
                "pca_lr"    = {temp <- predict(test_models[[i]], data.frame(x=pcaTestX), type="response")
                                temp <- temp > .5
                                temp[temp] <- 1
                                predicted[,i] <- factor(temp)},
                "pca_rf"    = predicted[,i] <- predict(test_models[[i]], pcaTestX, type="class"),                         
                "nnet"      = predicted[,i] <- predictNNET(test_models[[i]], testing), # custom fn defined above
                "rpart"     = predicted[,i] <- predict(test_models[[i]], newdata=data.frame(testing), type="class"),
                "adaboost"  = predicted[,i] <- predict.boost(test_models[[i]], newdata=data.frame(testing))$class,
                "pls"       = predicted[,i] <- predict.pls(test_models[[i]], newdata=testing)$class
        )
    }

    predicted <- predicted - 1
    
    # class probabilities by majority
    if(flagMult)
    {
      newclass <- factor(apply(predicted, 1, function(x) as.numeric(names(which.max(table(x))))), levels=EnsObject$levels)
      probabilities <- apply(predicted, 1, function(x) max(table(x))/M)
      res <- list()
      if(!is.null(y)){ # compute validation measures
        valM <- c("accuracy")
        acc <- accuracy(y, newclass)
        ensemblePerformance <- matrix(c(acc),1,1)
        colnames(ensemblePerformance) <- valM
        rownames(ensemblePerformance) <- "ensemble"
      }
    }
    else
    {
      newclass <- factor(apply(predicted, 1, function(x) ifelse(sum(x) > floor(M/2), 1, 0)), levels=EnsObject$levels)
      probabilities <- apply(predicted, 1, function(x) sum(x)/M)
      res <- list()
      if(!is.null(y)){ # compute validation measures
          valM <- c("accuracy", "sensitivity", "specificity", "auc","perf")
          acc <- accuracy(y, newclass)
          sens <- sensitivity(y, newclass)
          spec <- specificity(y, newclass)
          if(!is.null(graph))
            visible(graph)<-TRUE
          auc <- AUC(y, probabilities, plot=plot)
          pred <- prediction( probabilities, y)
          perf <- performance(pred,"auc")
          ensemblePerformance <- matrix(c(acc, sens, spec, auc,perf@y.values[[1]]),1,5)
          colnames(ensemblePerformance) <- valM
          rownames(ensemblePerformance) <- "ensemble"
      }
    }

    
    #############################################################################
    # predict using individual models
    #############################################################################
    if(length(EnsObject$indModels) > 0){  
       indPred <- matrix(0, nrow(newdata), length(EnsObject$indModels))
       indProb <- matrix(0, nrow(newdata), length(EnsObject$indModels))
       testing <- newdata
       
       for(i in 1:length(EnsObject$indModels)){
            # constract components for PLS and PCA
            if(attr(EnsObject$indModels[[i]], "algorithm") %in% c("pls_lda", "pls_qda", "pls_rf", "pls_lr")){
                R <- attr(EnsObject$indModels[[i]], "R")
                meanX <- attr(EnsObject$indModels[[i]], "meanX")
                plsTestX <- scale(testing, scale=FALSE, center=meanX)%*%R
            }
            if(attr(EnsObject$indModels[[i]], "algorithm") %in% c("pca_lda", "pca_qda", "pca_rf", "pca_lr")){
                pcaTestX <- prcomp(testing)$x[,1:as.numeric(attr(EnsObject$indModels[[i]], "ncomp"))]
            }      
            
            switch(attr(EnsObject$indModels[[i]], "algorithm"),            
                "svm"       = {temp <- predict(EnsObject$indModels[[i]], testing, prob=TRUE)
                                indPred[,i] <- temp
                                indProb[,i] <- attr(temp, "probabilities")[,which(colnames(attr(temp, "probabilities")) == "1")]},
                "rf"        = {indPred[,i] <- predict(EnsObject$indModels[[i]], testing, type="class")
                                temp <- predict(EnsObject$indModels[[i]], testing, type="prob")
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "lda"       = {indPred[,i] <- predict(EnsObject$indModels[[i]], testing)$class
                                temp <- predict(EnsObject$indModels[[i]], testing)$posterior
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "qda"       = {indPred[,i] <- predict(EnsObject$indModels[[i]], testing)$class
                                temp <- predict(EnsObject$indModels[[i]], testing)$posterior
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "lr"        = {indProb[,i] <- predict(EnsObject$indModels[[i]], data.frame(x=testing), type="response")
                                temp <- indProb[,i] > .5
                                temp[temp] <- 1
                                indPred[,i] <- factor(temp)},
                "plr"       = {indProb[,i] <- predict(EnsObject$indModels[[i]], testing,type="response")
                                temp <- indProb[,i] > .5
                                temp[temp] <- 1
                                indPred[,i] <- factor(temp)},
                "pls_lda"   = {indPred[,i] <- predict(EnsObject$indModels[[i]], plsTestX)$class
                                temp <- predict(EnsObject$indModels[[i]], plsTestX)$posterior
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "pls_qda"   = {indPred[,i] <- predict(EnsObject$indModels[[i]], plsTestX)$class
                                temp <- predict(EnsObject$indModels[[i]], plsTestX)$posterior
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "pls_lr"    = {indProb[,i] <- predict(EnsObject$indModels[[i]], data.frame(x=plsTestX), type="response")
                                temp <- indProb[,i] > .5
                                temp[temp] <- 1
                                indPred[,i] <- factor(temp)},
                "pls_rf"    = {indPred[,i] <- predict(EnsObject$indModels[[i]], plsTestX, type="class")
                                temp <- predict(EnsObject$indModels[[i]], plsTestX, type="prob")
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "pca_lda"   = {indPred[,i] <- predict(EnsObject$indModels[[i]], pcaTestX)$class
                                temp <- predict(EnsObject$indModels[[i]], pcaTestX)$posterior
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "pca_qda"   = {indPred[,i] <- predict(EnsObject$indModels[[i]], pcaTestX)$class
                                temp <- predict(EnsObject$indModels[[i]], pcaTestX)$posterior
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "pca_lr"    = {indProb[,i] <- predict(EnsObject$indModels[[i]], data.frame(x=pcaTestX), type="response")
                                temp <- indProb[,i] > .5
                                temp[temp] <- 1
                                indPred[,i] <- factor(temp)},
                "pca_rf"    = {indPred[,i] <- predict(EnsObject$indModels[[i]], pcaTestX, type="class")
                                temp <- predict(EnsObject$indModels[[i]], pcaTestX, type="prob")
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},                         
                "nnet"      = {indPred[,i] <- predictNNET(EnsObject$indModels[[i]], testing) # custom fn defined above
                                temp <- predict(EnsObject$indModels[[i]], testing)
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "rpart"     = {indPred[,i] <- predict(EnsObject$indModels[[i]], newdata=data.frame(testing), type="class")
                                temp <- predict(EnsObject$indModels[[i]], newdata=data.frame(testing))
                                indProb[,i] <- temp[,which(colnames(temp) == "1")]},
                "adaboost"  = {indPred[,i] <- predict.boost(EnsObject$indModels[[i]], newdata=data.frame(testing))$class
                                indProb[,i] <- predict.boost(EnsObject$indModels[[i]], newdata=data.frame(testing))$probs},
                "pls"       = {indPred[,i] <- predict.pls(EnsObject$indModels[[i]], newdata=testing)$class
                				indProb[,i] <- predict.pls(EnsObject$indModels[[i]], newdata=testing)$probs}                                        
               )
       }
        
        indPred <- indPred - 1
        
        if(!is.null(y)){
            indPerformance <- matrix(0, length(EnsObject$indModels), length(valM))
            rownames(indPerformance) <- unlist(sapply(EnsObject$indModels, FUN = function(x) attr(x, "algorithm")))
            colnames(indPerformance) <- valM

            if(flagMult)
              numfun=1
            else
              numfun=5
            truth <- y
            for(i in 1:length(EnsObject$indModels))
                for(j in 1:numfun)
                   indPerformance [i,j] <- switch(valM[j],
                    "accuracy" = accuracy(truth, factor(indPred[,i], levels=EnsObject$levels)),
                    "sensitivity" = sensitivity(truth, factor(indPred[,i], levels=EnsObject$levels)),
                    "specificity" = specificity(truth, factor(indPred[,i], levels=EnsObject$levels)),
                    "auc" = AUC(truth, indProb[,i]),
                    "perf"=performance(prediction(indProb[,i], truth),"auc")@y.values[[1]])
        }
    }
    
    if(is.null(y))
        res <- list(yhat=newclass, prob=probabilities, pred=predicted)
    else
        if(length(EnsObject$indModels) > 0)
            res <- list(yhat=newclass, prob=probabilities, pred=predicted, ensemblePerf=ensemblePerformance,
                    indPerf=indPerformance)
        else
            res <- list(yhat=newclass, prob=probabilities, pred=predicted, ensemblePerf=ensemblePerformance)
    class(res) <- "predictEnsemble"
    res
}
