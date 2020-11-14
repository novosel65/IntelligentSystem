proc_qda<-function(X,Y)
{
  var=unique(Y)
  fe<-sapply(var, function(z){apply(X[Y==z,],2,function(x) {out<-unique(x);
        res<-length(out)<=1;res})})
  res<-apply(fe,1,any)
  matr<-rep(TRUE,ncol(X))
  
  tmp <- sapply(var,function(z) {out=cor(X[Y==z,!res]);out[upper.tri(out)] <- 0;
              diag(out) <- 0;apply(out,2,function(x) any(abs(x) > 0.7))})
  res1<-apply(tmp,1,any)
  matr[!res]=res1
  data.new <- X[,!matr,drop=FALSE]
  cat(matr,"\n")
  #
  tryCatch(
    {
      out<-qda(data.new, Y)
    },
    error = function(e) {cat("Rank deficiency:\n",matr,"\n",table(Y))}
  )
  return(list(out=out,res=matr))
}

transformy_ver<-function(y)
{
  y<-as.numeric(y)
  K<-max(y)
  if (K>2)
  {
    Y<-matrix(0,length(y),K)
    for (k in 1:K)
    {
      Y[,k]<-as.numeric(y==k)
      Y[,k]<-Y[,k]-mean(Y[,k])
    }
  }
  else
  {
    Y<-matrix(y-mean(y),length(y),1)
  }
  
  Y
}
convertScores <- function(scores){
    scores <- t(scores)
    ranks <- matrix(0, nrow(scores), ncol(scores))
    weights <- ranks
    for(i in 1:nrow(scores)){
        ms <- sort(scores[i,], decr=TRUE, ind=TRUE)
        ranks[i,] <- colnames(scores)[ms$ix]
        weights[i,] <- ms$x
    }
    list(ranks = ranks, weights = weights)
}

predictNNET <- function(obj, newdata){
    # original didn't return class easily
    probs <- predict(obj, newdata)          
    probs <- probs[,which(colnames(probs) == "1")]
    probs <- probs > .5
    probs[probs] <- 1
    as.factor(probs)
}

predict.boost <- function (object, newdata, ...){ 
    # the original function form adabag package (predict.boosting) didn't return probabilities  
    vardep <- as.factor(object[[5]])
    mfinal <- length(object[[2]])
    n <- length(newdata[, 1])
    nclases <- nlevels(vardep)
    pesos <- rep(1/n, n)
    newdata <- data.frame(newdata, pesos)
    pond <- object[[3]]
    pred <- data.frame(rep(0, n))
    for (m in 1:mfinal) {
        if (m == 1) {
            pred <- predict(object[[2]][[m]], newdata, type = "class")
        }
        else {
            pred <- data.frame(pred, predict(object[[2]][[m]], 
                newdata, type = "class"))
        }
    }
    classfinal <- array(0, c(n, nlevels(vardep)))
    for (i in 1:nlevels(vardep)) {
        classfinal[, i] <- matrix(as.numeric(pred == levels(vardep)[i]), 
            nrow = n) %*% pond
    }
    predclass <- rep("O", n)
    for (i in 1:n) {
        predclass[i] <- as.character(levels(vardep)[(order(classfinal[i, 
            ], decreasing = TRUE)[1])])
    }
    probs <- classfinal[,which(levels(vardep) == "1")]/sum(classfinal[1,])
    output <- list(class = as.factor(predclass), probs = probs)
    output
}

predict.pls <- function(mod, newdata){
	# get the yhat from the pls regression in mod
	Yhat <- scale(newdata, scale = FALSE, center = mod$meanX) %*% mod$B
	probs <- as.vector((Yhat - min(Yhat))/(max(Yhat)-min(Yhat)))	
	class <- rep(0, length(Yhat))
	class[probs > .5] <- 1
	class <- as.factor(class)
	
	list(probs=probs, class=class)
}
	






