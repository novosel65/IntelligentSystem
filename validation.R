accuracy <- function(truth, predicted)
    if(length(truth) > 0)
        sum(truth==predicted)/length(truth) else    
        return(0)

sensitivity <- function(truth, predicted)
    # 1 means positive (present)
    if(sum(truth==1) > 0)
        sum(predicted[truth==1]==1)/sum(truth==1)   else
        return(0)

specificity <- function(truth, predicted)
    if(sum(truth==0) > 0)
        sum(predicted[truth==0]==0)/sum(truth==0)   else
        return(0)


AUC <- function(truth, probs, plot=FALSE){
    # probs - probability of class 1
    q <- seq(0, 1, .01)
    sens <- rep(0, length(q))
    spec <- rep(0, length(q))
    ly <- levels(truth)

    for(i in 1:length(q)){
        pred <- probs >= q[i]
        pred[pred] <- 1
        pred <- factor(pred, levels=ly)
        sens[i] <- sensitivity(truth, pred)
        spec[i] <- specificity(truth, pred)
    }

    # make sure it starts and ends at 0, 1
    sens <- c(1, sens, 0)
    spec <- c(0, spec, 1)
    
    trap.rule <- function(x,y) sum(diff(x)*(y[-1]+y[-length(y)]))/2 
    auc <- trap.rule(rev(1-spec), rev(sens))

    if(plot){
        plot(1-spec, sens, type="l", xlab="1-Specificity",
                ylab="Sensitivity", main="ROC Curve")
        legend("bottomright", legend=paste("AUC = ", round(auc, 3)), bty="n")   
    }
    auc
}
