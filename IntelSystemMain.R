chk_function <- function(filename,chk_dec,chk_sep)
{
  if(chk_dec)
    if (chk_sep)
      aa=read.table(filename,header=T,dec = ",",sep = "\t",row.names=1)
    else
      aa=read.table(filename,header=T,dec = ",")
  else
      if (chk_sep)
        aa=read.table(filename,header=T,dec = ".",sep = "\t",row.names=1)
      else
        aa=read.table(filename,header=T,dec = ".")
      aa
}

use_comma_for_decimal <- function()
{
  unname(Sys.localeconv()["decimal_point"] == ",")
}

use_tab_for_sep <- function()
{
  unname(Sys.localeconv()["thousands_sep"] == "\t")
}

makePlot <- function(ind,pca) {
  plot(PC2 ~ PC1, pca,cex=2, pch=16, col=as.numeric(dfr[,ncol(dfr)])+2)
  if(!missing(ind) && any(ind))
    points(PC2 ~ PC1, pca[ind,], cex=3, pch=1, col="red")
}

calcEnsemble<-function(data, M, algorithms, validation, boost, summary,
                       numfold,cbutton,pb,graph4,flagMult=FALSE)
{
  x<-as.matrix(data[,-ncol(data)])
  y <- factor(as.numeric(data[,ncol(data)]) - 1)
  out=createFolds(y, k = numfold, list = TRUE, returnTrain = T)
  Res <- list()
  Boost <- list()
  Ens_model<-list()
  start_time <- Sys.time()
  fold=1
  step<-(numfold*M)%/%100+1
  if(step==0)
    step<--(100%/%(numfold*M))
  
  bestAlg=rep(0,length(algorithms))
  names(bestAlg)=algorithms
  
  for(ff in out)
  {
    # while(gtkEventsPending()){
    #   gtkMainIteration()
      if(svalue(cbutton)=="Stopping..."){
        svalue(cbutton) <- "Stop simulation"
        cat("interrupted")
        #gtkMainQuit()
        return()
      }
    #}
    ens <- ensembleClassifier(x[ff,], y[ff], M=M, lambda1=10,boost=boost,
                              algorithms=algorithms,validation=validation,
                              progress=pb,fold=fold,step=step)
    
    # predict using the test data
    xtest <- x[-ff,]
    trueClass <- y[-ff]

    pred_model <- predictEns(ens, xtest, trueClass,graph=graph4,flagMult=flagMult)
    if(!is.null(boost))
    {
      pred_boost <- predictEns(ens, xtest, trueClass,FALSE)
      Boost=c(Boost,list(pred_boost))
    }

    Res=c(Res,list(pred_model))
    Ens_model=c(Ens_model,ens)
    vrem=table(ens$bestAlg)
    sapply(names(vrem),function(z) bestAlg[z]<<-bestAlg[z]+vrem[z])
    
    cat("The fold number=", fold, "\n")
    fold=fold+1
  }
  x <- Res[[1]]$indPerf
  lapply(seq_along(Res)[-1], function(i){
    x <<- x + Res[[i]]$indPerf
  })
  predAvgInd=x/length(Res)
  
  x <- Res[[1]]$ensemblePerf
  lapply(seq_along(Res)[-1], function(i){
    x <<- x + Res[[i]]$ensemblePerf
  })
  predAvgEns=x/length(Res)
  
  if(!is.null(boost))
  {
    x <- Boost[[1]]$ensemblePerf
    lapply(seq_along(Boost)[-1], function(i){
      x <<- x + Boost[[i]]$ensemblePerf
    })
    predAvgBoost=x/length(Boost)
  }
  if(flagMult)
    nfun<-1
  else
    nfun<-5
  end_time <- Sys.time()
  out=end_time - start_time
  cat(paste("Time elapsed",out))
  nalg=length(algorithms)
  t1=sapply(1:numfold, function(z) t(Res[[z]]$indPerf))
  dim(t1)=c(nfun,nalg,numfold)
  res1=do.call(rbind,lapply(1:nalg,function(z) apply(t1[,z,,drop=F],1,median)))
  colnames(res1)=colnames(Res[[1]]$indPerf)
  rownames(res1)=rownames(Res[[1]]$indPerf)
  
  t2=sapply(1:numfold, function(z) t(Res[[z]]$ensemblePerf))
  dim(t2)=c(nfun,numfold)
  res2=apply(t2,1,median)
  names(res2)=colnames(Res[[1]]$ensemblePerf)
  
  if(!is.null(boost))
  {
    t3=sapply(1:numfold, function(z) t(Boost[[z]]$ensemblePerf))
    dim(t3)=c(nfun,numfold)
    res3=apply(t3,1,median)
    names(res3)=colnames(Boost[[1]]$ensemblePerf)
  }
  ss=svalue(pb)
  if(ss<100)
    for(i in (ss+1):99) {svalue(pb) <- i}
  
  insert(summary,"Performance measures", font.attr=c(family="monospace"))
  insert(summary,"\nIndividual classifiers", font.attr=c(family="monospace"))
  insert(summary,paste(colnames(predAvgInd), collapse = "\t"),font.attr=c(family="monospace"))
  sapply(1:nalg,function(x) insert(summary,paste(round(predAvgInd[x,],4), collapse = "\t\t"),font.attr=c(family="monospace")))
  insert(summary,"\nEnsemble classifier", font.attr=c(family="monospace"))#, font.attr=c(color="red"))
  insert(summary,paste(colnames(predAvgEns), collapse = "\t"),font.attr=c(family="monospace"))
  insert(summary,paste(round(predAvgEns,4), collapse = "\t\t"),font.attr=c(family="monospace"))
  insert(summary,"\nNumber the classifiers used", font.attr=c(family="monospace"))
  insert(summary,paste(names(bestAlg), collapse = "\t\t"),font.attr=c(family="monospace"))
  insert(summary,paste(bestAlg, collapse = "\t\t"),font.attr=c(family="monospace"))
  ii=1
}

IntelSystemMain<-function()
{
  library(caret)
  library(gWidgets2)
  source("ensemble.R")
  
  boost=NULL
  numfold=3
  flag_data=""#DLBCL"
  flagMult=FALSE
  
  options(guiToolkit = 'RGtk2')
  
  window <- gwindow ("Intelligent data analysis", visible = F)
  #size(window)=c(300,600)
  notebook <- gnotebook (cont = window,tab.pos = 3,closebuttons = F,
                         font.attr=list(foreground="red"))
  group1 <- ggroup(cont = notebook, label = "Data set", horizontal=F)
  group2 <- ggroup(cont = notebook, label = "Feature Selection", horizontal=T)
  group3 <- ggroup(cont = notebook, label = "Clustering", horizontal=T)
  group4 <- ggroup(cont = notebook, label = "Classifier modeling", horizontal=F)
  window$set_rgtk2_font(window$children[[1]]$widget,list(background="yellow"))
#--------Dataset tab----------------------
  Tab1_Dataset=ggroup(container = group1,expand=T)
  Bigroup1=ggroup(container = Tab1_Dataset,horizontal=F,
                 font.attr=list(background="red"))

  grp_name <- ggroup(container = Bigroup1)
  lbl_data_frame_name <- glabel(
    "Variable to save data to: ",
    container = grp_name
  )
  txt_data_frame_name <- gedit("dfr", container = grp_name)

  grp_upload <- ggroup(container = Bigroup1,horizontal=F)
  gg<-ggroup(cont=grp_upload)
  
  chk_eurostyle <- gcheckbox(
    text      = "Use comma for decimal place",
    checked   = use_comma_for_decimal(),
    container = gg
  )
  chk_sep <- gcheckbox(
    text      = "Use tab for separation",
    checked   = use_tab_for_sep(),
    container = grp_upload
  )

  chk_data <- gcheckbox(
    text      = "Samples in columns",
    checked   = FALSE,
    container = grp_upload
  )

a <- gfilebrowse("Upload file",cont=gg,expand=T,
                   action  = "chk_function",
                   handler=function(h,...){
                     tryCatch(
                       {
                         data_frame_name <- make.names(svalue(txt_data_frame_name))
                         the_data <- do.call(h$action, list(svalue(a),svalue(chk_eurostyle),svalue(chk_sep)))
                         #the_data <- do.call(h$action, list(svalue(a)))
                         #the_data=read.table(svalue(a),dec = ",",sep = "\t")
                         if (svalue(chk_data))
                         {
                            yy=colnames(the_data)
                            tt=unlist(lapply(gregexpr(pattern = '\\.', yy), min))
                            names=sapply(1:length(yy),function(z){ifelse(tt[z] > 0, substr(yy[z],1,tt[z]-1), yy[z])})
                            vrem=t(the_data)
                            #colnames(vrem)=paste("V",1:ncol(vrem),sep="")
                            the_data=cbind(vrem,data.frame(class=names))
                            #the_data$class <- factor(the_data$class, levels=rev(levels(the_data$class)))
                            the_data$class <- factor(the_data$class)
                         }
                         dfNames <- names(the_data)
                         print(names(the_data)[1:10])
                         if(flag_data=="DLBCL")
                         {
                           lev=c("DLBCL","FL")
                           sapply(1:2, function(z) {ind=which(the_data$Result==z);the_data$Result[ind]<<-lev[z]})
                           the_data$Result=factor(the_data$Result,levels=rev(lev))
                         }
                         tbl[,] <- data.frame(the_data[,c(1:10,ncol(the_data))])#data.frame(variables=dfNames, stringsAsFactors=FALSE)#

                         assign(data_frame_name, the_data, envir = globalenv())
                         svalue(status_bar) <-
                           paste(nrow(the_data), "records saved to variable", data_frame_name)
                       },
                       error = function(e) svalue(status_bar) <- "Could not upload data"
                     )
                   },
                   filter = list(
                     "Tab delimited" = list(patterns = c("*.txt","*.dlm","*.tab")),
                     "All files" = list(patterns = c("*"))
                   )
  )

  blankDF = data.frame(variables=character(0), stringsAsFactors=FALSE)
  frmX <- gframe("Data set", container = Bigroup1)
  tbl <- gtable(blankDF,  container=frmX, filter.FUN="manual")
  size(tbl) <- c(500, 500)
  
  
  #gbutton("dismiss", container=Bigroup1, handler=function(h,...) dispose(window))
  
  graph1 <- ggraphics(container=Tab1_Dataset,expand=T)

  ID <- addHandlerChanged(graph1, handler=function(h,...) {
    x <- h$x; y <- h$y
    res<-prcomp(dfr[,-ncol(dfr)])
    pca=as.data.frame(res$x)
    print(summary(res)$importance[3,])
    print(as.numeric(dfr[,ncol(dfr)]))
    ind <- (pca$PC1 >= x[1]) & (pca$PC1 <= x[2]) &
      (pca$PC2 >= y[1]) & (pca$PC2 <= y[2])
    ## udpate graphic and data frame
    makePlot(ind,pca)
    if(any(ind))
      visible(tbl) <- ind
  })

  status_bar <- gstatusbar("", container = group1)
#-----------Dataset tab----------------------------------------  
  
#------------Classifier Tab---------------------------
  models <- c("svm", "rf", "lda", "qda", "lr", "nnet", "pls_lda", "pls_qda", "pls_rf", "pls_lr",
               "pca_lda", "pca_qda", "pca_rf", "pca_lr", "rpart", "adaboost", "plr", "pls")
  measures=c("accuracy", "sensitivity", "specificity")
  algorithms=NULL
  validation=NULL
  M=51
  
  Tab1_Classifier=ggroup(container = group4,expand=T)
  
  Bigroup4=ggroup(container = Tab1_Classifier,horizontal=F)
  size(Bigroup4)=c(300,1)
  
  frmM <- gframe("Models", container = Bigroup4)               
  # cmodels <- gcheckboxgroup(models, container=frmM, handler = function(h,...){
  #   algorithms <<- svalue(h$obj)})#,index=TRUE)}))
  cmodels <- gtable(models, container=frmM, multiple = T,expand=T)#,index=TRUE)}))
  size(cmodels)=c(1,200)
  frmV <- gframe("Validation", container = Bigroup4)
  cmeasures<-gcheckboxgroup(measures, container=frmV, handler = function(h,...){
    validation <<- svalue(h$obj)})#,index=TRUE)})
  
  bootstrap<-gslider(from=10,to=100,by=10, value=50, handler = function(h,...){
    M <<- svalue(h$obj)+1})
  frmS<-gframe("Bootstrap number", container=Bigroup4)
  add(frmS,bootstrap,expand=TRUE)
  
  b<-gbutton("Start", container=Bigroup4)
  cbutton <- gbutton("Stop simulation",cont=Bigroup4,handler=function(h,...) 
  {svalue(h$obj) <- "Stopping..."})
  
  pb <- gprogressbar(cont=group4,value=0)
  
  #status_bar <- gstatusbar("", container = group4)

  addHandlerClicked(b, handler=function(h,...) {
    tryCatch(
    {
      data<-get(svalue(txt_data_frame_name),envir=.GlobalEnv)
    },
      error = function(e) {svalue(status_bar) <- paste("Dataset ",svalue(txt_data_frame_name)," is empty");break}
    )
      algorithms=svalue(cmodels, drop=FALSE)
      if(is.null(algorithms)||is.null(validation))
      {
        svalue(status_bar) <-
          paste("Models or Validation are not defined")
      }
      else
      {
        svalue(textsummary)=""
        svalue(pb)=0
        svalue(cbutton) <- "Stop simulation"
        algorithms=algorithms[,1]
        calcEnsemble(data, M, algorithms, validation, boost,textsummary,
                     numfold,cbutton,pb,graph4,flagMult=flagMult)
        svalue(status_bar) <-"Ready"
      }
  })
  resgroup=ggroup(container=Tab1_Classifier,horizontal=F,expand=T)
  textsummary <- gtext("", container=resgroup,expand=T)
  frmgraph4 <- gframe("Plots", container = resgroup,expand=T)
  graph4 <- ggraphics("", container=frmgraph4,expand=T)
#------------Classifier Tab---------------------------
  
  
  svalue(notebook) <- 1
  visible(window) <- TRUE
}