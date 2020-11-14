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

calcFeature<-function(data=c(), nom.attr=c(),type=c(), disc="MDL", thresh=0.3,
                       incons=0,status=NULL,
                       flagMult=FALSE,flagTwo=FALSE,comboClass="")
{
  if(flagMult&&flagTwo)
  {
    y<-data[,ncol(data)]==svalue(comboClass)
    y[y]<-1
    y<-factor(y)
    data[,ncol(data)]=y
  }

  attrs.no=numeric()

  if(length(nom.attr)>0)
  {
      attrs.no<-numeric()
      for(i in 1:length(nom.attr)){
        attrs<-grep(nom.attr[i],colnames(data),fixed=T)[1]
        if(!is.na(attrs))
        {
          data[,attrs]<-as.factor(data[,attrs])
          attrs.no=c(attrs.no,attrs)
        }
      }
  }
    
  if(type=="FastFilter")
  {
      info.val<-select.fast.filter(data,disc,thresh,attrs.no)
  }
  if (type=="MRMR"){
    start_time<-Sys.time()
    info.val=select.mrmr.filter(data,disc,attrs.nominal,100)# more feature is time expensive
    end_time<-Sys.time()
    cat("Time=",end_time-start_time)
  }
  if(type=="auc")
  {
      
      index.auc=setdiff(1:(ncol(data)-1),attrs.no)
      auc.val<-compute.aucs(data[,c(index.auc,ncol(data)),drop=FALSE])
      
      aucs.all=rep(-1,ncol(data)-1)
      aucs.all[index.auc]=auc.val[,2]
      val <- sort(aucs.all,decreasing=T,index.return=TRUE)
      
      info.val <- data.frame(names(data)[val$ix[1:(ncol(data)-1)]],val$ix,val$x)
      names(info.val) <- c("Biomarker","Index","AUC")
  }
  if(type=="HUM")
  {
      indexF=1:(ncol(data)-1)
      indexClass=ncol(data)
      indexLabel=levels(data[,indexClass])
      
      index=setdiff(indexF,attrs.no)
      out=CalculateHUM_seq(data,indexF[index],indexClass,indexLabel)
      out.all=rep(-1,ncol(data)-1)
      out.all[index]=out$HUM
      val <- sort(out.all,decreasing=T,index.return=TRUE)
      info.val <- data.frame(names(data)[val$ix[1:(ncol(data)-1)]],val$ix,val$x)
      names(info.val) <- c("Biomarker","Index","HUM")
  }
  if(type=="CFS")
  {
      start_time<-Sys.time()
      info.val<-select.cfs(data)
      end_time<-Sys.time()
      cat("Time=",end_time-start_time)
  }
  if(type=="Relief")
  {
      start_time<-Sys.time()
      info.val<-select.relief(data)
      end_time<-Sys.time()
      cat("Time=",end_time-start_time)
  }
  if(type=="Forward search")
  {
      nn<-sapply(colnames(data),function(z) gsub('([-/])','_',z))
      nnold<-colnames(data)
      colnames(data)=nn
      start_time<-Sys.time()
      subset <- select.forward.wrapper(data)
      end_time<-Sys.time()
      cat("Time=",end_time-start_time)
      subset<-sapply(subset, function(z) which(names(data)==z))
      info.val <- data.frame(names(data)[subset],subset)
      colnames(data)=nnold
      names(info.val) <- c("Biomarker","Index")
  }
  if(type=="Chi2-algorithm")
  {
      start_time=Sys.time()
      out<-chi2.algorithm(data,attrs.no,incons)
      end_time <- Sys.time()
      cat("Time=",end_time-start_time)
      subset<-sapply(out$subset, function(z) which(names(data)==z))
      info.val <- data.frame(names(data)[subset],subset)
      names(info.val) <- c("Biomarker","Index")
  }
  if(type=="CorrSF")
  {
      start_time<-Sys.time()
      subset <- select.forward.Corr(data,disc,attrs.no)
      end_time <- Sys.time()
      cat("Time=",end_time-start_time)
      subset<-sapply(subset, function(z) which(names(data)==z))
      subset<-unlist(subset)
      if(length(subset)==0)
        info.val <- data.frame("","")
      else
        info.val <- data.frame(names(data)[subset],subset)
      names(info.val) <- c("Biomarker","Index")
  }
  if(type=="Information gain")
  {
      info.val<-select.inf.gain(data,disc,attrs.no)
  }
  if(type=="Symmetrical uncertainty")
  {
      info.val<-select.inf.symm(data,disc,attrs.no)
  }
  if(type=="Chi-square")
  {
      info.val<-select.inf.chi2(data,disc,attrs.no)
  }
  return(info.val)
}

calcEnsemble<-function(data, setFeature, numFeature,rfs,M, algorithms, validation, boost, summary,
                       numfold,cbutton,pb,status_bar,graphRoc,graphCross,flagMult=FALSE,flagTwo=FALSE,comboClass="")
{
  if(rfs==3)
  {
    x<-data[,setFeature[1:numFeature]]
    rfs<-1
  }
  else
  {
    x<-data[,-ncol(data)]
  }
  x<-as.matrix(x)
  xvrem<-x
  
  y<-data[,ncol(data)]
  if(flagMult&&flagTwo)
  {
      y<-y==svalue(comboClass)
      y[y]<-1
      y<-factor(y)
  }
  else
    y <- factor(as.numeric(data[,ncol(data)]) - 1)
  out=createFolds(y, k = numfold, list = TRUE, returnTrain = T)
  Res <- list()
  Boost <- list()
  Ens_model<-list()
  start_time <- Sys.time()
  fold=1
  step<-(numfold*M)%/%100
  if(step==0)
    step<--(100%/%(numfold*M))
  else
    step<-step+1
  
  bestAlg=rep(0,length(algorithms))
  names(bestAlg)=algorithms
  
  fit.individual=TRUE
  newpredEns<-rep(0,length(y))
  newprobEns<-rep(0,length(y))
  if(fit.individual)
  {
    newpredInd<-matrix(0,length(y),length(algorithms))
    newprobInd<-matrix(0,length(y),length(algorithms))
  }
  
  
  for(ff in out)
  {
    # # while(gtkEventsPending()){
    # #   gtkMainIteration()
    #   if(svalue(cbutton)=="Stopping..."){
    #     #cbutton$set_value("Stop simulation")
    #     blockHandler(cbutton)
    #     cat("interrupted")
    #     #gtkMainQuit()
    #     return()
    #   }
    # #}
    #------ select features---------rfs=2---------
    if(rfs==2)
    {
      res<-calcFeature(data=data[ff,], type="Information gain",status=status_bar,
                     flagMult=flagMult,flagTwo=flagTwo,comboClass=comboClass)
      x<-xvrem[,res[1:numFeature,"Biomarker"]]
    }
    ens <- ensembleClassifier(x[ff,], y[ff], M=M, fit.individual=fit.individual,
                              lambda1=10, boost=boost,
                              algorithms=algorithms, validation=validation,
                              progress=pb,status=status_bar,fold=fold,step=step,cbutton=cbutton,
                              flagMult=flagMult,flagTwo=flagTwo)
    if(svalue(cbutton)=="Stopping...")
    {
      blockHandler(cbutton)
      svalue(status_bar)<-"interrupted"
      return(FALSE)
    }
    if(is.null(ens))
    {
      return(FALSE)
    }
    
    # predict using the test data
    xtest <- x[-ff,]
    trueClass <- y[-ff]

    pred_model <- predictEns(ens, xtest, trueClass,flag=TRUE,graph=graphRoc,graph1=graphCross,flagMult=flagMult,flagTwo=flagTwo)
    if(!is.null(boost))
    {
      pred_boost <- predictEns(ens, xtest, trueClass,flag=FALSE,graph=graphRoc,graph1=graphCross,flagMult=flagMult,flagTwo=flagTwo)
      Boost=c(Boost,list(pred_boost))
    }

    Res=c(Res,list(pred_model))
    Ens_model=c(Ens_model,ens)
    vrem=table(ens$bestAlg)
    sapply(names(vrem),function(z) bestAlg[z]<<-bestAlg[z]+vrem[z])
    
    newpredEns[-ff]<-pred_model$yhat
    newprobEns[-ff]<-pred_model$prob
    if(length(ens$indModels) > 0)
    {
      newpredInd[-ff,]<-pred_model$indPred
      newprobInd[-ff,]<-pred_model$indProb
    }
    
    cat("The fold number=", fold, "\n")
    fold=fold+1
  }
  newpredEns<-newpredEns-1
  newpredEns<-factor(newpredEns,levels=levels(y))
#------ plot graphs-------------------
  if((flagMult&&flagTwo)||(!flagMult))
  {
  if(!is.null(graphRoc))
    visible(graphRoc)<-TRUE
  auc <- AUC(y, newprobEns, plot=TRUE)
  }
  
  if(!is.null(graphCross))
    visible(graphCross)<-TRUE
  ct <- crosstab(y, newpredEns,dnn = c("True", "Predicted"),drop.level=F)
  ct<-t(ct$tab)
  barp(ct,main="Crosstabulation",ylab="Value",
       names.arg=colnames(ct),col=2:(1+ncol(ct)),ylim=c(0,max(ct)+4))
  addtable2plot(x="topright",y=NULL,ct,bty="o",display.rownames=TRUE,hlines=TRUE,
                vlines=TRUE,title="The table",cex=1,xpad=0.3)
  
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
  if(flagMult&&!flagTwo)
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
  return(TRUE)
}

IntelSystemMain<-function()
{
#-- new for features-------------  
  library(pROC)
  library(RWeka)
  library(FSelector)
  library(arules)
  library(pamr)
  library(rgl)
  library(gtools)
  library(Rcpp)
  library(varrank)
  
  source("plotRoc.curves.R")
  source("ROC_Plot.R")
  
  #new
  source("select.feature.info.R")

  #new
  source("CalculateHUM_seq.R")
  source("CalculateHUM_ROC.R")
  source("CalculateHUM_Ex.R")
  sourceCpp("HUMseq.cpp")
  sourceCpp("CHI2.cpp")
  sourceCpp("CorrF.cpp")
#-- new for features------------- 
  
  library(caret)
  library(gWidgets2)
  #library(dprep)  # also for relief

  source("ensemble.R")

  boost=NULL
  flag_data="DLBCL"
  flagMult=FALSE
  flagFeature<-0
  
  options(guiToolkit = 'RGtk2')
  
  window <- gwindow ("Intelligent data analysis", visible = F)
  #size(window)=c(300,600)
  notebook <- gnotebook (cont = window,tab.pos = 3,closebuttons = F,
                         font.attr=list(foreground.colors="red"))
  group1 <- ggroup(cont = notebook, label = "Data set", horizontal=F)
  group2 <- ggroup(cont = notebook, label = "Feature Selection", horizontal=F)
  group3 <- ggroup(cont = notebook, label = "Clustering", horizontal=F)
  group4 <- ggroup(cont = notebook, label = "Classifier modeling", horizontal=F)
  window$set_rgtk2_font(window$children[[1]]$widget,list(background="yellow"))
#--------Dataset tab----------------------
  Tab1_Dataset=ggroup(container = group1,expand=T)

  Bigroup1=ggroup(container = Tab1_Dataset,horizontal=F,
                 font.attr=list(foreground.colors="red"))
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

a <- gfilebrowse("Upload file",cont=gg,
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
                         #-- the last column is the class column as factor---
                         the_data[,ncol(the_data)] <- factor(the_data[,ncol(the_data)])
                         tbl[,] <- data.frame(the_data[,c(1:10,ncol(the_data))])#data.frame(variables=dfNames, stringsAsFactors=FALSE)#
                         flagMult<<-length(table(the_data[,ncol(the_data)]))>2
                         if(flagMult)
                         {
                          enabled(frmClass)=TRUE
                          comboClass[]<-as.character(levels(the_data[,ncol(the_data)]))
                          svalue(comboClass,index=T)=1
                          svalue(checkClass)=F
                          #cat(comboClass[])
                         }
                         else
                         {
                           comboClass[]<-""
                           svalue(checkClass)=F
                           enabled(frmClass)=FALSE
                         }
                         assign(data_frame_name, the_data, envir = globalenv())
                         svalue(status_bar) <-
                           paste(nrow(the_data), "records saved to variable", data_frame_name)
                         enabled(buttonPlotBar)<-F
                         enabled(buttonDown)<-F
                         visible(frmAttrAUC)<-F
                         visible(frmAttrHUM)<-F
                         enabled(buttonPlotAUC)<-T
                         enabled(buttonPlotHUM)<-T
                       },
                       error = function(e) svalue(status_bar) <- "Could not upload data"
                     )
                   },
                   filter = list(
                     "Tab delimited" = list(patterns = c("*.txt","*.dlm","*.tab")),
                     "All files" = list(patterns = c("*"))
                   )
  )

  addSpace(Bigroup1, 10)
  frmClass <- gframe("Make two-class problem",container=Bigroup1)
  checkClass <- gcheckbox("Select case class", checked=FALSE,container=frmClass,handler=function(h,...)
              {
                cat("flagFeature=",flagFeature,"\n")
                cat("svalue=",svalue(h$obj),"\n")
                if((flagFeature==1&&svalue(h$obj))||(flagFeature==2&&!svalue(h$obj)))
                {
                  enabled(buttonPlotBar)<-F
                  if(visible(frmAttrAUC))
                    enabled(buttonPlotAUC)=F
                  if(visible(frmAttrHUM))
                    enabled(buttonPlotHUM)=F
                }
                else
                {
                  enabled(buttonPlotBar)<-T
                  if(visible(frmAttrAUC))
                    enabled(buttonPlotAUC)=T
                  if(visible(frmAttrHUM))
                    enabled(buttonPlotHUM)=T
                }
              })
  
  comboClass<-gcombobox("", container=frmClass,expand=TRUE)
  enabled(frmClass)=FALSE
  
  blankDF = data.frame(variables=character(0), stringsAsFactors=FALSE)
  frmX <- gframe("Data set", container = Bigroup1,expand=T)
  tbl <- gtable(blankDF,  container=frmX, filter.FUN="manual")
  size(tbl) <- c(1, 400)
  
  
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
#----------------- Feature tab---------------------------------
  Tab1_Feature=ggroup(container = group2,expand=T)
  types=c("Feature ranking", "Feature selection", "Wrapper algorithms")
  smethods=c(list(c("auc","HUM","Chi-square","Information gain",
             "Symmetrical uncertainty","Relief")),list(c("FastFilter","MRMR","CFS",
             "CorrSF","Chi2-algorithm")),list(c("Forward search")))
  names(smethods)=types
  dmethods=c("MDL","equal interval width","equal frequency")
  
  Bigroup2=ggroup(container = Tab1_Feature,horizontal=F)
  size(Bigroup2)=c(300,1)
  frmSelect <- gframe("Feature ranking and subset selection",
                      container = Bigroup2,horizontal=F)
  addSpace(frmSelect, 20)
  frmType <- gframe("Selection type", container = frmSelect)
  comboType<-gcombobox(types, container=frmType,expand=TRUE,handler=function(h,...)
             {
                comboMethod[]<-smethods[[svalue(h$obj)]]
                svalue(comboMethod,index=T)<-1
             })
  addSpace(frmSelect, 20)
  frmMethod <- gframe("Selection method",container = frmSelect)
  comboMethod<-gcombobox(smethods[[svalue(comboType)]], container=frmMethod,expand=TRUE,
                         handler=function(h,...)
                         {
                           if(svalue(h$obj)=="FastFilter")
                             visible(frmThre)=T
                           else
                             visible(frmThre)=F
                           if(svalue(h$obj)=="Chi2-algorithm")
                             visible(frmIncons)=T
                           else
                             visible(frmIncons)=F
                         })
  addSpace(frmSelect, 20)
  frmDisc <- gframe("Discretization method",container = frmSelect)
  comboDisc<-gcombobox(dmethods, container=frmDisc,expand=TRUE)
  
  addSpace(frmSelect, 20)
  frmThre <- ggroup(container = frmSelect,horizontal=F)
  visible(frmThre)=F
  glabel("Threshold for class mutual information",container=frmThre)
  vedit <- gedit("0.3", container=frmThre, coerce.with=as.numeric)
  
  addSpace(frmSelect, 20)
  frmIncons <- ggroup(container = frmSelect,horizontal=F)
  visible(frmIncons)=F
  glabel("Threshold for inconsistensy",container=frmIncons)
  vslider <- gslider(from=0,to=0.33,by=0.005, value=0,container=frmIncons)
  
  addSpace(Bigroup2, 100)
  fsubset<-gbutton("Feature subset", container=Bigroup2)
  fbutton <- gbutton("Stop simulation",cont=Bigroup2,handler=function(h,...) 
  {svalue(h$obj) <- "Stopping..."})
  
  addHandlerClicked(fsubset, handler=function(h,...) {
    if(!svalue(txt_data_frame_name)%in% ls(envir = .GlobalEnv))
    {
      svalue(status_bar) <- paste("Dataset ",
                    svalue(txt_data_frame_name)," is empty")
      return()
    }
    data<-get(svalue(txt_data_frame_name),envir=.GlobalEnv)
    cat(paste(svalue(txt_data_frame_name),"\n"))
    method=svalue(comboMethod)
    type<-svalue(comboType)
    disc<-svalue(comboDisc)
    thresh<-svalue(vedit)
    incons<-svalue(vslider)
    if((method=="FastFilter")&&(length(thresh)==0||thresh>1||thresh<0))
    {
      svalue(status_bar) <-
        paste("Threshold value should be between 0 and 1")
      return()
    }
    # if(method=="auc"&&flagMult&&!svalue(checkClass))
    # {
    #   svalue(status_bar) <-
    #     paste("The method cannot be applied for dataset with more than 2 classes")
    #   return()
    # }
        if(!is.null(graphAUC))
          visible(graphAUC)<-TRUE
        plot.new()
        # svalue(pb)=0
        svalue(chk_legend)<-FALSE
        svalue(chk_percent)<-FALSE
        svalue(chk_auc)<-FALSE
        svalue(status_bar) <-""
        tblRes[,] <- blankDF
        # svalue(cbutton) <- "Stop simulation"
        # unblockHandler(cbutton)
        res<-calcFeature(data=data, type=method, disc=disc, thresh=thresh,
                       incons=incons,status=status_bar,
                       flagMult=flagMult,flagTwo=svalue(checkClass),comboClass=comboClass)
        tblRes[,] <- data.frame(res)
        #cat(str(res),"\n")
        #cat(res[,"Biomarker"])
        svalue(status_bar) <-"Ready"
        enabled(buttonPlotBar)<-T
        enabled(buttonDown)<-T
        if(flagMult&&!svalue(checkClass))
        {
          visible(frmAttrAUC)=FALSE
          visible(frmAttrHUM)=TRUE
          enabled(buttonPlotHUM)=TRUE
          comboSelect[]<-as.character(tblRes[,"Biomarker"])
          svalue(comboSelect,index=T)=1
          flagFeature<<-1
        }
        else
        {
          visible(frmAttrAUC)=TRUE
          enabled(buttonPlotAUC)=TRUE
          visible(frmAttrHUM)=FALSE
          flagFeature<<-2
        }
  })
  grp_res<-ggroup(container = Tab1_Feature,horizontal=F,expand=T)
  grp_tbl<-ggroup(container = grp_res,expand=T)
  frmRes <- gframe("Feature selection results", container = grp_tbl,expand=T)
  tblRes <- gtable(blankDF,  container=frmRes, multiple = T, filter.FUN="manual")
  #size(tblRes)=c(300,200)
  frmAttrPlot <- ggroup(container = grp_tbl,horizontal=F)
  gButtons=gframe("",container=frmAttrPlot)
  buttonPlotBar<-gbutton("Bar Plot", container=gButtons)
  enabled(buttonPlotBar)<-F
  addSpace(gButtons, 50)
  buttonDown<-gbutton("Download", container=gButtons)
  enabled(buttonDown)<-F
  
  frmAttrAUCHUM <- ggroup(container = frmAttrPlot)
  frmAttrAUC <- ggroup(container = frmAttrAUCHUM,horizontal=F)
  visible(frmAttrAUC)=F
  glabel("AUC plot attributes",container=frmAttrAUC)
  chk_legend <- gcheckbox(text="Add Legends",checked=F,container=frmAttrAUC)
  chk_percent<-gcheckbox(text="Percent",checked=F,container=frmAttrAUC)
  chk_auc<-gcheckbox(text="Show the auc-values",checked=F,container=frmAttrAUC)
  addSpace(frmAttrAUC, 20)
  buttonPlotAUC<-gbutton("Plot AUC", container=frmAttrAUC)
  
  frmAttrHUM <- ggroup(container = frmAttrAUCHUM,horizontal=F,expand=T)
  visible(frmAttrHUM)=F
  #size(frmAttrHUM)=c(100,1)
  gHUM<-gframe("Select feature",container=frmAttrHUM)
  comboSelect <- gcombobox("",container=gHUM,expand=T)
  addSpace(frmAttrHUM, 20)
  buttonPlotHUM<-gbutton("Plot HUM", container=frmAttrHUM)
  
  frmgraph2 <- gframe("Plots", container = grp_res,expand=T)
  size(frmgraph2)=c(1,200)
  graphAUC <- ggraphics("", container=frmgraph2,expand=T)
  
  addHandlerClicked(buttonDown, handler=function(h,...) {
    res<-tblRes[,]
    cat(str(res))
    name<-paste(svalue(txt_data_frame_name),"feature",sep="_")
    if(flagMult)
    {
      if(flagFeature==2)  name<-paste(name,2,sep="_")
    }
    write.table(data.frame(res),file=name,sep="\t",row.names=F)
  })
  
  addHandlerClicked(buttonPlotBar, handler=function(h,...) {
    cat("Enter\n")
    data<-get(svalue(txt_data_frame_name),envir=.GlobalEnv)
    if(flagMult&&svalue(checkClass))
    {
      y<-data[,ncol(data)]==svalue(comboClass)
      y[y]<-1
      y<-factor(y)
      data[,ncol(data)]=y
    }
    class=data[,ncol(data)]
    label=levels(class)
    if(length(label)<2) return()
    
    index<-svalue(tblRes,index=T)
    if(length(index)<1) return()
    attrs<- levels(tblRes[,"Biomarker"])[tblRes[index,"Biomarker"]]
    attrs<-sapply(attrs,function(z){ if(regexpr("X",z)==1) out<-substr(z,2,nchar(z)) else out<-z})
    cat("attr=",attrs)
    
    attrs.no<-numeric()
    for(i in 1:length(attrs)){
      if(!is.factor(data[,attrs[i]]))
      {
        attrs.no<-c(attrs.no,grep(attrs[i],colnames(data),fixed=T)[1])
      }
    }
    if(length(attrs.no)<1) return()
    
    out=CalculateHUM_seq(data,attrs.no,ncol(data),label)
    HUM<-out$HUM
    barplot(HUM,main="AUC/HUM values-Bar Chart")
    y<-1/factorial(length(label))
    abline(h=y,col="red",lwd=2)
    text(1,y,paste("Threshold=",sprintf("%.3f",y),sep=""),col="red",adj = c(0, -.2))
  })
  
  addHandlerClicked(buttonPlotHUM, handler=function(h,...) {
    if(!svalue(txt_data_frame_name)%in% ls(envir = .GlobalEnv))
    {
      svalue(status_bar) <- paste("Dataset ",
                                  svalue(txt_data_frame_name)," is empty")
      return()
    }
    data<-get(svalue(txt_data_frame_name),envir=.GlobalEnv)
    class=data[,ncol(data)]
    label=levels(class)
    cat(label)
    if(length(label)!=3) return()
    attrs<-svalue(comboSelect)
    cat("attrs=",attrs)
    if(is.factor(data[,attrs])) return()
    out=CalculateHUM_seq(data,attrs,ncol(data),label)
    seq<-out$seq
    HUM<-out$HUM
    out=CalculateHUM_ROC(data,attrs,ncol(data),label,seq)
    print(Calculate3D(attrs,out$Sn,out$Sp,out$S3,out$optSn,out$optSp,out$optS3,out$thresholds,HUM,label[seq]))
  })
  
  addHandlerClicked(buttonPlotAUC, handler=function(h,...) {
    if(!svalue(txt_data_frame_name)%in% ls(envir = .GlobalEnv))
    {
      svalue(status_bar) <- paste("Dataset ",
                                  svalue(txt_data_frame_name)," is empty")
      return()
    }
    cat("Enter\n")
    data<-get(svalue(txt_data_frame_name),envir=.GlobalEnv)
    if(flagMult&&svalue(checkClass))
    {
      y<-data[,ncol(data)]==svalue(comboClass)
      y[y]<-1
      y<-factor(y)
      data[,ncol(data)]=y
    }
    class=data[,ncol(data)]
    label=levels(class)
    cat(label)
    if(length(label)>2) return()
    
    index<-svalue(tblRes,index=T)
    if(length(index)<1) return()
    
    cat("index=",index,"\n")
    #cat("tblRes=",levels(tblRes[,"Biomarker"])[tblRes[index,"Biomarker"]],"\n")
    attrs<- levels(tblRes[,"Biomarker"])[tblRes[index,"Biomarker"]]
    attrs<-sapply(attrs,function(z){ if(regexpr("X",z)==1) out<-substr(z,2,nchar(z)) else out<-z})
    cat("attr=",attrs)
    
    attrs.no<-numeric()
    for(i in 1:length(attrs)){
      if(!is.factor(data[,attrs[i]]))
      {
        attrs.no<-c(attrs.no,grep(attrs[i],colnames(data),fixed=T)[1])
      }
    }
    cat("attrs.no=",attrs.no)
    if(length(attrs.no)<1) return()
    
    add.legend<-svalue(chk_legend)
    is.percent<-svalue(chk_percent)
    include.auc<-F
    if(add.legend){
      include.auc<-svalue(chk_auc)
    }
    roc.curve<-plotRoc.curves(data[,c(attrs.no,ncol(data))],add.legend=add.legend,
                              include.auc=include.auc,ispercent=is.percent)
    print(roc.curve)
  })
#----------------- Feature tab------------------------  
#------------Classifier Tab---------------------------
  models <- c("svm", "rf", "lda", "qda", "lr", "nnet", "pls_lda", "pls_qda", "pls_rf", "pls_lr",
               "pca_lda", "pca_qda", "pca_rf", "pca_lr", "rpart", "adaboost", "plr", "pls")
  measures=c("accuracy", "sensitivity", "specificity")
  algorithms=NULL
  validation=NULL
  numfold=3
  M=51
  setFeature=NULL

  Tab1_Classifier=ggroup(container = group4,expand=T)
  
  Bigroup4=ggroup(container = Tab1_Classifier,horizontal=F)
  size(Bigroup4)=c(300,1)
  
  # frmClass <- gframe("Make two-class problem",container=Bigroup4)
  # checkClass <- gcheckbox("Select case class", checked=FALSE,container=frmClass)
  # 
  # comboClass<-gcombobox("", container=frmClass,expand=TRUE)
  # enabled(frmClass)=FALSE
  frmSel <- gframe("Feature selection", container = Bigroup4)
  grb <- gradio(c("Without selection","Embedded selection","Load from file"), container=frmSel,horizontal = FALSE)
  gbtn<-ggroup(container=frmSel,horizontal=FALSE,expand=T)
  gnumFeature<-gslider(from=50,to=500,by=50, value=100,container=gbtn)
  
  buttonLoad<-gbutton("Load",container=gbtn)
  addHandlerClicked(buttonLoad, handler=function(h,...) {
    if(!svalue(txt_data_frame_name)%in% ls(envir = .GlobalEnv))
    {
      svalue(status_bar) <- paste("Dataset ",
                                  svalue(txt_data_frame_name)," is empty")
      return()
    }
    tryCatch(
      {
        data_file_name <- make.names(svalue(txt_data_frame_name))
        data_file_name<-paste(data_file_name,"feature",sep="_")
        if(flagMult&&svalue(checkClass))
        {
          data_file_name<-paste(data_file_name,2,sep="_")
        }
        res<-read.table(data_file_name,header=T,sep = "\t")
        data<-get(svalue(txt_data_frame_name),envir=.GlobalEnv)
      },
      error = function(e) {svalue(status_bar) <- paste("Dataset ",data_file_name," is not exist");break}
    )
    setFeature<<-res[,"Biomarker"]
  })

  frmM <- gframe("Models", container = Bigroup4)               
  # cmodels <- gcheckboxgroup(models, container=frmM, handler = function(h,...){
  #   algorithms <<- svalue(h$obj)})#,index=TRUE)}))
  cmodels <- gtable(models, container=frmM, multiple = T,expand=T)#,index=TRUE)}))
  size(cmodels)=c(1,200)
  frmV <- gframe("Validation", container = Bigroup4)
  cmeasures<-gcheckboxgroup(measures, container=frmV, handler = function(h,...){
    validation <<- svalue(h$obj)})#,index=TRUE)})
  foldstrap<-gslider(from=3,to=10,by=1, value=3, handler = function(h,...){
    numfold <<- svalue(h$obj)})
  frmF<-gframe("Crossvalidation folds", container=Bigroup4)
  add(frmF,foldstrap,expand=TRUE)
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
    #cat(paste(svalue(txt_data_frame_name),"\n"))
      algorithms=svalue(cmodels, drop=FALSE)
      cat("Validation: ",validation, length(validation),"\n")
      if(length(algorithms)==0||length(validation)==0)
      {
        svalue(status_bar) <-
          paste("Models or Validation are not defined")
      }
      else
      {
        if(flagMult&&!svalue(checkClass)&&any(validation%in%c("sensitivity","specificity")))
        {
          svalue(status_bar) <-
            paste("For datasets with more than 2 classes validation can be only accuracy")
        }
        else
        {
          rfs<-svalue(grb,index=T)
          if(rfs==3)
          {
            if(is.null(setFeature))
            {
              svalue(status_bar) <- paste("Feature set is empty")
              return()
            }
          }
          if(!is.null(graphRoc))
            visible(graphRoc)<-TRUE
          plot.new()
          if(!is.null(graphCross))
            visible(graphCross)<-TRUE
          plot.new()
          svalue(textsummary)=""
          svalue(pb)=0
          svalue(status_bar) <-""
          svalue(cbutton) <- "Stop simulation"
          unblockHandler(cbutton)
          algorithms=algorithms[,1]
          if(calcEnsemble(data, setFeature, svalue(gnumFeature),rfs, M, algorithms, validation, boost,textsummary,
                       numfold,cbutton,pb,status_bar,graphRoc,graphCross,flagMult=flagMult,flagTwo=svalue(checkClass),comboClass=comboClass))
            svalue(status_bar) <-"Ready"
        }
      }
  })
  resgroup=ggroup(container=Tab1_Classifier,horizontal=F,expand=T)
  textsummary <- gtext("", container=resgroup,expand=T)
  frmgraph4 <- gframe("Plots", container = resgroup,expand=T)
  size(frmgraph4)=c(1,200)
  graphRoc <- ggraphics("", container=frmgraph4,expand=T)
  graphCross <- ggraphics("", container=frmgraph4,expand=T)
#------------Classifier Tab---------------------------
  
  
  svalue(notebook) <- 1
  visible(window) <- TRUE
}