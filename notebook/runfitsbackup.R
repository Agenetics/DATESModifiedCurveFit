runfits<-function(directorylist=NULL,pvalue=0.001,foldernameformat="NoFormat",outputfilename=NULL){

  
  out1=data.frame()
  len=length(directorylist)
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                      max = len, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar

  for (i in 1:len[[1]]){

      b=gregexpr(pattern ='/',directorylist[i])
      c=length(b[[1]])
      str=substring(directorylist[i],b[[1]][c]+1)
      filelist=list.files(directorylist[i], ".pdf", full.names=FALSE)
      c=nchar(filelist[[1]])
      folder=paste(directorylist[i],"/",substring(filelist[[1]],1,c-3),"out:")
      folder=gsub(" ","",folder)
      headers=strsplit(str,"_")
    
      ans=trysinglefit(folder,0)
    
    if (foldernameformat!="NoFormat"){
      for (j in 1:4){
        out1[i,j]=headers[[1]][j]}
    }    
    else {
      out1[i,1]=str
      out1[i,2:4]=""
    }
    
    out1[i,5]=ans$meant1
    out1[i,6]=ans$t1SE
    out1[i,7]=ans$shapiro$p.value
    out1[i,8:12]=""
    
    if (ans$shapiro$p.value < pvalue){
      ans1=try100twofit(folder,0)
      if (!is.null(ans1)){
          if (ans1$meant1>ans1$meant2){
            out1[i,8]=ans1$meant1
            out1[i,9]=ans1$meant2
            out1[i,10]=ans1$t1SE
            out1[i,11]=ans1$t2SE
            out1[i,12]=ans1$shapiro$p.value
          }
        else {
          out1[i,8]=ans1$meant2
          out1[i,9]=ans1$meant1
          out1[i,10]=ans1$t2SE
          out1[i,11]=ans1$t1SE
          out1[i,12]=ans1$shapiro$p.value
        }
        }
    }
    setTxtProgressBar(pb, i)
  } 
  close(pb)
  
  if (foldernameformat!="NoFormat"){
  colnames(out1)=c("T1","T2","adm%1","adm%2","T est", "T SE", "SWilk p-val","T1 Est","T2 Est","T1 SE","T2 SE","SWilk P-val2")
  }
  else {
    colnames(out1)=c("FolderName","NA","NA","NA","T est", "T SE", "SWilk p-val","T1 Est","T2 Est","T1 SE","T2 SE","SWilk P-val2")
    
  }
  if (!is.null(outputfilename)){
    write_tsv(out1,outputfilename)
  }
  return(out1)
}

twofit<-function(prefix,aa=0.05,bb=0.04, tt=125,t22=10){
  
  twofit<-data.frame()
  spgt<-list()
  for (i in 1:22){
    strr=paste(prefix,i)
    strr=gsub(" ", "", strr)
    spgt[[i]]=read.table(strr)
    spgt[[i]]=spgt[[i]][,-c(2,4,5)]
    spgt[[i]]=spgt[[i]][5:300,]
  }
  
  nonlin2 <- function(d, a,b, t,t2) { a * exp(d/100*(-t+1)) + b * exp(d/100*(-t2+1))}
  start=list(a=aa,b=bb,t=tt,t2=t22)
  for (i in 1:22){
    ##c.0<-min(spgt[[i]]$V3)-0.0000001
    nlc <- nls.control(maxiter = 1000)
    nlsfit <- nls(V3 ~ nonlin2(V1,a,b,t,t2), data=spgt[[i]], start=start,control = nlc)
    for (j in 1:4){
      twofit[i,j]=coef(nlsfit)[j]
    }
  }
  out1<-list()
  out1$meant1=mean(twofit$V3)
  out1$meant2=mean(twofit$V4)
  out1$alpha=mean(twofit$V1)
  out1$beta=mean(twofit$V2)
  out1$t1SE <- sqrt(21/22 * sum((twofit[,3]-out1$meant1)^2))
  out1$t2SE <- sqrt(21/22 * sum((twofit[,4]-out1$meant2)^2))
  
  strr=gsub(".out:",".fit",prefix)
  hello=read.table(strr,header=FALSE)
  hello=hello[,-c(3,4)]
  
  nlsfit <- nls(V2 ~ nonlin2(V1,a,b,t,t2), data=hello, start=start,control = nlc)
  hello$fit=predict(nlsfit)
  hello$residual=residuals(nlsfit)
  out1$directfit=coef(nlsfit)
  
  hello=hello[c(1:296),]
  out1$shapiro=shapiro.test(hello$residual)
  out1$LL=-296/2*(log(2*pi)+1-log(296)+log(sum(hello$residual^2)))
  
  out1$model=nlsfit
  out1$outdata=hello
  
  rm(twofit,spgt,i,j,nlsfit,nonlin2,strr)
  
  return(out1)
}

singlefit<-function(prefix,aa=0.03,tt=10){
  
  singlefit<-data.frame()
  spgt<-list()
  for (i in 1:22){
    strr=paste(prefix,i)
    strr=gsub(" ", "", strr)
    spgt[[i]]=read.table(strr)
    spgt[[i]]=spgt[[i]][,-c(2,4,5)]
    spgt[[i]]=spgt[[i]][5:300,]
  }
  ## 
  
  nonlin <- function(d, a, t) {  a * exp(d/100*(-t+1))}
  start=list(a=aa,t=tt)
  
  for (i in 1:22){
    ##c.0 <- min(spgt[[i]]$V3)-0.0000001
    ##model.0 <- lm(log(V3 - c.0) ~ V1, data=spgt[[i]])
    ##start <- list(a=exp(coef(model.0)[1]), t=1-100*coef(model.0)[2], c=c.0)
    nlc <- nls.control(maxiter = 1000)
    nlsfit <- nls(V3 ~ nonlin(V1,a,t), data=spgt[[i]], start=start,control = nlc)
    for (j in 1:2){
      singlefit[i,j]=coef(nlsfit)[j]
    }
  }
  out1<-list()
  out1$meant1=mean(singlefit$V2)
  out1$t1SE <- sqrt(21/22 * sum((singlefit[,2]-out1$meant1)^2))
  out1$alpha<- mean(singlefit$V1)
  ##out1$constantC<-mean(singlefit$V3)
  
  strr=gsub(".out:",".fit",prefix)
  hello=read.table(strr,header=FALSE)
  hello=hello[,-c(3,4)]
  
  nlsfit <- nls(V2 ~ nonlin(V1,a,t), data=hello, start=start,control = nlc)
  hello$fit=predict(nlsfit)
  hello$residual=residuals(nlsfit)
  out1$directfit=coef(nlsfit)
  hello=hello[c(1:296),]
  out1$shapiro=shapiro.test(hello$residual)
  out1$ks=ks.test(hello$residual,rnorm(length(hello$residual),0,sd(hello$residual)))
  out1$LL=-296/2*(log(2*pi)+1-log(296)+log(sum(hello$residual^2)))
  out1$model=nlsfit
  out1$outdata=hello
  
  
  
  rm(singlefit,spgt,i,j,nlsfit,nonlin,strr)
  
  return(out1)
}

trytwofit<-function (folder,i=0){
  if (i==50){
    return(NULL)
  }      
  aa1=runif(1,min=-0.02,max=0.1)
  bb1=runif(1,min=-0.02,max=0.1)
  year1=runif(1,min=1,max=300)
  year2=runif(1,min=1,max=300)
  
  ans1=try(twofit(folder,aa1[[1]],bb1[[1]],year1[[1]],year2[[1]]),silent = TRUE)
  if("try-error" %in% class(ans1)) ans1=trytwofit(folder,i+1)
  return(ans1)
  
}

trysinglefit<-function(folder,i=0){
  if (i==50){
    return(NULL)
  }      
  aa1=runif(1,min=-0.02,max=0.1)
  year1=runif(1,min=1,max=300)
  
  ans=try(singlefit(folder,aa=aa1[[1]],tt=year1[[1]]),silent=TRUE)
  if("try-error" %in% class(ans)) ans=trysinglefit(folder,i+1)
  return(ans)
  
}

try100twofit<-function(folder,i=0,jj=50){
  
  z=0
  ans2=NULL
  ans1=NULL
  for (j in 1:jj){
    ans2=trytwofit(folder,i=0)
    if (!is.null(ans2)){
      if ((ans2$meant1/ans2$t1SE + ans2$meant2/ans2$t2SE)>z){
        ans1=ans2 
        z=ans2$meant1/ans2$t1SE + ans2$meant2/ans2$t2SE
      }
    }
    
  }
  return(ans1)
}