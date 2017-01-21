Addvec1<-function(X){
  return (cbind(rep(1,nrow(X)),X))
}

Hentr<-function(x){
  y=-x*log(x)-(1-x)*log(1-x)
  y[is.nan(y)==TRUE]=0
  return (y)

}



Signoid<-function(x){
  return (1/(1+exp(-x)))
}

Lexp<-function(x){
  return (log(1+exp(-x)))
}

F<-function(w,X1,X2=NA,t1,lambda1=0,lambda2=0,alpha){

  if(is.na(X2[1])==FALSE){
    return (sum(Lexp(t1*(X1%*%w)))
            + lambda1*( (1-alpha)/2*sum(w^2) + alpha*sum(abs(w)))
            + lambda2*sum(Hentr(Signoid(X2%*%w))))

  }else{
    return (sum(laply(t1*(X1%*%w),Lexp))
            +lambda1*( (1-alpha)/2*sum(w^2) + alpha*sum(abs(w))))


  }
}





GetPredict<-function(w,Xt){
  return (sign(sign(Xt%*%w)+.5))
}

MissErr<-function(tpred,t){
  if(length(tpred)!=length(t))stop("lengths are different in MissErr!")
  return (sum(abs(matrix(tpred,ncol=1)-matrix(t,ncol=1))/(2*length(t))) )
}




getlambda1seq<-function(w0,X1,X2,Total,t1,lambda2,
                        sparsete,lambda1.min.ratio,nlambda1,te=te,wn=FALSE){
  for(i in 4:6){
    tempresult<-PALM(w0,X1,X2,t1,lambda1=10^(i/2),lambda2=lambda2,alpha=1,te=te)

    if(sum(abs(tempresult)/Total<sparsete)){
      lambda1max =10^(i/2)
    }
  }
  if(lambda1max ==10^(i/2)&&wn==TRUE)warning("lamdba1 reaches the maximum!")
  return (exp(seq(from=log(lambda1max),
                  to=log(lambda1max*lambda1.min.ratio),length=nlambda1)))
}

writetxt<-function(X1,X2,X3=NULL,n1,n2,dim,t1,t3,
                   svmlight=FALSE,svm_light_dir,svmlin=FALSE,svmlin_dir,universvm=FALSE){
  twd=getwd()
  n3=length(t3)
  if(svmlight==TRUE){
    setwd(svm_light_dir)
    sink("train.txt")
    for(i in 1: n1){
      cat(c(t1[i]," "))
      for(j in 1:length(X1[1,])){
        cat(paste0(j,":",X1[i,j]," "))
      }
      cat("\n")
    }
    for(i in 1:n2){
      cat(0," ")
      for(j in 1:length(X1[1,])){
        cat(paste0(j,":",X2[i,j]," "))
      }
      cat("\n")
    }
    sink()
    sink("test.txt")
    for(i in 1:n3){
      cat(t3[i]," ")
      for(j in 1:length(X1[1,])){
        cat(paste0(j,":",X3[i,j]," "))
      }
      cat("\n")
    }
    sink()
  }
  if(svmlin==TRUE){
    setwd(svmlin_dir)
    sink("training_examples")
    for(i in 1: n1){
      for(j in 1:(length(X1[1,])-1)){
        cat(paste0(j,":",X1[i,j+1]," "))
      }
      cat("\n")
    }
    for(i in 1:n2){
      for(j in 1:(length(X1[1,])-1)){
        cat(paste0(j,":",X2[i,j+1]," "))
      }
      cat("\n")
    }
    sink()
    sink("training_labels")
    for(i in 1: n1){
      cat(t1[i],"\n")
    }
    for(i in 1:n2){
      cat(0,"\n")
    }
    sink()
    sink("test_examples")
    for(i in 1:n3){
      for(j in 1:(length(X1[1,])-1)){
        cat(paste0(j,":",X3[i,j+1]," "))
      }
      cat("\n")
    }
    sink()
    sink("test_labels")
    for(i in 1:n3){
      cat(t3[i],"\n")
    }
    sink()
  }
  if(universvm==TRUE){
    sink("train_file")
    for(i in 1: n1){
      cat(c(t1[i]," "))
      for(j in 1:(length(X1[1,])-1)){
        cat(paste0(j,":",X1[i,j+1]," "))
      }
      cat("\n")
    }
    for(i in 1:n2){
      cat(-3," ")
      for(j in 1:(length(X1[1,])-1)){
        cat(paste0(j,":",X2[i,j+1]," "))
      }
      cat("\n")
    }
    sink()
    sink("test_file")
    for(i in 1:n3){
      cat(t3[i]," ")
      for(j in 1:(length(X1[1,])-1)){
        cat(paste0(j,":",X3[i,j+1]," "))
      }
      cat("\n")
    }
    sink()
  }
  setwd(twd)
}



testsepa<-function(import=FALSE,X=NULL,X1=NULL,X2=NULL,X3=NULL,t1=NULL,t2=NULL,t3=NULL,dim=NULL,w0=NULL,standardize=FALSE,
                   n1,n2,n3,Total,pnr=.5,transductive=FALSE,lambda1,lambda2,alpha,
                   alpha_1=0.2,beta_1=0.8,seed=0,te,
                   goodlambda=FALSE,useunlabel=FALSE,sa,showcoefcmp=FALSE,detail=FALSE,psparse=.1,wt=NA,
                   showlbd1=FALSE,crv=FALSE,type.measure="c",lambda1.min.ratio=.0001,nlambda1=20,lambda1seq=NA,
                   warmstart=FALSE,warminit=FALSE,writeout=FALSE,svmlight=FALSE,svm_light_dir="/Users/user/Downloads/svm_light/",svmlin=FALSE,
                   svmlin.mode=NULL,svmlin_dir="/Users/user/Downloads/svmlin-v1.0/",
                   universvm=FALSE,universvm.mode=1,os="windows",sparsete=.001,palmquiet=TRUE,nlambda2=8,acccomp=FALSE,rho=0.4){

  set.seed(seed)
  te=te*(n1+n2)
  if(import==FALSE){print("Doing simulation")
    Total=n1+n2+n3
    Beta<-rep(0,dim+1)
    Beta[sample(c(2:(dim+1)),psparse*dim,replace=FALSE)]=c(1+runif(psparse*dim/2),-1-runif(psparse*dim/2))
    Edcrm=rho^( abs( matrix(1:dim,nrow=dim,ncol=dim,byrow=TRUE) - matrix(1:dim,nrow=dim,ncol=dim,byrow=FALSE) ))
    X<-mvrnorm(n=n1*pnr,mu=-sa*Beta[-1],Sigma=Edcrm)

    X<-rbind(X,mvrnorm(n=n1*(1-pnr),mu=sa*Beta[-1],Sigma=Edcrm))
    if(n2>0){
      X<-rbind(X,mvrnorm(n=n2*pnr,mu=-sa*Beta[-1],Sigma=Edcrm))
      X<-rbind(X,mvrnorm(n=n2*(1-pnr),mu=sa*Beta[-1],Sigma=Edcrm))
    }
    if(n3>0){
      X<-rbind(X,mvrnorm(n=n3*pnr,mu=-sa*Beta[-1],Sigma=Edcrm))
      X<-rbind(X,mvrnorm(n=n3*(1-pnr),mu=sa*Beta[-1],Sigma=Edcrm))
    }
    X=Addvec1(X)
    X1<-X[1:n1,]
    t<-2*rbinom(Total,1,1/(1+exp(-X%*%Beta)))-1
    t1<-t[1:n1]
    if(n2>0){
      X2<-X[(n1+1):(n1+n2),]
      t2<-t[(n1+1):(n1+n2)]
      if(transductive==TRUE){
        X3=X2
        t3=t2
      }else{
        X3=X[(n1+n2+1):Total,]
        t3=t[(n1+n2+1):Total]
      }
    }
    if(standardize==TRUE){
      X1=scale(X1[,-1])
      X2=scale(X2[,-1])
      X3=scale(X3[,-1])

    }else{
      X1=X1[,-1]
      X2=X2[,-1]
      X3=X3[,-1]
    }
  }else{
    arry=sample(x=1:(n1+n2+n3),n1+n2+n3,replace=FALSE)
    t1= Label[ arry[1:n1]  ]
    t2= Label[ arry[(n1+1):(n1+n2)] ]
    t3= Label[ arry[(n1+n2+1):(n1+n2+n3)]  ]
    X1= t( fdata[ ,arry[1:n1] ] )
    X2= t( fdata[ ,arry[(n1+1):(n1+n2)] ] )
    X3= t( fdata[ ,arry[(n1+n2+1):(n1+n2+n3)] ] )
    if(standardize==TRUE){
      X1=scale(X1)
      X2=scale(X2)
      X3=scale(X3)

    }
  }


  X1=Addvec1(X1)
  X2=Addvec1(X2)
  X3=Addvec1(X3)


  if(writeout==TRUE&n2>0){
    writetxt(X1,X2,X3,n1,n2,dim,t1,t3,svmlight,svm_light_dir,svmlin,svmlin_dir,universvm)
  }
  svm_light_e<-NA
  svmlin_l2_e<-NA
  svmlin_DA_e<-NA
  if(os=="linux"){
    if(svmlight==TRUE){
      if(writeout==FALSE)stop("writeout must be TRUE if svmlight is TRUE!")
      setwd(svm_light_dir)
      #system(paste0("cd ",svm_light_dir))
      system("./svm_learn train.txt model.txt")
      system("./svm_classify test.txt model.txt prediction.txt")
      svm_light_pred<-read.table(file="prediction.txt")
      print(svm_light_pred)
      svm_light_e<-MissErr(sign(svm_light_pred[,1]),t3)
      print("svmlight miser:")
      print(svm_light_e)
    }
    if(svmlin==TRUE){
      setwd(svmlin_dir)
      if("l2hatloss"%in%svmlin.mode){
        system(paste0("./svmlin -A ",2," -W 1 -U 1 -R 0.5 training_examples training_labels"))
        system("./svmlin -f  training_examples.weights  test_examples test_labels")
        svmlin_pred<-read.table(file="test_examples.outputs")
        svmlin_l2_e<-MissErr(sign(svmlin_pred[,1]),t3)
        print("svmlin l2hatloss miser:")
        print(svmlin_l2_e)
      }
      if("DA"%in%svmlin.mode){
        system(paste0("./svmlin -A ",3," -W 1 -U 1 -R 0.5 training_examples training_labels"))
        system("./svmlin -f  training_examples.weights  test_examples test_labels")
        svmlin_pred<-read.table(file="test_examples.outputs")
        svmlin_DA_e<-MissErr(sign(svmlin_pred[,1]),t3)
        print("svmlin DA miser:")
        print(svmlin_DA_e)
      }
    }
    if(universvm==TRUE){
      system("./universvm -o 1 -T test_file train_file")
    }
  }else{#os=="windows"
    if(svmlight==TRUE){
      if(writeout==FALSE)stop("writeout must be TRUE if svmlight is TRUE!")
      system("svm_learn train.txt model.txt")
      system("svm_classify test.txt model.txt prediction.txt")
      svm_light_pred<-read.table(file="prediction.txt")
      svm_light_e<-MissErr(sign(svm_light_pred),t3)
      print("svmlight miser:")
      print(svm_light_e)
    }
    if(svmlin==TRUE){
      system(paste0("svmlin -A ",svmlin.mode," -W 1 -U 1 -R 0.5 training_examples training_labels"))
      system("svmlin -f  training_examples.weights  test_examples test_labels")
      svmlin_pred<-read.table(file="test_examples.outputs")
      svmline<-MissErr(sign(svmlin_pred[,1]),t3)
      print("svmlin miser:")
      print(svmline)
    }
    if(universvm==TRUE){
      system("universvm -o 1 -T test_file train_file")
    }
  }
  if(is.null(w0)){
    w0<-rep(0,dim+1)
    if(length(w0)!=ncol(X1)){print("dimension error!")}
  }else{
    w0=w0
  }




  lbd1<-NULL
  if(goodlambda==TRUE){
    cvfit<-cv.glmnet(x=X1[,-1],y=as.factor(t1),family="binomial",standardize=FALSE,
                     type.measure = "class",nfolds=5)
    lbd1=n1*cvfit$lambda.1se
    if(showlbd1==TRUE){
      print("lbd1=")
      print(lbd1)
      plot(cvfit)
      print("\n")
    }
  }
  else {
    if(crv==FALSE){
      lbd1=lambda1
      lbd2=lambda2
      alpha=alpha
    }else{
      if(warmstart==FALSE){
        if(is.na(lambda1seq)){
          lambda1seq=getlambda1seq(w0,X1,X2,Total,t1,lambda2=1,sparsete,lambda1.min.ratio,nlambda1,te=te)
        }
        in1<-sample(1:n1)
        in2<-sample(1:n2)
        l1<-list(0,0,0,0,0)
        l2<-list(0,0,0,0,0)
        for(i in 1:5){
          l1[[i]]<-in1[floor((i-1)*n1/5+1):floor(i*n1/5)]
          l2[[i]]<-in2[floor((i-1)*n2/5+1):floor(i*n2/5)]
        }

        Minl1<-matrix(nrow=nlambda2,ncol=3)
        startw=NULL
        lcount<-0
        for(k in 0:nlambda2){
          Fseq<-matrix(ncol=6,nrow=length(lambda1seq))
          for(j in 1:length(lambda1seq)){
            for(al in 0:5){
              Fval<-0
              for(i in 1:5){
                resu<-PALM(w0,X1[-l1[[i]],],X2[-l2[[i]],],t1[-l1[[i]]],lambda1=lambda1seq[j],lambda2=0.25*k,alpha=al*0.2,te=te)
                if(type.measure=="c"){
                  #    print(paste0( "missclassification error on jth lambda1 with ",length(l1[[i]])," samples for ith cv: "))
                  #    print(MissErr(GetPredict(resu[[1]],X1[l1[[i]],]),t1[l1[[i]]]))
                  Fval=Fval+MissErr(GetPredict(resu,X1[l1[[i]],]),t1[l1[[i]]])
                }else{
                  tempProb=F(resu,X1[l1[[i]],],NA,t1[l1[[i]]],lambda1=0,lambda2=0,alpha=1)
                  #   print(paste0( "neg-log likelihood on jth lambda1 with ",length(l1[[i]])," samples for ith cv: "))
                  #   print(tempProb)
                  Fval=Fval+tempProb
                  startw=c(startw,resu)
                }

              }#for i
              print((al+1+(j-1)*6 + length(lambda1seq)*6*(k)  )/(6*length(lambda1seq)*(nlambda2+1)))
              Fseq[j,al+1]<-Fval/5#! al+1 not al !
            }#for al
          }#for j
          that=which(Fseq==min(Fseq))[1]
          Minl1[k,]=c(min(Fseq),lambda1seq[ nrow(Fseq)-(-that%%nrow(Fseq)) ],ceiling(that/nrow(Fseq)))

        }# for k
        lbd2=0.25 *  which(Minl1[,1]==min(Minl1[,1]))[1]
        lbd1=Minl1[which(Minl1[,1]==min(Minl1[,1]))[1],2]
        alpha=Minl1[which(Minl1[,1]==min(Minl1[,1]))[1],3]*0.2-0.2
        print("the lbd1,lbd2 and alpha by cv is")
        print(c(lbd1,lbd2,alpha))


      }else{#warmstart==TRUE
        if(acccomp==FALSE){
          l1seq=c(10,8,6,4,3,2,1.5,1.2,1,0.9,0.8,0.75,0.7,0.65)
          w0s=w0
          for(i in 1:length(l1seq)){
            resu<-iPALM(w0s,X1,X2,t1,lambda1=l1seq[i],lambda2=1,alpha=alpha,te=te,palmquiet=FALSE,maxtimes=500)
            tmyu<-GetPredict(resu,X3)
            myeu<-MissErr(tmyu,t3)
            if(warminit==TRUE){
              w0s=resu
            }
            l1seq[i]=myeu
          }
          print("warstart")
          l2seq=seq(from=0,to=2,by=0.1)
          w0s=w0
          for(i in 1:length(l2seq)){
            resu<-iPALM(w0s,X1,X2,t1,lambda1=0.68,lambda2=l2seq[i],alpha=alpha,te=te,palmquiet=FALSE,maxtimes=500)
            tmyu<-GetPredict(resu,X3)
            myeu<-MissErr(tmyu,t3)
            if(warminit==TRUE){
              w0s=resu
            }
            l2seq[i]=myeu
          }
          return(list(l1seq,l2seq))
        }else{
          objvalue=NULL
          objvalue=rep(F(w0,X1,X2,t1,lambda1=lambda1,lambda2=lambda2,alpha=alpha),2)
          for(i in 1:10){
            resu1<-PALM(w0,X1,X2,t1,lambda1=lambda1,lambda2=lambda2,alpha=alpha,te=te,palmquiet=FALSE,maxtimes=i*20)
            resu2<-iPALM(w0,X1,X2,t1,lambda1=lambda1,lambda2=1,alpha=alpha,te=te,alpha_1=alpha_1,beta_1=beta_1,palmquiet=FALSE,maxtimes=i*20)
            objvalue=c(objvalue,F(resu1,X1,X2,t1,lambda1=lambda1,lambda2=lambda2,alpha=alpha),
                       F(resu2,X1,X2,t1,lambda1=lambda1,lambda2=lambda2,alpha=alpha))
          }
          return (matrix(objvalue,nrow=2))
        }
      }
    }
  }
  if(n2>0){
    if(useunlabel==FALSE){
      mycoef<-PALM(w0_r=w0,X1_r=X1,X2_r=as.matrix(NA),t1_r=t1,lambda1=lbd1,lambda2=lbd2,alpha=alpha,te=te)
      glmfit<-glmnet(x=X1[,-1],y=as.factor(t1),family="binomial",
                     standardize=FALSE,lambda=lbd1/nrow(X1),alpha=alpha)
      glmcoef<-coef( glmfit)
      print("Beta=")
      print(Beta)
      print("\n")
      print(glmfit)
      print("coeff:")
      return(list(mycoef,"\t",glmcoef,F(mycoef,X1,NA,t1,lambda1=lambda1,lambda2=0,alpha=alpha)))
    }else{
      resu<-PALM(w0,X1,X2,t1,lambda1=lbd1,lambda2=lbd2,alpha=alpha,te=te,palmquiet=FALSE)
      tmyu<-GetPredict(resu,X3)
      myeu<-MissErr(tmyu,t3)
      # s1=proc.time()
      cvmin=matrix(ncol=2,nrow=6)
      foldid=sample(1:5,size=n1,replace=TRUE)
      if(type.measure=="c"){
        typec="class"
      }else{
        typec="deviance"
      }
      for(al in 0:5){
        cvgela=cv.glmnet(x=X1[,-1],y=as.factor(t1),family="binomial",
                         type.measure=typec,foldid=foldid,alpha=al*0.2)
        cvmin[al+1,]=c( min(cvgela$cvm), cvgela$lambda.min)
      }
      that=which(cvmin[,1]==min(cvmin[,1]))
      lcvfit=cv.glmnet(x=X1[,-1],y=as.factor(t1),family="binomial",
                       type.measure=typec,foldid=foldid,alpha=(that-1)*0.2)
      tmyl=as.numeric(predict(lcvfit,newx=X3[,-1],s="lambda.min",type="class"))
      myel<-MissErr(as.vector(tmyl),t3)
      print(F(resu,X1,X2,t1,lambda1,lambda2,alpha))
      if(detail==TRUE){
        print("t3=\n")
        print("l1 norm entropy:\n")
        print("lasso only using labeled:\n ")
        print(matrix(c(t3,tmyu,tmyl),ncol=3,byrow=FALSE))
        print("the alpha by glmnet is")
        print((that-1)*0.2)
        print("the l1 by glmnet is")
        print(lcvfit$lambda.min * n1)
        print("miserror of l1entr,l1,l2entr")
      }
      # print(proc.time()-s1)

      return (c(myeu,myel,svm_light_e,svmlin_l2_e,svmlin_DA_e))

    }
  }

}
