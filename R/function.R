

if(require(MASS)==FALSE){
  install.packages("MASS")
  library(MASS)
}


if(require(glmnet)==FALSE){
  install.packages("glmnet")
  library(glmnet)
}

if(require(plyr)==FALSE){
  install.packages("plyr")
  library(plyr)
}

if(require(inline)==FALSE){
  install.packages("inline")
  library(inline)
}
if(require(Rcpp)==FALSE){
  install.packages("Rcpp")
  library(Rcpp)
}

if(require(RcppArmadillo)==FALSE){
  install.packages("RcppArmadillo")
  library(RcppArmadillo)
}
if(require(devtools)==FALSE){
  install.packages("devtools")
  library(devtools)
}
if(require(crayon)==FALSE){
  install.packages("crayon")
  library(crayon)
}
#.sourceCpp_1_DLLInfo` <- dyn.load('../cs.so')
#PALM <- Rcpp:::sourceCppFunction(function(w0_r, X1_r, X2_r, t1_r, lambda1, lambda2 = 0, alpha = 1, te = 0.0001, maxtimes = 0L, palmquiet = TRUE) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_PALM')
#iPALM <- Rcpp:::sourceCppFunction(function(w0_r, X1_r, X2_r, t1_r, lambda1, lambda2 = 0, alpha = 1, te = 0.0001, maxtimes = 500L, alpha_1 = 0.2, beta_1 = 0.8, palmquiet = TRUE) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_iPALM')

sourceCpp("testarcpp.cpp")

#' Addvec1
#'
#' This function allows you to add constant 1s to the features so the intercept of the model could be absorded.
#' @param X data matrix with nrow the number of instances and ncol the number of features.
#' @keywords add a constant feature 1 for any instances to treat intercept or bias in a uniform way with all other features
#' @export
#' @examples
#' Addvec1(X)
Addvec1<-function(X){
  return (cbind(rep(1,nrow(X)),X))
}

#' Hentr
#'
#' This function allows you to compute entropy value of a bernoullie random variable.
#' @param x a vector of probabilities between 0 and 1.0 or 1 will give 0 output.Values outside the range will generate errors.
#' @keywords entropy
#' @export
#' @examples
#' Hentr(X)
Hentr<-function(x){
  y=-x*log(x)-(1-x)*log(1-x)
  y[is.nan(y)==TRUE]=0
  return (y)
  
}

#' Signoid
#'
#' This function allows you to compute sigmoid function value of a real number vectorwisely.Input of positive infinity will return 1 
#' and negative infinity 0 
#' @param x vector.
#' @keywords sigmoid
#' @export
#' @examples
#' Signoid(x)
Signoid<-function(x){
  return (1/(1+exp(-x)))
}

#' Lexp
#'
#'This function allows you to compute log sigmoid function value of a real number vectorwisely.
#'Input of positive infinity will return 0 
#' @param x vector.
#' @keywords log-sigmoid
#' @export
#' @examples
#' Lexp(x)
Lexp<-function(x){
  return (log(1+exp(-x)))
}

#' F
#'
#' This function allows you to compute the value of the objective function of EELR.
#' @param w coefficient vector
#' @param X1 labeled data matrix(each row representing features for an instance).
#' @param X2 unlabeled data matrix(each row representing features for an instance).NA for supervised learning
#' @param t1 instance labels of labeled data, only taking -1  and 1.
#' @param lambda1 regularization parameter for naive elasticnet penalty.
#' @param lambda2 regularization parameter for unlabled data contribution.
#' @param alpha regularization parameter for tradeoff between l1 and l2 penalty.
#' @keywords Objective function (with penalty) value
#' @export
#' @examples
#' F(w,X1,X2=NA,t1,lambda1=0,lambda2=0,alpha)
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
#' GetPredict
#'
#'This function allows you to computing predicted value by inner product of w coefficients and data matrix.1,-1
#'will be the typical label output.0 can also be the output if the prediction ties between two classes.
#' @param w0.coefficients vector in the model
#' @param Xt.data matrix
#' @keywords computing predicted value from w coefficients and data matrix
#' @export
#' @examples
#' GetPredict(w,Xt)
GetPredict<-function(w,Xt){
  return (sign(sign(Xt%*%w)+.5))
}

#' MissErr
#'
#'This function allows you to compute classification error between predicted labels and known labels(with 0-1 loss).Note that
#'the predicted labels can take 0 and the loss is 0.5.
#' @param tpred prediction label vector(taking 1 and -1).
#' @param t known data label vector.
#' @keywords log-sigmoid
#' @export
#' @examples
#' MissErr(tpred,t)
MissErr<-function(tpred,t){
  if(length(tpred)!=length(t))stop(paste0("lengths are different in MissErr with first",length(tpred),"second",length(t)))
  return (sum(abs(matrix(tpred,ncol=1)-matrix(t,ncol=1))/(2*length(t))) )
}

#' getlambda1seq
#'
#'This function allows you to computing sequence of lambda1 for cross validation by first computing a maximum value and then
#'generate an exponentially arithmetic dropping sequence.
#' @param w0. as in PALM
#' @param X1. as in PALM
#' @param X2. as in PALM 
#' @param t1. as in PALM
#' @param te. as in PALM
#' @param lambda2. 
#' @param sparsete.tolerance for searching maximum lambda1, when the l-1 norm of the parameter 
#' vector is below this threshold or lamdda1 = 10^(6) the maximum lambda1 is reached
#' @param lambda1.min.ratio.ratio of min and max lambda1 searched
#' @param nlambda1.number of lambda1 searched
#' @keywords sequence of lambda1 in the grid search
#' @export
#' @examples
#' getlambda1seq(w0,X1,X2,t1,lambda2,
#' sparsete,lambda1.min.ratio,nlambda1,te=te,wn=FALSE)
getlambda1seq<-function(w0,X1,X2,t1,lambda2,
                        sparsete,lambda1.min.ratio,nlambda1,te=te,wn=FALSE){
  for(i in 4:6){
    tempresult<-PALM(w0,X1,X2,t1,lambda1=10^(i/2),lambda2=lambda2,alpha=1,te=te)
    if(mean(abs(tempresult))<sparsete||i==6){
      lambda1max =10^(i/2)
      break
    }
  }
  if(lambda1max ==10^(i/2)&&wn==TRUE)warning("lamdba1 reaches the maximum!")
  return (exp(seq(from=log(lambda1max),
                  to=log(lambda1max*lambda1.min.ratio),length=nlambda1)))
}

#' lambda2tr
#'
#'This function generates the exponentially rising sequence of lambda2 in the grid search.
#' @param k.transformed to be as.numeric(k!=0)*2^(k/2-3.5)
#' @keywords exponentially arithmetic increasing sequence
#' @export
#' @examples
#' lambda2tr(k)

lambda2tr<-function(k){
  return ( as.numeric(k!=0)*2^(k/2-3.5) )
}

#' writesvmdata2
#'
#'This function writes txt data for software svmlight and svmlin.
#' @param X.data matrix
#' @param t.label
#' @param filename1.name for files written for software svmlight
#' @param filename2.name for files written for software svmlin
#' @param svm_light_dir.directory for svm_light
#' @param svmlin_dir. directory for svm_lin

#' @keywords write txt data for software svmlight and svmlin
#' @export
#' @examples
#' writesvmdata2(X,t,filename1,filename2,svm_light_dir,svmlin_dir)



#filename1 for light filename2 for lin




#filename1 for light filename2 for lin
writesvmdata2<-function(X,t,filename1,filename2,svm_light_dir,svmlin_dir){
  twd=getwd()
  if(NA%in%X | nrow(X) == 0 | ncol(X) == 0)stop("invalid input for makesvmdata")
  
  nr = nrow(X)
  nc = ncol(X)
  #t3=proc.time()
  outmatrix = matrix(nrow=nr,ncol=4*nc)
  outmatrix[,seq(2,4*nc,by=4)]=':'
  outmatrix[,seq(4,4*nc,by=4)]=' '
  outmatrix[,seq(1,4*nc,by=4)]=matrix(1:nc,nrow=nr,ncol=nc,byrow=TRUE)
  outmatrix[,seq(3,4*nc,by=4)]=X
  #t4=t4 + proc.time()-t3
  if(length(t)!=nr){
    print(nr)
    print("\n")
    print(length(t))
    print("\n")
    print(X)
    print("\n")
    print(t)
    stop("length of t must be equal to nrow(X)!")
  }
  #t3=proc.time()
  setwd(svmlin_dir)
  write.table(outmatrix,file=filename2,sep='',quote=FALSE,row.names=FALSE,col.names=FALSE )
  setwd(twd)
  setwd(svm_light_dir)
  outmatrix=cbind(t,rep(' ',length(t)),outmatrix)
  write.table(outmatrix,file=filename1,sep='',quote=FALSE,row.names=FALSE,col.names=FALSE )
  setwd(twd)
  #t5=t5 + proc.time()-t3
}

#' writetxt2
#'
#'This function generate and then write data for software svmlight and svmlin..
#' @param jobname.the same of in testsepa
#' @param X1. the same as in testsepa.
#' @keywords write txt data for software svmlight and svmlin
#' @export
#' @examples
#' writetxt2(jobname,X1,X2,X3,t1,t3,
#' svm_light_dir,svmlin_dir)




writetxt2<-function(jobname,X1,X2,X3,t1,t3,
                    svm_light_dir,svmlin_dir){
  n1=nrow(X1)
  n2=nrow(X2)
  twd=getwd()
  n3=length(t3)
  if(length(t1)!=n1|length(t3)!=nrow(X3))stop("inconsistent input with X1,t1 or X3,t3 in writetxt")
  
  setwd(svm_light_dir)
  writesvmdata2( rbind(X1[,-1],X2[,-1]) ,filename1=(paste0(jobname,"train.txt")) ,
                 filename2=paste0(jobname,"training_examples") ,t=c(t1, rep(0,n2)),
                 svm_light_dir=svm_light_dir,svmlin_dir=svmlin_dir )
  writesvmdata2( X3[,-1], filename1=paste0(jobname,"test.txt"), filename2=paste0(jobname,"test_examples"), t=rep(0,n3),
                 svm_light_dir=svm_light_dir,svmlin_dir=svmlin_dir )
  writesvmdata2( X2[,-1], filename1=paste0(jobname,"trainunlabel.txt"), 
                 filename2=paste0(jobname,"trainingunlabel_examples"), t=rep(0,n2),
                 svm_light_dir=svm_light_dir,svmlin_dir=svmlin_dir )
  setwd(svmlin_dir)
  write.table(matrix( c(t1,rep(0,n2)),ncol=1   ), file=paste0(jobname,"training_labels"),row.names=FALSE,col.names=FALSE )
  setwd(twd)
}

#' svmtran
#'
#'This function efficiently utilize saved results for doing resductive learning with svm software.
#' @param jobname.See parameters of svmtest
#' @param t3. label of test data set
#' @keywords transductive learning
#' @export
#' @examples
#' svmtran(jobname=jobname,t3,svmlight,svm_light_dir,svmlin,
#' svmlin.mode,svmlin_dir,
#' twd)
svmtran<-function(jobname=jobname,t3,svmlight,svm_light_dir,svmlin,
                  svmlin.mode,svmlin_dir,
                  twd){
  svm_light_e<-NA
  svmlin_l2_e<-NA
  svmlin_DA_e<-NA
  
  
  
  if(svmlight==TRUE){
    
    setwd(twd)
    setwd(svm_light_dir)
    #system(paste0("cd ",svm_light_dir))
    #system(paste0( "./svm_learn -b 1 -c ", svmlightC, " ",jobname,"train.txt ",jobname,"model.txt")  )
    # no biase since it is included in the data
    system(paste0("./svm_classify ",jobname,"trainunlabel.txt ",jobname,"model.txt ",jobname,"unlabelprediction.txt"))
    svm_light_pred<-read.table(file=paste0(jobname,"unlabelprediction.txt"))
    # print(svm_light_pred)
    svm_light_e<-MissErr(sign(svm_light_pred[,1]),t3)
    print("svmlight miser:")
    print(svm_light_e)
    
  }
  if(svmlin==TRUE){
    setwd(twd)
    setwd(svmlin_dir)###note that the training_examples contain labeled data soin transductive learning modethe new training
    ###_example file is required to be created to contain only the unlabeled data.
    if("l2hatloss"%in%svmlin.mode){
      #system(paste0("./svmlin -A ",2," -W 1 -U 1 -R 0.5 training_examples training_labels"))
      #system(paste0("./svmlin -A ",2," -W ",svmlinWv[1]/n1," -U ",svmlinUv[1]*n2/n1," -R 0.5 ",jobname,"training_examples ",jobname,"training_labels"))
      system(paste0("./svmlin -f  ",jobname,"training_examples_ms.weights  ",jobname,"trainingunlabel_examples "))
      svmlin_pred<-read.table(file=paste0(jobname,"trainingunlabel_examples.outputs"))
      svmlin_l2_e<-MissErr(sign(svmlin_pred[,1]),t3)
      print("svmlin l2hatloss miser:")
      print(svmlin_l2_e)
      
      
    }
    if("DA"%in%svmlin.mode){
      #system(paste0("./svmlin -A ",3," -W 1 -U 1 -R 0.5 training_examples training_labels"))
      #system(paste0("./svmlin -A ",3," -W ",svmlinWv[2]/n1," -U ",svmlinUv[2]*n2/n1," -R 0.5 ",jobname,"training_examples ",jobname,"training_labels"))
      system(paste0("./svmlin -f  ",jobname,"training_examples_da.weights  ",jobname,"trainingunlabel_examples "))
      svmlin_pred<-read.table(file=paste0(jobname,"trainingunlabel_examples.outputs"))
      svmlin_DA_e<-MissErr(sign(svmlin_pred[,1]),t3)
      print("svmlin DA miser:")
      print(svmlin_DA_e)
      
    }
  }
  
  
  setwd(twd)
  return ( c(svm_light_e,svmlin_l2_e,svmlin_DA_e) )
}

#' svmtest
#'
#'This function does classification with svmlight and svmlin.
#' @param jobname.the same as parameter jobname for testsepa
#' @param X1. X1 and the remaining parameters are the same as in testsepa
#' @param svmlightC. in svmlight,the parameter for svmlight controlling tradeoff between l2norm and margin
#' @param svmlinWv.default:c(1,1) in svmlin, the parameter controlling trade off between l2norm
#' @param svmlinUv.default:c(1,1)
#' @param svmlin.mode. "l2hatloss" quandratic loss for labeled data with multiple switch for transductive learning,
#' "DA" global optimization method with deterministic annealing
#' @param cpsvmlin. a parameter to improve efficiency
#' @keywords svmlight svmlin
#' @export
#' @examples
#' svmtest(jobname,X1,X2,X3,t1,t3,writeout,svmlight,svm_light_dir,svmlin,
#' svmlin.mode,svmlin_dir,
#' twd,svmlightC,svmlinWv=c(1,1),svmlinUv=c(1,1),cpsvmlin=FALSE)
#' 
#' 
#' 
svmtest<-function(jobname,X1,X2,X3,t1,t3,writeout,svmlight,svm_light_dir,svmlin,
                  svmlin.mode,svmlin_dir,
                  twd,svmlightC,svmlinWv=c(1,1),svmlinUv=c(1,1),cpsvmlin=FALSE){
  
  svm_light_e<-NA
  svmlin_l2_e<-NA
  svmlin_DA_e<-NA
  
  if(writeout==TRUE){
    writetxt2(jobname=jobname,X1,X2,X3,t1,t3,svm_light_dir=svm_light_dir,svmlin_dir=svmlin_dir)
  }
  n1=nrow(X1)
  n2=ncol(X2)
  if(svmlight==TRUE){
    
    setwd(twd)
    setwd(svm_light_dir)
    #system(paste0("cd ",svm_light_dir))
    system(paste0( "./svm_learn -b 1 -c ", 1/svmlightC, " ",jobname,"train.txt ",jobname,"model.txt")  )
    # no biase since it is included in the data
    system(paste0("./svm_classify ",jobname,"test.txt ",jobname,"model.txt ",jobname,"prediction.txt"))
    svm_light_pred<-read.table(file=paste0(jobname,"prediction.txt"))
    # print(svm_light_pred)
    svm_light_e<-MissErr(sign(svm_light_pred[,1]),t3)
    print("svmlight miser:")
    print(svm_light_e)
    
  }
  if(svmlin==TRUE){
    setwd(twd)
    setwd(svmlin_dir)
    if("l2hatloss"%in%svmlin.mode){
      #system(paste0("./svmlin -A ",2," -W 1 -U 1 -R 0.5 training_examples training_labels"))
      system(paste0("./svmlin -A ",2," -W ",svmlinWv[1]/n1," -U ",svmlinUv[1]*n2/n1," -R 0.5 ",jobname,"training_examples ",jobname,"training_labels"))
      system(paste0("./svmlin -f  ",jobname,"training_examples.weights  ",jobname,"test_examples "))
      svmlin_pred<-read.table(file=paste0(jobname,"test_examples.outputs"))
      svmlin_l2_e<-MissErr(sign(svmlin_pred[,1]),t3)
      print("svmlin l2hatloss miser:")
      print(svmlin_l2_e)
      if(cpsvmlin==TRUE){# cp training_example.weights to store coefficients
        system(paste0("cp ",jobname,"training_examples.weights ",jobname,"training_examples_ms.weights"))
      }
      
    }
    if("DA"%in%svmlin.mode){
      #system(paste0("./svmlin -A ",3," -W 1 -U 1 -R 0.5 training_examples training_labels"))
      system(paste0("./svmlin -A ",3," -W ",svmlinWv[2]/n1," -U ",svmlinUv[2]*n2/n1," -R 0.5 ",jobname,"training_examples ",jobname,"training_labels"))
      system(paste0("./svmlin -f  ",jobname,"training_examples.weights  ",jobname,"test_examples "))
      svmlin_pred<-read.table(file=paste0(jobname,"test_examples.outputs"))
      svmlin_DA_e<-MissErr(sign(svmlin_pred[,1]),t3)
      print("svmlin DA miser:")
      print(svmlin_DA_e)
      if(cpsvmlin==TRUE){
        system(paste0("cp ",jobname,"training_examples.weights ",jobname,"training_examples_da.weights"))
      }
    }
  }
  
  
  setwd(twd)
  return ( c(svm_light_e,svmlin_l2_e,svmlin_DA_e) )
}



#' testsepa
#'
#'This function does inductive and transductive learning performance comparison with supervised learning by 
#'glmnet, semi-sueprvsied learnin svmlight and svmlin for simulated dataset
#'or assigned dataset or 3 microarray dataset provided.It returns a list error1 containing the classification 
#'error rates in the first element and the coefficients learned in the next elements.if transductive
#'=TRUE,error1[[1]] contains CER of RSLR(tr),RSLR(in),glmnet(tr),glmnet(in),svmlight(tr),
#'svmlight(in),svmlin_ms(tr),svmlin_ms(in),svmlin_da(tr),svmlin_da(in) where "tr" is short for
#'transductive learning mode and "in" inductive learning mode.If transductive=FALSE,only inductive
#'mode CER are in error1[[1]].error[[2]] ~ error[[7]] contains the coefficients learned by the
#'5 methods in the order mentioned above.the lambda1 parameter is for the elastic net penalty.
#'And the lambda2 parameter is for the entropy regularizer for unlabeleda data.Empirically,
#'lambda2 can be fixed to 1 by setting lambda2=1,and crv = FALSE.
#'If svmlight(http://svmlight.joachims.org/) or svmlin(http://vikas.sindhwani.org/svmlin.html)
#' method is used for making comparisons of methods,please install them in the parent directory
#' of current working directory and assign the path by parameter svm_light_dir and svmlin_dir
#'
#' @param jobname.Name of this simulation job.It must be assigned if import is 0.default:1.
#' @param import. 0 in simulation mode, or 1 if using given data set with dataset
#' =1(for leukemia data set),2(for colon data set) or 3(for prostate data set), or 1 with dataset != 1,2 or 3 if using assigned data
#' @param X1. labeled data matrix
#' @param X2. unlabeled data matrix
#' @param X3. test data matrix
#' @param t1. labeled data set labels
#'@param t2. unlabeled data set labels
#'@param t3. test data set labels
#' @param dim. dimensions for simulation(import=1)
#' @param w0. initial value of the w coefficients, deault zero vector
#'@param standardize. TRUE to scale the features, except for the later added constant feature of 1 in simulation mode
#'@param n1. labeled data size in simulation mode
#'@param n2. unlabeled data size in simulation mode
#'@param n3. test data size in simulation mode
#'@param pnr. ratio of sizes of the two classes generated in simulation mode
#'@param transductive. TRUE to compute both inductive and transductive learning classification error rates;FALSE only to compute 
#'inductive learning classification error rates
#'@param lambda1. assigned value of lambda1 in simulation mode with crv = FALSE
#'@param lambda2. assigned value of lambda2 in simulation mode with crv = FALSE
#'@param alpha. assigned value of alpha in simulation mode with crv = FALSE
#'@param seed. seed for set.seed() to make repreducible results
#'@param te. tolerance for converging criteria ,in simulation experiment 0.001 is typically for cross validation and 0.0001 for
#'training on the whole data set
#'@param sa. parameter mu to control distance of the two classes mean in simulation mode
#'@param psparse. number of effective or nonzero coefficients in w in simulation mode
#'@param crv. TRUE to do 5fold cross validation to select the lambda1,lambda2 and alpha; FALSE to use assigned values
#'@param type.measure. "c" currently only 0-1 loss classification error implemented
#'@param lambda1.min.ratio. parameter for function getlambda1seq
#'@param nlambda1. number of lambda1 searched in a sequence
#'@param lambda1seq. default= NA to generate by the program, or to use assigned sequence of lambda1 to search by cross validation
#'@param writeout. parameter for function svmtest
#'@param svmlight. TRUE to turn svmlight on
#'@param svm_light_dir. path of directory of svm_light, parameter for svmtest and svmtran
#'@param svmlin. TRUE to turn svmlin on
#'@param svmlin_dir. path of directory of svmlin-v1.0, parameter for svmtest and svmtran
#'@param svmlin.mode. parameter of function svmtest and function svmtran
#'@param sparsete. parameter for function getlambda1seq,see help file for getlambda1seq for details
#' @keywords the function for simulation(import=0), test on microarray dataset(import=1 and dataset
#' =1 , 2 or 3) or user assigned data(import = 1 and dataset != 1,2 or 3)
#' @export
#' @examples
#' testsepa(
#' jobname=1,import=0,X=NULL,X1=NULL,X2=NULL,X3=NULL,t1=NULL,t2=NULL,t3=NULL,dim=NULL,w0=NULL,standardize=FALSE,
#' n1,n2,n3,pnr=.5,transductive=FALSE,lambda1,lambda2,alpha,
#' seed=0,te,
#' sa,psparse=.1
#' ,crv=FALSE,type.measure="c",lambda1.min.ratio=.0001,nlambda1=20,lambda1seq=NA,
#' writeout=FALSE,svmlight=FALSE,svm_light_dir,svmlin=FALSE,
#' svmlin.mode=NULL,svmlin_dir,
#' sparsete=.001,palmquiet=TRUE,nlambda2=8,rho=0.4,testbratio=0.5,dataset=NULL)
testsepa<-function(jobname=1,import=0,X1=NULL,X2=NULL,X3=NULL,t1=NULL,t2=NULL,t3=NULL,dim=NULL,w0=NULL,standardize=FALSE,
                   n1,n2,n3,Total,pnr=.5,transductive=FALSE,lambda1,lambda2,alpha,
                   seed=0,te,
                   sa,psparse=20
                   ,crv=FALSE,type.measure="c",lambda1.min.ratio=.0001,nlambda1=20,lambda1seq=NA,
                   writeout=FALSE,svmlight=FALSE,svm_light_dir,svmlin=FALSE,
                   svmlin.mode=NULL,svmlin_dir,
                   sparsete=.001,palmquiet=TRUE,nlambda2=8,rho=0.5,testbratio=0.5,dataset=NULL)
{
  twd=getwd()
  set.seed(seed)
  finalresult=vector("list",14)#[[1]]for cer others for storing coefficients of the methods
  if(import==0){
    print("Doing simulation")
    Total=n1+n2+n3
    Beta<-rep(0,dim)
    #Beta[sample(c(2:(dim+1)),psparse*dim,replace=FALSE)]=c(1+runif(psparse*dim/2),-1-runif(psparse*dim/2))
    # Beta[c(  (1+dim*0.3-psparse/2):(dim*0.3),(dim*0.7+1):(dim*0.7+psparse/2))]=c(4+runif(psparse/2),-4-runif(psparse/2))
    Beta[c(  (1+dim*0.3-psparse/2):(dim*0.3),(dim*0.7+1):(dim*0.7+psparse/2))]=c(rep(4,psparse/2),rep(-4,psparse/2))
    
    finalresult[[7]]=Beta
    
    Edcrm=rho^( abs( matrix(1:dim,nrow=dim,ncol=dim,byrow=TRUE) - matrix(1:dim,nrow=dim,ncol=dim,byrow=FALSE) ))
    X<-mvrnorm(n=n1*pnr,mu=-sa*Beta,Sigma=Edcrm)
    
    X<-rbind(X,mvrnorm(n=n1*(1-pnr),mu=sa*Beta,Sigma=Edcrm))
    if(n2>0){
      X<-rbind(X,mvrnorm(n=n2*pnr+(0.5-testbratio)*n2,mu=-sa*Beta,Sigma=Edcrm))
      X<-rbind(X,mvrnorm(n=n2*(1-pnr)-(0.5-testbratio)*n2,mu=sa*Beta,Sigma=Edcrm))
    }
    if(n3>0){
      X<-rbind(X,mvrnorm(n=n3*pnr,mu=-sa*Beta,Sigma=Edcrm))
      X<-rbind(X,mvrnorm(n=n3*(1-pnr),mu=sa*Beta,Sigma=Edcrm))
    }
    
    X1<-X[1:n1,]
    #t<-2*rbinom(Total,1,1/(1+exp(-X%*%Beta)))-1
    #generate t according to gaussian distribution
    t=c(rep(1,n1*pnr),rep(-1,n1*(1-pnr)),rep(1,n2*pnr+(0.5-testbratio)*n2),
        rep(-1,n2*(1-pnr)-(0.5-testbratio)*n2),
        rep(1,n3*pnr),rep(-1,n3*(1-pnr)))
    t1<-t[1:n1]
    if(n2>0){
      X2<-X[(n1+1):(n1+n2),]
      t2<-t[(n1+1):(n1+n2)]
      X3=X[(n1+n2+1):Total,]
      t3=t[(n1+n2+1):Total]
      
    }
    
  }else{
    if(dataset==1){
      load("leukemia.rda")
      ylen=length(ly.original)
      ss=sample(1:ylen)
      tset=ss[61:ylen]
      
      uset=ss[27:60]
      lset=ss[1:26]
      
      X1=lx.original[lset,]
      t1=ly.original[lset]
      X2=lx.original[uset,]
      t2=ly.original[uset]
      X3=lx.original[tset,]
      t3=ly.original[tset]
      
      t1=2*t1-1##turn 0 to -1
      t3=2*t3-1
      t2=2*t2-1
      
      n1=nrow(X1)
      n2=nrow(X2)
      dim=ncol(X1)
    }
    if(dataset==2){
      load("colon.rda")
      ylen=length(colon.y)
      ss=sample(1:ylen)
      tset=ss[53:ylen]
      
      uset=ss[27:52]
      lset=ss[1:26]
      
      X1=colon.x[lset,]
      t1=colon.y[lset]
      X2=colon.x[uset,]
      t2=colon.y[uset]
      X3=colon.x[tset,]
      t3=colon.y[tset]
      
      
      t1=2*t1-1
      t3=2*t3-1
      t2=2*t2-1
      
      n1=nrow(X1)
      n2=nrow(X2)
      dim=ncol(X1)
    }
    if(dataset==3){
      load("prostate.rda")
      
      ylen=length(prostate.y)
      ss=sample(1:ylen)
      tset=ss[79:ylen]
      
      uset=ss[27:78]
      lset=ss[1:26]
      
      X1=prostate.x[lset,]
      t1=prostate.y[lset]
      X2=prostate.x[uset,]
      t2=prostate.y[uset]
      X3=prostate.x[tset,]
      t3=prostate.y[tset]
      
      t1=2*t1-1
      t3=2*t3-1
      t2=2*t2-1
      
      n1=nrow(X1)
      n2=nrow(X2)
      dim=ncol(X1)
    }

    
    
  }
  if(standardize==TRUE){######this scale is before adding the w0
    X1=scale(X1)
    X2=scale(X2)
    X3=scale(X3)
    
  }
  # finalresult[[8]]=X1
  # finalresult[[9]]=X2
  # finalresult[[10]]=X3
  # finalresult[[11]]=t1
  # finalresult[[12]]=t2
  # finalresult[[13]]=t3
  
  X1=Addvec1(X1)##adding const for w0
  X2=Addvec1(X2)
  X3=Addvec1(X3)
  
  
  if( FALSE==(setequal ( unique(t1) ,c(1,-1)  )  ) ){
    stop("the category indicator must only takes 1 or -1!")
  }
  
  #initialize the parameters to be computed
  if(is.null(w0)){
    w0<-rep(0,dim+1)
    if(length(w0)!=ncol(X1)){print("dimension error!")}
  }else{
    w0=w0
  }
  
  setwd(twd)
  
  
  lbd1<-NULL
  #cv or not cv
  if(crv==FALSE){
    lbd1=lambda1
    lbd2=lambda2
    alpha=alpha
    svmlightC=1
    svmlinUl2=1
    svmlinWl2=1
    svmlinUDA=1
    svmlinWDA=1
  }else{
    
    
    if(is.na(lambda1seq)){
      lambda1seq=getlambda1seq(w0=w0,X1=X1,X2=X2,t1=t1,lambda2=1,sparsete=sparsete,
                               lambda1.min.ratio=lambda1.min.ratio,nlambda1=nlambda1,te=te)
    }
    in1<-sample(1:n1)
    in2<-sample(1:n2)
    l1<-list(0,0,0,0,0)
    l2<-list(0,0,0,0,0)
    for(i in 1:5){
      l1[[i]]<-in1[floor((i-1)*n1/5+1):floor(i*n1/5)]
      l2[[i]]<-in2[floor((i-1)*n2/5+1):floor(i*n2/5)]
    }
    #writing the txtfiles in preparation for training in cross validation on svms
    jobnamecv = paste(jobname,"cv",1:5,sep="")#not paste0 here! note for parameter collapse
    for(i in 1:5){
      writetxt2(jobname=jobnamecv[i],X1=X1[-l1[[i]],],X2=X2[-l2[[i]],],X3=X1[l1[[i]],],t1=t1[-l1[[i]]],t3=t1[l1[[i]]],
                svm_light_dir,svmlin_dir)
      
      
    }
    
    Minl1<-matrix(nrow=nlambda2,ncol=3) #to store each mincv error for lambda2(col1),col2=best lambda1 for the lambda2 and col3=best al
    Svmlinl2mat<-matrix(nrow=nlambda2,ncol=length(lambda1seq))
    SvmlinDAmat<-matrix(nrow=nlambda2,ncol=length(lambda1seq))
    Svmlightvec<-rep(NA,nlambda2)# svmlight with c in svmlight equal to our lambda2
    
    for(k in 1:nlambda2 ){
      
      #######svmlight cv start
      Svmval=rep(0,3)###store 3 svmtest error in one cross validation
      for(i in 1:5){
        
        Svmval=Svmval+svmtest(jobname=jobnamecv[i],X1=X1[-l1[[i]],],X2=X2[-l2[[i]],],X3=X1[l1[[i]],],t1=t1[-l1[[i]]],t3=t1[l1[[i]]],writeout=FALSE,
                              svmlight=svmlight,svm_light_dir=svm_light_dir,svmlin=FALSE,
                              svmlin.mode=svmlin.mode,svmlin_dir=svmlin_dir,
                              twd=twd,svmlightC=lambda2tr(k),svmlinWv=rep(NA,2),svmlinUv=rep(NA,2)  )
        
      }
      Svmlightvec[k]=Svmval[1]/5   #k for hyperparameter lambda2 j for lambda1
      
      
      #######svmlight cv end
      
      
      
      #Fseq<-matrix(ncol=6,nrow=length(lambda1seq))
      Fseq<-matrix(ncol=11,nrow=length(lambda1seq))  #store crossvalidation error for each lambda1 and al in 0:10, with lambda2 fixed
      for(j in 1:length(lambda1seq)){
        
        ######svmlin cv
        Svmval=rep(0,3)###store 3 svmtest error in one cross validation
        for(i in 1:5){
          
          Svmval=Svmval+svmtest(jobname=jobnamecv[i],X1=X1[-l1[[i]],],X2=X2[-l2[[i]],],X3=X1[l1[[i]],],t1=t1[-l1[[i]]],t3=t1[l1[[i]]],writeout=FALSE,
                                svmlight=FALSE,svm_light_dir=svm_light_dir,svmlin=svmlin,
                                svmlin.mode=svmlin.mode,svmlin_dir=svmlin_dir,
                                
                                twd=twd,svmlightC=NULL,svmlinWv=rep(lambda1seq[j],2),svmlinUv=rep(lambda2tr(k),2) )
          
        }
        Svmlinl2mat[k,j]=Svmval[2]/5   #k for hyperparameter lambda2 j for lambda1
        SvmlinDAmat[k,j]=Svmval[3]/5
        
        ######svmlin cv
        
        
        for(al in 0:10){
          Fval<-0
          
          
          
          
          for(i in 1:5){
            #resu<-PALM(w0,X1[-l1[[i]],],X2[-l2[[i]],],t1[-l1[[i]]],lambda1=lambda1seq[j],lambda2=0.25*k,alpha=al*0.1,te=te)
            resu<-PALM(w0,X1[-l1[[i]],],X2[-l2[[i]],],t1[-l1[[i]]],lambda1=lambda1seq[j],lambda2=lambda2tr(k),alpha=al*0.1,te=te,palmquiet=palmquiet)
            if(type.measure=="c"){
              Fval=Fval+MissErr(GetPredict(resu,X1[l1[[i]],]),t1[l1[[i]]])
            }else{
              tempProb=F(resu,X1[l1[[i]],],NA,t1[l1[[i]]],lambda1=0,lambda2=0,alpha=1)
              Fval=Fval+tempProb
            }
            
          }#for i
          #print((al+1+(j-1)*6 + length(lambda1seq)*6*(k)  )/(6*length(lambda1seq)*(nlambda2+1)))
          print((al+1+(j-1)*11 + length(lambda1seq)*11*(k)  )/(11*length(lambda1seq)*(nlambda2+1)))
          Fseq[j,al+1]<-Fval/5#! al+1 not al !
        }#for al
      }#for j
      that=which(Fseq==min(Fseq),arr.ind=TRUE) #Fseq is :Fseq[i,j] cv error when lbd1=lambda1seq[i] and alpha =al[j]
      Minl1[k,]=c(min(Fseq),lambda1seq[ that[1,1] ],that[1,2])
      #the cv result for fixed lbd2 as k,col2 is the corresponding lamdba1,col3 is the corresponding al
      
    }# for k
    
    lbd2=lambda2tr( which(Minl1[,1]==min(Minl1[,1]))[1] )
    lbd1=Minl1[which(Minl1[,1]==min(Minl1[,1]))[1],2]
    alpha=Minl1[which(Minl1[,1]==min(Minl1[,1]))[1],3]*0.1-0.1#col3 in Minl1 is only the index so subtracted by 0.1
    
    lightloc=which( Svmlightvec==min(Svmlightvec) )[1]
    svmlightC=lambda2tr(lightloc)######svmlight hyperparameter
    
    linl2loc=which( Svmlinl2mat==min(Svmlinl2mat) ,arr.ind=TRUE )
    linDAloc=which( SvmlinDAmat==min(SvmlinDAmat) ,arr.ind=TRUE )
    svmlinUl2=lambda2tr(linl2loc[1,1])######svmlin hyperparameter
    svmlinWl2=lambda1seq[linl2loc[1,2]]
    svmlinUDA=lambda2tr(linDAloc[1,1])
    svmlinWDA=lambda1seq[linDAloc[1,2]]
    
    print("the lbd1,lbd2 and alpha by cv is")
    print(c(lbd1,lbd2,alpha))
    #print(Minl1)#########
    #print(Fseq)#########
    
    
  }
  
  
  
  #eelr predict
  resu<-PALM(w0,X1,X2,t1,lambda1=lbd1,lambda2=lbd2,alpha=alpha,te=0.1*te,palmquiet=FALSE)
  tmyu<-GetPredict(resu,X3)
  myeu<-MissErr(tmyu,t3)
  #glmnet predict supervised
  # s1=proc.time()
  cvmin=matrix(ncol=2,nrow=11)
  if(crv==TRUE){
    foldid=rep(1,n1)
    for(i in 2:5){
      foldid[l1[[i]]]=i
    }
  }else{
    foldid=c( sample(1:5,size=5,replace=FALSE),sample(1:5,size=n1-5,replace=TRUE))
  }
  
  if(type.measure=="c"){
    typec="class"
  }else{
    typec="deviance"
  }
  for(al in 0:10){
    cvgela=cv.glmnet(x=X1[,-1],y=as.factor(t1),family="binomial",
                     type.measure=typec,foldid=foldid,alpha=al*0.1)
    cvmin[al+1,]=c( min(cvgela$cvm), cvgela$lambda.min)
  }
  that=which(cvmin[,1]==min(cvmin[,1]))[1]
  lcvfit=cv.glmnet(x=X1[,-1],y=as.factor(t1),family="binomial",
                   type.measure=typec,foldid=foldid,alpha=(that-1)*0.1)
  tmyl=as.numeric(predict(lcvfit,newx=X3[,-1],s="lambda.min",type="class"))
  myel<-MissErr(as.vector(tmyl),t3)
  #svms predict
  
  svm_inductive<-NA
  svm_transductive<-NA
  svm_inductive=svmtest(jobname=jobname,X1=X1,X2=X2,X3=X3,t1=t1,t3=t3,writeout=writeout,svmlight=svmlight,svm_light_dir=svm_light_dir,svmlin=svmlin,
                        svmlin.mode=svmlin.mode,svmlin_dir=svmlin_dir,
                        twd=twd,svmlightC=svmlightC,
                        svmlinWv=c(svmlinWl2,svmlinWDA),svmlinUv=c(svmlinUl2,svmlinUDA),cpsvmlin=TRUE)
  
  
  ##to store coefficients of the methods
  finalresult[[2]]=resu
  finalresult[[3]]=as.vector(coef(lcvfit,s="lambda.min"))
  system(paste0("perl ",svm_light_dir,"svm2weight.pl ",svm_light_dir,jobname,"model.txt > ",svm_light_dir,jobname,"svmlightw8s.out"))
  finalresult[[4]]=read.table(paste0(svm_light_dir,jobname,"svmlightw8s.out"))
  finalresult[[5]]= read.table(paste0(svmlin_dir,jobname,"training_examples_ms.weights"))
  #in the final inductive svmtest()the cpsvmlin has been set to be TRUE
  finalresult[[6]]= read.table(paste0(svmlin_dir,jobname,"training_examples_da.weights"))
  
  #handling transductive learning
  
  
  
  if(transductive==TRUE){
    #transductive version of our methods
    tmyu_tran<-GetPredict(resu,X2)
    myeu_tran<-MissErr(tmyu_tran,t2)
    #transductive version of svm
    svm_transductive=svmtran(jobname=jobname,t2,svmlight,svm_light_dir,svmlin,
                             svmlin.mode,svmlin_dir,
                             twd)
    # transductive prediction of glmnet
    tmyl_tran=as.numeric(predict(lcvfit,newx=X2[,-1],s="lambda.min",type="class"))
    myel_tran<-MissErr(as.vector(tmyl_tran),t2)
    
  }
  
  
  
  
  
  
  if(transductive==FALSE){
    finalresult[[1]]=c(myeu,myel,svm_inductive)
    
  }else{
    finalresult[[1]]=c(myeu,myeu_tran,myel,myel_tran,svm_inductive,svm_transductive)
  }
  
  return (finalresult)
  
  
}    

#'runjob
#'
#'this funciont conducts experiment and saves the result in current working directory with name paste0(
#'"error",jobname,".RData"), the result is a list named error1 containing element as each 
#'round of experiment with different seed numbers.It is assumed that svm_light and svmlin
#'directories are in the parent directory of the current working directory
#'@param startseed the seed for the first experiment, by default startseed=endseed and only
#'one round of experiment is conducted
#'@param endseed the seed for the last experiment
#'@param jobname the name in the middle of the file name to save the result. default:1
#'@param l_size the size of labeled dataset in simulation mode. default:20
#'@param dim the dimension of the (labeled) dataset in simulation mode. default:300
#'@param testsize the size of the test dataset in simulation mode.default:1200
#'@param sa the parameter for the distance between two gaussian means in simulation mode. default:0.6
#'@param psparse dimensions of effective features in simulation mode
#'@param testbratio proportion of class 1 in the whole unlabeled data set
#'@param u_size size of unlabeled dataset in simulation mode
#'@param rho parameter for covariance matrix, see the paper
#'@param import 0 simulation mode;1 use real data with "dataset" must be assigned
#'@param dataset if "import" parameter is 1,dataset must take 1,2 or 3 indicating the 
#'"leukemia","colon",or "prostate" microarray data set.The data has been preprocessed.
#' @keywords run jobs of simulation experiment or on microarray dataset
#' @export
#' @examples
runjob<-function(startseed=100,
              endseed=100,
              jobname=1,
              l_size=20,
              dim=300,
              testsize=1200,
              sa=0.6,
              psparse=20,
              testbratio=0.5,
              u_size=100,
              rho=0.5,
              import=0,
              dataset = 1){

error1=list()


for(s in startseed:endseed){
  cerr=testsepa(jobname = jobname, import=import,dim=dim,n1=l_size,n2=u_size,n3=testsize,Total=440,lambda1=2,lambda2=1.2,alpha=0.8,standardize=TRUE,
                seed=s,te=10^(-3),
                sa=l_sa,psparse=l_psparse ,
                crv=TRUE,type.measure="c",lambda1.min.ratio=.0001,nlambda1=9,
                sparsete=.02,
                writeout=TRUE,svmlight=TRUE,svm_light_dir="../svm_light/",svmlin=TRUE,svmlin_dir="../svmlin-v1.0/",svmlin.mode=c("DA","l2hatloss"),
                nlambda2=9,palmquiet=TRUE,,rho=rho,transductive=TRUE,testbratio=testbratio,dataset=dataset)
  error1[[length(error1) + 1]] = cerr
}

save(error1,file=paste0("error",jobname,"i.RData"))
}

