
Update=function(m,y,X,Psi,thetas,alpha,Design_List,Distance_List,Covtype){
  # Design_List=Updates$Design_List;Distance_List=Updates$Distance_List
  # This function updates model parameters based on the new Psi, y, and
  # parameter values

  # for(k in 1:length(Design_List)){
  #   Design_List[[k]]=Psi[,k]
  #   Distance_List[[k]]=as.matrix(dist(Design_List[[k]]),method="manhattan")
  # }

  #-----------------------------------------------------------------------------
  Corr_List=list()
  for(i in 1:(length(Distance_List))){
    Corr_List[[i]]=corr_compute(thetas[i],Distance_List[[i]],Covtype=Covtype)
  }
  Corr=alpha*Reduce("*", Corr_List)+(1-alpha)*diag(rep(1,ncol(Distance_List[[1]])))
  outisigbet=sigma2_compute(Corr=Corr,y=y,X=X,m=m)
  v=outisigbet$sigma2_v
  beta=outisigbet$beta
  InvCholCorr=solve(outisigbet$CholCorr)
  sigma2=alpha*v
  #-----------------------------------------------------------------------------
  return(list("v"=v,"beta"=beta,"sigma2"=sigma2,"Design_List"=Design_List,"Distance_List"=Distance_List,"InvCholCorr"=InvCholCorr))
}


ReInterpolation=function(m,X,y,Psi,thetas,alpha,beta,InvCholCorr,Design_List,Distance_List,Noise_Type,Covtype){
  # This function computes 1) the data for the re-interpolating model
  #                        2) the parameters for the re-interpolating model

  yhat=apply(Psi,1,FUN=posterior_mean,Y=y,XX=X,thetas=thetas,alpha=alpha,beta=beta,InvCholCorr=InvCholCorr,Design_List=Design_List,Noise_Type=Noise_Type,Covtype=Covtype)
  #-----------------------------------------------------------------------------
  Corr_ri_List=list()
  for(i in 1:(length(Distance_List))){
    Corr_ri_List[[i]]=corr_compute(thetas[i],Distance_List[[i]],Covtype=Covtype)
  }
  Corr_ri=Reduce("*", Corr_ri_List)
  #+diag(rep(-0.00005,25))


  outisigbet_ri=sigma2_compute(as.matrix(Corr_ri),y=as.vector(yhat),X=X,m=m)
  sigma2_ri=outisigbet_ri$sigma2_v
  beta_ri=outisigbet_ri$beta
  InvCholCorr_ri=solve(outisigbet_ri$CholCorr)
  Y_max_ri=max(yhat)


  index_max=which(yhat==Y_max_ri)
  x_max_ri=Psi[index_max,]
  if(length(index_max)>1)x_max_ri=x_max_ri[1,]

  return(list("yhat"=yhat
              ,"sigma2_ri"=sigma2_ri
              ,"beta_ri"=beta_ri
              ,"InvCholCorr_ri"=InvCholCorr_ri
              ,"Y_max_ri"=Y_max_ri
              ,"x_max_ri"=x_max_ri))
}

#' Obtaining posterior mean of Gaussian process.
#' @description function fits the posterior mean.
#' @param GP_Object an object returned from either DesignFit or SequenceFit
#' @export
PostMean=function(x,GP_Object){

X=model.matrix(~1,data=data.frame(GP_Object$Updates$Design_List[[1]]))
Y=GP_Object$y
thetas=GP_Object$thetas
Design_List=GP_Object$Updates$Design_List
alpha=GP_Object$Updates$sigma2/GP_Object$Updates$v
beta=GP_Object$Updates$beta
InvCholCorr=GP_Object$Updates$InvCholCorr
Covtype=GP_Object$GP_Part$Covtype

  Design_New_List=list()
  for(k in 1:length(Design_List)){
    Design_New_List[[k]]=x[k]
  }
  diagy=offdiag_compute(Design_List=Design_List,
                        Design_New_List=Design_New_List,
                        thetas=thetas,
                        Covtype=Covtype)

  corr_off_diag=diagy$corr_off_diag
  mu_new=t(diagy$x_new)%*%beta+alpha*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)

  return(mu_new)
}

#' A Gaussian Process Functions
#' @description This function fits Gaussian Process parameters to an initial set of design points.
#' @param Data dataset used for analysis
#' @param Treat_M_List list of treatment models
#' @param Outcome_M_List list of outcome models
#' @param Outcome_Var outcome variable
#' @param Treat_vars treatment variables
#' @param PatID patient identified
#' @param G_List decision rule list
#' @param Psi points to explore in the parameter space
#' @param DR boolean indicating whether double robust estimator should be used
#' @param Bayes boolean indicating whether Bayesian methods will be used
#' @param Bayes_Seed seed to control random number generation
#' @param Numbr_Samp number of random starts in the maximum likelihood
#' @param IthetasU upper limits in the bounded BFGS algorithm
#' @param IthetasL lower limits in the bounded BFGS algorithm
#' @param IalphaU  lower limits in the bounded BFGS algorithm
#' @param IalphaL  lower limits in the bounded BFGS algorithm
#' @param Likelihood_Limits limits for where to plot the parameters in the Gaussian Process Likelihood
#' @param Covtype 1 is for Matern 3/2; 2 is for Matern 5/2
#' @param Prior_List thsia
#' @param Prior_Der_List those
#' @return This function returns the input parameters required for further fitting of the Gaussian Process when more samples are added
#'
#'         as well as the current values of the parameters
#'
#' This function allows you to fit a Bayesian MSM
#' @export
DesignFit=function(PatID
                  ,Data
                  ,Treat_M_List
                  ,Outcome_M_List=NULL
                  ,Outcome_Var
                  ,Treat_Vars
                  ,G_List
                  ,Psi
                  ,Normalized=TRUE
                  ,DR=FALSE
                  ,...
                  #,Bayes=FALSE
                  #,Bayes_Seed=1


                  ,Covtype
                  ,Numbr_Samp
                  ,IthetasU
                  ,IthetasL
                  ,IalphaU
                  ,IalphaL
                  ,Likelihood_Limits=NA
                  ,Prior_List=NULL
                  ,Prior_Der_List=NULL
                  ){


  if(!hasArg(Bayes)){
    Bayes=FALSE
    Bayes_Seed=1
  }else{
    Bayes=TRUE
    Bayes_Seed=loop_index_global
  }
  # print(Bayes)

  start_num=nrow(Psi)
  m=start_num
  Noise_Type="Homosk"
  type_simul=2

  if(!is.list(Prior_List)){
    Prior_List=list()
    Prior_Der_List=list()
    for(i in 1:ncol(Psi)){
      Prior_List[[i]]="0"
      Prior_Der_List[[i]]="0"
    }
  }
  y=c(BayesMSM(PatID=PatID
               ,Outcome_Var=Outcome_Var
               ,Treat_Vars=Treat_Vars
               ,Treat_M_List=Treat_M_List
               ,Outcome_M_List=Outcome_M_List
               ,Data=BayesDat
               ,G_List=G_List
               ,Psi=Psi
               ,DR=DR
               ,Normalized=Normalized
               ,Bayes=Bayes
               ,B=1
               ,Bayes_Seed=Bayes_Seed))

  Design_List=list()
  Distance_List=list()
  for(k in 1:length(G_List)){
    Design_List[[k]]=Psi[,k]
    Distance_List[[k]]=as.matrix(dist(Design_List[[k]]),method="manhattan")
  }

  X=model.matrix(~1,data=data.frame(Design_List[[1]]))


  for(i in 1:ncol(Psi)){
    U=runif(Numbr_Samp,0.01,IthetasU[i])
    if(i==1){DD=U
    }else{DD=cbind(DD,U)}
  }
  DD=cbind(DD,runif(Numbr_Samp,0.01,0.99))
  pari=multi_eval(DD=DD
                  ,y=y
                  ,XX=X
                  ,IthetasL=IthetasL
                  ,IthetasU=IthetasU
                  ,IalphaL=IalphaL
                  ,IalphaU=IalphaU
                  ,Noise_Type=Noise_Type
                  ,Distance_List=Distance_List
                  ,Prior_List=Prior_List
                  ,Prior_Der_List=Prior_Der_List
                  ,m=m
                  ,Covtype=Covtype)

  thetas=pari[1:ncol(Psi)]
  alpha=pari[ncol(Psi)+1]

  if(Bayes==FALSE & sum(is.na(Likelihood_Limits))==0){
    for(k in 1:length(thetas)){
      plot_lik(
        thetas=thetas,
        alpha=alpha,
        numb_varv=k,
        abliny=thetas[k],
        xvals=Likelihood_Limits[[k]],
        counter=0,
        Distance_List=Distance_List,
        y=y,
        X=X,
        Prior_List=Prior_List,
        Prior_Der_List=Prior_Der_List,
        Noise_Type=Noise_Type,
        m=m,
        Covtype=Covtype)

    }
  }

  Updates=Update(m=m,y=y,X=X,Psi=Psi,thetas=thetas,alpha=alpha,Design_List=Design_List,Distance_List=Distance_List,Covtype=Covtype)
  #plotlik

  #######################################
  # reinterpolation parameter calculation
  #######################################

  ReInter=ReInterpolation(m=m
                          ,X=X
                          ,y=y
                          ,Psi=Psi
                          ,thetas=thetas
                          ,alpha=alpha
                          ,beta=Updates$beta
                          ,InvCholCorr=Updates$InvCholCorr
                          ,Design_List=Updates$Design_List
                          ,Distance_List=Updates$Distance_List
                          ,Noise_Type=Noise_Type
                          ,Covtype=Covtype)


  return(list(
    "BayesMSM"=list("Data"=Data
                    ,"Treat_M_List"=Treat_M_List
                    ,"Outcome_M_List"=Outcome_M_List
                    ,"Outcome_Var"=Outcome_Var
                    ,"Treat_Vars"=Treat_Vars
                    ,"PatID"=PatID
                    ,"G_List"=G_List
                    ,"Psi"=Psi
                    ,"Normalized"=Normalized
                    ,"DR"=DR
                    ,"Bayes"=Bayes
    )
    ,"GP_Part"=list("Numbr_Samp"=Numbr_Samp
              ,"IthetasU"=IthetasU
              ,"IthetasL"=IthetasL
              ,"IalphaU"=IalphaU
              ,"IalphaL"=IalphaL
              ,"Covtype"=Covtype
              ,"Prior_List"=Prior_List
              ,"Prior_Der_List"=Prior_Der_List
             )
    ,"Updates"=Updates
    ,"ReInter"=ReInter
    ,"Bayes_Seed"=Bayes_Seed
    ,"thetas"=thetas
    ,"alpha"=alpha
    ,"y"=y
    ,"Likelihood_Limits"=Likelihood_Limits
    ))
}


#' Fitting Gaussian Process parameters for sequential samples
#' @description This function fits Gaussian Process parameters to an initial set of design points in addition to those sampled via the Expected
#'              improvement criterion
#' @param Previous_Fit an object that has resulted from DesignFit or SequenceFit
#' @param Additional_Samp Number of points to sample sequentially
#' @param Bayes_Seed allows user to set seed
#' @param Control_Genoud
#'
#' @return This function returns the input parameters required for further fitting of the Gaussian Process when more samples are added
#'         as well as the current values of the parameters
#'
#' @export
SequenceFit=function(Previous_Fit,Additional_Samp,Control_Genoud=list(),Psi_new,N){
                                                                        #these last two parameters need not be used by the user N and Psi_New
  options(warn=1)

  Data=Previous_Fit$BayesMSM$Data
  Treat_M_List=Previous_Fit$BayesMSM$Treat_M_List
  Outcome_M_List=Previous_Fit$BayesMSM$Outcome_M_List
  Outcome_Var=Previous_Fit$BayesMSM$Outcome_Var
  Treat_Vars=Previous_Fit$BayesMSM$Treat_Vars
  PatID=Previous_Fit$BayesMSM$PatID
  G_List=Previous_Fit$BayesMSM$G_List
  Psi=Previous_Fit$BayesMSM$Psi
  DR=Previous_Fit$BayesMSM$DR
  Normalized=Previous_Fit$BayesMSM$Normalized
  Bayes=Previous_Fit$BayesMSM$Bayes
  Bayes_Seed=Previous_Fit$Bayes_Seed




  Numbr_Samp=Previous_Fit$GP_Part$Numbr_Samp
  IthetasU=Previous_Fit$GP_Part$IthetasU
  IthetasL=Previous_Fit$GP_Part$IthetasL
  IalphaU=Previous_Fit$GP_Part$IalphaU
  IalphaL=Previous_Fit$GP_Part$IalphaL
  Covtype=Previous_Fit$GP_Part$Covtype
  Prior_List=Previous_Fit$GP_Part$Prior_List
  Prior_Der_List=Previous_Fit$GP_Part$Prior_Der_List


  Updates=Previous_Fit$Updates
  ReInter=Previous_Fit$ReInter


  start_num=nrow(Psi)
  m=start_num
  Noise_Type="Homosk"
  type_simul=2

  thetas=Previous_Fit$thetas
  alpha=Previous_Fit$alpha
  y=Previous_Fit$y
  Likelihood_Limits=Previous_Fit$Likelihood_Limits
  X=model.matrix(~1,data=data.frame(Updates$Design_List[[1]]))
  EI_hist=vector()
  Posterior_Draws=NA

  if(Bayes==TRUE){Design_New_List=list(Psi_new[,1],Psi_new[,2])}

  for(counter in 1:Additional_Samp){
        print(paste0("Sequential Sample Number ", counter))
        m=m+1
        if(Bayes==TRUE){
          set.seed((Bayes_Seed-1)*Additional_Samp+counter)
        }

        #setting up parameters for genoud function
        Control_Genoud_Internal=list(nvars=length(G_List)
                                     ,pop.size=500
                                     ,BFGS=TRUE
                                     ,hard.generation.limit=FALSE
                                     ,max=TRUE
                                     ,boundary.enforcement=2
                                     ,print.level=0)
        Control_Genoud=c(Control_Genoud,Control_Genoud_Internal)
        keep=!duplicated(names(Control_Genoud))
        Control_Genoud=Control_Genoud[keep]

        Sub_Control_Genoud=list(fn=EI
                                ,ymax=ReInter$Y_max_ri
                                ,Y=ReInter$yhat
                                ,X=X
                                ,Design_List=Updates$Design_List
                                ,thetas=thetas
                                ,sigma2=ReInter$sigma2_ri
                                ,beta=ReInter$beta_ri
                                ,InvCholCorr=ReInter$InvCholCorr_ri
                                ,Covtype=Covtype
        )
        All_Control_Genoud=c(Control_Genoud,Sub_Control_Genoud)

        tryCatch({
          outi=do.call(genoud
                      ,All_Control_Genoud)
        },
        error = function(c){
          print("genoud function error")
          print(c)
        })

        x_new=outi$par;names(x_new)=colnames(Psi)
        EI_hist[counter]=outi$value
        Psi=rbind(Psi,x_new)


        for(k in 1:length(Updates$Design_List)){
          Updates$Design_List[[k]]=Psi[,k]
          Updates$Distance_List[[k]]=as.matrix(dist(Updates$Design_List[[k]]),method="manhattan")
        }

        Y_new=c(BayesMSM(PatID=PatID
                         ,Outcome_Var=Outcome_Var
                         ,Treat_Vars=Treat_Vars
                         ,Treat_M_List=Treat_M_List
                         ,Outcome_M_List=Outcome_M_List
                         ,Data=BayesDat
                         ,G_List=G_List
                         ,Psi=x_new
                         ,DR=DR
                         ,Normalized=Normalized
                         ,Bayes=Bayes
                         ,B=1
                         ,Bayes_Seed=Bayes_Seed))

        X=model.matrix(~1,data=data.frame(Updates$Design_List[[1]]))
        y=c(y,Y_new)


        #generating the different random starting positions
        for(i in 1:ncol(Psi)){
          U=runif(Numbr_Samp,0.01,IthetasU[i])
          if(i==1){DD=U
          }else{DD=cbind(DD,U)}
        }
        DD=cbind(DD,runif(Numbr_Samp,0.01,0.99))
        pari=multi_eval(DD=DD
                        ,m=m
                        ,y=y
                        ,XX=X
                        ,IthetasL=IthetasL
                        ,IthetasU=IthetasU
                        ,IalphaL=IalphaL
                        ,IalphaU=IalphaU
                        ,Noise_Type=Noise_Type
                        ,Covtype=Covtype
                        ,Distance_List=Updates$Distance_List
                        ,Prior_List=Prior_List
                        ,Prior_Der_List=Prior_Der_List)


        thetas=pari[1:ncol(Psi)]
        alpha=pari[ncol(Psi)+1]

        if(Bayes==FALSE & sum(is.na(Likelihood_Limits))==0){
          for(k in 1:length(thetas)){
            plot_lik(
                   thetas=thetas,
                   alpha=alpha,
                   numb_varv=k,
                   abliny=thetas[k],
                   xvals=Likelihood_Limits[[k]],
                   counter=counter,
                   Distance_List=Updates$Distance_List,
                   y=y,
                   X=X,
                   Prior_List=Prior_List,
                   Prior_Der_List=Prior_Der_List,
                   Noise_Type=Noise_Type,
                   m=m,
                   Covtype=Covtype)

          }
        }


        Updates=Update(Psi=Psi,thetas=thetas,alpha=alpha,Design_List=Updates$Design_List,Distance_List=Updates$Distance_List,m=m,y=y,X=X,Covtype=Covtype)

        if(counter==Additional_Samp & Bayes==TRUE){
          Posterior_Draws=CIdraw_posterior(Y=y
                                        ,X=X
                                        ,N=N
                                        ,Design_List=Updates$Design_List
                                        ,Design_New_List=Design_New_List
                                        ,InvCholCorr=Updates$InvCholCorr
                                        ,thetas=thetas
                                        ,alpha=alpha
                                        ,beta=as.numeric(Updates$beta)
                                        ,v=as.numeric(Updates$v)
                                        ,sigma2=as.numeric(Updates$sigma2)
                                        ,Covtype=Covtype)
        }
        #######################################
        # reinterpolation parameter calculation
        #######################################
        ReInter=ReInterpolation(m=m
                                ,X=X
                                ,y=y
                                ,Psi=Psi
                                ,thetas=thetas
                                ,alpha=alpha
                                ,beta=as.numeric(Updates$beta)
                                ,InvCholCorr=Updates$InvCholCorr
                                ,Design_List=Updates$Design_List
                                ,Distance_List=Updates$Distance_List
                                ,Noise_Type=Noise_Type
                                ,Covtype=Covtype)

 }
  return(list(
    "BayesMSM"=list("Data"=Data
                    ,"Treat_M_List"=Treat_M_List
                    ,"Outcomoe_M_List"=Outcome_M_List
                    ,"Outcome_Var"=Outcome_Var
                    ,"Treat_Vars"=Treat_Vars
                    ,"PatID"=PatID
                    ,"G_List"=G_List
                    ,"Psi"=Psi
                    ,"DR"=DR
                    ,"Normalize"="Normalize"
                    ,"Bayes"=Bayes)
    ,"GP_Part"=list("Numbr_Samp"=Numbr_Samp
                    ,"IthetasU"=IthetasU
                    ,"IthetasL"=IthetasL
                    ,"IalphaU"=IalphaU
                    ,"IalphaL"=IalphaL
                    ,"Covtype"=Covtype)
    ,Control_Genoud=Control_Genoud
    ,"Updates"=Updates
    ,"ReInter"=ReInter
    ,"Posterior_Draws"=Posterior_Draws
    ,"Bayes_Seed"=Bayes_Seed
    ,"thetas"=thetas
    ,"alpha"=alpha
    ,"y"=y
    ,"EI_hist"=EI_hist
    ,"Likelihood_Limits"=Likelihood_Limits
  ))
}


#' Credible Intervals for Optimizing Regime
#'
#' @description This function allows the user to obtain credible intervals for the optimizing regimes
#' @param Sequential_Object object returned from SequenceFit
#' @param Boot Number of Bayesian bootstrap samples
#' @param Psi_new grid of points where posterior paths will be evaluated
#' @param N number of posterior samples
#' @param Location Location where to save results for bootstrapped samples, in case process is interrupted
#' @param Control_Genoud parameters to be fed into the genoud function
#' @return this function returns a matrix with the full posterior draws
#' @export
FitInfer=function(Design_Object,Boot_Start,Boot_End,Psi_new,N,Location=NA,Additional_Samp,Control_Genoud=list()){
  # yhist=vector()
  # xhist=matrix(NA,nrow=Boot,ncol=length(Design_Object$thetas))
  Boot=(Boot_End-Boot_Start+1)
  Coverage=matrix(NA,nrow=N*Boot,ncol=length(Design_Object$thetas)+1)

  for(loop_index in Boot_Start:Boot_End){
    print(paste("Bootstrap Number", loop_index))
    assign("loop_index_global", loop_index, envir = .GlobalEnv)
    start_fit=DesignFit(Data=Design_Object$Data
                       ,Treat_M_List=Design_Object$BayesMSM$Treat_M_List
                       ,Outcome_Var=Design_Object$BayesMSM$Outcome_Var
                       ,Treat_Vars=Design_Object$BayesMSM$Treat_Vars
                       ,PatID=Design_Object$BayesMSM$PatID
                       ,G_List=Design_Object$BayesMSM$G_List
                       ,Psi=Design_Object$BayesMSM$Psi
                       ,Bayes=TRUE # this values actually does not get used except to check if it was passed
                       ,Bayes_Seed=loop_index #this value actually does not get used except to check if it was passed

                       ,Numbr_Samp=Design_Object$GP_Part$Numbr_Samp
                       ,IthetasU=Design_Object$GP_Part$IthetasU
                       ,IthetasL=Design_Object$GP_Part$IthetasL
                       ,IalphaU=Design_Object$GP_Part$IalphaU
                       ,IalphaL=Design_Object$GP_Part$IalphaL
                       ,Covtype=Design_Object$GP_Part$Covtype)

    second_fit=SequenceFit(Previous_Fit=start_fit
                           ,Additional_Samp=Additional_Samp
                           ,N=N
                           ,Psi_new=Psi_new
                           ,Control_Genoud=Control_Genoud)
    start=(1+N*(loop_index-Boot_Start+1-1))
    stop=N*(loop_index-Boot_Start+1)
    Coverage[start:stop,]=second_fit$Posterior_Draws

    print(paste0("start:",start,"stop:",stop,"lengthdraw:",nrow(second_fit$Posterior_Draws)))
    if(!is.na(Location)){
      write.table(second_fit$Posterior_Draws,sep=",",file=Location,append=TRUE,col.names=FALSE,row.names=FALSE)
    }
  }
  # return(list("hists"=cbind(xhist,yhist), "draws"=Coverage))
  return(list("draws"=Coverage))
}


