#' A BayesMSM Function
#'
#' This function allows you to fit a Bayesian MSM
#' @keywords Bayes
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
#' @param MSM_Model model for continuous MSM
#' @param Normalized boolean indicating whether to use normalized or non-normalized estimators
#' @param B number of Bayesian bootstrap samples
#' @export
#' @importFrom gtools rdirichlet
BayesMSM=function(Data,PatID,Treat_Vars,Outcome_Var,Treat_M_List,Outcome_M_List,MSM_Model=NULL,G_List=NULL,Psi,DR=FALSE,Bayes=TRUE,Normalized=FALSE,B=100,Bayes_Seed=1){
  nstages=length(Treat_M_List)
  Coef_Return=NULL

  #checkin that regime indices match with psi parameters
  g_string=paste(G_List,collapse=" ");
  parameter_checks=sapply(colnames(Psi),FUN=function(X){grepl(X, g_string, fixed = TRUE)})
  if(sum(parameter_checks==FALSE)>0){
    print("Warning: parameters in G_List may not match parameter names in Psi matrix")
  }

  if(Bayes==FALSE){
    B=1
  }

  if(!is.null(G_List)&!is.null(MSM_Model)){
    #If requesting a Dynamic MSMs
    Aug_Data=Build_Augmented(Data=Data,Treat_Vars=Treat_Vars,Psi=Psi,G_List=G_List)
  }
  for(b in 1:B){
    if(Bayes==TRUE){
      # old <- .Random.seed #setting seed locally
      set.seed(b+Bayes_Seed) #bayes_seed is set to the bootstrap counter from the GP function
      Dirich_Vect=t(rdirichlet(1,rep(1,nrow(Data))))
      # .Random.seed <<- old
    }else{
      Dirich_Vect=rep(1,nrow(Data))
    }
    Data_Analysis=Data
    Data_Analysis$Dirich=Dirich_Vect
    Prob_Matrix=cbind(Data_Analysis[,PatID,drop=FALSE])
    for(i in 1:nstages){
      model_treat=glm(Treat_M_List[[i]], data=Data_Analysis, family="quasibinomial",weights=Dirich)
      prob_treat=predict(model_treat, type="response")
      prob_obs_treat=Data_Analysis[,Treat_Vars[i]]*prob_treat+(1-Data_Analysis[,Treat_Vars[i]])*(1-prob_treat)
      Prob_Matrix[,paste0("ProbTreat",i)]=prob_obs_treat
    }
    Data_Analysis=data.frame(Data_Analysis,Prob_Matrix)

    if(DR==TRUE){  #DR Grid
      if(is.null(nrow(Psi))){ #if just one regime
        Coef=DR_Est(Psi=Psi,Data=Data_Analysis,Prob_Matrix=Prob_Matrix,nstages=nstages,Outcome_Var=Outcome_Var,Bayes=Bayes,Normalized=Normalized)
      }else{                    #if a grid of regimes
        Coef=apply(X=Psi,MARGIN=1,FUN=DR_Est,Data=Data_Analysis,Prob_Matrix=Prob_Matrix,nstages=nstages,Outcome_Var=Outcome_Var,Bayes=Bayes,Normalized=Normalized)
      }
    }else if(is.null(MSM_Model)){ #IPW Grid
      if(is.null(nrow(Psi))){ #if just one regime
        Coef=IPW_Est(Psi=Psi,Data=Data_Analysis,Treat_Vars=Treat_Vars,nstages=nstages,G_List=G_List,Outcome_Var=Outcome_Var,Bayes=Bayes,Normalized=Normalized)
      }else{                    #if a grid of regimes
        Coef=apply(X=Psi,MARGIN=1,FUN=IPW_Est,Data=Data_Analysis,Treat_Vars=Treat_Vars,nstages=nstages,G_List=G_List,Outcome_Var=Outcome_Var,Bayes=Bayes,Normalized=Normalized)
      }
    }else{
      Prob_Matrix_w_Dirich=cbind(Prob_Matrix,"Dirich"=Dirich_Vect)
      Model_Data=merge(Aug_Data,Prob_Matrix_w_Dirich, by = PatID)#merging probabilities which depend on Dirichlet draw
      Coef=IPW_Marginal_Model(MSM_Model=MSM_Model
                              ,Data=Model_Data
                              ,Treat_Vars=Treat_Vars
                              ,nstages=nstages)
    }
    if(is.null(nrow(Psi))){
      Coef_Return=c(Coef_Return,Coef)
    }else{
      Coef_Return=rbind(Coef_Return,Coef)
    }
  }
  return(Coef_Return)
}

IPW_Marginal_Model=function(MSM_Model,Data,Treat_Vars,nstages,G_List){
  for(i in 1:nstages){
    if(i==1){
      Cum_Proba=Data[,paste0("ProbTreat",i)]
      Cum_Adhere=Data[,paste0("Adhere",i)]
    }else{
      Cum_Proba=Cum_Proba*Data[,paste0("ProbTreat",i)]
      Cum_Adhere=Cum_Adhere*Data[,paste0("Adhere",i)]
    }
  }
  Data$Cum_Proba=Cum_Proba
  Data$Cum_Adhere=Cum_Adhere
  MarginalModel=lm(MSM_Model,data=Data,weights= Cum_Adhere*Dirich/Cum_Proba)
  return(coef(MarginalModel))
}


IPW_Est=function(Psi,Data,Treat_Vars,nstages,G_List,Outcome_Var,Bayes,Normalized){
  Data=g_compliant(X=Psi,Data=Data,Treat_Vars=Treat_Vars,G_List=G_List)
  for(i in 1:nstages){
    if(i==1){
      Cum_Proba=Data[,paste0("ProbTreat",i)]
      Cum_Adhere=Data[,paste0("Adhere",i)]
    }else{
      Cum_Proba=Cum_Proba*Data[,paste0("ProbTreat",i)]
      Cum_Adhere=Cum_Adhere*Data[,paste0("Adhere",i)]
    }
  }
  Data$Cum_Proba=Cum_Proba
  Data$Cum_Adhere=Cum_Adhere
  Y=Data[,Outcome_Var]
  if(Normalized==TRUE){ #normalized either bayes or not
    returny=with(Data,sum(Cum_Adhere*Dirich*Y/Cum_Proba)/sum(Cum_Adhere*Dirich/Cum_Proba))
  }else if(Bayes==FALSE){ #not normalized not bayes
    returny=with(Data,mean(Cum_Adhere*Dirich*Y/Cum_Proba))
  }else if(Bayes==TRUE){ #not normalized bayes
    returny=with(Data,sum(Cum_Adhere*Dirich*Y/Cum_Proba))
  }
    return(returny)
}

g_compliant=function(X,Data,Treat_Vars,G_List){
  Data[,names(X)]=t(replicate(nrow(Data),X)) # creates new columns
  for(i in 1:length(Treat_Vars)){
    g_rule=with(Data,eval(G_List[[i]]))
    Z=Data[,Treat_Vars[i]]
    Data[,paste0("Adhere",i)]=Z*g_rule+(1-Z)*!g_rule
  }
  Data$Regime=paste(X,collapse="_")
  return(Data)
}

Build_Augmented=function(Data,Treat_Vars,Psi,G_List){
    Aug_Data=apply(X=Psi,1,FUN=g_compliant,Data=Data,Treat_Vars=Treat_Vars,G_List=G_List)
    Aug_Data=do.call(rbind,Aug_Data)
  return(Aug_Data)
}

DR_Est=function(Psi,Data,Prob_Matrix,nstages,Outcome_Var,Normalized,Bayes){
  Phi_Adhere_Matrix=Phi_Adhere_Create(Psi=Psi,Outcome_M_List=Outcome_M_List,Data=Data,G_List=G_List,Treat_Vars=Treat_Vars,nstages=nstages) # note that this is regime specific
  SF=Phi_Adhere_Matrix


  if(Normalized==TRUE){
    if(Bayes==TRUE){
      AugY=Data$Dirich*SF[,"Phi1"]
    }else{
      AugY=Data$Dirich*SF[,"Phi1"]/nrow(SF)
    }
  }else{
    AugY=Data$Dirich*SF[,"Phi1"]
  }
  # if(Bayes==TRUE){#works for both normalized and unormalized as there is no extra w term
  #   AugY=Data$Dirich*SF[,"Phi1"]
  # }else{#works for both normalized and unormalized as there is no extra w term
  #   AugY=SF[,"Phi1"]/nrow(SF)
  # }
  Adhere_Product=rep(1,nrow(SF))
  Prob_Product=rep(1,nrow(SF))
  for(i in 1:nstages){
    Adhere_Product=Adhere_Product*SF[,paste0("Adhere",i)]
    Prob_Product=Prob_Product*Prob_Matrix[,paste0("ProbTreat",i)]
    if(i<nstages){
      if(Normalized==TRUE){
        AugY=AugY+Data$Dirich*(Adhere_Product*(SF[,paste0("Phi",i+1)]-SF[,paste0("Phi",i)])/Prob_Product)/sum(Adhere_Product*Data$Dirich/Prob_Product)
      }else{
        AugY=AugY+Data$Dirich*Adhere_Product*(SF[,paste0("Phi",i+1)]-SF[,paste0("Phi",i)])/Prob_Product
      }
    }else{
      Y=Data[,Outcome_Var]
      if(Normalized==TRUE){
        AugY=AugY+Data$Dirich*(Adhere_Product*(Y-SF[,paste0("Phi",i)])/Prob_Product)/sum(Adhere_Product*Data$Dirich/Prob_Product)
      }else{
        AugY=AugY+Data$Dirich*Adhere_Product*(Y-SF[,paste0("Phi",i)])/Prob_Product
      }
    }
  }

  if(Normalized==TRUE){
    return(sum(AugY))
  }else{
    if(Bayes==TRUE){
      return(sum(AugY))
    }else{
      return(mean(AugY)) #the frequentist needs the extra 1/n factor
    }
  }
}


Phi_Adhere_Create=function(Psi,Outcome_M_List,Data,G_List,Treat_Vars,nstages){
  #does one regime at a time
  Data[,names(Psi)]=t(replicate(nrow(Data),Psi))
  Data_g=Data
  for(i in nstages:1){
    Temp_Data=Data
    g_rule=with(Data,eval(G_List[[i]]))+0
    Temp_Data[,Treat_Vars[i]]=g_rule                         # follow regime in present and forward
    Data_g[,Treat_Vars[i]]=g_rule                            # follow regime for entire history
    model_outcome=lm(Outcome_M_List[[i]],data=Data, weights = Dirich)
    Data$Pseudo_Outcome=predict(model_outcome, newdata=Temp_Data) #this is the pseudo outcome for one stage prior

    Z=Data[,Treat_Vars[i]]
    if(i==nstages){
      Phi_Matrix=data.frame(predict(model_outcome, newdata=Data_g));colnames(Phi_Matrix)=paste0("Phi",i)
      Adhere_Matrix=data.frame(Z*g_rule+(1-Z)*!g_rule);colnames(Adhere_Matrix)=paste0("Adhere",i)
    }else{
      Phi_Matrix[,paste0("Phi",i)]=predict(model_outcome, newdata=Data_g)
      Adhere_Matrix[,paste0("Adhere",i)]=Z*g_rule+(1-Z)*!g_rule
    }
  }
  return(data.frame(Phi_Matrix,Adhere_Matrix))
}



