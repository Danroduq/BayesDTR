#rm(list=ls(all=TRUE))  # removes all objects from the current workspace
library(mvtnorm)
library(ggplot2)
library(plotly)#used for plotting
library(rgenoud)
library(reshape2)
library(matrixcalc)
library(BayesDTR)
set.seed(1)
#-------------------------------------------------------------------------------
Treat_M_List=list(formula1="z1~karnof+race+gender+symptom+str2+cd4.0+wtkg"
                  ,formula2="z2~karnof+race+gender+symptom+str2+cd4.20+wtkg+z1")

Outcome_Var="cd4.outcome"
Treat_Vars=c("z1","z2")
PatID="pidnum"

G_List=list(g1=expression(cd4.0>=theta1)
          ,g2=expression(cd4.20>=theta2))

sequence1=seq(200,500,100)
sequence2=seq(200,500,100)
theta=as.matrix(expand.grid(sequence1,sequence2));colnames(theta)=c("theta1","theta2")

#-------------------------------------------------------------------------------
checki=exists("argies")
start_num=16
m=start_num

Noise_Type="Homosk"
type_simul=2
Trigger=0
numbr_samp=45

#For normalized IPW
IthetasU=c(600,600);IthetasL=c(0.01,0.01) # for 2 vs 1 comparison
IalphaU=0.9999999; IalphaL=0.1

expressi1="0";expressi1_der="0"
expressi2="0";expressi2_der="0"

genoud_limits=matrix(c(200,200,600,600),ncol=2)

thetas_lag=c(80,100)
alpha_lag=0.99

#-------------------------------------------------------------------------------

y=c(BayesMSM(PatID=PatID
         ,Outcome_Var=Outcome_Var
         ,Treat_Vars=Treat_Vars
         ,Treat_M_List=Treat_M_List
         ,Data=BayesDat
         ,G_List=G_List
         ,Theta=theta
         ,Normalized=TRUE
         ,Bayes=FALSE))

design1=theta[,1]
design2=theta[,2]

distance1=as.matrix(dist(design1),method="manhattan")
distance2=as.matrix(dist(design2),method="manhattan")
X=model.matrix(~1,data=data.frame(design1))




#------Generalized Addition-----------#
# thetas
Design_List=list(design1,design2)
Distance_List=list(distance1,distance2)
Expressi_List=list(expressi1,expressi2)
Expressi_Der_List=list(expressi1_der,expressi2_der)

#---------------------------------------------------------------------------------------


U1=runif(numbr_samp,0.01,IthetasU[1])
U2=runif(numbr_samp,0.01,IthetasU[2])
U3=runif(numbr_samp,0.01,0.99)
DD=cbind(U1,U2,U3)
pari=multi_eval(DD=DD,
                y=y,
                IthetasL=IthetasL,
                IthetasU=IthetasU,
                IalphaL=IalphaL,
                IalphaU=IalphaU,
                Noise_Type=Noise_Type,
                thetas_lag=thetas_lag,
                alpha_lag=alpha_lag
)



thetas=pari[1:2]
alpha=pari[3]
thetas;alpha
#------------

Corr_List=list()
for(i in 1:(length(Distance_List))){
  Corr_List[[i]]=corr_compute(thetas[i],Distance_List[[i]])
}
Corr=alpha*Reduce("*", Corr_List)+(1-alpha)*diag(rep(1,ncol(Distance_List[[1]])))
outisigbet=sigma2_compute(Corr=Corr,y=y,X=X);if(Trigger==1){next}
v=outisigbet$sigma2_v
beta=outisigbet$beta
InvCholCorr=solve(outisigbet$CholCorr)

#variance
sig2=alpha*v;sig2
#nugget
v-sig2

plot_lik(thetas=thetas,
         alpha=alpha,
         numb_varv=1,
         abliny=thetas[1],
         xvals=seq(0.01,500,1),
         counter=0,
         Distance_List=Distance_List,
         y=y,
         Expressi_List=Expressi_List,
         Expressi_Der_List=Expressi_Der_List,
         Noise_Type=Noise_Type)

plot_lik(thetas=thetas,
         alpha=alpha,
         abliny=thetas[2],
         numb_varv=2,
         xvals=seq(0.01,500,1),
         counter=0,
         Distance_List=Distance_List,
         y=y,
         Expressi_List=Expressi_List,
         Expressi_Der_List=Expressi_Der_List,
         Noise_Type=Noise_Type)

#######################################
# reinterpolation parameter calculation
#######################################

yhat=apply(theta,1,FUN=posterior_mean,XX=as.matrix(X),Y=y,thetas=thetas,alpha=alpha,beta=beta,InvCholCorr=InvCholCorr,Design_List=Design_List,Noise_Type=Noise_Type)

Corr_ri_List=list()
for(i in 1:(length(Distance_List))){
  Corr_ri_List[[i]]=corr_compute(thetas[i],Distance_List[[i]])
}
Corr_ri=Reduce("*", Corr_ri_List)

outisigbet_ri=sigma2_compute(Corr_ri,y=yhat,X=X);if(Trigger==1){next}
sigma2_ri=outisigbet_ri$sigma2_v
beta_ri=outisigbet_ri$beta
InvCholCorr_ri=solve(outisigbet_ri$CholCorr)


Y_true_max_init=max(yhat)
index_max=which(yhat==Y_true_max_init)
x_true_max_init=theta[index_max,]
if(length(index_max)>1)x_true_max_init=x_true_max_init[1,]
Y_true_max_ri=Y_true_max_init


beta;beta_ri
sig2;sigma2_ri


########################
# Itirative Selection
########################
breaky=0
counter=0
xhist=matrix(nrow=25,ncol=length(Distance_List))
yhist=rep(NA,25)
Pari_Hist_Matrix=matrix(NA,nrow=25,ncol=length(Distance_List)+1)


EI_hist=vector()
alphahist=rep(NA,25)

Distx_Hist_List=list()
MaxxCI_Hist_List=list()
LBx_Hist_List=list()
UBx_Hist_List=list()
for(k in 1:length(Distance_List)){
  Distx_Hist_List[[k]]=c(NA,NA,NA,NA,NA,NA)
  MaxxCI_Hist_List[[k]]=c(NA,NA,NA,NA,NA,NA)
  LBx_Hist_List[[k]]=c(NA,NA,NA,NA,NA,NA)
  UBx_Hist_List[[k]]=c(NA,NA,NA,NA,NA,NA)
}

while(counter<15){
  counter=counter+1
  m=m+1

  TryAgain=0
  outi=tryCatch({
    genoud(EI, nvars=2,pop.size=1000,BFGS=TRUE,max=TRUE,ymax=Y_true_max_ri,
           Y=yhat,Design_List=Design_List,thetas=thetas,sigma2=sigma2_ri,beta=beta,InvCholCorr=InvCholCorr_ri,
           Domains=genoud_limits, boundary.enforcement = 2)
  },
  error = function(c){
      assign("Trigger", 1, env=globalenv())
  } )
  if(Trigger==1){break}

  # EI(x=c(370,340),ymax=Y_true_max_ri,
  #        Y=yhat,Design_List=Design_List,thetas=thetas,sigma2=sigma2_ri,beta=beta,InvCholCorr=InvCholCorr_ri)

  EI_hist[counter]=outi$value
  x_new=outi$par;names(x_new)=colnames(theta)
  theta=rbind(theta,x_new)


  for(k in 1:length(Design_List)){
    Design_List[[k]]=theta[,k]
    Distance_List[[k]]=as.matrix(dist(Design_List[[k]]),method="manhattan")
  }

  X=model.matrix(~1,data=data.frame(Design_List[[1]]))
  Y_new=c(BayesMSM(PatID=PatID
               ,Outcome_Var=Outcome_Var
               ,Treat_Vars=Treat_Vars
               ,Treat_M_List=Treat_M_List
               ,Data=BayesDat
               ,G_List=G_List
               ,Theta=x_new
               ,Normalized=TRUE
               ,Bayes=FALSE))



  y=c(y,Y_new)


    U1=c(runif(numbr_samp,0.01,IthetasU[1]),IthetasL[1],IthetasL[1],IthetasU[1],IthetasU[1])
    U2=c(runif(numbr_samp,0.01,IthetasU[2]),IthetasL[2],IthetasU[2],IthetasL[2],IthetasU[2])
    U3=c(runif(numbr_samp,0.01,0.999),0.99,.99,0.99,0.99)
    DD=cbind(U1,U2,U3)
    pari=multi_eval(DD=DD,
                    y=y,
                    IthetasL=IthetasL,
                    IthetasU=IthetasU,
                    IalphaL=IalphaL,
                    IalphaU=IalphaU,
                    Noise_Type=Noise_Type,thetas_lag=thetas,alpha_lag=alpha)


  thetas=pari[1:2]
  alpha=pari[3]

  Pari_Hist_Matrix[counter,]=pari


  for(i in 1:(length(Distance_List))){
    Corr_List[[i]]=corr_compute(thetas[i],Distance_List[[i]])
  }
  Corr=alpha*Reduce("*", Corr_List)+(1-alpha)*diag(rep(1,ncol(Distance_List[[1]])))

  outisigbet=sigma2_compute(Corr=Corr,y=y,X=X);if(Trigger==1){break}
  v=outisigbet$sigma2_v
  beta=outisigbet$beta
  InvCholCorr=solve(outisigbet$CholCorr)

  #variance
  sig2=alpha*v;sig2
  #nugget
  v-sig2

  par(mfrow=c(1,2))
  plot_lik(theta=thetas,
           alpha=alpha,
           numb_varv=1,
           abliny=thetas[1],
           xvals=seq(0.01,thetas[1]+250,1),
           counter=counter,
           Distance_List=Distance_List,
           y=y,
           Expressi_List=Expressi_List,
           Expressi_Der_List=Expressi_Der_List,
           Noise_Type=Noise_Type)

  plot_lik(thetas=thetas,
           alpha=alpha,
           numb_varv=2,
           abliny=thetas[2],
           xvals=seq(0.01,thetas[2]+250,1),
           counter=counter,
           Distance_List=Distance_List,
           y=y,
           Expressi_List=Expressi_List,
           Expressi_Der_List=Expressi_Der_List,
           Noise_Type=Noise_Type)


  #######################################
  # reinterpolation parameter calculation
  #######################################

  yhat=apply(theta,1,FUN=posterior_mean,XX=X,Y=y,thetas=thetas,alpha=alpha,beta=beta,InvCholCorr=InvCholCorr,Design_List=Design_List,Noise_Type=Noise_Type)

  for(i in 1:(length(Distance_List))){
    Corr_ri_List[[i]]=corr_compute(thetas[i],Distance_List[[i]])
  }
  Corr_ri=Reduce("*", Corr_ri_List)


  # y=yhat
  Corr=Corr_ri

  outisigbet_ri=sigma2_compute(as.matrix(Corr_ri),y=as.vector(yhat),X=X);if(Trigger==1){break}
  sigma2_ri=outisigbet_ri$sigma2_v
  beta_ri=outisigbet_ri$beta
  InvCholCorr_ri=solve(outisigbet_ri$CholCorr)


  Y_true_max_ri=max(yhat)




  yhist[counter]=Y_true_max_ri
  index_max=which(yhat==Y_true_max_ri)
  x_true_max=theta[index_max,]
  if(length(index_max)>1)x_true_max=x_true_max[1,]
  xhist[counter,]=x_true_max


  beta;beta_ri
  sig2;sigma2_ri
}
#
# ReInterpolation=function(X,y,thetas,alpha,beta,InvCholCorr,Design_List,Noise_Type){
#   yhat=apply(theta,1,FUN=posterior_mean,XX=X,Y=y,thetas=thetas,alpha=alpha,beta=beta,InvCholCorr=InvCholCorr,Design_List=Design_List,Noise_Type=Noise_Type)
#
#   for(i in 1:(length(Distance_List))){
#     Corr_ri_List[[i]]=corr_compute(thetas[i],Distance_List[[i]])
#   }
#   Corr_ri=Reduce("*", Corr_ri_List)
#
#   y_save=y
#   Corr_save=Corr
#   y=yhat
#   Corr=Corr_ri
#
#   outisigbet_ri=sigma2_compute(as.matrix(Corr_ri),y=as.vector(yhat),X=X);if(Trigger==1){break}
#   sigma2_ri=outisigbet_ri$sigma2_v
#   beta_ri=outisigbet_ri$beta
#   InvCholCorr_ri=solve(outisigbet_ri$CholCorr)
#   Y_true_max_ri=max(yhat)
#
#   return(list(outisigbet_ri=outisigbet_ri
#               ,sigma2_ri=sigma2_ri
#               ,beta_ri=beta_ri
#               ,InvCholCorr_ri=InvCholCorr_ri
#               ,Y_true_max_ri))
# }
