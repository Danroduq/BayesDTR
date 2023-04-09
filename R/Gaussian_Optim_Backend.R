#'@import mvtnorm plotly rgenoud matrixcalc
NULL

covi=function(h,theta,Covtype){
    #matern 5/2 covariance
    if(Covtype==1){
      to_return=(1+sqrt(3)*abs(h)/theta)*exp(-sqrt(3)*abs(h)/theta)
    }else if(Covtype==2){
      to_return=(1+sqrt(5)*abs(h)/theta + 5*h^2/(3*theta^2))*exp(-sqrt(5)*abs(h)/theta)
    }
  return(to_return)
}
dercovi=function(h,theta,Covtype){
  if(Covtype==1){
    to_return=3*h^2/theta^3*exp(-sqrt(3)*abs(h)/theta) # it should be abs(h)
  }else if(Covtype==2){
    to_return=5*h^2*(theta+sqrt(5)*abs(h))*exp(-sqrt(5)*abs(h)/theta)/(3*theta^4)
  }
  return(to_return)
}

corr_compute=function(theta,distance,Covtype){
  Corr=structure(vapply(distance, covi, numeric(1),theta=theta,Covtype=Covtype), dim=dim(distance))
  return(Corr)
}

dif_func=function(x,design,theta,Covtype){
  distance_between=x-design
  to_return=covi(distance_between,theta=theta,Covtype=Covtype)
  return(to_return)
}


sigma2_compute=function(Corr,y,X,m){
  tryCatch({
    cholcorr= chol(Corr)
    term1= backsolve(t(cholcorr),y,upper.tri=FALSE)

    #Step 4:  backsolve t(T)X=F (steps refer to the likelihood optimization paper followed)
    term2= backsolve(t(cholcorr),X, upper.tri=FALSE)
    #Step 5:
    mo=lm(term1~-1+term2)
    z=mo$residuals
    sigma2=(1/m)*t(z) %*%z #depending on which case we are in, this can be either be sigma2 or v
    beta=mo$coefficients

    return(list(sigma2_v=sigma2,beta=beta, CholCorr=cholcorr))
  },
  error = function(c){
    print("Error Updating Model Parameters")
    print(c)
  })
}


GP_likelihood=function(params,Distance_List,y,XX,indic=0,varyv,thetasv,alphav,numb_varv,Prior_List,Prior_Der_List,Noise_Type,m,Covtype){
  X=XX
  if(indic==0){
    params_to_name=sub("\\_.*","", names(Prior_List))
    for(i in 1:length(Distance_List)){
      assign(params_to_name[i],params[i])
    }
    if(Noise_Type=="Homosk")alpha=params[3]
  }else{
    params=thetasv
    params[numb_varv]=varyv

    params_to_name=sub("\\_.*","", names(Prior_List))
    for(i in 1:length(Distance_List)){
      assign(params_to_name[i],params[i])
    }
    if(Noise_Type=="Homosk"){params=c(params,alphav);alpha=alphav}
  }

  Corr_List=list()
  for(i in 1:(length(Distance_List))){
    Corr_List[[i]]=corr_compute(theta=params[i],distance=Distance_List[[i]],Covtype=Covtype)
  }
  #Step 1:
  if(Noise_Type=="Interpol"){
    Corr=Reduce("*", Corr_List)
  }else if(Noise_Type=="Homosk"){
    Corr=alpha*Reduce("*", Corr_List)+(1-alpha)*diag(rep(1,ncol(Distance_List[[1]])))
  }
  #Step 2:
  cholcorr=chol(Corr)
  #Step 3:
  term1= backsolve(r=t(cholcorr),x=y,upper.tri=FALSE)
  #Step 4:  backsolve t(T)X=F
  term2= backsolve(t(cholcorr),X, upper.tri=FALSE)

  #Step 5:
  mo=lm(term1~-1+term2)
  z=mo$residuals
  #step 6
  sigma2_v=(1/m)*t(z) %*%z
  #paste(expressi_der)
  logdetcholcorr=2*sum(log(diag(cholcorr)))

  out=(m)*log(sigma2_v) +logdetcholcorr
  for(i in 1:(length(Distance_List))){
    out=out-2* eval(parse(text=Prior_List[[i]]))
  }
  return(out)
}


derloglik=function(params,Distance_List,y,XX,indic=0, thetasv,alphav,numb_varv,Prior_List,Prior_Der_List,Noise_Type,m,Covtype){
  X=XX
  if(indic==0){
     params_to_name=sub("\\_.*","", names(Prior_List))
    for(i in 1:length(Distance_List)){
      assign(params_to_name[i],params[i])
    }
    if(Noise_Type=="Homosk")alpha=params[3]
  }else{
    params=thetasv
    params[numb_varv]=varyv
    params_to_name=sub("\\_.*","", names(Prior_List))
    for(i in 1:length(Distance_List)){
      assign(params_to_name[i],params[i])
    }
    if(Noise_Type=="Homosk"){params=c(params,alphav);alpha=alphav}
  }
  #-----------------------------#
  # Derivative of Log Likelihood
  #-----------------------------#
  Corr_List=list()
  for(i in 1:(length(Distance_List))){
    Corr_List[[i]]=corr_compute(theta=params[i],distance=Distance_List[[i]],Covtype=Covtype)
  }
  #Step 1:
  if(Noise_Type=="Interpol"){
    Corr=Reduce("*", Corr_List)
  }else if(Noise_Type=="Homosk"){
    R=Reduce("*", Corr_List)
    Corr=alpha*R+(1-alpha)*diag(rep(1,ncol(Distance_List[[1]])))
  }

  cholcorr= chol(Corr)
  term1= backsolve(t(cholcorr),y,upper.tri=FALSE)
  #Step 4:  backsolve t(T)X=F
  term2= backsolve(t(cholcorr),X, upper.tri=FALSE)
  #Step 5:
  mo=lm(term1~-1+term2)
  z=mo$residuals
  sigma2_v=as.numeric((1/m)*t(z) %*%z)

  #computing gradient
  #step 1
  gterm1=backsolve(cholcorr,z)
  #step 2
  Corr_inv=chol2inv(cholcorr)
  #step 3
  gterm2_List=list()

    for(k in 1:(length(Distance_List))){
      gterm2_List[[k]]=structure(vapply(Distance_List[[k]], dercovi, numeric(1),theta=params[k],Covtype=Covtype),dim=dim(Distance_List[[k]]))*Reduce("*", Corr_List[-k])
    }
    if(Noise_Type=="Homosk")gterm2_alpha=R-diag(nrow(R))


  #step 4
  gterm3_vect=rep(0,length(Distance_List))
  # gterm3_1=0
  # gterm3_2=0
  if(Noise_Type=="Homosk")gterm3_alpha=0

  #step 5
  gterm4_vect=rep(0,length(Distance_List))
  # gterm4_1=0
  # gterm4_2=0
  if(Noise_Type=="Homosk")gterm4_alpha=0


  if(Noise_Type=="Interpol"){
    for(i in 1:(m-1)){for(j in (i+1):m){
      #step 4
      producty=gterm1[i]*gterm1[j]
      for(k in 1:(length(Distance_List))){
        gterm3_vect[k]=gterm3_vect[k]+producty*gterm2_List[[k]][i,j]
        #step 5
        gterm4_vect[k]=gterm4_vect[k]+Corr_inv[i,j]*gterm2_List[[k]][i,j]
      }
    }}
  }else if(Noise_Type=="Homosk"){
    for(i in 1:(m-1)){for(j in (i+1):m){
      producty=gterm1[i]*gterm1[j]
      for(k in 1:(length(Distance_List))){
        #step 4
        gterm3_vect[k]=gterm3_vect[k]+producty*gterm2_List[[k]][i,j]
        #step 5
        gterm4_vect[k]=gterm4_vect[k]+Corr_inv[i,j]*gterm2_List[[k]][i,j]
      }
      gterm3_alpha=gterm3_alpha+gterm1[i]*gterm1[j]*gterm2_alpha[i,j]
      gterm4_alpha=gterm4_alpha+Corr_inv[i,j]*gterm2_alpha[i,j]
    }}
  }

  #step 4
  gterm3_vect=(2/sigma2_v)*gterm3_vect
  if(Noise_Type=="Homosk")gterm3_alpha=(2/sigma2_v)*gterm3_alpha

  #step5
  gterm4_vect=2*gterm4_vect
  if(Noise_Type=="Homosk")gterm4_alpha=2*gterm4_alpha

  # paste(expressi)# need to use this parameter otherwise optim is unhappy
  #step6
  dertheta_vect=vector()
  if(Noise_Type=="Interpol"){
    for(k in 1:(length(Distance_List))){
      dertheta_vect[k]=(-gterm3_vect[k]+gterm4_vect[k])-2*eval(parse(text=Prior_Der_List[k]))
    }
    return(c(as.numeric(dertheta_vect)))
  }else if(Noise_Type=="Homosk"){
    for(k in 1:(length(Distance_List))){
      dertheta_vect[k]=alpha*(-gterm3_vect[k]+gterm4_vect[k])-2*eval(parse(text=Prior_Der_List[k]))
    }
    deralpha=(-gterm3_alpha+gterm4_alpha)
    return(c(dertheta_vect, deralpha))
  }
}



offdiag_compute=function(Design_List,Design_New_List,thetas,Covtype){
  x_new=matrix(model.matrix(~1, data=data.frame(Design_New_List[[1]])))

  R_New_List=list()
  Corr_Off_Diag_List=list()
  for(k in 1:length(Design_List)){
    distance_new=as.matrix(dist(Design_New_List[[k]]),method="manhattan")
    # distance2_new=as.matrix(dist(design2_new),method="manhattan")

    R_New_List[[k]]=structure(vapply(distance_new,covi,numeric(1),theta=thetas[k],Covtype=Covtype),dim=dim(distance_new))

    Corr_Off_Diag_List[[k]]=sapply(Design_New_List[[k]],FUN=dif_func,design=Design_List[[k]],theta=thetas[k],Covtype=Covtype)

  }
  R_new=Reduce("*",R_New_List)
  corr_off_diag=Reduce("*",Corr_Off_Diag_List)

  # R_new=R1_new*R2_new
  # corr_off_diag=corr_off_diag1*corr_off_diag2
  return(list("x_new"=x_new,"corr_off_diag"=corr_off_diag,"R_new"=R_new))
}

EI=function(x,Y,ymax,Design_List,thetas,sigma2,beta,InvCholCorr,X,Covtype){

  Design_New_List=list()
  for(k in 1:length(Design_List)){
    Design_New_List[[k]]=x[k]
  }

  beta=matrix(beta) #making sure its column wise vector
  diagy=offdiag_compute(Design_List=Design_List,Design_New_List=Design_New_List,thetas=thetas,Covtype=Covtype)
  corr_off_diag=diagy$corr_off_diag
  mu_new=t(diagy$x_new)%*%beta+t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)
  sigma_new=sigma2*((diagy$R_new-t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%corr_off_diag))


  if(sigma_new>10e-12){
    argi=(mu_new-ymax)/sqrt(sigma_new)
    to_return=(mu_new-ymax)*pnorm(argi)+sqrt(sigma_new)*dnorm(argi)
  }else{
    to_return=0
  }
  return(to_return)
}

CIdraw_posterior<- function(N,
                            Y,
                            X,
                            Design_List,
                            Design_New_List,
                            InvCholCorr,
                            thetas,
                            sigma2,
                            alpha,
                            beta,
                            v,
                            Noise_Type
                            ,Covtype){
  # Y=y
  # X=X
  # N=100
  # Design_List=Updates$Design_List
  # Design_New_List=Design_New_List
  # InvCholCorr=Updates$InvCholCorr
  # thetas=thetas
  # alpha=alpha
  # beta=as.numeric(Updates$beta)
  # v=as.numeric(Updates$v)
  # sig2=as.numeric(Updates$sig2)
  ###########################################
  # Finding the maximizer of the max function
  ###########################################
  beta=matrix(beta) #making sure its column wise vector
  diagy=offdiag_compute(Design_List=Design_List,Design_New_List=Design_New_List,thetas=thetas,Covtype=Covtype)
  corr_off_diag=diagy$corr_off_diag


  mu_new=diagy$x_new%*%beta+alpha*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)
  sigma_new=sigma2*(diagy$R_new-sigma2*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%corr_off_diag/v)


  sigma_new[lower.tri(sigma_new)] = t(sigma_new)[lower.tri(sigma_new)]

  index=which(mu_new==max(mu_new))

  max_x_mu_vect=vector()
  for(k in 1:length(Design_List)){
    max_x_mu_vect[k]=Design_New_List[[k]][index]
  }

  checki2=is.positive.definite(sigma_new)

  tryCatch({
    Drawy=CIdraw_posterior_sub1(N=N,mu_new=mu_new,sigma_new=sigma_new, Design_New_List=Design_New_List,checki2)
  },
  error = function(c){
    print("Error drawing posterior paths")
    print(c)
  })


    return(matrix(unlist(Drawy), ncol =(length(Design_List)+1), nrow =N))
}

CIdraw_posterior_sub1=function(N=N,mu_new,sigma_new, Design_New_List,checki2){
  if(checki2==FALSE){
    #sampling full normal vector
    Y_post <- matrix(NA, nrow = length(Design_New_List[[1]]), ncol = N)

    Maxx_List=list()
    maxy=vector()
    for (n in 1:N) {
      Y_post[, n] <-rmvnorm(1,mean=mu_new,sigma=sigma_new)
      indexio=which(Y_post[, n]==max(Y_post[, n]))
      for(k in 1:length(Design_New_List)){
        if(n==1){Maxx_List[[k]]=vector()}
        Maxx_List[[k]][n]=Design_New_List[[k]][indexio]
      }
      maxy[n]=Y_post[indexio,n]
    }
    Maxx_List
  }else if(checki2==TRUE){
    #Sampling Normal vector in thirds
    #setting the vector limits
    mu=nrow(sigma_new)
    mu1=floor(mu/3)
    mu2=mu1+floor(mu/3)

    # Splitting up variance for for y1 conditional on y2,y3
    sigmac23_11=sigma_new[1:mu1,1:mu1]
    sigmac23_22=sigma_new[(mu1+1):mu,(mu1+1):mu]
    sigmac23_12=sigma_new[1:mu1,(mu1+1):mu]
    sigmac23_21=sigma_new[(mu1+1):mu,1:mu1]

    #Computing Cholesky decomposition for more precise inverse
    cholc23_22=chol(sigmac23_22)
    InvCholCorrc23_22=solve(cholc23_22)

    to_sub=sigmac23_12%*%InvCholCorrc23_22%*%t(InvCholCorrc23_22)%*%sigmac23_21
    sigma_c23=sigmac23_11-to_sub
    sigma_c23[lower.tri(sigma_c23)] = t(sigma_c23)[lower.tri(sigma_c23)]
    # is.positive.semi.definite(sigma_c23)
    #computing the conditional covariance

    #Splitting up covariance for y2 conditional on y3
    sigmac3_11=sigma_new[(mu1+1):mu2,(mu1+1):mu2]
    sigmac3_22=sigma_new[(mu2+1):mu,(mu2+1):mu]
    sigmac3_12=sigma_new[(mu1+1):mu2,(mu2+1):mu]
    sigmac3_21=sigma_new[(mu2+1):mu,(mu1+1):mu2]

    # computing cholesky decomposition for more precise inverse computation
    cholc3_22=chol(sigmac3_22)
    InvCholCorrc3_22=solve(cholc3_22)
    # Invc3_22=InvCholCorrc3_22%*%t(InvCholCorrc3_22)

    #computing the conditional variance
    sigma_c3=sigmac3_11-sigmac3_12%*%InvCholCorrc3_22%*%t(InvCholCorrc3_22)%*%sigmac3_21
    sigma_c3[lower.tri(sigma_c3)] = t(sigma_c3)[lower.tri(sigma_c3)]

    Y_post <- matrix(NA, nrow = length(Design_New_List[[1]]), ncol = N)

    maxy=vector()
    Maxx_List=list()
    for(n in 1:N){
      #sampling y3
      Y_post[(mu2+1):mu, n]=rmvnorm(1,mean=mu_new[(mu2+1):mu],sigma=sigmac3_22)

      #getting ready to sample y2/y3
      mu_c3=mu_new[(mu1+1):mu2]+sigmac3_12%*%InvCholCorrc3_22%*%t(InvCholCorrc3_22)%*%(Y_post[(mu2+1):mu,n]-mu_new[(mu2+1):mu])
      Y_post[(mu1+1):mu2,n]=rmvnorm(1,mean=mu_c3,sigma=sigma_c3)

      #now we have a sample from (y2,y3)
      #sampling from y1/y2y3
      mu_c23=mu_new[1:mu1]+sigmac23_12%*%InvCholCorrc23_22%*%t(InvCholCorrc23_22)%*%(Y_post[(mu1+1):mu, n]-mu_new[(mu1+1):mu])

      Y_post[1:mu1, n]=rmvnorm(1,mean=mu_c23,sigma=sigma_c23)
      indexio=which(Y_post[, n]==max(Y_post[, n]))
      for(k in 1:length(Design_New_List)){
        if(n==1){Maxx_List[[k]]=vector()}
        Maxx_List[[k]][n]=Design_New_List[[k]][indexio]
      }
      maxy[n]=Y_post[indexio,n]
    }
  }
  Maxx_List[[k+1]]=maxy
  return(Maxx_List)
}

# This function allows you to fit a Bayesian MSM
posterior_mean=function(x,XX,Y
                        ,thetas
                        ,Design_List
                        ,sigma2 #Note that this is only needed for the Interpolating model
                        ,alpha  #This is needed for the homoskedastic model
                        ,beta,InvCholCorr,Noise_Type,Covtype){
  X=XX

  Design_New_List=list()
  for(k in 1:length(Design_List)){
    Design_New_List[[k]]=x[k]
  }
  diagy=offdiag_compute(Design_List=Design_List,
                        Design_New_List=Design_New_List,
                        thetas=thetas,
                        Covtype=Covtype)

  corr_off_diag=diagy$corr_off_diag

  # print(length(corr_off_diag))
  # print(dim(InvCholCorr))
  # print(length(Y))
  # print(length(X))
  #note this alpha term is new and was missing before
    mu_new=t(diagy$x_new)%*%beta+alpha*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)

  return(mu_new)
}

multi_eval=function(
  DD
  ,y
  ,XX
  ,IthetasL
  ,IthetasU
  ,IalphaL=NULL
  ,IalphaU=NULL
  ,Noise_Type
  ,Distance_List
  ,Prior_List
  ,Prior_Der_List
  ,m
  ,Covtype){


  returny=apply(DD,1,FUN=multi_eval_sub1
                ,Noise_Type=Noise_Type
                ,y=y
                ,XX=XX
                ,IthetasL=IthetasL
                ,IthetasU=IthetasU
                ,IalphaL=IalphaL
                ,IalphaU=IalphaU
                ,Distance_List=Distance_List
                ,Prior_List=Prior_List
                ,Prior_Der_List=Prior_Der_List
                ,m=m
                ,Covtype=Covtype)
  NA_count=is.na(returny[1,])
  returny=returny[,!NA_count]
  valis=apply(returny,2, FUN=GP_likelihood,
              y=y,
              XX=XX,
              Distance_List=Distance_List,
              Prior_List=Prior_List,
              Prior_Der_List=Prior_Der_List,
              Noise_Type=Noise_Type,
              m=m,
              Covtype=Covtype)
  pari=returny[,which(valis==min(valis))[1]]
  return(pari)
}

multi_eval_sub1=function(pari,Noise_Type,y,XX,IthetasL,IthetasU,IalphaL=NULL,IalphaU=NULL
                         ,Distance_List, Prior_List,Prior_Der_List,m,Covtype){
  tryCatch({
    coin=rbinom(1,1,0.5)
    if(coin==1){
      pari=optim(par=pari
                 ,fn=GP_likelihood
                 ,gr=derloglik
                 ,method = "L-BFGS-B"
                 ,lower=c(IthetasL,IalphaL)
                 ,upper=c(IthetasU,IalphaU)
                 ,Distance_List=Distance_List
                 ,y=y
                 ,XX=XX
                 ,Prior_List=Prior_List
                 ,Prior_Der_List=Prior_Der_List
                 ,Noise_Type=Noise_Type
                 ,m=m
                 ,Covtype=Covtype
                 # ,control=list(parscale=c(thetas_lag,alpha_lag))
      )
    }else{
      pari=nlminb(start=pari
                  ,objective=GP_likelihood
                  ,gradient=derloglik
                  ,lower=c(IthetasL,IalphaL)
                  ,upper=c(IthetasU,IalphaU)
                  ,Distance_List=Distance_List
                  ,y=y
                  ,XX=XX
                  ,Prior_List=Prior_List
                  ,Prior_Der_List=Prior_Der_List
                  ,m=m
                  ,Noise_Type=Noise_Type
                  ,Covtype=Covtype
      )
    }
    to_return=pari$par
  },
  error = function(c){
    to_return=pari*NA
  })
  return(to_return)
}


plot_lik=function(thetas,
                  numb_varv,
                  alpha=NULL,
                  abliny,
                  xvals,
                  counter=0,
                  Distance_List,
                  y,
                  X,
                  Prior_List,
                  Prior_Der_List,
                  Noise_Type,
                  m,
                  Covtype){

  xaxis=xvals
  GP_likelihoodv=Vectorize(FUN=GP_likelihood, vectorize.args = c("varyv"), SIMPLIFY = TRUE)
  yout=GP_likelihoodv(varyv=xaxis,
                      thetasv=thetas,
                      alphav=alpha,
                      numb_varv=numb_varv,
                      indic=1,
                      Distance_List=Distance_List,
                      y=y,
                      X=X,
                      Prior_List=Prior_List,
                      Prior_Der_List=Prior_Der_List,
                      Noise_Type=Noise_Type,
                      m=m,
                      Covtype=Covtype)
  # print(yout)
  plot(xaxis,yout,main=counter, pch=20,xlab=expression(theta),ylab="Negative Log Likelihood")
  abline(v=abliny)
}
