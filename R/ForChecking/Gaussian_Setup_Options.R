
covi=function(h,theta){
    (1+sqrt(5)*abs(h)/theta + 5*h^2/(3*theta^2))*exp(-sqrt(5)*abs(h)/theta)
}
dercovi=function(h,theta){

    5*h^2*(theta+sqrt(5)*abs(h))*exp(-sqrt(5)*abs(h)/theta)/(3*theta^4)

}

corr_compute=function(theta,distance){
  Corr=structure(vapply(distance, covi, numeric(1),theta=theta), dim=dim(distance))
  return(Corr)
}

dif_func=function(x,design,theta){
  distance_between=x-design
  to_return=covi(distance_between,theta=theta)
  return(to_return)
}


sigma2_compute=function(Corr,y,X){
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

      assign("Trigger", 1, env=globalenv())

  })
}
# thetas
# Design_List=list(design1,design2)
# Design_New_list=list(design1_new,design2_new)
# Distance_List=list(distance1,distance2)
# Expressi_List=list(expressi1,expressi2)
# Expressi_Der_List=list(Expressi1_der,Expressi2_der)
# Reduce("*", Distance_List)

GP_likelihood=function(params,Distance_List,y,indic=0, varyv,thetasv,alphav,numb_varv,Expressi_List,Expressi_Der_List,Noise_Type){
  # print(params)
  #params=c(0.5,0.5,0.5)
  if(indic==0){
    # theta1=params[1]; theta2=params[2]
    if(Noise_Type=="Homosk")alpha=params[3]
  }else{
    params=thetasv
    # theta1=params[1]; theta2=params[2]
    params[numb_varv]=varyv
    if(Noise_Type=="Homosk"){params=c(params,alphav);alpha=alphav}
  }

  Corr_List=list()
  for(i in 1:(length(Distance_List))){
    Corr_List[[i]]=corr_compute(params[i],Distance_List[[i]])
  }
  #Step 1:
  if(Noise_Type=="Interpol"){
    # Corr=corr_compute(theta1,covtype,distance1)*corr_compute(theta2,covtype,distance2)
    Corr=Reduce("*", Corr_List)
  }else if(Noise_Type=="Homosk"){
    # Corr=alpha*corr_compute(theta1,covtype,distance1)*corr_compute(theta2,covtype,distance2)+(1-alpha)*diag(rep(1,ncol(distance1)))
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
    out=out-2* eval(parse(text=Expressi_List[[i]]))
  }
  return(out)
}


derloglik=function(params,Distance_List,y,indic=0, thetasv,alphav,numb_varv,Expressi_List,Expressi_Der_List,Noise_Type){
  if(indic==0){
    # theta1=params[1]; theta2=params[2]
    if(Noise_Type=="Homosk")alpha=params[3]
  }else{
    params=thetasv
    # theta1=params[1]; theta2=params[2]
    params[numb_varv]=varyv
    if(Noise_Type=="Homosk"){params=c(params,alphav);alpha=alphav}
  }
  #-----------------------------#
  # Derivative of Log Likelihood
  #-----------------------------#
  Corr_List=list()
  for(i in 1:(length(Distance_List))){
    Corr_List[[i]]=corr_compute(params[i],Distance_List[[i]])
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
      gterm2_List[[k]]=structure(vapply(Distance_List[[k]], dercovi, numeric(1),theta=params[k]),dim=dim(Distance_List[[k]]))*Reduce("*", Corr_List[-k])
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
      dertheta_vect[k]=(-gterm3_vect[k]+gterm4_vect[k])-2*eval(parse(text=Expressi_Der_List[k]))
    }
    return(c(as.numeric(dertheta_vect)))
  }else if(Noise_Type=="Homosk"){
    for(k in 1:(length(Distance_List))){
      dertheta_vect[k]=alpha*(-gterm3_vect[k]+gterm4_vect[k])-2*eval(parse(text=Expressi_Der_List[k]))
    }
    deralpha=(-gterm3_alpha+gterm4_alpha)
    return(c(dertheta_vect, deralpha))
  }
}



offdiag_compute=function(Design_List,Design_New_List,thetas){
  x_new=matrix(model.matrix(~1, data=data.frame(Design_New_List[[1]])))

  R_New_List=list()
  Corr_Off_Diag_List=list()
  for(k in 1:length(Design_List)){
    distance_new=as.matrix(dist(Design_New_List[[k]]),method="manhattan")
    # distance2_new=as.matrix(dist(design2_new),method="manhattan")

    R_New_List[[k]]=structure(vapply(distance_new,covi,numeric(1),theta=thetas[k]),dim=dim(distance_new))
    # R1_new=structure(vapply(distance1_new, covi, numeric(1),theta=theta1,  covtype=covtype), dim=dim(distance1_new))
    # R2_new=structure(vapply(distance2_new, covi, numeric(1),theta=theta2,  covtype=covtype), dim=dim(distance2_new))

    Corr_Off_Diag_List[[k]]=sapply(Design_New_List[[k]],FUN=dif_func,design=Design_List[[k]],theta=thetas[k])
    # corr_off_diag1=sapply(design1_new,FUN=dif_func, design=design1, theta=theta1,covtype=covtype)
    # corr_off_diag2=sapply(design2_new,FUN=dif_func, design=design2, theta=theta2,covtype=covtype)
  }
  R_new=Reduce("*",R_New_List)
  corr_off_diag=Reduce("*",Corr_Off_Diag_List)

  # R_new=R1_new*R2_new
  # corr_off_diag=corr_off_diag1*corr_off_diag2
  return(list("x_new"=x_new,"corr_off_diag"=corr_off_diag,"R_new"=R_new))
}

EI=function(x,Y,ymax,Design_List,thetas,sigma2,beta,InvCholCorr){

  Design_New_List=list()
  for(k in 1:length(Design_List)){
    Design_New_List[[k]]=x[k]
  }
  # Design_New_Lect=x
  # design1_new=x[1]
  # design2_new=x[2]

  beta=matrix(beta) #making sure its column wise vector
  diagy=offdiag_compute(Design_List=Design_List,Design_New_List=Design_New_List,thetas=thetas)
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



CIdraw_posterior<- function(N,Y,
                            Design_List,
                            InvCholCorr,
                            theta,
                            sig2,
                            alpha,
                            v,
                            # Design_New_List,
                            checki=FALSE,Noise_Type){
  #Y=y;N=100;design1=design1; design2=design2;InvCholCorr=InvCholCorr;theta1=theta1;theta2=theta2;sig2=as.numeric(sigma2)
  # design1_new=seq(50,100,4)
  # design2_new=seq(200,600,7.5)
  ###########################################
  # Finding the maximizer of the max function
  ###########################################

  thetay=as.matrix(expand.grid(design1_new,design2_new))
  Design_New_List=list()
  for(k in 1:length(Design_List)){
    Design_New_List[[k]]=thetay[,k]
  }
  # design1_new=thetay[,1]
  # design2_new=thetay[,2]

  beta=matrix(beta) #making sure its column wise vector
  diagy=offdiag_compute(Design_List=Design_List,Design_New_List=Design_New_List,thetas=thetas)
  corr_off_diag=diagy$corr_off_diag

  if(Noise_Type=="Interpol"){
    mu_new=diagy$x_new%*%beta+t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)
    sigma_new=sig2*(diagy$R_new-t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%corr_off_diag)
  }else if(Noise_Type=="Homosk"){
    mu_new=diagy$x_new%*%beta+alpha*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)
    sigma_new=sig2*(diagy$R_new-sig2*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%corr_off_diag/v)
  }else if(Noise_Type=="Heterosk"){
    InvCholVar=InvCholCorr
    mu_new=x_new%*%beta+sig2*t(corr_off_diag)%*%InvCholVar%*%t(InvCholVar)%*%(matrix(Y)-X%*%beta)
    sigma_new=sig2*(diagy$R_new-sig2*t(corr_off_diag)%*%InvCholVar%*%t(InvCholVar)%*%corr_off_diag)
  }

  sigma_new[lower.tri(sigma_new)] = t(sigma_new)[lower.tri(sigma_new)]

  index=which(mu_new==max(mu_new))
  max1_mu=design1_new[index]
  max2_mu=design2_new[index]

  checki2=is.positive.definite(sigma_new)
  Drawy=CIdraw_posterior_sub1(mu_new,sigma_new, design1_new, design2_new,checki2)
  if(checki==FALSE){
    quantix1=quantile(Drawy$maxx1, p=c(0.025,0.5,0.975))
    quantix2=quantile(Drawy$maxx2, p=c(0.025,0.5,0.975))
    quantiy=quantile(Drawy$maxy, p=c(0.025,0.5,0.975))

    LBx1=quantix1[1];Midx1=quantix1[2];UBx1=quantix1[3]
    LBx2=quantix2[1];Midx2=quantix2[2];UBx2=quantix2[3]
    LBy=quantiy[1];Midy=quantiy[2];UBy=quantiy[3]

    return(c(LBx1,Midx1,UBx1,max1_mu,LBx2,Midx2,UBx2,max2_mu,LBy,Midy,UBy,max(mu_new)))
  }else{
    return(list("distx1"=Drawy$maxx1,"distx2"=Drawy$maxx2,"disty"=Drawy$maxy))
  }
}

CIdraw_posterior_sub1=function(mu_new,sigma_new, Design_New_List,checki2){
  if(checki2==FALSE){
    #sampling full normal vector
    Y_post <- matrix(NA, nrow = length(Design_New_List[[k]]), ncol = N)

    Maxx_List=List()
    maxx1=vector()
    maxx2=vector()
    maxy=vector()
    for (n in 1:N) {
      Y_post[, n] <-rmvnorm(1,mean=mu_new,sigma=sigma_new)
      indexio=which(Y_post[, n]==max(Y_post[, n]))
      for(k in 1:length(Design_New_List)){
        Maxx_List[[k]][n]=Design_New_List[[k]][indexio]
      }
      maxy[n]=Y_post[indexio,n]
    }
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
    maxx1=vector()
    maxx2=vector()
    maxy=vector()
    for(n in 1:N){
      #sampling y3
      Y_post[(mu2+1):mu, n]=rmvnorm(1,mean=mu_new[(mu2+1):mu],sigma=sigmac3_22)

      #getting ready to sample y2/y3
      mu_c3=mu_new[(mu1+1):mu2]+sigmac3_12%*%InvCholCorrc3_22%*%t(InvCholCorrc3_22)%*%(Y_post[(mu2+1):mu, n]-mu_new[(mu2+1):mu])
      Y_post[(mu1+1):mu2, n]=rmvnorm(1, mean=mu_c3, sigma=sigma_c3)

      #now we have a sample from (y2,y3)
      #sampling from y1/y2y3
      mu_c23=mu_new[1:mu1]+sigmac23_12%*%InvCholCorrc23_22%*%t(InvCholCorrc23_22)%*%(Y_post[(mu1+1):mu, n]-mu_new[(mu1+1):mu])

      Y_post[1:mu1, n]=rmvnorm(1,mean=mu_c23,sigma=sigma_c23)
      Maxx_List=List()
      indexio=which(Y_post[, n]==max(Y_post[, n]))
      for(k in 1:length(Design_New_List)){
        Maxx_List[[k]][n]=Design_New_List[[k]][indexio]
      }
      # maxx1[n]=design1_new[indexio]
      # maxx2[n]=design2_new[indexio]
      maxy[n]=Y_post[indexio,n]
    }
  }
  Maxx_List[[k+1]]=maxy
  return(Maxx_List)
}



posterior_mean=function(x,XX,Y
                        ,thetas
                        ,Design_List
                        ,sigma2 #Note that this is only needed for the Interpolating model
                        ,alpha  #This is needed for the homoskedastic model
                        ,beta,InvCholCorr,Noise_Type){
  #X=as.matrix(X[1:m,]);Y=y[1:m];theta1=theta1_converge;theta2=theta2_converge;alpha=alpha_converge;beta=beta;InvCholCorr=InvCholCorr;design1=design1[1:m];design2= design2[1:m]
  X=XX
  # design1_new=x[1]
  # design2_new=x[2]

  Design_New_List=list()
  for(k in 1:length(Design_List)){
    Design_New_List[[k]]=x[k]
  }
  diagy=offdiag_compute(
                        Design_List=Design_List,
                        Design_New_List=Design_New_List,
                        thetas=thetas)

  corr_off_diag=diagy$corr_off_diag

  #note this alpha term is new and was missing before
  if(Noise_Type=="Interpol"){
    mu_new=t(diagy$x_new)%*%beta+t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(Y-X%*%beta)
  }else if(Noise_Type=="Homosk"){
    mu_new=t(diagy$x_new)%*%beta+alpha*t(corr_off_diag)%*%InvCholCorr%*%t(InvCholCorr)%*%(matrix(Y)-X%*%beta)
  }else if(Noise_Type=="Heterosk"){
    InvCholVar=InvCholCorr
    mu_new=t(diagy$x_new)%*%beta+sig2*t(corr_off_diag)%*%InvCholVar%*%t(InvCholVar)%*%(matrix(Y)-X%*%beta)
  }
  return(mu_new)
}
#
# derloglik(c(5.400346,0.364818)
#           ,covtype=covtype
#           ,Distance_List=Distance_List
#           ,y=y
#           ,Expressi_List=Expressi_List
#           ,Expressi_Der_List=Expressi_Der_List
#           ,Noise_Type=Noise_Type)

multi_eval_sub1=function(pari,Noise_Type,y,IthetasL,IthetasU,IalphaL=NULL,IalphaU=NULL,thetas_lag=NULL,alpha_lag=NULL){
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
                 ,Expressi_List=Expressi_List
                 ,Expressi_Der_List=Expressi_Der_List
                 ,Noise_Type=Noise_Type
                 ,control=list(parscale=c(thetas_lag,alpha_lag))
      )
    }else{
      pari=nlminb(start=pari
                  ,objective=GP_likelihood
                  ,gradient=derloglik
                  ,lower=c(IthetasL,IalphaL)
                  ,upper=c(IthetasU,IalphaU)
                  ,Distance_List=Distance_List
                  ,y=y
                  ,Expressi_List=Expressi_List
                  ,Expressi_Der_List=Expressi_Der_List
                  ,Noise_Type=Noise_Type
      )
    }
    to_return=pari$par
  },
  error = function(c){
    to_return=pari*NA
  })
  return(to_return)
}



multi_eval=function(
  # U1,U2,U3=NULL
  DD,
  y,
  IthetasL,
  IthetasU,
  IalphaL=NULL,
  IalphaU=NULL,
  Noise_Type,
  thetas_lag,
  alpha_lag=NULL){
  # U1=runif(15,0.01,Itheta1U);
  # U2=runif(15,0.01,Itheta2U);
  # U3=runif(15,0.01,0.999);
  # y=z;
  # IthetaL=IthetaL;
  # Itheta1U=Itheta1U;
  # Itheta2U=Itheta2U;
  # IalphaL=IalphaL;
  # IalphaU=IalphaU;
  # Noise_Type="Homosk";
  # theta1_lag=theta1_z;theta2_lag=theta2_z;alpha_lag=alpha_z

  # DD=cbind(U1,U2,U3) # when U3 is Null nothing extra is linked
  returny=apply(DD,1,FUN=multi_eval_sub1,
                Noise_Type=Noise_Type,
                y=y,
                IthetasL=IthetasL,
                IthetasU=IthetasU,
                IalphaL=IalphaL,
                IalphaU=IalphaU,
                thetas_lag=thetas_lag,
                alpha_lag=alpha_lag)
  NA_count=is.na(returny[1,])
  print(sum(NA_count))
  returny=returny[,!NA_count]
  valis=apply(returny,2, FUN=GP_likelihood,
              Distance_List=Distance_List,
              y=y,
              Expressi_List=Expressi_List,
              Expressi_Der_List=Expressi_Der_List,
              Noise_Type=Noise_Type)

  pari=returny[,which(valis==min(valis))[1]]
  return(pari)
}

# xvals1=seq(0.01,10,0.1)
# xvals2=seq(thetas[2],thetas[2])
# xvals_List=list(xvals1,xvals2)
plot_lik=function(thetas,
                  numb_varv,
                  alpha=NULL,
                  abliny,
                  xvals,
                  counter=0,
                  Distance_List,
                  y,
                  Expressi_List,
                  Expressi_Der_List,
                  Noise_Type){

  xaxis=xvals
  # domain=expand.grid(xvals_List)
  #
  # length_checker=lapply(xvals_List,length)
  # length_checker_index=which(length_checker>1)
  #
  # xaxis=domain[,length_checker_index]
  # if(length(xvals1)>1){xaxis=domain$theta1
  # }else{xaxis=domain$theta2}


  # varyv=paste0("theta",numb_varv,"v")
  GP_likelihoodv=Vectorize(FUN=GP_likelihood, vectorize.args = c("varyv"), SIMPLIFY = TRUE)
  yout=GP_likelihoodv(varyv=xaxis,
                      thetasv=thetas,
                      alphav=alpha,
                      numb_varv=numb_varv,
                      indic=1,
                      Distance_List=Distance_List,
                      y=y,
                      Expressi_List=Expressi_List,
                      Expressi_Der_List=Expressi_Der_List,
                      Noise_Type=Noise_Type)
  # print(yout)
  plot(xaxis,yout,main=counter, pch=20,xlab=expression(theta),ylab="Negative Log Likelihood")
  abline(v=abliny)
}

# # ---------------------
# #3-D plot
# #----------------------
# xvals1=seq(0.01,0.5,0.01);xvals2=seq(0.01,0.2,0.01);counter=0;
#
# domain=expand.grid(theta1=xvals1, theta2=xvals2)
# domain=cbind(domain,alpha)
# yy=apply(domain,1
#         ,FUN=GP_likelihood
#         ,covtype=covtype
#         ,Distance_List=Distance_List
#         ,y=y
#         ,Expressi_List=Expressi_List
#         ,Expressi_Der_List=Expressi_Der_List
#         ,Noise_Type=Noise_Type)
#
# data_grid2=data.frame(theta1=domain[,1],
#                       theta2=domain[,2],
#                       value=yy)
# plot_matrix2 <- t(acast(data_grid2, theta1~theta2, value.var="value"))
#
# plot_ly(
#   x = as.numeric(colnames(plot_matrix2)),
#   y = as.numeric(rownames(plot_matrix2)),
#   z = plot_matrix2) %>%
#   add_surface()%>%
#   layout(title="Value of Regime for Varying Cutoffs \n IPW",
#          scene=list(xaxis = list(title="Psi_1"),
#                     yaxis = list(title="Psi_2"),
#                     zaxis =  list(title="Y")))
#
#
# ####################################
# # 1-D
# # First version
# ####################################
# xaxis=seq(0.01,0.1,0.001)
# GP_likelihoodv=Vectorize(FUN=GP_likelihood, vectorize.args = c("varyv"), SIMPLIFY = TRUE)
# yout=GP_likelihoodv(varyv=xaxis,
#                     thetasv=thetas,
#                     alphav=alpha,
#                     numb_varv=2,
#                     indic=1,
#                     covtype=covtype,
#                     Distance_List=Distance_List,
#                     y=y,
#                     Expressi_List=Expressi_List,
#                     Expressi_Der_List=Expressi_Der_List,
#                     Noise_Type=Noise_Type)
# thetas
# par(mfrow=c(1,2))
# plot(xaxis,yout,main=counter, pch=20,xlab=expression(theta),ylab="Negative Log Likelihood")
# abline(v=thetas)
####################################
#1-D
#2nd version
####################################
# xvals=seq(0.01,0.3,0.001);counter=0;
#
# domain=cbind(thetas[1],xvals,alpha)
# yy=apply(domain,1
#          ,FUN=derloglik #GP_likelihood
#          ,covtype=covtype
#          ,Distance_List=Distance_List
#          ,y=y
#          ,Expressi_List=Expressi_List
#          ,Expressi_Der_List=Expressi_Der_List
#          ,Noise_Type=Noise_Type)
# plot(xvals,yy[2,])
#
#
# derloglik(params=c(0.3390568,0.1353130,0.9999999)
#       ,covtype=covtype
#       ,Distance_List=Distance_List
#       ,y=y
#       ,Expressi_List=Expressi_List
#       ,Expressi_Der_List=Expressi_Der_List
#       ,Noise_Type=Noise_Type)

# ###
# returny=apply(DD,1,FUN=multi_eval_sub1,
#               Noise_Type=Noise_Type,
#               y=y,
#               IthetasL=IthetasL,
#               IthetasU=IthetasU,
#               IalphaL=IalphaL,
#               IalphaU=IalphaU,
#               thetas_lag=thetas_lag,
#               alpha_lag=alpha_lag)
# NA_count=is.na(returny[1,])
# print(sum(NA_count))
# returny=returny[,!NA_count]
# valis=apply(returny,2, FUN=GP_likelihood,
#             covtype=covtype,
#             Distance_List=Distance_List,
#             y=y,
#             Expressi_List=Expressi_List,
#             Expressi_Der_List=Expressi_Der_List,
#             Noise_Type=Noise_Type)
#
# returny[,which(valis==min(valis))[1]]


