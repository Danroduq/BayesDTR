rm(list=ls(all=TRUE))
set.seed(1)


library(LongCART)
library(dplyr)
library(reshape2)
#---------------------------------------------------------
# Processing LongCART data.
#---------------------------------------------------------
data("ACTG175")
group1=c(2) #this is the group that patients switch into
group2=c(1)
ACTG175_wide=reshape(ACTG175,idvar = c("pidnum"), timevar = "time", direction = "wide", v.names=c("cd4"))
ACTG175_wide=ACTG175_wide%>%filter(arms%in%c(group1,group2))
ACTG175_wide$z1=(ACTG175_wide$arms==2)+0
ACTG175_wide$z2=rbinom(1046,1,0.5)

ACTG175_wide$cd4.outcome=with(ACTG175_wide,0.2*(5*cd4.0+ 6*gender+ (-3000+9*cd4.0)*z1+(-3000+9*cd4.20)*z2))
sum(ACTG175_wide$cd4.outcome<0)
ACTG175_wide$cd4.outcome[ACTG175_wide$cd4.outcome<0]=0

BayesDat=ACTG175_wide
usethis::use_data(BayesDat,compress="xz")
