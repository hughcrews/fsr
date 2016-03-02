#####################################################################################################
#This function implements the Fast FSR methodology for handling interactions  
#
#x: matrix containing covariates
#   Note: Only main effects need to be included as function handles interactions internally
#
#y: vector containing response variable
#
#method: "main" = main effects only (if quad=T supplied then quadratic terms also considered)
#        "strong" = strong hierarchy
#        "weak" = weak hierarchy
#        "weakadj" = weak hierarchy with adjusted pvalues
#        "nohier" = no hierarchy
#        "nohieradj" = no hierarchy with adjusted pvalues
#
#gamma0: desired false selection rate (set to 0.05 by default)
#
#center: by default covariates are centered (can be turned off using center=F)
#
#quad: by default quadratic terms are not included but if quad=T provided then quadratic 
#      terms are handled based on the particular hierarchy approach specified
#
#maxsteps: maximum number of steps in Forward addition sequence to complete
#          by default the sequence continues until either all terms are selected or the df are spent
#
#print: option of whether to print default FSR results (turned on by default)
#
#digits: number of digits for rounding of printed results (set to 4 by default)
#####################################################################################################




fsr.int<-function(x,y,method="strong",gamma0=.05,center=T,quad=F,maxsteps=kt,print=T,digits=4){

require(MASS)

#Check to see if a proper method was specified
mchk<- method %in% c("main","strong","weak","weakadj","nohier","nohieradj")

if(mchk==F) 

{
print("Please specify a method on the list: main, strong, weak, weakadj, nohier, nohieradj")
return(NULL)
}


ok<-complete.cases(x,y)
x<-x[ok,]                            # get rid of na's
y<-y[ok]                             
p<-ncol(x)
n<-nrow(x)
x<-as.matrix(x)                      # in case x is not a matrix
one <- rep(1, n)

#Center covariates by default
if(center) 
{
meanx <- drop(one %*% x)/n
x <- scale(x, meanx, FALSE)
}

xs<-x 

xsone<-cbind(one,xs)




#Code for main effects (and quadratic terms) method
if(method=="main")

{
if(quad) kt=2*p else kt=p


maxstep<-min(n-2,kt,maxsteps)

muhat<-mean(y)


#Set up vector to keep track of main effects in model
maininmodel<-matrix(0,p,1)
if(quad) intinmodel<-matrix(0,p,1)


#Set up vectors to keep the forward addition sequence info
var.step=matrix("NA",maxstep+1,1)
pval.step=matrix(NA,maxstep+1,1)
maxpval.step=matrix(NA,maxstep+1,1)
ghat.step=matrix(NA,maxstep+1,1)
ghigh.step=matrix(NA,maxstep+1,1)
totcand.step=matrix(NA,maxstep+1,1)
alphahat.step=matrix(NA,maxstep+1,1)
rsquare.step=matrix(NA,maxstep+1,1)



#Initialize values for Step 0 (Intercept)
var.step[1]="Intercept"
pval.step[1]=0
maxpval.step[1]=0
ghat.step[1]=0
totcand.step[1]=kt
alphahat.step[1]=0
rsquare.step[1]=0




i=1

#Begin loop to get each step of forward addition sequence

while(i<=maxstep) {




#Figure out which main effects to consider
mainind<-1-maininmodel
maincand<-sum(mainind)

if(quad)
{
intind<-1-intinmodel
intcand<-sum(intind)

}


j<-1
mcnt<-0
while(mcnt<maincand) {

if(mainind[j]==1 & mcnt==0)  {
 xnew<-xs[,j]
 candvars<-paste("x",deparse(j),sep="")
 mainspot<-j
}
if(mainind[j]==1 & mcnt >0)  {
 xnew<-cbind(xnew,xs[,j])
      candvars<-cbind(candvars,paste("x",deparse(j),sep=""))
      mainspot<-cbind(mainspot,j)
}
if(mainind[j]==1) mcnt<-mcnt+1

j<-j+1
}


#This code handles quadratic terms if requested

if(quad==T)


{

k<-1
intcnt<-0

while(k<=p)

{


  if(intinmodel[k]==1) intind[k]=0

  if( intinmodel[k]==0){
   if(maincand==0 & intcnt==0) xnew<-xs[,k]*xs[,k] else xnew<-cbind(xnew,xs[,k]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(k),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(k),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-k else intspot<-rbind(intspot,k)
   intcnt<-intcnt+1
  }


  k<-k+1
 }

}

 

#Count number of candidate terms
if(quad) totcand<-maincand+intcand 
 else totcand<-maincand

#Get new residuals
rstep<-y-muhat

#Get correlations to pick term to enter
if(i>1) rxnew<-residuals(lm(xnew~xcand)) else rxnew<-xnew


scorr<-abs(cor(rstep,rxnew))
select<-which.max(scorr)
selectvar<-candvars[select]
xnew<-as.matrix(xnew)
if (i<kt) selectx<-xnew[,select] else selectx<-xnew

if(i==1) xcand<-cbind(one,selectx) else xcand<-cbind(xcand,selectx)


if(select<=maincand) maininmodel[mainspot[select]]<-1  
else intinmodel[intspot[select-maincand]]<-1



invxpx<-ginv(crossprod(xcand))
tcandy<-crossprod(xcand,y)
beta<-invxpx%*%tcandy
muhat<-xcand%*%beta
r<-y-muhat
sigmahat<-sum(r*r)/(n-i-1)
sxx<-invxpx*sigmahat
t<-beta[i+1]/(sqrt(sxx[i+1,i+1]))
pvalue<-2*(1-pt(abs(t),n-i-1))



var.step[i+1]<-selectvar
pval.step[i+1]<-pvalue
ok<-complete.cases(pval.step)
maxpval.step[i+1]<-max(pval.step[ok],pval.step[i+1])
if (quad) totcand.step[i+1]<-2*p-sum(maininmodel)-sum(intinmodel)
 else totcand.step[i+1]<-p-sum(maininmodel)

if(quad) totinmodel<-sum(maininmodel)+sum(intinmodel) else totinmodel<-sum(maininmodel)

ghat.step[i+1]<-(totcand.step[i+1]*maxpval.step[i+1])/(1+totinmodel)
ghigh.step[i]<-(totcand.step[i]*maxpval.step[i+1])/(totinmodel)
if(i==maxstep) ghigh.step[i+1]<-(totcand.step[i+1]*1)/(1+totinmodel)

alphahat.step[i+1]<-((1+totinmodel)*gamma0)/(totcand.step[i+1])

rsquare.step[i+1]<-1-sum(r*r)/sum((y-mean(y))*(y-mean(y)))


if(i==1) beta.step<-rbind(mean(y),beta) else beta.step<-rbind(beta.step,beta)


i<-i+1 
}

}





#Code for strong hierarchy
if(method=="strong")

{
if(quad) kt=2*p+choose(p,2) else kt=p+choose(p,2)


maxstep<-min(n-2,kt,maxsteps)

muhat<-mean(y)


#Set up vector to keep track of main effects in model
maininmodel<-matrix(0,p,1)


#Set up matrices to keep track of which interactions are candidates and which are in the model
intind<-matrix(0,p,p)
intinmodel<-matrix(0,p,p)


#Set up vectors to keep the forward addition sequence info
var.step=matrix("NA",maxstep+1,1)
pval.step=matrix(NA,maxstep+1,1)
maxpval.step=matrix(NA,maxstep+1,1)
ghat.step=matrix(NA,maxstep+1,1)
ghigh.step=matrix(NA,maxstep+1,1)
totcand.step=matrix(NA,maxstep+1,1)
alphahat.step=matrix(NA,maxstep+1,1)
rsquare.step=matrix(NA,maxstep+1,1)



#Initialize values for Step 0 (Intercept)
var.step[1]="Intercept"
pval.step[1]=0
maxpval.step[1]=0
ghat.step[1]=0
totcand.step[1]=p
alphahat.step[1]=0
rsquare.step[1]=0




i=1

#Begin loop to get each step of forward addition sequence

while(i<=maxstep) {





#Figure out which main effects to consider
mainind<-1-maininmodel
maincand<-sum(mainind)


j<-1
mcnt<-0
while(mcnt<maincand) {

if(mainind[j]==1 & mcnt==0)  {
 xnew<-xs[,j]
 candvars<-paste("x",deparse(j),sep="")
 mainspot<-j
}
if(mainind[j]==1 & mcnt >0)  {
 xnew<-cbind(xnew,xs[,j])
      candvars<-cbind(candvars,paste("x",deparse(j),sep=""))
      mainspot<-cbind(mainspot,j)
}
if(mainind[j]==1) mcnt<-mcnt+1

j<-j+1
}

 

#Figure out which interactions/quadratic terms to consider
#Code depends on whether quadratic terms are allowed using quad call in function


#This code only allows interactions

if(quad==F)
{
j<-1
k<-j+1
intcnt<-0


while(j<p)
{
 while(k<=p){


  if(maininmodel[j]==1 & maininmodel[k]==1 & intinmodel[j,k]==0){
   intind[j,k]=1
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j+1
}

}


#This code handles quadratic terms with interactions if requested

if(quad==T)


{
j<-1
k<-j
intcnt<-0

while(j<=p)
{
 while(k<=p){


  if(maininmodel[j]==1 & maininmodel[k]==1 & intinmodel[j,k]==0){
   intind[j,k]=1
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j
}
}

#Count number of candidate terms
intcand<-sum(intind)
totcand<-maincand+intcand

#Get new residuals
rstep<-y-muhat

#Get correlations to pick term to enter
if(i>1) rxnew<-residuals(lm(xnew~xcand)) else rxnew<-xnew

scorr<-abs(cor(rstep,rxnew))
select<-which.max(scorr)
selectvar<-candvars[select]
xnew<-as.matrix(xnew)
if (i<kt) selectx<-xnew[,select] else selectx<-xnew

if(i==1) xcand<-cbind(one,selectx) else xcand<-cbind(xcand,selectx)


if(select<=maincand) maininmodel[mainspot[select]]<-1  
else intinmodel[intspot[select-maincand,1],intspot[select-maincand,2]]<-1




invxpx<-ginv(crossprod(xcand))
tcandy<-crossprod(xcand,y)
beta<-invxpx%*%tcandy
muhat<-xcand%*%beta
r<-y-muhat
sigmahat<-sum(r*r)/(n-i-1)
sxx<-invxpx*sigmahat
t<-beta[i+1]/(sqrt(sxx[i+1,i+1]))
pvalue<-2*(1-pt(abs(t),n-i-1))



var.step[i+1]<-selectvar
pval.step[i+1]<-pvalue
ok<-complete.cases(pval.step)
maxpval.step[i+1]<-max(pval.step[ok],pval.step[i+1])
totcand.step[i+1]<-p+intcand-sum(maininmodel)-sum(intinmodel)


ghat.step[i+1]<-(totcand.step[i+1]*maxpval.step[i+1])/(1+sum(maininmodel)+sum(intinmodel))
ghigh.step[i]<-(totcand.step[i]*maxpval.step[i+1])/(1+sum(maininmodel)+sum(intinmodel)-1)
if(i==maxstep) ghigh.step[i+1]<-(totcand.step[i+1]*1)/(1+sum(maininmodel)+sum(intinmodel))

alphahat.step[i+1]<-((1+sum(maininmodel)+sum(intinmodel))*gamma0)/(totcand.step[i+1])

rsquare.step[i+1]<-1-sum(r*r)/sum((y-mean(y))*(y-mean(y)))


if(i==1) beta.step<-rbind(mean(y),beta) else beta.step<-rbind(beta.step,beta)


i<-i+1 
}


}



#Code for weak hierarchy
if(method=="weak")

{
if(quad) kt=2*p+choose(p,2) else kt=p+choose(p,2)


maxstep<-min(n-2,kt,maxsteps)

muhat<-mean(y)


#Set up vector to keep track of main effects in model
maininmodel<-matrix(0,p,1)


#Set up matrices to keep track of which interactions are candidates and which are in the model
intind<-matrix(0,p,p)
intinmodel<-matrix(0,p,p)


#Set up vectors to keep the forward addition sequence info
var.step=matrix("NA",maxstep+1,1)
pval.step=matrix(NA,maxstep+1,1)
maxpval.step=matrix(NA,maxstep+1,1)
ghat.step=matrix(NA,maxstep+1,1)
ghigh.step=matrix(NA,maxstep+1,1)
totcand.step=matrix(NA,maxstep+1,1)
alphahat.step=matrix(NA,maxstep+1,1)
rsquare.step=matrix(NA,maxstep+1,1)



#Initialize values for Step 0 (Intercept)
var.step[1]="Intercept"
pval.step[1]=0
maxpval.step[1]=0
ghat.step[1]=0
totcand.step[1]=p
alphahat.step[1]=0
rsquare.step[1]=0




i=1

#Begin loop to get each step of forward addition sequence

while(i<=maxstep) {





#Figure out which main effects to consider
mainind<-1-maininmodel
maincand<-sum(mainind)


j<-1
mcnt<-0
while(mcnt<maincand) {

if(mainind[j]==1 & mcnt==0)  {
 xnew<-xs[,j]
 candvars<-paste("x",deparse(j),sep="")
 mainspot<-j
}
if(mainind[j]==1 & mcnt >0)  {
 xnew<-cbind(xnew,xs[,j])
      candvars<-cbind(candvars,paste("x",deparse(j),sep=""))
      mainspot<-cbind(mainspot,j)
}
if(mainind[j]==1) mcnt<-mcnt+1

j<-j+1
}

 

#Figure out which interactions/quadratic terms to consider
#Code depends on whether quadratic terms are allowed using quad call in function


#This code only allows interactions

if(quad==F)
{
j<-1
k<-j+1
intcnt<-0


while(j<p)
{
 while(k<=p){


  if((maininmodel[j]==1 | maininmodel[k]==1) & intinmodel[j,k]==0){
   intind[j,k]=1
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j+1
}

}


#This code handles quadratic terms with interactions if requested

if(quad==T)


{
j<-1
k<-j
intcnt<-0

while(j<=p)
{
 while(k<=p){


  if((maininmodel[j]==1 | maininmodel[k]==1) & intinmodel[j,k]==0){
   intind[j,k]=1
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j
}
}

#Count number of candidate terms
intcand<-sum(intind)
totcand<-maincand+intcand

#Get new residuals
rstep<-y-muhat

#Get correlations to pick term to enter
if(i>1) rxnew<-residuals(lm(xnew~xcand)) else rxnew<-xnew

scorr<-abs(cor(rstep,rxnew))
select<-which.max(scorr)
selectvar<-candvars[select]
xnew<-as.matrix(xnew)
if (i<kt) selectx<-xnew[,select] else selectx<-xnew

if(i==1) xcand<-cbind(one,selectx) else xcand<-cbind(xcand,selectx)


if(select<=maincand) maininmodel[mainspot[select]]<-1  
else intinmodel[intspot[select-maincand,1],intspot[select-maincand,2]]<-1


invxpx<-ginv(crossprod(xcand))
tcandy<-crossprod(xcand,y)
beta<-invxpx%*%tcandy
muhat<-xcand%*%beta
r<-y-muhat
sigmahat<-sum(r*r)/(n-i-1)
sxx<-invxpx*sigmahat
t<-beta[i+1]/(sqrt(sxx[i+1,i+1]))
pvalue<-2*(1-pt(abs(t),n-i-1))


var.step[i+1]<-selectvar
pval.step[i+1]<-pvalue
ok<-complete.cases(pval.step)
maxpval.step[i+1]<-max(pval.step[ok],pval.step[i+1])
totcand.step[i+1]<-p+intcand-sum(maininmodel)-sum(intinmodel)


ghat.step[i+1]<-(totcand.step[i+1]*maxpval.step[i+1])/(1+sum(maininmodel)+sum(intinmodel))
ghigh.step[i]<-(totcand.step[i]*maxpval.step[i+1])/(1+sum(maininmodel)+sum(intinmodel)-1)
if(i==maxstep) ghigh.step[i+1]<-(totcand.step[i+1]*1)/(1+sum(maininmodel)+sum(intinmodel))

alphahat.step[i+1]<-((1+sum(maininmodel)+sum(intinmodel))*gamma0)/(totcand.step[i+1])

rsquare.step[i+1]<-1-sum(r*r)/sum((y-mean(y))*(y-mean(y)))


if(i==1) beta.step<-rbind(mean(y),beta) else beta.step<-rbind(beta.step,beta)


i<-i+1 
}


}





#Code for adjusted weak hierarchy
if(method=="weakadj")

{
if(quad) kt=2*p+choose(p,2) else kt=p+choose(p,2)


maxstep<-min(n-2,kt,maxsteps)

muhat<-mean(y)


#Set up vector to keep track of main effects in model
maininmodel<-matrix(0,p,1)


#Set up matrices to keep track of which interactions are candidates and which are in the model
intind<-matrix(0,p,p)
intinmodel<-matrix(0,p,p)


#Set up vectors to keep the forward addition sequence info
var.step=matrix("NA",maxstep+1,1)
pval.step=matrix(NA,maxstep+1,1)
maxpval.step=matrix(NA,maxstep+1,1)
ghat.step=matrix(NA,maxstep+1,1)
ghigh.step=matrix(NA,maxstep+1,1)
totcand.step=matrix(NA,maxstep+1,1)
alphahat.step=matrix(NA,maxstep+1,1)
rsquare.step=matrix(NA,maxstep+1,1)
cadj.step=matrix(NA,maxstep+1,1)
maininmod.step=matrix(NA,maxstep+1,1)
intinmod.step=matrix(NA,maxstep+1,1)


#Initialize values for Step 0 (Intercept)
var.step[1]="Intercept"
pval.step[1]=0
maxpval.step[1]=0
ghat.step[1]=0
totcand.step[1]=p
alphahat.step[1]=0
rsquare.step[1]=0
cadj.step[1]=1
maininmod.step[1]=0
intinmod.step[1]=0



i=1

#Begin loop to get each step of forward addition sequence

while(i<=maxstep) {





#Figure out which main effects to consider
mainind<-1-maininmodel
maincand<-sum(mainind)


j<-1
mcnt<-0
while(mcnt<maincand) {

if(mainind[j]==1 & mcnt==0)  {
 xnew<-xs[,j]
 candvars<-paste("x",deparse(j),sep="")
 mainspot<-j
}
if(mainind[j]==1 & mcnt >0)  {
 xnew<-cbind(xnew,xs[,j])
      candvars<-cbind(candvars,paste("x",deparse(j),sep=""))
      mainspot<-cbind(mainspot,j)
}
if(mainind[j]==1) mcnt<-mcnt+1

j<-j+1
}

 

#Figure out which interactions/quadratic terms to consider
#Code depends on whether quadratic terms are allowed using quad call in function


#This code only allows interactions

if(quad==F)
{
j<-1
k<-j+1
intcnt<-0


while(j<p)
{
 while(k<=p){


  if((maininmodel[j]==1 | maininmodel[k]==1) & intinmodel[j,k]==0){
   intind[j,k]=1
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j+1
}

}


#This code handles quadratic terms with interactions if requested

if(quad==T)


{
j<-1
k<-j
intcnt<-0

while(j<=p)
{
 while(k<=p){


  if((maininmodel[j]==1 | maininmodel[k]==1) & intinmodel[j,k]==0){
   intind[j,k]=1
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j
}
}

#Count number of candidate terms
intcand<-sum(intind)
totcand<-maincand+intcand

#Get new residuals
rstep<-y-muhat

#Get correlations to pick term to enter
if(i>1) rxnew<-residuals(lm(xnew~xcand)) else rxnew<-xnew

scorr<-abs(cor(rstep,rxnew))

xnew<-as.matrix(xnew)


if (maincand>0) 
{
mainscorr<-scorr[1:maincand]
select1<-which.max(mainscorr)
selectvar1<-candvars[select1]


if (i<kt) selectx1<-xnew[,select1] else selectx1<-xnew
if(i==1) xcand1<-cbind(one,selectx1) else xcand1<-cbind(xcand,selectx1)



invxpx1<-ginv(crossprod(xcand1))
tcandy1<-crossprod(xcand1,y)
beta1<-invxpx1%*%tcandy1
muhat1<-xcand1%*%beta1
r1<-y-muhat1
sigmahat1<-sum(r1*r1)/(n-i-1)
sxx1<-invxpx1*sigmahat1
t1<-beta1[i+1]/(sqrt(sxx1[i+1,i+1]))
pvalue1<-2*(1-pt(abs(t1),n-i-1))

}


if (intcand>0) 
{
lb<-maincand+1
intscorr<-scorr[lb:totcand]
select2<-which.max(intscorr)


ub<-select2+maincand
selectvar2<-candvars[ub]
if (i<kt) selectx2<-xnew[,ub] else selectx2<-xnew
if(i==1) xcand2<-cbind(one,selectx2) else xcand2<-cbind(xcand,selectx2)



invxpx2<-ginv(crossprod(xcand2))
tcandy2<-crossprod(xcand2,y)
beta2<-invxpx2%*%tcandy2
muhat2<-xcand2%*%beta2
r2<-y-muhat2
sigmahat2<-sum(r2*r2)/(n-i-1)
sxx2<-invxpx2*sigmahat2
t2<-beta2[i+1]/(sqrt(sxx2[i+1,i+1]))

cstep<-min((1+intcand-sum(intinmodel))/(1+p-sum(maininmodel)),2*(1+kt-p)/(1+p))


pvalue2<-2*(1-pt(abs(t2),n-i-1))*cstep
}






if(maincand>0 & intcand>0 )
{
if(pvalue1<=pvalue2)
{
pvalue<-pvalue1
r<-r1
beta<-beta1
select<-select1
xcand<-xcand1
selectvar<-selectvar1
}


if(pvalue1>pvalue2)
{
pvalue<-pvalue2
r<-r2
beta<-beta2
select<-select2+maincand
xcand<-xcand2
selectvar<-selectvar2
}


}


if(maincand==0)

{
pvalue<-pvalue2
r<-r2
beta<-beta2
select<-select2
xcand<-xcand2
selectvar<-selectvar2
}


if(intcand==0)

{
pvalue<-pvalue1
r<-r1
beta<-beta1
select<-select1
xcand<-xcand1
selectvar<-selectvar1
}







if(select<=maincand) maininmodel[mainspot[select]]<-1  
else intinmodel[intspot[select-maincand,1],intspot[select-maincand,2]]<-1

maininmod.step[i+1]=sum(maininmodel)
intinmod.step[i+1]=sum(intinmodel)

var.step[i+1]<-selectvar
pval.step[i+1]<-pvalue
ok<-complete.cases(pval.step)
maxpval.step[i+1]<-max(pval.step[ok],pval.step[i+1])
totcand.step[i+1]<-p+intcand-sum(maininmodel)-sum(intinmodel)

cadj.step[i+1]<-min((1+intcand-sum(intinmodel))/(1+p-sum(maininmodel)),2*(1+kt-p)/(1+p))

ghat.step[i+1]<-(((p-sum(maininmodel))*maxpval.step[i+1]) +((intcand-sum(intinmodel))*maxpval.step[i+1]/cadj.step[i+1]))/(1+sum(maininmodel)+sum(intinmodel))

ghigh.step[i]<-(((p-maininmod.step[i])*maxpval.step[i+1])+((intcand-intinmod.step[i])*maxpval.step[i+1]/cadj.step[i]))/(1+sum(maininmodel)+sum(intinmodel)-1)

if(i==maxstep) ghigh.step[i+1]<-(((p-sum(maininmodel))*1) +((intcand-sum(intinmodel))*1/cadj.step[i+1]))/(1+sum(maininmodel)+sum(intinmodel))


alphahat.step[i+1]<-((1+sum(maininmodel)+sum(intinmodel))*gamma0)/((p-sum(maininmodel))+((intcand-sum(intinmodel))/cadj.step[i+1]))


rsquare.step[i+1]<-1-sum(r*r)/sum((y-mean(y))*(y-mean(y)))


if(i==1) beta.step<-rbind(mean(y),beta) else beta.step<-rbind(beta.step,beta)


i<-i+1 
}


}




#Code for no hierarchy
if(method=="nohier")

{
if(quad) kt=2*p+choose(p,2) else kt=p+choose(p,2)


maxstep<-min(n-2,kt,maxsteps)

muhat<-mean(y)


#Set up vector to keep track of main effects in model
maininmodel<-matrix(0,p,1)


#Set up matrices to keep track of which interactions are candidates and which are in the model
intind<-matrix(1,p,p)
intinmodel<-matrix(0,p,p)


#Set up vectors to keep the forward addition sequence info
var.step=matrix("NA",maxstep+1,1)
pval.step=matrix(NA,maxstep+1,1)
maxpval.step=matrix(NA,maxstep+1,1)
ghat.step=matrix(NA,maxstep+1,1)
ghigh.step=matrix(NA,maxstep+1,1)
totcand.step=matrix(NA,maxstep+1,1)
alphahat.step=matrix(NA,maxstep+1,1)
rsquare.step=matrix(NA,maxstep+1,1)



#Initialize values for Step 0 (Intercept)
var.step[1]="Intercept"
pval.step[1]=0
maxpval.step[1]=0
ghat.step[1]=0
totcand.step[1]=kt
alphahat.step[1]=0
rsquare.step[1]=0




i=1

#Begin loop to get each step of forward addition sequence

while(i<=maxstep) {





#Figure out which main effects to consider
mainind<-1-maininmodel
maincand<-sum(mainind)


j<-1
mcnt<-0
while(mcnt<maincand) {

if(mainind[j]==1 & mcnt==0)  {
 xnew<-xs[,j]
 candvars<-paste("x",deparse(j),sep="")
 mainspot<-j
}
if(mainind[j]==1 & mcnt >0)  {
 xnew<-cbind(xnew,xs[,j])
      candvars<-cbind(candvars,paste("x",deparse(j),sep=""))
      mainspot<-cbind(mainspot,j)
}
if(mainind[j]==1) mcnt<-mcnt+1

j<-j+1
}

 

#Figure out which interactions/quadratic terms to consider
#Code depends on whether quadratic terms are allowed using quad call in function


#This code only allows interactions

if(quad==F)
{
j<-1
k<-j+1
intcnt<-0


while(j<p)
{
 while(k<=p){

 intind[j,k]=1

 intind[k,j]=0
 intind[j,j]=0
 intind[k,k]=0

 if(intinmodel[j,k]==1) intind[j,k]=0


  if( intinmodel[j,k]==0){

   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j+1
}

}


#This code handles quadratic terms with interactions if requested

if(quad==T)


{
j<-1
k<-j
intcnt<-0

while(j<=p)
{
 while(k<=p){

  intind[j,k]=1

  if (k>j)  intind[k,j]=0

  if(intinmodel[j,k]==1) intind[j,k]=0

  if( intinmodel[j,k]==0){
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j
}
}

#Count number of candidate terms
intcand<-sum(intind)
totcand<-maincand+intcand

#Get new residuals
rstep<-y-muhat

#Get correlations to pick term to enter
if(i>1) rxnew<-residuals(lm(xnew~xcand)) else rxnew<-xnew

scorr<-abs(cor(rstep,rxnew))
select<-which.max(scorr)
selectvar<-candvars[select]
xnew<-as.matrix(xnew)
if (i<kt) selectx<-xnew[,select] else selectx<-xnew

if(i==1) xcand<-cbind(one,selectx) else xcand<-cbind(xcand,selectx)


if(select<=maincand) maininmodel[mainspot[select]]<-1  
else intinmodel[intspot[select-maincand,1],intspot[select-maincand,2]]<-1


invxpx<-ginv(crossprod(xcand))
tcandy<-crossprod(xcand,y)
beta<-invxpx%*%tcandy
muhat<-xcand%*%beta
r<-y-muhat
sigmahat<-sum(r*r)/(n-i-1)
sxx<-invxpx*sigmahat
t<-beta[i+1]/(sqrt(sxx[i+1,i+1]))
pvalue<-2*(1-pt(abs(t),n-i-1))


var.step[i+1]<-selectvar
pval.step[i+1]<-pvalue
ok<-complete.cases(pval.step)
maxpval.step[i+1]<-max(pval.step[ok],pval.step[i+1])
totcand.step[i+1]<-totcand-1


ghat.step[i+1]<-(totcand.step[i+1]*maxpval.step[i+1])/(1+sum(maininmodel)+sum(intinmodel))
ghigh.step[i]<-(totcand.step[i]*maxpval.step[i+1])/(1+sum(maininmodel)+sum(intinmodel)-1)
if(i==maxstep) ghigh.step[i+1]<-(totcand.step[i+1]*1)/(1+sum(maininmodel)+sum(intinmodel))


alphahat.step[i+1]<-((1+sum(maininmodel)+sum(intinmodel))*gamma0)/(totcand.step[i+1])

rsquare.step[i+1]<-1-sum(r*r)/sum((y-mean(y))*(y-mean(y)))


if(i==1) beta.step<-rbind(mean(y),beta) else beta.step<-rbind(beta.step,beta)


i<-i+1 
}


}




#Code for adjusted no hierarchy
if(method=="nohieradj")

{
if(quad) kt=2*p+choose(p,2) else kt=p+choose(p,2)


maxstep<-min(n-2,kt,maxsteps)

muhat<-mean(y)


#Set up vector to keep track of main effects in model
maininmodel<-matrix(0,p,1)


#Set up matrices to keep track of which interactions are candidates and which are in the model
intind<-matrix(1,p,p)
intinmodel<-matrix(0,p,p)


#Set up vectors to keep the forward addition sequence info
var.step=matrix("NA",maxstep+1,1)
pval.step=matrix(NA,maxstep+1,1)
maxpval.step=matrix(NA,maxstep+1,1)
ghat.step=matrix(NA,maxstep+1,1)
ghigh.step=matrix(NA,maxstep+1,1)
totcand.step=matrix(NA,maxstep+1,1)
alphahat.step=matrix(NA,maxstep+1,1)
rsquare.step=matrix(NA,maxstep+1,1)
cadj.step=matrix(NA,maxstep+1,1)
maininmod.step=matrix(NA,maxstep+1,1)
intinmod.step=matrix(NA,maxstep+1,1)



#Initialize values for Step 0 (Intercept)
var.step[1]="Intercept"
pval.step[1]=0
maxpval.step[1]=0
ghat.step[1]=0
totcand.step[1]=kt
alphahat.step[1]=0
rsquare.step[1]=0
cadj.step[1]=(1+kt-p)/(1+p)
maininmod.step[1]=0
intinmod.step[1]=0



i=1

#Begin loop to get each step of forward addition sequence


while(i<=maxstep) {



#Figure out which main effects to consider
mainind<-1-maininmodel
maincand<-sum(mainind)


j<-1
mcnt<-0
while(mcnt<maincand) {

if(mainind[j]==1 & mcnt==0)  {
 xnew<-xs[,j]
 candvars<-paste("x",deparse(j),sep="")
 mainspot<-j
}
if(mainind[j]==1 & mcnt >0)  {
 xnew<-cbind(xnew,xs[,j])
      candvars<-cbind(candvars,paste("x",deparse(j),sep=""))
      mainspot<-cbind(mainspot,j)
}
if(mainind[j]==1) mcnt<-mcnt+1

j<-j+1
}

 

#Figure out which interactions/quadratic terms to consider
#Code depends on whether quadratic terms are allowed using quad call in function


#This code only allows interactions

if(quad==F)
{
j<-1
k<-j+1
intcnt<-0


while(j<p)
{
 while(k<=p){

 intind[j,k]=1

 intind[k,j]=0
 intind[j,j]=0
 intind[k,k]=0

 if(intinmodel[j,k]==1) intind[j,k]=0


  if( intinmodel[j,k]==0){

   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j+1
}

}


#This code handles quadratic terms with interactions if requested

if(quad==T)


{
j<-1
k<-j
intcnt<-0

while(j<=p)
{
 while(k<=p){

  intind[j,k]=1

  if (k>j)  intind[k,j]=0

  if(intinmodel[j,k]==1) intind[j,k]=0


  if( intinmodel[j,k]==0){
   if(maincand==0 & intcnt==0) xnew<-xs[,j]*xs[,k] else xnew<-cbind(xnew,xs[,j]*xs[,k])
   if(maincand==0 & intcnt==0) candvars<-paste("x",deparse(j),"*x",deparse(k),sep="") 
     else candvars<-cbind(candvars,paste("x",deparse(j),"*x",deparse(k),sep="") )
   if(intcnt==0) intspot<-cbind(j,k) else intspot<-rbind(intspot,cbind(j,k))
   intcnt<-intcnt+1
  }


  k<-k+1
 }
 
 j<-j+1
      k<-j
}
}

#Count number of candidate terms
intcand<-sum(intind)
totcand<-maincand+intcand

#Get new residuals
rstep<-y-muhat

#Get correlations to pick term to enter
if(i>1) rxnew<-residuals(lm(xnew~xcand)) else rxnew<-xnew

scorr<-abs(cor(rstep,rxnew))

xnew<-as.matrix(xnew)


if (maincand>0) 
{
mainscorr<-scorr[1:maincand]
select1<-which.max(mainscorr)
selectvar1<-candvars[select1]


if (i<kt) selectx1<-xnew[,select1] else selectx1<-xnew
if(i==1) xcand1<-cbind(one,selectx1) else xcand1<-cbind(xcand,selectx1)



invxpx1<-ginv(crossprod(xcand1))
tcandy1<-crossprod(xcand1,y)
beta1<-invxpx1%*%tcandy1
muhat1<-xcand1%*%beta1
r1<-y-muhat1
sigmahat1<-sum(r1*r1)/(n-i-1)
sxx1<-invxpx1*sigmahat1
t1<-beta1[i+1]/(sqrt(sxx1[i+1,i+1]))
pvalue1<-2*(1-pt(abs(t1),n-i-1))

}


if (intcand>0) 
{
lb<-maincand+1
intscorr<-scorr[lb:totcand]
select2<-which.max(intscorr)


ub<-select2+maincand
selectvar2<-candvars[ub]
if (i<kt) selectx2<-xnew[,ub] else selectx2<-xnew
if(i==1) xcand2<-cbind(one,selectx2) else xcand2<-cbind(xcand,selectx2)



invxpx2<-ginv(crossprod(xcand2))
tcandy2<-crossprod(xcand2,y)
beta2<-invxpx2%*%tcandy2
muhat2<-xcand2%*%beta2
r2<-y-muhat2
sigmahat2<-sum(r2*r2)/(n-i-1)
sxx2<-invxpx2*sigmahat2
t2<-beta2[i+1]/(sqrt(sxx2[i+1,i+1]))

cstep<-min((1+kt-p-sum(intinmodel))/(1+p-sum(maininmodel)),2*(1+kt-p)/(1+p))

pvalue2<-2*(1-pt(abs(t2),n-i-1))*cstep
}






if(maincand>0 & intcand>0 & pvalue1<=pvalue2)

{
pvalue<-pvalue1
r<-r1
beta<-beta1
select<-select1
xcand<-xcand1
selectvar<-selectvar1
}


if(maincand>0 & intcand>0 & pvalue1>pvalue2)

{
pvalue<-pvalue2
r<-r2
beta<-beta2
select<-select2+maincand
xcand<-xcand2
selectvar<-selectvar2
}


if(maincand==0)

{
pvalue<-pvalue2
r<-r2
beta<-beta2
select<-select2
xcand<-xcand2
selectvar<-selectvar2
}


if(intcand==0)

{
pvalue<-pvalue1
r<-r1
beta<-beta1
select<-select1
xcand<-xcand1
selectvar<-selectvar1
}






if(select<=maincand) maininmodel[mainspot[select]]<-1  
else intinmodel[intspot[select-maincand,1],intspot[select-maincand,2]]<-1

maininmod.step[i+1]=sum(maininmodel)
intinmod.step[i+1]=sum(intinmodel)


var.step[i+1]<-selectvar
pval.step[i+1]<-pvalue
ok<-complete.cases(pval.step)
maxpval.step[i+1]<-max(pval.step[ok],pval.step[i+1])
totcand.step[i+1]<-totcand-1

cadj.step[i+1]<-min((1+kt-p-sum(intinmodel))/(1+p-sum(maininmodel)),2*(1+kt-p)/(1+p))


ghat.step[i+1]<-(((p-sum(maininmodel))*maxpval.step[i+1]) +((kt-p-sum(intinmodel))*maxpval.step[i+1]/cadj.step[i+1]))/(1+sum(maininmodel)+sum(intinmodel))

ghigh.step[i]<-(((p-maininmod.step[i])*maxpval.step[i+1])+((kt-p-intinmod.step[i])*maxpval.step[i+1]/cadj.step[i]))/(1+sum(maininmodel)+sum(intinmodel)-1)
if(i==maxstep) ghigh.step[i+1]<-(((p-sum(maininmodel))*1) +((kt-p-sum(intinmodel))*1/cadj.step[i+1]))/(1+sum(maininmodel)+sum(intinmodel))

alphahat.step[i+1]<-((1+sum(maininmodel)+sum(intinmodel))*gamma0)/((p-sum(maininmodel))+((kt-p-sum(intinmodel))/cadj.step[i+1]))

rsquare.step[i+1]<-1-sum(r*r)/sum((y-mean(y))*(y-mean(y)))


if(i==1) beta.step<-rbind(mean(y),beta) else beta.step<-rbind(beta.step,beta)


i<-i+1 
}



}



i<-1
mtest<-matrix(T,maxstep+1,1)
while(i<=maxstep)

{
mtest[i]=(maxpval.step[i]<maxpval.step[i+1])
i<-i+1
}


#This is code common to all programs that just prints out the fsr results


stepmax<-which.max(ghigh.step)
#stepmax<-which.max(ghat.step)

alphamax<-maxpval.step[stepmax]

#Select last model where ghat<=gamma0 and maxpval<alpha where ghat is maximized
fsrselect<-which.max((ghat.step<=gamma0 & maxpval.step<=alphamax & mtest)*row(ghat.step))


alpha.fsr<-alphahat.step[fsrselect]
var.fsr<-matrix("NA",fsrselect,1)




var.fsr<-var.step[1:fsrselect]

rsquare.fsr<-rsquare.step[1:fsrselect]
rsquare.fsr<-max(rsquare.fsr)



lospot<-1+((fsrselect-1)*fsrselect/2)
hispot<-fsrselect*(fsrselect+1)/2

if(fsrselect==1) beta.fsr<-beta.step[1] else beta.fsr<-beta.step[lospot:hispot]




model.fsr<-matrix(NA,fsrselect,2)
model.fsr[,1]<-var.fsr
model.fsr[,2]<-round(beta.fsr,digits)

colnames(model.fsr)<-c("Term entered","Coefficient")

if(print){print(paste("Estimated alpha:", round(alpha.fsr,digits))) 
      print(paste("Number of terms in model:", fsrselect))
	#print(paste("FSR model term (coefficient):",mod.fsr) )
	print(paste("R-Square of FSR model:",round(rsquare.fsr,digits)))
      print("FSR model")
	print(model.fsr)
      }

info.step<-cbind(var.step,pval.step,maxpval.step,ghat.step,ghigh.step, rsquare.step)

if(method=="nohieradj" | method=="weakadj") 
colnames(info.step)<-c("Term entered","Adjusted p-value","Max adjusted p-value","Ghat low","Ghat high","Cumulative R-square")
else colnames(info.step)<-c("Term entered","p-value","Max p-value","Ghat low","Ghat high","Cumulative R-square")



return(list(info.step=info.step,alpha.fsr=alpha.fsr,rsquare.fsr=rsquare.fsr,model.fsr=model.fsr,size=fsrselect))







}







