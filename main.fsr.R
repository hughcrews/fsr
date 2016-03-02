#This function implements the Fast FSR methodology using main effects only.  



main.fsr<-function(x,y,gamma0=.05,center=T,maxsteps=kt,quad=F,print=T,digits=4){
require(MASS)
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


i<-1
mtest<-matrix(T,maxstep+1,1)
while(i<=maxstep)

{
mtest[i]=(maxpval.step[i]<maxpval.step[i+1])
i<-i+1
}


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
colnames(info.step)<-c("Term entered","p-value","Max p-value","Ghat low","Ghat high","Cumulative R-square")



return(list(info.step=info.step,alpha.fsr=alpha.fsr,rsquare.fsr=rsquare.fsr,model.fsr=model.fsr,size=fsrselect))






}

