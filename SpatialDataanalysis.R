#--------------------------------#
#       require library          #
#--------------------------------#
library(readxl)
library(sgeostat)
library(gstat)
library(geospt)
library(sp)
library(geoR)
library(scatterplot3d)
####################################
#########reading & introducing data#
####################################
data=read_excel(file.choose())
str(data)
x=data$X
y=data$Y
temp=data$avg_Temp
summary(temp)
#--------------------------------#
#          EDA                   #
#--------------------------------#
### plots
plot(x, y, type="n",lab=c(16,23,1), xlab="X",ylab="Y")    
text(x,y,format(round(temp,2)))     
title("Average Of Temprature ",cex=1,col=2)
plot(x,y)
plot(x,temp)
plot(y,temp)
####################
######ravand zodoyi#
model=lm(temp~y)
summary(model)
ravand=63.50872-1.48717*y
temp1=temp-ravand 
plot(y,temp1)
plot(x,temp1)
####################
####normality#######
hist(temp,breaks=7,xlab=c(8,30))
shapiro.test(temp)
############################
####normality ravand zodode#
hist(temp1,breaks=7,xlab=c(8,30))
shapiro.test(temp1)
##########################
########## h scatter plot
df<- data.frame(x,y,temp)
library(sp)
coordinates(df) <- c("x", "y") 
hscat(temp ~ 1, data=df, breaks=c(0,7,10,13,16))
###########################################
########## h scatter plot ravand zodode####
df1<- data.frame(x,y,temp1)
library(sp)
coordinates(df1) <- c("x", "y") 
hscat(temp1 ~ 1, data=df1, breaks=c(0,7,10,13,16))
########################
###  haininig for outlier
boxplot(temp)
max(temp)
min(temp)
Q3<-quantile(temp,0.75)    
Ql<-quantile(temp,0.25)
h1<-(Ql-1.5*(Q3-Ql))
h2<-(Q3+1.5*(Q3-Ql))
outlier<-function(temp){
  a<-vector(mode="numeric")
  b<-vector(mode="numeric")
  for(i in 1:length(temp)){
    if(temp[i]>h2 | temp[i]<h1){
      a[i]<-temp[i]
      b<-c(b,i)}
    i=i+1
  }
  return(list("outlier"=a,"number of outlier"=b))
}
outlier(temp)
###############################
### variogram cloud###########
f<- as.geodata(df)
variogcl<-variog(f, bin.cloud=TRUE,,estimator.type="classic")   
variogcl
plot(variogcl, bin.cloud=TRUE)
############################
#### isotropic##############

f<- as.geodata(df)
var4<-variog4(f)     
par(mfrow=c(1,1)) 
plot(var4,lwd=2,type='b' ,option = c(1,2,100),pch=1,type="l")

####find k 

H<-dist(cbind(x,y),method = "euclidean")
k<-floor(1/2*max(H))+1
k

#--------------------------------#
#      fitting variogram         #
#--------------------------------#
### emprical variogram

df1<- data.frame(x,y,temp1)
f1<- as.geodata(df1)

varigc<-variog(f1,estimator.type="classical",max.dist=k)
varigc
variogr<-variog(f1,estimator.type="modulus",max.dist=k)
variogr

par(mfrow=c(1,2)) 
plot(variog(f1,estimator.type="classical",max.dist=k),type= "b", main="Classic")
plot(variog(f1,estimator.type="modulus",max.dist=k),type= "l", main="Robust")

#####################################
#### fiting model for variogram######

va<-variog(f1,max.dist=k,estimator.type="modulus")

ini.valsM <-expand.grid(seq(11,12,l=50),seq(5,6,l=50))
wlsM <- variofit(va, ini=ini.valsM, fix.nug=FALSE,cov.model="matern",nugget = 3 )

ini.valsP.E<-expand.grid(seq(3,7,l=50), seq(1,2,l=50))
wlsP.E <- variofit(va, ini=ini.valsP.E, fix.nug=FALSE,cov.model="powered.exponential",nugget = 1 ,kappa= 2.69)

ini.valsE<-expand.grid(seq(10,11,l=50), seq(7,8,l=50))
wlsE <- variofit(va, ini=ini.valsE, fix.nug=FALSE,cov.model="exponential",nugget=3)

ini.valsSPER<-expand.grid(seq(7,10,l=50), seq(2,3,l=50))
wlsS<- variofit(va, ini=ini.valsSPER, fix.nug=FALSE,cov.model="spherical",nugget =  3)

ini.valsC<- expand.grid(seq(8,9,l=50), seq(10,11,l=50))
wlsC<- variofit(va, ini=ini.valsC, fix.nug=FALSE,cov.model="cauchy",kappa=0.45,nugget = 3)

ini.valsG <- expand.grid(seq(10,11,l=50), seq(4,5,l=50))
wlsG<- variofit(va, ini=ini.valsG, fix.nug=FALSE,cov.model="gaussian",nugget = 3.5)

ini.valsw <- expand.grid(seq(11,12,l=50), seq(4,5,l=50))
wlsw<- variofit(va, ini=ini.valsw, fix.nug=FALSE,cov.model="wave",nugget = 3.5 )

wlsM$value
wlsP.E$value
wlsE$value 
wlsS$value 
wlsG$value # selected 
wlsC$value
wlsw$value

##################

ols <- variofit(va,ini=ini.valsG,fix.nug=FALSE,cov.model="gaussian",nugget = 3, wei="equal")
summary(ols)
wls <- variofit(va,ini=ini.valsG,fix.nug=FALSE,cov.model="gaussian",nugget = 3)
summary(wls)
plot(va)
lines(wls)
lines(ols,col="red")


#------------------------------------------#
# kiriging & pish bini noghat jadid        # 
#------------------------------------------#
df1<- data.frame(x,y,temp1)  
f1=as.geodata(df1)
kriging1<-krige.control(type.krige="ok",cov.model="gaussian", cov.pars=c(11.05,1.03), nugget=0.78)
z1=c(33.73, 55.92)
ok1<-krige.conv(f1, coords=f1$coords,locations=z1,data=f1$data,krige=kriging1)
ravand1=63.50872-1.48717*z1[2]             
pr1=ok1$predict+ravand1
cI1=c(pr1-1.96*sqrt(ok1$krige.var),pr1+1.96*sqrt(ok1$krige.var))
z2=c(55.92,36.66)
ok2<-krige.conv(f1, coords=f1$coords,locations=z2,data=f1$data,krige=kriging1)             
ravand2=63.50872-1.48717*z2[2]             
pr2=ok2$predict+ravand2
cI2=c(pr2-1.96*sqrt(ok2$krige.var),pr2+1.96*sqrt(ok2$krige.var))
z3=c(33.98,48.38)
ok3<-krige.conv(f1, coords=f1$coords,locations=z3,data=f1$data,krige=kriging1)             
ravand3=63.50872-1.48717*z3[2]             
pr3=ok3$predict+ravand3
cI3=c(pr3-1.96*sqrt(ok3$krige.var),pr3+1.96*sqrt(ok3$krige.var))
z4=c(30.90,56.44)
ok4<-krige.conv(f1, coords=f1$coords,locations=z4,data=f1$data,krige=kriging1)             
ravand4=63.50872-1.48717*z4[2]             
pr4=ok4$predict+ravand4
cI4=c(pr4-1.96*sqrt(ok4$krige.var),pr4+1.96*sqrt(ok4$krige.var))
####################################