
#As model 1 from Arribas and Romo (2014)

source('LFunctions_to_calculate_real_archetypes_with_swap2.R')
source('LstepArchetypesMod.R')
source('LstepArchetypoids.R')
source('fhdr.R')
source('robout.R')


library(archetypes)
library(Anthropometry)
library( roahd )
library(rainbow)
library(mrfDepth)
library(fda)
library(ks)
library(cluster)
library(fda.usc)


# Auxiliary rep.col and rep.row functions (http://www.r-bloggers.com/a-quick-way-to-do-row-repeat-and-col-repeat-rep-row-rep-col/)
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

# Covariance functions      
cov.fun=function(d,k,c,mu){ #Modelo Hydmann y Shang 2010
        k*exp(-c*d^mu)
}
cov.fun2=function(d,k,c){   #Modelo Fraiman y Svarc 2013
        k*exp(-c*d)
}

# Contamination models
f.out1=function(t,m){            
	u<-runif(1,0,1) 
	if (u<0.5){
		f<-4*t+exp(-(t-m)^2/0.02)/sqrt(2*pi*0.01)-1.8
	}
	else{f<-4*t-exp(-(t-m)^2/0.02)/sqrt(2*pi*0.01)+1.8}
	return(f)
	}
f.out2=function(t,th){4*t+2*sin(4*(t+th)*pi)}
f.out5=function(t){30*(1-t)*t^(3/2)}  



########################################################
#Table simulated data in Section 3.1
########################################################

out.rate=.05 #0, 0.15


VP<-c()
FP<-c()
VPr<-c()
FPr<-c()

p=50
t=seq(0,1,len=p)
d=dist(t,upper=T,diag=T)
d.matrix=as.matrix(d)

#Model 1
t.cov=cov.fun2(d.matrix,0.3,1/0.3)  #covariance function in time
L=chol(t.cov)  # Cholesky Decomposition
mu=30*t*(1-t)^(3/2)



th=0.5
K=3
n=100
nrun=100
 norep=20


VP.or<-c()
FP.or<-c()
VPr.or<-c()
FPr.or<-c()
  
for (ite in 1:nrun) {
		cat(ite, '\n')
		set.seed(ite)
### Generate Data
e=matrix(rnorm(n*p),p,n)
y=mu+t(L)%*%e    

### Generate Outliers
trueout<-sort(sample(1:n,n*out.rate))
t.cov=cov.fun(d.matrix,0.5,1,1)  #covariance function in time
L=chol(t.cov) # Cholesky Decomposition     

for (i in trueout){
	e=rnorm(p)
    y[,i]= f.out5(t)+t(L)%*%e                     #Outliers for Model 1               
}

p=length(t)
n=ncol(y)

#################################### ADA

 X=t(y)
 set.seed(1234)
 lass10d <- stepLArchetypoids3(data=X, k=3, norep)

a3=t(y[,lass10d[[1]][[1]]])
data=t(y)
huge=200
k=3 #

  n <- ncol(t(data))
  x_gvv <- rbind(t(data), rep(huge, n))
  zs <- rbind(t(a3), rep(huge, k))  
  zs <- as.matrix(zs)
  alphascce <- matrix(0, nrow = k, ncol = n)
  
  for (j in 1 : n){
   alphascce[, j] = coef(nnls(zs, x_gvv[,j]))
  }

wm= which.min(c(sum(alphascce[1,]>th[1]),sum(alphascce[2,]>th[1]),sum(alphascce[3,]>th[1])))
out= which(alphascce[wm,]>th[1])

#outliers by threshold
aa.faltan<-length(setdiff(trueout,out))
aa.sobran<-length(setdiff(out,trueout))


#outliers by multivariate: Mah
sco=t(alphascce)
rownames(sco) = as.numeric(1:n)
s = cbind(sco, rep(1, n))

if( (out.rate == 0.15) & (ite == 35)){outr = robout(s, 1, "mve") # error with mcd Error in solve.default(cov, ...) 
}else{
outr = robout(s, 1, "mcd")}



aa.faltanr<-length(setdiff(trueout,outr$outliers))
aa.sobranr<-length(setdiff(outr$outliers,trueout))



################ OutlierGram (Arribas and Romo, 2014)

fD = fData( t, t(y) )
oo=outliergram( fD, display = F )
sout.OG.faltan<-length(setdiff(trueout,oo$ID_outliers))
sout.OG.sobran<-length(setdiff(oo$ID_outliers,trueout))


################ Robust Mahalanobis Distance (See Hydmann and Shang 2010)  

yf=fds(t,y)
 colnames(yf$y)=1:n

FO.robMah=foutliers(yf,method = "robMah") 
sout.robMah.faltan<-length(setdiff(trueout,FO.robMah$outliers))
sout.robMah.sobran<-length(setdiff(FO.robMah$outliers,trueout))

################ Integrated Square Error method (Hyndman and Ullah (2007); See Hydmann and Shang 2010)        
FO.HU<-foutliers(data=yf, method = "HUoutliers",order = 2, lambda = 3.29)
soutHU.faltan<-length(setdiff(trueout,FO.HU$outliers))
soutHU.sobran<-length(setdiff(FO.HU$outliers,trueout))


#################################################################
x=t(y)
rownames(x)<-as.character(1:n)

################ Package fda.usc          
Datosf<-fdata(x,t)

############################LRT (see Febrero et al. 2007)
flo=outliers.lrt(Datosf,nb=200,smo=0.05,dfunc=depth.mode,trim=out.rate)
#flo=foutliers(yf,method = "lrt",trim=out.rate) 
soutFLO.faltan<-length(setdiff(trueout,flo$outliers))
soutFLO.sobran<-length(setdiff(flo$outliers,trueout))


################ Depth Based Trimming (Febrero, Galeano, Gonzalez-Manteiga, 2008)  
#dt=foutliers(yf,method ="depth.trim",trim=out.rate)
    
dt<-outliers.depth.trim(Datosf,nb=200,smo=0.05,trim=out.rate,dfunc=depth.mode)
soutDT.faltan<-length(setdiff(trueout,dt$outliers))
soutDT.sobran<-length(setdiff(dt$outliers,trueout))


################ Depth Based Weighting (Febrero, Galeano, Gonzalez-Manteiga, 2008)  
#dp=foutliers(yf,method = "depth.pond",trim=out.rate)
    
dp<-outliers.depth.pond(Datosf,nb=200,smo=0.05,dfunc=depth.mode)
soutDP.faltan<-length(setdiff(trueout,dp$outliers))
soutDP.sobran<-length(setdiff(dp$outliers,trueout))


######################FOM (see Hubert et. al 2015)

yt=array(y,dim=c(dim(y),1))
Result <- fOutl(yt,  type = "fAO", diagnostic = TRUE);

# The user may opt to draw a cut off line separating the outliers.
# which will be plotted in red
ffom=fom(Result, cutoff = TRUE) 

soutFOM.faltan<-length(setdiff(trueout,which(ffom$data[,4]=='red')))
soutFOM.sobran<-length(setdiff(which(ffom$data[,4]=='red'),trueout))

######################Functional Boxplots (Sun and Genton 2011)
out.fb=fda::fbplot(y,plot=F)$outpoint 
sout.fb.faltan<-length(setdiff(trueout,out.fb))
sout.fb.sobran<-length(setdiff(out.fb,trueout))

#########################High density regions (Hydmann and Shang 2010) 
Datos=fds(t, y)
colnames(Datos$y)<-seq(1,n,1)
#fboxplot(data = Datos, plot.type = "functional", type = "hdr", alpha = c(out.rate,0.5), projmethod="PCAproj",ylab="x(t)",xlab="",xpd=NA)
hdrout<-fhdr(Datos, alpha = c(out.rate,0.5),projmethod="PCAproj",xlab='',ylab='',plotlegend=F)
sout.hdr.faltan<-length(setdiff(trueout,hdrout))
sout.hdr.sobran<-length(setdiff(hdrout,trueout))



################ SUMARY
VPr.or<-rbind(VPr.or,1-c(aa.faltan,aa.faltanr,sout.OG.faltan,sout.robMah.faltan,soutHU.faltan,soutFLO.faltan,soutDT.faltan,soutDP.faltan,soutFOM.faltan,sout.fb.faltan,sout.hdr.faltan)/length(trueout))
FPr.or<-rbind(FPr.or,c(aa.sobran,aa.sobranr,sout.OG.sobran,sout.robMah.sobran,soutHU.sobran,soutFLO.sobran,soutDT.sobran,soutDP.sobran,soutFOM.sobran,sout.fb.sobran,sout.hdr.sobran)/(n-length(trueout)))
VP.or<-rbind(VP.or,length(trueout)-c(aa.faltan,aa.faltanr,sout.OG.faltan,sout.robMah.faltan,soutHU.faltan,soutFLO.faltan,soutDT.faltan,soutDP.faltan,soutFOM.faltan,sout.fb.faltan,sout.hdr.faltan))
FP.or<-rbind(FP.or,c(aa.sobran,aa.sobranr,sout.OG.sobran,sout.robMah.sobran,soutHU.sobran,soutFLO.sobran,soutDT.sobran,soutDP.sobran,soutFOM.sobran,sout.fb.sobran,sout.hdr.sobran))

}


apply(VPr.or,2,mean)
apply(VPr.or,2,sd)
apply(FPr.or,2,mean)
apply(FPr.or,2,sd)




##############################################
#Figure example: Sect. 2.4.3 Functional outlier detection by archetype analysis
#############################################

ite=1
K=3
n=100
out.rate=.02

set.seed(ite)
### Generate Data
e=matrix(rnorm(n*p),p,n)
y=mu+t(L)%*%e    

### Generate Outliers
trueout<-sort(sample(1:n,n*out.rate))
t.cov=cov.fun(d.matrix,0.5,1,1)  #covariance function in time
L=chol(t.cov) # Cholesky Decomposition     

for (i in trueout){
	e=rnorm(p)
    y[,i]= f.out5(t)+t(L)%*%e                     #Outliers for Model 1                    
}

p=length(t)
n=ncol(y)


 X=t(y)
 set.seed(1234)
 norep=20
 lass10d <- stepLArchetypoids3(data=X, k=3, norep)



a3=t(y[,lass10d[[1]][[1]]])
data=t(y)
huge=200
k=3 #

  n <- ncol(t(data))
  x_gvv <- rbind(t(data), rep(huge, n))
  zs <- rbind(t(a3), rep(huge, k))  
  zs <- as.matrix(zs)
  alphascce <- matrix(0, nrow = k, ncol = n)
  
  for (j in 1 : n){
   alphascce[, j] = coef(nnls(zs, x_gvv[,j]))
  }

library(vcd)
ternaryplot(t(alphascce),grid=T,dimnames=c(1,2,3),col=1,cex=.6,main='')


 matplot(t, y, type = "l",lty=1,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,col=8) #data
lines(t, y[,trueout[1]], type = "l",lty=1,xlab="",ylab="",col=1) #outliers
lines(t, y[,trueout[2]], type = "l",lty=1,xlab="",ylab="",col=1) #outliers


matplot(t, y, type = "l",lty=1,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,col=8) #data
lines(t, a3[1,], type = "l",lty=1,xlab="",ylab="",col=2,lwd=3) #archetypoids
lines(t, a3[2,], type = "l",lty=2,xlab="",ylab="",col=3,lwd=3) #archetypoids
lines(t, a3[3,], type = "l",lty=3,xlab="",ylab="",col=4,lwd=3) #archetypoids

  

#############################
#with out.rate=0-> to see ternaryplot without holes


K=3
n=100
out.rate=0
ite=1

set.seed(ite)
### Generate Data
e=matrix(rnorm(n*p),p,n)
y=mu+t(L)%*%e    

### Generate Outliers
trueout<-sort(sample(1:n,n*out.rate))
t.cov=cov.fun(d.matrix,0.5,1,1)  #covariance function in time
L=chol(t.cov) # Cholesky Decomposition     

for (i in trueout){
	e=rnorm(p)
    y[,i]= f.out5(t)+t(L)%*%e                     #Outliers for Model 1
                         
}

p=length(t)
n=ncol(y)


 X=t(y)
 set.seed(1234)
 norep=20
 lass10d <- stepLArchetypoids3(data=X, k=3, norep)



a3=t(y[,lass10d[[1]][[1]]])
data=t(y)
huge=200
k=3 #

  n <- ncol(t(data))
  x_gvv <- rbind(t(data), rep(huge, n))
  zs <- rbind(t(a3), rep(huge, k))  
  zs <- as.matrix(zs)
  alphascce <- matrix(0, nrow = k, ncol = n)
  
  for (j in 1 : n){
   alphascce[, j] = coef(nnls(zs, x_gvv[,j]))
  }

ternaryplot(t(alphascce),grid=T,dimnames=c(1,2,3),col=1,cex=.6,main='')



matplot(t, y, type = "l",lty=1,xlab="",ylab="",cex.lab=1.5,cex.axis=1.5,col=8) #data
lines(t, a3[1,], type = "l",lty=1,xlab="",ylab="",col=2,lwd=3) #archetypoids
lines(t, a3[2,], type = "l",lty=2,xlab="",ylab="",col=3,lwd=3) #archetypoids
lines(t, a3[3,], type = "l",lty=3,xlab="",ylab="",col=4,lwd=3) #archetypoids

