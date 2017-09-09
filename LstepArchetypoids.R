library(FNN) #for the function knn.

intpm<-function(ri){
(t(ri)%*%ri)}

norm2.Frnormfn <- function(m) {


  return(sum(apply(m,2,intpm))) # 
}

LstepArchetypoids <- function(i,nearest,data,ArchObj,sequ,aux){


  
  N = dim(data)[1]
  
  if(sequ){
   ai <- bestModel(ArchObj[[i]])
  }else{
   ai <- bestModel(ArchObj[[i-aux]])  
  }
  
  if(is.null(parameters(ai))){
    stop("No archetypes computed")  
  }
  
    
  if(nearest){
  
    	dime=dim(parameters(ai))
  	ar=parameters(ai)
  	R=matrix(0,nrow=dime[1], ncol=N)
  	for(di in 1:dime[1]){ #dime[1] is the number of archetypes
  
  		#Raux=matrix(ar[di,],byrow=T,nrow=N,ncol=dime[2])-data
  Raux=matrix(ar[di,],byrow=F,ncol=N,nrow=dime[2])-t(data)
  
  
  
  		
  R[di,]= apply(Raux,2,intpm)

    	}

  ini_arch <- apply(R, 1, which.min)

 
    
    if( any(duplicated(ini_arch)) == TRUE){ #for very rare ocasions, it can be improved
     # ini_arch=unique(ini_arch)      

	 k=1
       neig=apply(R, 1, rank)     
      indices1 <- neig[1,]
      ini_arch <- indices1[,k]
      
      while(any(duplicated(ini_arch))){
        k=k+1  
       
        indicesk <- neig[k,]
               
        dupl <- anyDuplicated(indices1)
        ini_arch <- c(indices1[-dupl],indicesk[dupl])
       
      }
    }
    
  }else{
    ini_arch <- apply(coef(ai, "alphas"), 2, which.max) 
  }
  
  res <- LrealArchetypes(data, huge = 200, i, ini_arch) #as.vector(ini_arch)
  cat("Done!") 
  return(list(res[[1]],res[[2]],ini_arch))
}




stepLArchetypoidsn <- function (data, k, norep=3){
  mycall <- match.call()
  as <- list()
  for (i in 1:length(k)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypoids"
    lass <- LstepArchetypesMod(data,k[i],verbose=FALSE,nrep=norep,saveHistory=F) 
          as[[i]]<- LstepArchetypoids(k[i],TRUE,data,lass,FALSE,k[i]-1)
     
      
    
  }
  return(structure(as, class = "stepArchetypoids", call = mycall))
}

screeplot.stepArchetypoids <- function(as,nar=NULL){
xx=c()
l=length(as)

for (i in 1:l){
xx[i]=as[[i]][[2]]
}
if (is.null(nar)){
plot(1:l,xx,type='b',xlab='Archetypoids', ylab='RSS')}
else{
plot(nar,xx,type='b',xlab='Archetypoids', ylab='RSS')}
}


stepLArchetypoidsw <- function (data, k, norep=3){
  mycall <- match.call()
  as <- list()
  for (i in 1:length(k)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypoids"
    lass <- LstepArchetypesMod(data,k[i],verbose=FALSE,nrep=norep,saveHistory=F) 
          as[[i]]<- LstepArchetypoids(k[i],FALSE,data,lass,FALSE,k[i]-1)
     
      
    
  }
  return(structure(as, class = "stepArchetypoids", call = mycall))
}


stepLArchetypoidsb <- function (data, k, norep=3){
  mycall <- match.call()
  as <- list()
  for (i in 1:length(k)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypoids"
    lass <- LstepArchetypesMod(data,k[i],verbose=FALSE,nrep=norep,saveHistory=F) 
    ai <- bestModel(lass[[1]]) #archetypes
    ini_arch=c()
 for (j in 1:k[i]){
 ini_arch[j]=which.max(ai$betas[j,])
 }
 
as[[i]] <- LrealArchetypes(data,huge = 200, k[i], ini_arch) #which beta archetypoid
      
  }
  return(structure(as, class = "stepArchetypoids", call = mycall))
}




stepLArchetypoids3 <- function (data, k, norep=3){
  mycall <- match.call()
  as <- list()
   asn=list()
   asw=list()
   asb=list()
  for (i in 1:length(k)) {
    asb[[i]] <- list()
    asn[[i]]=list()
    asw[[i]]=list()
    as[[i]]=list()
    class(asn[[i]]) <- "repArchetypoids"
    lass <- LstepArchetypesMod(data,k[i],verbose=FALSE,nrep=norep,saveHistory=F) 
    ai <- bestModel(lass[[1]]) #archetypes
    ini_arch=c()
 for (j in 1:k[i]){
 ini_arch[j]=which.max(ai$betas[j,])
 }
 
#from cand_beta
asb[[i]] <- LrealArchetypes(data,huge = 200, k[i], ini_arch) 

#from nearest
asn[[i]]<- LstepArchetypoids(k[i],TRUE,data,lass,FALSE,k[i]-1)
#from cand_alpha
asw[[i]]<- LstepArchetypoids(k[i],FALSE,data,lass,FALSE,k[i]-1)

aux=c(asn[[i]][[2]],asw[[i]][[2]],asb[[i]][[2]])
wmi=which.min(aux)
auxa=list()
auxa[[1]]=asn[[i]]
auxa[[2]]=asw[[i]]
auxa[[3]]=asb[[i]]
as[[i]]= auxa[[wmi]]
      
  }
  return(structure(as, class = "stepArchetypoids", call = mycall))
}

