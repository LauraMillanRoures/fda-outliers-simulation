#Loading the archetypes package:
library(archetypes)

#Function to calculate real archetypes (archetypoids):
LrealArchetypes <- function(data, huge=200, k, nearest, ...){

  
  archs_ini <- nearest
  
  n <- ncol(t(data))
  x_gvv <- rbind(t(data), rep(huge, n))
  
  zs <- x_gvv[, archs_ini] 
  zs <- as.matrix(zs)
  
  alphas <- matrix(0, nrow = k, ncol = n)
  for (j in 1 : n){
    alphas[, j] = coef(nnls(zs, x_gvv[,j]))
  }
  
  #resid <- zs %*% alphas - x_gvv
  resid <- zs[1: (nrow(zs)-1),] %*% alphas - x_gvv[1: (nrow(x_gvv)-1),]
  
  
  rss_ini<-norm2.Frnormfn(resid)/n
  
   
  res_def <- Lswap(archs_ini, rss_ini, huge, k, x_gvv, n)
  
  return(res_def)
}


#Function to improve the set of archetypes that uses as initial vector of archetypes, the vector of the nearest 
#individuals to the archetypes calculated by the algorithm of Cutler & Breiman:
Lswap <- function(vect_arch_ini, rss_arch_ini, huge=200, k, x_gvv, n){

 vect_arch_end <- vect_arch_ini
 rss <- rss_arch_ini

 for(l in 1 : k){

  rss1 <- c()

  setpossibles <- setdiff(1:n,vect_arch_ini)

   for(i in setpossibles){
    zs <- x_gvv[,c(i,vect_arch_ini[-l])] 
    zs <- as.matrix(zs)
    alphas <- matrix(0, nrow = k, ncol = n)
     for (j in 1 : n){
      alphas[, j] = coef(nnls(zs, x_gvv[,j]))
     }

    #resid <- zs %*% alphas - x_gvv
  resid <- zs[1: (nrow(zs)-1),] %*% alphas - x_gvv[1: (nrow(x_gvv)-1),]
  
    rss1[i] <- norm2.Frnormfn(resid)/n

    if(rss1[i] < rss){
     rss <- rss1[i]
     vect_arch_end = c(i,vect_arch_ini[-l])
    }
   }
 }

 if(k==1){
  result <- Lswap2_k1(vect_arch_end, vect_arch_ini, rss, huge, k, x_gvv, n)
 }else{
  result <- Lswap2(vect_arch_end, vect_arch_ini, rss, huge, k, x_gvv, n)
 } 

 return(result)
}


Lswap2 <- function(vect_arch_end, vect_arch_ini, rss, huge=200, k, x_gvv, n){
  
  vect_arch_ini_aux <- vect_arch_ini
  vect_arch_end_aux <- vect_arch_end
  
  rss_aux <- rss
  vect_arch_end2 <- c()
  
  while(any(sort(vect_arch_end_aux) != sort(vect_arch_ini_aux))){ #this loop is executed while both vectors 
    #are different in at least one element. Since we have not found a R function that tell us if two vectors 
    #are equal, we have done the following: first, we have ordered them (because both vectors can be the same, 
    #although their elements are in a different order, for example, c(1,2,3) y c(3,1,2). Then we have checked if 
    #any element do not match with the element that is placed in the same position in the other vector.
    
    se <- setdiff(vect_arch_end_aux,vect_arch_ini_aux) #this function looks for the distinct element between 
    #the initial vector (in the first iteration, the initial vector is either the nearest vector or the which 
    #vector, while in the second iteration, the initial vector is the final vector returned by the swap step 
    #and so on and so forth) and the final vector (the final vector in the first iteration is that one returned 
    #by the swap step but in the following iterations, it is the vector returned by the swap2).
    
    se1 <- setdiff(vect_arch_end_aux,se) #the elements different than the distinct one of the former 
    #setdiff function.
    
    for(l in 1 : length(se1)){
      
      rss1 <- c()
      
      comp <- c(se,se1[-l]) #the vector of the distinct element with the no distincts without one.
      setpossibles <- setdiff(1:n,vect_arch_end_aux)
      
      for(i in setpossibles){
        zs <- x_gvv[,c(i,comp)] 
        zs <- as.matrix(zs)
        alphas <- matrix(0, nrow = k, ncol = n)
        for (j in 1 : n){
          alphas[, j] = coef(nnls(zs, x_gvv[,j]))
        }
        
        #resid <- zs %*% alphas - x_gvv
  resid <- zs[1: (nrow(zs)-1),] %*% alphas - x_gvv[1: (nrow(x_gvv)-1),]
        
        rss1[i] <- norm2.Frnormfn(resid)/n
        
        if(rss1[i] < rss_aux){
          rss_aux <- rss1[i]
          vect_arch_end2 <- c(i,comp)
        }
      }
    }
    
    if(is.null(vect_arch_end2)){ #if vect_arch_end2 is NULL, it means that any vector improve the final vector 
      #of the swap function (called vect_arch_end_aux). Therefore, the vect_arch_end_aux that it is going to be 
      #returned is just the vect_arch_end_aux that it is the final vector of the first swap. If we don't add this 
      #vector, the function displays an error because it might happen that the vect_arch_end returned by the 
      #swap function already is the best vector (this happens with the nba2d database) and therefore the swap2 
      #function is not going to be able to improve it. 
      vect_arch_end_aux <- vect_arch_end_aux
      vect_arch_ini_aux <- vect_arch_end_aux #In addition, we also have to fix the initial vector as the final 
      #vector in order to the while loop stops.
    }else{
      vect_arch_ini_aux <- vect_arch_end_aux #the initial vector of the following iteration must be the final 
      #vector of the previous iteration to compare it with the final vector returned by the swap2 in the following
      #iteration.
      vect_arch_end_aux <- vect_arch_end2 #final vector returned by the swap2.
    }
  }
  
  return(list(archet=vect_arch_end_aux,rss=rss_aux))
  
}


Lswap2_k1 <- function(vect_arch_end, vect_arch_ini, rss, huge=200, k, x_gvv, n){
  
  vect_arch_ini_aux <- vect_arch_ini
  vect_arch_end_aux <- vect_arch_end
  
  rss_aux <- rss
  vect_arch_end2 <- c()
  
  while(vect_arch_end_aux != vect_arch_ini_aux){
    
    rss1 <- c()
    
    setpossibles <- setdiff(1:n,vect_arch_end_aux)
    
    for(i in setpossibles){
      zs <- x_gvv[,i] 
      zs <- as.matrix(zs)
      alphas <- matrix(0, nrow = k, ncol = n)
      for (j in 1 : n){
        alphas[, j] = coef(nnls(zs, x_gvv[,j]))
      }
      
      #resid <- zs %*% alphas - x_gvv
  resid <- zs[1: (nrow(zs)-1),] %*% alphas - x_gvv[1: (nrow(x_gvv)-1),]
      
      rss1[i] <- norm2.Frnormfn(resid)/n
      
      if(rss1[i] < rss_aux){
        rss_aux <- rss1[i]
        vect_arch_end2 = i
      }
    }
    
    if(is.null(vect_arch_end2)){ 
      vect_arch_end_aux <- vect_arch_end_aux
      vect_arch_ini_aux <- vect_arch_end_aux 
    }else{
      vect_arch_ini_aux <- vect_arch_end_aux 
      vect_arch_end_aux <- vect_arch_end2 
    }
  }
  return(list(archet=vect_arch_end_aux,rss=rss_aux))
}


