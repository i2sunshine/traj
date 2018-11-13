#To calculate derivative of a curve, i.e., to estimate slope of a curve
## Input:
## a2: a s-by-1 vector, where s is the number of simulations
## b: a s-by-K matrix, where K is the number of knots
## Output: a s-by-1 vector includes slopes of s simulated curves

DerivativeCurve <-function(t, knots, a2, b){

  #no. simulations
  s= dim(b)[1]
  
  #no. of knots
  K= length(knots)
  
  if(K!=dim(b)[2]) print("Error in dim(b)!")else{
     
     precess = matrix((t-knots)^2, nrow=K, ncol=1)
      
     ##find the index of first knot such that knot is larger than given follow-up time (t)
     midx = min(which(knots > t))
     
     #out is a vector with s*1
     out = a2
     
     
     if(is.infinite(midx)==TRUE){ ##i.e., all knot < t
          out = out + 3* b %*% precess
     }else{#i.e.,  exits a knot such that knot > t
         if(midx < 1)print(paste("Error in ", ppgid, "!", sep="")) else{
            ###i.e., midx >=1
           
            ##for knots <= t
            if(midx>1) out = out + 3 * matrix(b[, 1:(midx-1)], nrow=s, ncol=midx-1) %*% precess[1:(midx-1),1]
            
            ##for knots > t
            out = out - 3* matrix(b[, (midx:K)], nrow=s, ncol=(K-midx+1) ) %*% precess[midx:K, 1]   
         }
        
       
     }
    
  }
  
  return (out)
}