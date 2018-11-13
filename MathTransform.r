##Aim: to do necessary Math transform so that simulation could be easier.
#Ref: Liang et al. 2012 and Crainiceanu et al. 2005.


# Input:
## t: a vector of time intervals from patient first visit (by year)
## knots: a vector of knots that allows us to model a curve by thin-plate spline
# Output:
## a list including 3 items:
### a matrix (n * 2)
### a matrix (n * K), 
### a matrix (K * K): SVD decomposition
## where n is the number of visits and K is the number of knots. 
MathTransform <-function(t, knots){

   # no. of internal knots
   K=length(knots)
   
   ## J: #followup times
   J = length(t)
  
   ##
   X = cbind(rep(1,J), t) 

   
   #
   Z = (abs(outer(t, knots, "-"))^3)
 
   ##
   Omga = (abs(outer(knots, knots,"-")))^3


   ## SVD decomposition of Omga
   svd.Omga<-svd(Omga)
   msqrt<-t(svd.Omga$v %*%(t(svd.Omga$u)*sqrt(svd.Omga$d)))
   Y = t(solve(msqrt,t(Z)))
   msqrt_inv = Inverse(msqrt)

   
   out=list(X, Y, msqrt_inv)
   return(out)

}


#### ----- End: Preparation: Neccesary Math transform so that simulation could be easier ------------------
