

# prepare JAGs Model

modelstring = "    
      model{
   
         ## for each row, i.e., each follow-up
         for(j in 1:n){
	        # note tau_eps = 1/sigma2_eps is precision
            eGFR[j] ~ dnorm(m[j], tau_eps)  
        
            # note mean of eGFR
            m[j]  <-  mfe[j] + mre[j]
            mfe[j] <- inprod(X[j,], a[ ]) 
            mre[j] <- inprod(Y[j, ], u[ ])
        }

        ## ---------Begin:  prior distributions -------------------------------
        sigma_eps ~ dnorm(0, precision_eps)
	    sigma_b ~ dnorm(0, precision_b)
	
	    tau_eps <- pow(sigma_eps, -2) 
	    tau_b <- pow(sigma_b, -2)

        a[1] ~ dnorm(mu_a[1], precision_a[1])
        a[2] ~ dnorm(mu_a[2], precision_a[2])

        for(i in 1:K){     
           u[i] ~ dnorm(0, tau_b) 
        }

        ## ----------End:  prior distribution ------------

 
	
        ## ----------Begin: determinative values ------------
   
         # Retrieve original random effect parameters b from u
         for(i in 1:K){
           b[i]  <- inprod(msqrt_inv[i, ], u[ ])
         }
        ## -------- End: determinative values ---------------- 
				
		
     }"

	 # close quote for modelstring

# Write nags model to a file:
model.file.name=paste("modelfile_prior_normal.txt", sep="")
writeLines(modelstring, con= model.file.name)