
##plot posterior density as well as given prior density
PosteriorPlot<-function(mcmc.combi, out.dir){
  
  #number of simulations
  s = dim(mcmc.combi)[1]
  
  ## for sigma_b, i.e., standard deviation of b
  pdf.file.sigma_b = paste(out.dir, "/jags_prior_posterior_SD_b.pdf", sep="")
  pdf(pdf.file.sigma_b)
  plot(density(rnorm(s, 0, sd=sqrt(variance_b))),main="sigma_b (black=prior, red=mcmc)")
  lines(density(mcmc.combi[,"sigma_b"]),  col="red")
  dev.off();
  
  ##for sigma_eps, i.e., SD for regression error eps
  pdf.file.sigma_eps = paste(out.dir, "/jags_prior_posterior_SD_eps.pdf", sep="")
  pdf(pdf.file.sigma_eps)
  plot(density(rnorm(s, 0, sd=sqrt(variance_eps))),main="sigma_eps (black=prior, red=mcmc)")
  lines(density(mcmc.combi[,"sigma_eps"]),  col="red")
  dev.off();
  
  
  ## for a
  pdf.file.a1 = paste(out.dir, "/jags_prior_posterior_a1.pdf", sep="")
  pdf(pdf.file.a1)
  plot(density(rnorm(s, mu_a[1], sd= sqrt(variance_a[1]))), main="a[1] (black=prior, red = mcmc)")
  lines(density(mcmc.combi[,"a[1]"]), col="red")
  dev.off();
  
  pdf.file.a2 = paste(out.dir, "/jags_prior_posterior_a2.pdf", sep="") 
  pdf(pdf.file.a2)
  plot(density(rnorm(s, mu_a[2], sd= sqrt(variance_a[2]))), main="a[1] (black=prior, red = mcmc)")
  lines(density(mcmc.combi[,"a[2]"]), col="red")
  dev.off();
  
}