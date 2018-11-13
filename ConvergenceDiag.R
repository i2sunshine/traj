## Convergence diagnosis

ConvergenceDiag <-function(bayes.mod.fit.mcmc, out.dir){
  
  pdf.file=paste(out.dir, "/jags_fit_mcmc_convergdiag.pdf", sep="")
  pdf(pdf.file)
  plot(bayes.mod.fit.mcmc)
  dev.off()
  
  pdf.file.autocorr= paste(out.dir, "/jags_fit_mcmc_convergdiag_autocorr.pdf", sep="")
  pdf(pdf.file.autocorr)
  autocorr.plot(bayes.mod.fit.mcmc)
  dev.off()
  
  ### Gelman-Rubin convergence diagnostic. Need at least two chains
  pdf.file.gelman=paste(out.dir, "/jags_fit_mcmc_convergdiag_gelman.pdf", sep="")
  pdf(pdf.file.gelman)
  gelman.plot(bayes.mod.fit.mcmc)
  dev.off();
  
  ####
  pdf.file.geweke = paste(out.dir, "/jags_fit_mcmc_convergdiag_geweke.pdf", sep="")
  pdf(pdf.file.geweke)
  geweke.plot(bayes.mod.fit.mcmc)
  dev.off();
  
  
  #effective sample size
  #effectiveSize(bayes.mod.fit.mcmc)
  
  #gelman.diag(bayes.mod.fit.mcmc)
  ##PSRF = 1 if the chains describe perfectly identical distributions. 
  #L&eGFR (2010) recommend that we want PSRF < 1.1 and closer to 1 is better.  
  
  
  #END  
  
}