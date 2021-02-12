## Current lme4 versions have a problem with false
## convergence warnings that are not relevant.
## This function implements a quick test suggested
## by Ben Bolker (package maintainer) and another
## test I use.
glmerConverge <- function(mod){
  maxgrad <- max(abs(with(mod@optinfo$derivs, solve(Hessian,gradient))))
  nmod <- refit(mod)
  coefdiff <- all.equal(summary(mod)$coef[, 3],
                        summary(nmod)$coef[,3])

  return(list(maxgrad=maxgrad, coefdiff=coefdiff))
}


## function to test for overdispersion in a model
overDispTest <- function(mod)
{
  chisq <- sum(resid(mod, type='pearson')^2)
  df.resid <- diff(getME(mod, 'devcomp')$dims[c('nth', 'nmp')])
  phi <- chisq/df.resid
  p <- pchisq(chisq, df.resid, lower.tail=FALSE)
  return(c(chisq=chisq, df.resid=df.resid, phi=phi, p=p))
}

## function to simulate Chisq
ChiSim <- function(mod1, mod0, mod0sim, use.u, max.try=100) {
    converge <- 0 ## error trap to deal with non-convergence issues
    ## try to simulate a max
    while(converge < max.try)
        {
            newY <- simulate(mod0sim, 1, use.u=use.u) ## simulate new Y
            ## refit H1 model to null data
            newmod1 <- try(refit(mod1, newY), silent=TRUE)
            ## rerun if mod1 did not fit properly
            if(class(newmod1) == 'try-error')
              {
                converge <- converge + 1
                cat(converge)
                next
            }
            else
               if(length(newmod1@optinfo$warnings) != 0)
                {        converge <- converge+1
                    cat(converge)
                    next
                }

            ## refit H0 model to null data
            newmod0 <- try(refit(mod0, newY), silent=TRUE)
            ## rerun if mod0 did not fit properly
            if(class(newmod0) == 'try-error')
              {
                converge <- converge + 1
                cat(converge)
                next
                }
            else
              if(length(newmod1@optinfo$warnings) != 0)
              {        
                converge <- converge+1
                cat(converge)
                next
              }
            ## accept iteration if model fits converged
            else
                if(!(class(newmod1)=='try-error' |
                     class(newmod0)=='try-error'))
                    {

                        converge <- 1000
                    }
        }

    ## calculate the logLik ratio for simulated data
    simChi <- 2*(logLik(newmod1)[1] - logLik(newmod0)[1])
    ##return(list(simChi=simChi, nmod0=newmod0, nmod1=newmod1))
    return(simChi)
}

## function to run a parametric bootstrap on a GLMM, using parallel
## processing. Experimental so use with care.
## handles errors silently unless the model just cannot be bootstrapped.
simAnova <- function(mod1, mod0, nsim, REMLsim=TRUE, ncore=1, use.u=FALSE,
                     max.try=100){
  require(parallel)

  cl <-  makeCluster(ncore)
  RNGkind <- "L'Ecuyer-CMRG"
  clusterSetRNGStream(cl)

  on.exit({stopCluster(cl); print('clusters closed on exit')})

  if(REMLsim)
    mod0sim <- mod0

  if(isREML(mod1)){
    mod1 <- update(mod1, REML=FALSE)
    mod0 <- update(mod0, REML=FALSE)
  }

  if(!REMLsim)
    mod0sim <-  mod0

  obsChi <- 2*(logLik(mod1)[1] - logLik(mod0)[1])

  clusterExport(cl=cl, varlist=list('mod0', 'mod1', 'mod0sim', 'use.u', 'ChiSim'),
                envir=environment())
  clusterEvalQ(cl, library(lme4))

  sims <- parSapply(cl, 1:nsim, function(i, mod1, mod0, mod0sim, use.u)
                    {
                        ChiSim(mod1=mod1, mod0=mod0, mod0sim=mod0sim, use.u=use.u)
                    }, mod1=mod1, mod0=mod0, mod0sim=mod0sim, use.u=use.u,
                    simplify=TRUE)
  p <-  (sum(sims > obsChi)+1)/(nsim +1)
  df <- as.numeric(mod1@devcomp$dims['p'] - mod0@devcomp$dims['p'])
  return(list(ChiSq=obsChi, df=df,  p=p, simChi=sims))}

