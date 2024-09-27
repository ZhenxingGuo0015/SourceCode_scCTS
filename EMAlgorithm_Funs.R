library(doParallel)
### first round EM for each gene
### just a case with two mixture normal distribution with same variance
### in one normal distribution, the mean is fixed at 0

### INPUT, delta, init for logfold change (m_g in Eq 3 in manuscript)
### INPUT, tau, init for variance between samples
### INPUT, q, init for probability of CTS show DE signal
### INPUT, Y, gene expression
### INPUT, maxiter, maximum iteration number
### INPUT, tol, EM stop control
### OUPUT, estimates of delta  tau and q for each gene
EM_var1_g.indep <- function(delta, tau, q, sigma, Y, maxiter=1000, tol=1e-5){
  
  diff <- 1
  iter <- 0 

  
  while (diff>tol & iter<maxiter){
    
    ### E step:
    p0 <- dnorm(Y, mean = 0,     sd = sqrt(tau + sigma),log = T)
    p1 <- dnorm(Y, mean = delta, sd = sqrt(tau + sigma),log = T)
    
    denom <- 1/sigma + 1/tau
    
    if(q == 1){
      E.X <- rep(1, length(Y))
    }else if(q==0){
      E.X <- rep(0, length(Y))
    }else{    
      E.X <- q / (q + (1-q) * exp(p0 - p1)  )
    }

    E.delta <- E.X * (Y/sigma + delta/tau)/denom + (1-E.X) * (Y/sigma)/denom # Eq 11 in Supplementary file
    
    E.Xdelta <- E.X * (Y/sigma + delta/tau)/denom # Eq 12 in Supplementary file
    
    E.delta2 <- E.X * (((Y/sigma + delta/tau)/denom)^2 + 1/denom) +  # Eq 13 in supplementary file
      (1-E.X) * (((Y/sigma)/denom)^2 + 1/denom)
    
    ### M step:
    q.new <- mean(E.X) # Eq 17 in supplementary file
    delta.new <- sum(E.Xdelta) / sum(E.X) # Eq 15 in supplementary file
    tau.new <- (sum(E.delta2) - 2*delta.new*sum(E.Xdelta) +  delta.new^2 * sum(E.X) )/length(Y)  # Eq 16 in supplementary file
    
    if(q.new < 1e-300){ # sometimes it trapped, just end it, another naive way for robustness
      q.new <- 0
      break
    }
    ### check:
    diff=sqrt(sum((q.new-q)^2 + (delta.new-delta)^2 + (tau.new-tau)^2))
    q=q.new;
    delta=delta.new;
    tau=tau.new;
    iter=iter+1;
    cat("Iter", iter, ": delta=", delta.new, "; tau=",tau.new,  "; q=", q.new, "; diff=", diff,'\n')

  }
  return(c(delta.new, tau.new, q.new))
}


### Second step to estimate \pi and provide posterior prob
### INPUT: dat.info.sum, result from marker_summarization, summary statistics for each sample
### INPUT: keep.ix, index for genes with postive log fold change
### INPUT: min.cutoff, remove samples with extreme small log fold change for robust estimtation 
### INPUT: max.cutoff, remove samples with extreme small log fold change for robust estimtation 
### INPUT: celltype.ix, cell type index
### INPUT: delta.init, tau.init, q.init, estimated value from EM_var1_g.indep, fixed in this round
### INPUT: maxiter, tol, control for EM algorithm
### OUTPUT: posterior prob, estimates for parameters

EM_est.param_moment <- function(dat.info.sum, keep.ix, min.cutoff, max.cutoff,
                                celltype.ix=2, pi_init, delta.init, tau.init, q.init, 
                                maxiter=1000,tol=10^-10 ){
  ### thetas: delta_g, tau2_g, q_g, and p
  ###         for g = 1, 2, ..., G
  nGene <- dim(dat.info.sum$effect_info[keep.ix,,])[[1]]
  nInd <- dim(dat.info.sum$effect_info[keep.ix,,])[[3]]
  gene.names <- names(keep.ix)
  
  ### initial
  Y <- dat.info.sum$effect_info[keep.ix,celltype.ix,]
  sigma2 <- (dat.info.sum$effect_se_info^2)[keep.ix,celltype.ix,]
  
  sample.ix <- list()
  for(i in 1:nGene){
    sigma2[i,sigma2[i,]<10^-10] <- min(sigma2[i,sigma2[i,]>10^-10])
    
    min.thres <- quantile(Y[i,], min.cutoff)
    max.thres <- quantile(Y[i,], max.cutoff)
    #remain.ix <- which(Y[i,] > min.thres & Y[i,] < max.thres) # commented on August 28, 2023
    remain.ix <- which(Y[i,] >= min.thres & Y[i,] <= max.thres)
    
    sample.ix[[i]] <- remain.ix
  }
  
  delta <- delta.init
  tau2 <- tau.init
  q <- q.init
  p <- pi_init
  
  iter.num <- 0
  diff <- 1
  
  
  while(iter.num < maxiter & diff > tol  ){
    
    ### E-step
    p0 <- t(sapply(1:nGene,  function(x) dnorm(Y[x,sample.ix[[x]]], mean= 0, 
                                               sd = sqrt(tau2[x] + sigma2[x,sample.ix[[x]]]), log = T) ))
    p1 <- t(sapply(1:nGene,  function(x) dnorm(Y[x,sample.ix[[x]]], mean= delta[x], 
                                               sd = sqrt(tau2[x] + sigma2[x,sample.ix[[x]]]), log = T) ))
    
    p.y.c.d1 <- sapply(1:nGene, function(x) sum( Raddlog(p1[x,]+log(q[x]), p0[x,]+log(1-q[x]) ) ) )
    p.y.c.d0 <- sapply(1:nGene, function(x) sum(( p0[x,] ) ) )
    
    if(p == 0){
      exp.d.c <- rep(0,nrow(Y))
    }else{
      exp.d.c <- 1 / (1 + exp(p.y.c.d0 - p.y.c.d1) * (1-p)/p )
    }
    
    exp.dz.c <- t(sapply( 1:nGene, function(x) exp((p1[x,]+log(q[x])) - 
                                                     Raddlog(p1[x,]+log(q[x]), 
                                                             p0[x,]+log(1-q[x]))) )) * exp.d.c
    
    
    
    exp.delta.c <- t(sapply(1:nGene, function(x) exp.dz.c[x,] * (Y[x,sample.ix[[x]]]/sigma2[x,sample.ix[[x]]] + delta[x]/tau2[x])/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x]) + 
                              (1-exp.dz.c[x,]) * (Y[x,sample.ix[[x]]]/sigma2[x,sample.ix[[x]]] )/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x])  ))
    
    exp.dz.delta.c <- t(sapply(1:nGene, function(x)  exp.dz.c[x,] * (Y[x,sample.ix[[x]]]/sigma2[x,sample.ix[[x]]] + delta[x]/tau2[x])/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x]) ))
    
    exp.delta2.c <-  t(sapply(1:nGene, function(x)  exp.dz.c[x,] * ( ((Y[x,sample.ix[[x]]]/sigma2[x,sample.ix[[x]]] + delta[x]/tau2[x])/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x]))^2 + 1/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x])) + 
                                (1-exp.dz.c[x,]) * ( ((Y[x,sample.ix[[x]]]/sigma2[x,sample.ix[[x]]] )/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x]))^2 + 1/(1/sigma2[x,sample.ix[[x]]] + 1/tau2[x]) ) ))  
    
    ### M-step
    #  delta.new <- rowSums(exp.dz.delta.c) / rowSums(exp.dz.c)
    #  delta.new[is.na(delta.new)] <- 0
    #  delta.new <- delta
    #  tau2.new <- rowMeans(exp.delta2.c) - 2*delta.new*rowMeans(exp.dz.delta.c) + (delta.new^2)*rowMeans(exp.dz.c)
    #  tau2.new <- tau2
    #  q.new <- rowSums(exp.dz.c) / (exp.d.c*nInd)
    #  q.new[q.new>1] <- 1
    #  q.new[q.new<0] <- 0
    #  q.new[is.na(q.new)] <- 0
    p.new <- mean(exp.d.c) * pi_init
    
    ### converge check
    
    iter.num <- iter.num + 1
    diff <-  sqrt( ( ( mean( (p.new - p)^2)   ) ))
    # cat(iter.num,',',p,', ',q[1],', ',delta[1],', ',tau2[1],', ', diff,'\n ')
    # delta <- delta.new
    # tau2 <- tau2.new
    # q <- q.new
    p <- p.new
    
  }
  
  param.res <- cbind(delta,tau2, q)
  rownames(param.res) = names(exp.d.c) = rownames(exp.dz.c)<- gene.names

  return(list('param'=param.res, 'pi'=p, 'pp.d1'=exp.d.c, 'pp.z1' = exp.dz.c))
  
}

### INPUT: data.info.sum, output of function marker_summarization
###        contains summarized statistics of each sample
### INPUT: effect_thres, threshold for filtering genes with negative mean (m_gk in Eq 3 in manuscirpt), 
### INPUT: maxiter, maximum iteration number 
### INPUT: tol, EM stop control
### INPUT: numCores, number of cores for parallel computation
### INPUT: min.cutoff, remove samples with extreme small log fold change for robust estimtation 
###                    quantile 0.05 
### INPUT: max.cutoff, remove samples with extreme large log fold change for robust estimtation 
###                    quantile 0.95
### PURPOSE: EM algorithm main function
### OUTPUT:  posterior probability and estimates of parameters
csmarker_moment <- function( dat.info.sum, effect_thres= 0,#0.01, 
                             maxiter=1000, tol=10^-3,numCores=1,
                             min.cutoff = 0.05, max.cutoff = 0.95,
                             min.q = 0
                             ){
  
  registerDoParallel(numCores) # register cores number for parallel computation
  
  nGene <- dim(dat.info.sum$effect_info)[[1]] # gene number
  name.gene <- dimnames(dat.info.sum$effect_info)[[1]] # gene name
  nCell <- dim(dat.info.sum$effect_info)[[2]] # cell type number
  name.cell <- dimnames(dat.info.sum$effect_info)[[2]] # cell type name
  nInd <- dim(dat.info.sum$effect_info)[[3]] # sample/subjects number
  name.ind <- dimnames(dat.info.sum$effect_info)[[3]] # sample/subject name
  
  res.new <- list() # used to store all result
  
  res.all <- foreach( cell.ix = 1:nCell) %dopar%{ # loop for each cell type
    
    res.store <- matrix(NA,nrow=nGene,ncol=3) # init for storing result
    rownames(res.store) <- name.gene # assign gene name, row is gene
    cat(cell.ix, '\n') 
    
    for(i in 1:nGene){ # loop for each gene,  
      X <- dat.info.sum$effect_info[i,cell.ix,] # log2 fold change (a.k.a Y in Eq3 first line)
      s2 <- dat.info.sum$effect_se_info[i,cell.ix,]^2 # sigma^2 in Eq 3 first line
      s2[s2==0] <- min(s2[s2>0]) # for robust calculation purpose, replace 0 value with minimum non-zero value
      cat(i,',')
      de.ix <- dat.info.sum$wilcox.fdr_info[i, cell.ix, ] < 0.05 #  CTS genes identified by wilcoxon rank sum test
      
      min.thres <- quantile(X , min.cutoff) # thres to remove samples with extreme small value
      max.thres <- quantile(X, max.cutoff) # thres to remvoe samples with extreme big value
      

      remain.ix <- which(X >= min.thres & X <= max.thres) # remove samples 

      X <- X[remain.ix]
      de.ix <- de.ix[remain.ix]
      s2 <- s2[remain.ix]
      
      # estimate initial values for EM. 
      # the estimation is not very robust, since the normal distribution assumption does not hold well sometimes,  
      # especially with small sample size
      # if there more than five samples, use them to estimate the mean and variance of the samples showing DE signal
      if(sum(de.ix, na.rm = TRUE) > 5 ){     
        var.1 <- var(X[de.ix])
        if(var.1 < mean(s2[de.ix]) ){ # sometimes, the total variance is smaller than variance within sample
          var.1 <- var.1/10000        # we make it very small to start, indicates the total variance is close to 0
        }else{                        # if not, then the variance between samples = total variance - within sample variance
          var.1 <- var.1 -  mean(s2[de.ix])
        }
        
        tryCatch({ # may still exist error, if so, re-assign the delta with quantile 0.8
          res.store[i,]<-  EM_var1_g.indep(delta=mean(X[de.ix]),tau=rep(var.1,length(X)),q=0.5,
                                           sigma = s2, Y = X, maxiter = maxiter, tol=tol)
        },
        error=function(e){cat("ERROR:",conditionMessage(e),"\n")}) 
        
        if(sum(is.na(res.store[i,])) == 3){ # this condition represents error for above EM, re-run with another initial value 
          var.1 <- var(X)
          if(var.1 <= mean(s2) ){
            var.1 <- var.1/10000
          }else{
            var.1 <- var.1 -  mean(s2)
          }
          res.store[i,] <-  EM_var1_g.indep(delta=quantile(X, 0.8),tau=rep(var.1,length(X)),q=0.5,
                                           sigma = s2, Y = X, maxiter = maxiter, tol=tol)
        }
        
      }else{ # if there are no more than five samples showing DE, then use all samples to estimate initial values
        var.1 <- var(X)
        if(var.1 <= mean(s2) ){
          var.1 <- var.1/10000
        }else{
          var.1 <- var.1 -  mean(s2)
        }
        res.store[i,] <- EM_var1_g.indep(delta = quantile(X, 0.8),tau=rep(var.1,length(X)),q=0.5,
                                         sigma = s2, Y = X, maxiter = maxiter, tol=tol)
      }
    }
    
    
    ### commented on May 26, 2024
    # keep.ix <- which(res.store[,1] > effect_thres  &
    #                    res.store[,3] > 1/nInd &  ## with high prevalence
    #                    res.store[,2] > 0) ## remove negative CTS gene
    # 
    # pi_init <- mean(res.store[,1] > effect_thres  &
    #                   res.store[,3] > 1/nInd  &
    #                   res.store[,2] > 0) ## estimate initial value for pi
    ###
    
    keep.ix <- which(res.store[,1] > effect_thres  &
                       res.store[,3] > min.q &  ## with high prevalence
                       res.store[,2] > 0) ## remove negative CTS gene
    
    pi_init <- mean(res.store[,1] > effect_thres  &
                      res.store[,3] > min.q &
                      res.store[,2] > 0) ## estimate initial value for pi
    
    
    
    
    
    
    # ####### check what happens to the missed markers
    # rm.ix.effect = sum(res.store[,1] < effect_thres)
    # sum(res.store[,1] < effect_thres & 
    #       res.sim$FullParas$D_gk[, cell.ix] == 1)
    # rm.ix.preva = sum(res.store[,3] <= 1/nInd)
    # 
    # ix.marker.lowEstiq = which(res.store[,3] <= 1/nInd & 
    #                              res.sim$FullParas$D_gk[, cell.ix] == 1)
    # 
    # ## Seems there is an underestimate for the prevalence of majority of markers
    # 
    # pdf(file = paste0("./RawPattern_ID14/Examine_preva_CT", cell.ix, ".pdf"),
    #     height = 4, width = 4)
    # par(mfrow = c(1,1), mai = c(0.8, 0.8, 0.4, 0.2))
    # plot(res.sim$FullParas$q_gk[ix.marker.lowEstiq, cell.ix], 
    #      res.store[ix.marker.lowEstiq,3], 
    #      xlim = c(0, 0.4), ylim = c(0, 0.4),
    #      ylab = "Estimated prevalence",
    #      xlab = "True prevalence", pch = 16, col = "#00000060",
    #      main = paste0("Cell type ", cell.ix, ": N.miss = ", length(ix.marker.lowEstiq)))
    # abline(0,1, col = "red")
    # dev.off()
    # 
    # pdf(file = paste0("./RawPattern_ID14/Examine_logFC_CT", cell.ix, ".pdf"),
    #     height = 4, width = 4)
    # par(mfrow = c(1,1), mai = c(0.8, 0.8, 0.4, 0.2))
    # plot(m_gk[ix.marker.lowEstiq, cell.ix], 
    #      res.store[ix.marker.lowEstiq, 1], 
    #      xlab = "True logFC",
    #      ylab = "Estimated logFC",
    #      pch = 16, col = "#00000060",
    #      main = paste0("Cell type ", cell.ix, ": N.miss = ", length(ix.marker.lowEstiq)))
    # abline(0,1, col = "red")
    # dev.off()
    # #######
    # 
    
    
    
    # ##### update on August 28, 2023
    # keep.ix <- which(res.store[,1] > effect_thres  & 
    #                    res.store[,3] > min.q &  ## with high prevalence
    #                    res.store[,2] > 0) 
    # 
    # pi_init <- mean(res.store[,1] > effect_thres  & 
    #                   res.store[,3] > min.q  & 
    #                   res.store[,2] > 0) 
    # ######
    
    
    ### this function used to estimate pi and calculate posterior probability
    res.tmp <- EM_est.param_moment(dat.info.sum=dat.info.sum, 
                                   keep.ix = keep.ix,
                                   celltype.ix = cell.ix,  
                                   pi_init = pi_init, 
                                   min.cutoff = min.cutoff,
                                   max.cutoff = max.cutoff,
                                   delta.init = res.store[keep.ix,1], 
                                   tau.init = res.store[keep.ix,2],
                                   q.init = res.store[keep.ix,3], 
                                   maxiter = maxiter, tol = tol)
    
    res.tmp$remain.subj = remain.ix ### added on August 28, 2023
    res.tmp
  }
  res.new <- res.all
  names(res.new) <- name.cell
  return(res.new)
}



####### created on May 28, 2024
csmarker_moment.2 <- function( dat.info.sum, effect_thres = 0,#0.01, 
                             maxiter=1000, tol=10^-3,numCores=1,
                             min.cutoff = 0.05, max.cutoff = 0.95,
                             min.q = 0,
                             delta.ini = 0.01,
                             q.ini = 0.5,
                             pi.ini = 0.5,
                             NeedInitial = FALSE,
                             InitialDelta.Bydata = FALSE,
                             InitialPI.Bydata = TRUE
                             ){
  
  registerDoParallel(numCores) # register cores number for parallel computation
  
  nGene <- dim(dat.info.sum$effect_info)[[1]] # gene number
  name.gene <- dimnames(dat.info.sum$effect_info)[[1]] # gene name
  nCell <- dim(dat.info.sum$effect_info)[[2]] # cell type number
  name.cell <- dimnames(dat.info.sum$effect_info)[[2]] # cell type name
  nInd <- dim(dat.info.sum$effect_info)[[3]] # sample/subjects number
  name.ind <- dimnames(dat.info.sum$effect_info)[[3]] # sample/subject name
  
  res.new <- list() # used to store all result
  
  res.all <- foreach( cell.ix = 1:nCell) %dopar%{ # loop for each cell type
    
    res.store <- matrix(NA,nrow=nGene,ncol=3) # init for storing result
    rownames(res.store) <- name.gene # assign gene name, row is gene
    cat(cell.ix, '\n') 
    
    for(i in 1:nGene){ # loop for each gene,  
      X <- dat.info.sum$effect_info[i,cell.ix,] # log2 fold change (a.k.a Y in Eq3 first line)
      s2 <- dat.info.sum$effect_se_info[i,cell.ix,]^2 # sigma^2 in Eq 3 first line
      s2[s2==0] <- min(s2[s2>0]) # for robust calculation purpose, replace 0 value with minimum non-zero value
      cat(i,',')
      de.ix <- dat.info.sum$wilcox.fdr_info[i, cell.ix, ] < 0.05 #  CTS genes identified by wilcoxon rank sum test
      
      min.thres <- quantile(X , min.cutoff) # thres to remove samples with extreme small value
      max.thres <- quantile(X, max.cutoff) # thres to remvoe samples with extreme big value
      
      
      remain.ix <- which(X >= min.thres & X <= max.thres) # remove samples 
      
      X <- X[remain.ix]
      de.ix <- de.ix[remain.ix]
      s2 <- s2[remain.ix]
      
      if(!NeedInitial){
        # estimate initial values for EM. 
        # the estimation is not very robust, since the normal distribution assumption does not hold well sometimes,  
        # especially with small sample size
        # if there more than five samples, use them to estimate the mean and variance of the samples showing DE signal
        if(sum(de.ix, na.rm = TRUE) > 5 ){     
          var.1 <- var(X[de.ix])
          if(var.1 < mean(s2[de.ix]) ){ # sometimes, the total variance is smaller than variance within sample
            var.1 <- var.1/10000        # we make it very small to start, indicates the total variance is close to 0
          }else{                        # if not, then the variance between samples = total variance - within sample variance
            var.1 <- var.1 -  mean(s2[de.ix])
          }
          
          tryCatch({ # may still exist error, if so, re-assign the delta with quantile 0.8
            res.store[i,]<-  EM_var1_g.indep(delta=mean(X[de.ix]),tau=rep(var.1,length(X)),q=0.5,
                                             sigma = s2, Y = X, maxiter = maxiter, tol=tol)
          },
          error=function(e){cat("ERROR:",conditionMessage(e),"\n")}) 
          
          if(sum(is.na(res.store[i,])) == 3){ # this condition represents error for above EM, re-run with another initial value 
            var.1 <- var(X)
            if(var.1 <= mean(s2) ){
              var.1 <- var.1/10000
            }else{
              var.1 <- var.1 -  mean(s2)
            }
            res.store[i,] <-  EM_var1_g.indep(delta=quantile(X, 0.8),tau=rep(var.1,length(X)),q=0.5,
                                              sigma = s2, Y = X, maxiter = maxiter, tol=tol)
          }
          
        }else{ # if there are no more than five samples showing DE, then use all samples to estimate initial values
          var.1 <- var(X)
          if(var.1 <= mean(s2) ){
            var.1 <- var.1/10000
          }else{
            var.1 <- var.1 -  mean(s2)
          }
          res.store[i,] <- EM_var1_g.indep(delta = quantile(X, 0.8),tau=rep(var.1,length(X)),q=0.5,
                                           sigma = s2, Y = X, maxiter = maxiter, tol=tol)
        }
        
      }else if(NeedInitial){
        var.1 <- var(X)
        if(var.1 <= mean(s2) ){
          var.1 <- var.1/10000
        }else{
          var.1 <- var.1 -  mean(s2)
        }
        tau.ini = var.1
      
        if(InitialDelta.Bydata){
          delta.ini = quantile(X, 0.8)
        }
        res.store[i,] <- EM_var1_g.indep(delta = delta.ini, tau=rep(tau.ini,length(X)), q=q.ini, 
                                         sigma = s2, Y = X, maxiter = maxiter, tol=tol)
      }
      
    }
    
    
    keep.ix <- which(res.store[,1] > effect_thres  &
                       res.store[,3] > min.q &  ## with high prevalence
                       res.store[,2] > 0) ## remove negative CTS gene
    
    if(InitialPI.Bydata){
      pi.ini <- mean(res.store[,1] > effect_thres  &
                        res.store[,3] > min.q &
                        res.store[,2] > 0) ## estimate initial value for pi
      
    }
    

    ### this function used to estimate pi and calculate posterior probability
    res.tmp <- EM_est.param_moment(dat.info.sum=dat.info.sum, 
                                   keep.ix = keep.ix,
                                   celltype.ix = cell.ix,  
                                   pi_init = pi.ini, 
                                   min.cutoff = min.cutoff,
                                   max.cutoff = max.cutoff,
                                   delta.init = res.store[keep.ix,1], 
                                   tau.init = res.store[keep.ix,2],
                                   q.init = res.store[keep.ix,3], 
                                   maxiter = maxiter, tol = tol)
    
    res.tmp$remain.subj = remain.ix ### added on August 28, 2023
    res.tmp
  }
  res.new <- res.all
  names(res.new) <- name.cell
  return(res.new)
}

#######


### IGNORE! calculate posterior probability, but not used in current manuscript
### this is used in previous version, so just IGNORE
marker.pp <- function(dat.info, info.pop){
  nGene <- dim(dat.info$effect)[[1]]
  name.gene <-  dimnames(dat.info$effect)[[1]]
  nCell <- dim(dat.info$effect)[[2]]
  name.cell <-  dimnames(dat.info$effect)[[2]]
  ct.names <- names(info.pop)
  pp <- list()
  for( i in 1:nCell ){
    gene.keep <- rownames(info.pop[[i]]$param)
    
    # gene.keep <- gene.keep[gene.keep %in% names(which(p.adjust(dat.info$pval[,i],'fdr') < 0.05))]
    
    Y <- dat.info$effect[gene.keep,i]
    s2 <- dat.info$effect_se[gene.keep,i]^2
    for( g.ix in 1:length(gene.keep)){
      s2[ s2 < 10^-20] <- min(s2[ s2> 10^-20])  
    }
    
    p1 <- dnorm(Y, mean= info.pop[[i]]$param[gene.keep,1], sd= sqrt(info.pop[[i]]$param[gene.keep,2] + s2) )
    p0 <- dnorm(Y, mean= 0,                                sd= sqrt(info.pop[[i]]$param[gene.keep,2] + s2) )
    
    
    prior.d.p <- info.pop[[i]]$pp.d1 * 0.999
    prior.z.p <- info.pop[[i]]$param[gene.keep,3] * 0.999 
    pp[[ct.names[i]]] <- p1* prior.z.p * prior.d.p * prior.z.p / (p1* prior.z.p * prior.d.p + 
                                            p0 * (1-prior.z.p) * prior.d.p + 
                                            p0 * (1-prior.d.p) )
    names(pp[[ct.names[i]]]) <- gene.keep
    
  }
  
  return(pp)
  
}


### two functions for calculating sum of log P
Raddlog <- function(a, b) {
  result <- rep(0, length(a))
  idx1 <- a>b+200
  result[idx1] <- a[idx1]
  
  idx2 <- b>a+200
  result[idx2] <- b[idx2]
  
  idx0 <- !(idx1|idx2)
  result[idx0] <- a[idx0] + log1p(exp(b[idx0]-a[idx0]))
  result
}
Rsumlog <- function(a) {
  s <- a[1]
  for(i in 2:length(a))
    s <- Raddlog(s, a[i])
  s
}
