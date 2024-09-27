SimPara <- function(nCells, nGenes,
                    N = 16, K = 4,
                    Subj.Prop = rep(1/N, N),
                    CT.Prop = rep(1/K, K), 
                    pi_k = rep(0.05, K),
                    multiCT.Marker = TRUE,
                    seed = 12345){
  
  set.seed(seed) 
  n.cell.bySub = as.vector(rmultinom(1, size = nCells,
                                     prob = Subj.Prop
                                     )
                           )
  
  # In this case, maker may not be detected as maker by Wilcoxon
  NC = Sim.NC(n.cell.bySub = n.cell.bySub,
              Da_uniform = CT.Prop,
              Da_diff =  CT.Prop)
  
  ### simulate D_gk
  if(!multiCT.Marker){
    ## B. one gene can only be the marker of one celltype
    set.seed(seed)
    markerG = matrix(FALSE, nrow = nGenes, ncol = K)
    mark.idx = NULL
    D_gk = matrix(0, nrow = nGenes, ncol = K)
    for(ik in 1:K){
      if(ik==1){
        this.pi = rep(pi_k[ik],nGenes)
        d_g = rbinom(length(this.pi), size = 1, prob = this.pi)
        markerG[d_g == 1, ik] = TRUE
        mark.idx = which(d_g==1)
      }else{
        remainG = setdiff(1:nGenes, mark.idx) # which(!markerG[, ik-1])
        this.pi = rep(pi_k[ik], length(remainG))
        d_g = rbinom(length(this.pi), size = 1, prob = this.pi)
        markerG[remainG, ik] = as.logical(d_g)
        mark.idx = union(mark.idx, remainG[d_g==1])
      }
    }
    D_gk[markerG] = 1
    cat( "Adjusted pi_k is: ", colMeans(D_gk),  sep = "\t")
  }else{
    set.seed(125)
    D_gk = matrix(0, nrow = nGenes, ncol = K)
    for (ik in 1:K) {
      tmp.flag = rbinom(nGenes, size = 1, prob = pi_k[ik])
      D_gk[tmp.flag == 1, ik] = 1
    }
  }
  
  colnames(D_gk) = paste0("Celltype_", 1:4)
  rownames(D_gk)=paste0("gene_",1:nrow(D_gk))
  
  res = list(NC = NC, D_gk = D_gk)
}

SimMeanLFC <- function( D_gk, CTs, M_gk = NULL, Delta_Gk = NULL,
                        Log2.FC = NULL, seed = 12345){
  
  ### generate matrix of m_gk, either based on estimated Delta, 
  ###                          or by LFC directly from data
  if(length(Delta_Gk) > 0){
    Delta = Delta_Gk[, CTs]
    set.seed(seed)
    m_gk = matrix(0, nrow = nrow(D_gk), ncol = ncol(D_gk))
    for (ct.ix in 1:ncol(D_gk)) {
      id.marker = which(D_gk[, ct.ix] ==1)
      m_gk[id.marker, ct.ix] = sample(Delta[Delta[, ct.ix] > 0.01, ct.ix],
                                      length(id.marker),
                                      replace = FALSE)

      L = range(m_gk[id.marker, ct.ix])[1]
      U = range(m_gk[id.marker, ct.ix])[2]
      m_gk[-id.marker, ct.ix] = runif((nrow(D_gk)-length(id.marker)), L, U)
    }
  }else if(length(M_gk) > 0){
    set.seed(seed)
    m_gk = M_gk[, CTs]
    m_gk[m_gk <= 0] = runif(sum(m_gk <= 0), 0.05, 0.08) # make sure postive m_gk
  }else if(length(Log2.FC) > 0){
    set.seed(seed)
    m_gk = Log2.FC[, CTs]
    m_gk[m_gk <= 0] = runif(sum(m_gk <= 0), 0.05, 0.08) # make sure postive m_gk
  }
  
  colnames(m_gk) = paste0("Celltype_", 1:ncol(D_gk))
  rownames(m_gk)=paste0("gene_",1:nrow(D_gk))
  return(m_gk)
}

SimPrevalenceByLFC <- function(m_gk, ctypes, Downrate = 1){
  #### Simulate prevalence based on LFC
  
  # Load estimated relationship between prevalence and LFC in real lupus data
  load("./data/Coefs_log2qOdds_vs_log2LFC_lupus_onlyCtrl.rda")
  
  # generate prevalence based on estimated relationship
  CT.Coefs = Coefs[ctypes, ]
  logit.q = matrix(NA, ncol = ncol(m_gk), nrow = nrow(m_gk))
  for (ct.ix in 1:ncol(m_gk)) {
    logit.q[, ct.ix] = Downrate*CT.Coefs[ct.ix, 1] + log2(m_gk[, ct.ix])* CT.Coefs[ct.ix, 2]
  }
  q_gk = 2^logit.q/(1+2^logit.q)
  q_gk[q_gk==0] = 0.0001
  q_gk[q_gk==1] = 0.9999
  return(q_gk)
}


SimData <- function(NC, alpha, D_gk, m_gk, tau_gk, q_gk,
                    multiCT.Marker,
                    seed = 12345){
  
  N = nrow(NC); K = ncol(NC)
  nGenes = nrow(D_gk)
  
  set.seed(seed)
  ### 1. first fix k and g, simulate Z_gik for all i
  Z_gi = vector("list", length = K)
  Z.P = D_gk*q_gk
  for (k in 1: K) {
    this.P = as.vector(Z.P[, k])
    tmp = matrix(rep(this.P, each = N), nrow = length(this.P), byrow = TRUE)
    Z_gi[[k]] = matrix(rbinom(length(tmp), size = 1, prob = as.vector(tmp)), 
                       nrow = nrow(tmp), byrow = FALSE)
  }
  
  
  ### 2. simulate by subject
  set.seed(seed)
  beta = vector("list", length = N)
  mu_g = disp_g = Counts = sf = vector("list", length = N)
  for (i in 1:N) {
    mu_g[[i]] = disp_g[[i]] = 
      Counts[[i]] = sf[[i]] = vector("list", length = K)
    beta[[i]] = matrix(NA, nrow = nGenes, ncol = K)
    for (k in 1:K) {
      for (g in 1:nGenes) {
        if(D_gk[g, k] == 1 & Z_gi[[k]][g,i] == 1){
          beta[[i]][g, k] = rnorm(1, mean = 10*m_gk[g, k], sd = tau_gk[g, k]) ### log2FC
          if(beta[[i]][g, k] < 0)
            beta[[i]][g, k] = m_gk[g, k]
        }else if((D_gk[g, k] == 1 & Z_gi[[k]][g, i] == 0)){
          beta[[i]][g, k]  =  0 # markers not show DE
        }else if(D_gk[g, k] == 0){
          ### nonmarkers
          beta[[i]][g, k]  = rnorm(1, 0, 0.01) # 0
        }
      } 
    }
    
    if(multiCT.Marker){
      ### adjust value of beta to make sure multiple celltype markers hold
      tmp.beta = beta[[i]]
      idx = which(rowSums(tmp.beta > 0) >= 2)
      for (ig in idx) {
        i.k = which(tmp.beta[ig, ]>0)
        tmp.beta[ig, i.k] = mean(tmp.beta[ig, ])
      }
      beta[[i]] = tmp.beta
    }
    
    
    for (k in 1:K) {
      mu_g[[i]][[k]] = 2^{alpha[, i] + as.vector(beta[[i]][, k])}
     # mu_g[[i]][[k]] = (2^{alpha[, i]}-1)*2^{as.vector(beta[[i]][, k])} # modified on Sept 6
      
      
      disp_g[[i]][[k]] = exp(rnorm( nGenes,
                                   mean = log(0.01/(mu_g[[i]][[k]] + 0.1) + 0.05),
                                   sd = 0.001))/10
      
      mu = matrix(rep(mu_g[[i]][[k]], each = NC[i, k]),
                  nrow = nGenes, byrow = TRUE)
      disp = matrix(rep(disp_g[[i]][[k]], each = NC[i, k]),
                    nrow = nGenes, byrow = TRUE)
      
      ### size factors
      sf[[i]][[k]] = runif(NC[i, k], min = 0.5, max = 5)
      inflate.mu = sweep(mu, 2, sf[[i]][[k]], FUN = "*")
      ###
      
      Counts[[i]][[k]] = matrix(rnbinom(length(mu),
                                        mu = as.vector(inflate.mu),
                                        size = as.vector(1/disp)),
                                nrow = nGenes, byrow = FALSE)
    }
    
  }
  
  #### reorganize the data structure
  nInds = nrow(NC)
  nCT = ncol(NC)
  
  array.tmp <- array(NA, dim=c(nGenes, N,  K))
  dimnames(array.tmp) <- list(paste0("gene_", 1:nGenes),
                              rownames(NC),
                              colnames(NC))
  
  mu = phi = Z = Beta = array.tmp
  for(sub.ix in 1: N){
    mu[, sub.ix, ] = matrix(unlist(mu_g[[sub.ix]]), 
                            ncol = length(mu_g[[sub.ix]]), byrow = FALSE)
    phi[, sub.ix, ] = matrix(unlist(disp_g[[sub.ix]]), 
                             ncol = length(disp_g[[sub.ix]]), byrow = FALSE)
    Beta[, sub.ix, ] = matrix(unlist(beta[[sub.ix]]), 
                              ncol = length(beta[[sub.ix]]), byrow = FALSE)
  }
  
  
  for(ct.ix in 1:K){
    Z[, , ct.ix ] = Z_gi[[ct.ix]]
  }
  
  
  ### put all count from all subjects and all celltypes together  
  library(SingleCellExperiment)
  Count.wide = NULL
  for (i in 1:N) {
    tmp.i = NULL
    sf.i = NULL
    for (k in 1:K) {
      if(k==1){
        tmp.i = Counts[[i]][[1]]
        sf.i = sf[[i]][[1]]
      }else{
        tmp.i = cbind(tmp.i, Counts[[i]][[k]])
        sf.i = c(sf.i, sf[[i]][[k]])
      }
    }
    if(i == 1){
      cellType.all = rep(paste0("CellType_", 1:K), NC[i, ])
      subject.all = rep(paste0("Subject_", i), sum(NC[i, ]))
      Count.wide = tmp.i
      sf.wide = sf.i
    }else{
      cellType.all = c(cellType.all, rep(paste0("CellType_", 1:K), NC[i, ]))
      subject.all =  c(subject.all, rep(paste0("Subject_", i), sum(NC[i, ])))
      Count.wide = cbind(Count.wide, tmp.i)
      sf.wide = c(sf.wide, sf.i)
    }
  }
  
  rownames(Count.wide) = paste0("gene_", 1:nrow(Count.wide))
  rownames(q_gk) = rownames(Count.wide)
  rownames(D_gk) = rownames(Count.wide)
  
  
  coldata <- DataFrame(
    subject = subject.all,
    celltype= cellType.all,
    sf = sf.wide)
  
  rowdata <- DataFrame(
    Genes = rownames(Count.wide))
  
  
  
  library(SummarizedExperiment)
  sce = SummarizedExperiment(assays=list(counts=Count.wide),
                             colData =coldata,
                             rowData =rowdata )
  
  
  res.raw = list(Counts = Counts, mu = mu, 
                 phi = phi, D_gk = D_gk, 
                 q_gk = q_gk, Z = Z,
                 beta = Beta)
  
  return(list(FullParas = res.raw, Counts.wide = sce))
  
}



Sim.NC <- function(n.cell.bySub,
                   Da_uniform = c(2, 2, 2, 2) , 
                   Da_diff = c(0.5, 0.5, 4, 2.5), 
                   nTimes = 100,
                   seed = 12345){
  
  set.seed(seed)
  K = length(Da_uniform)
  N = length(n.cell.bySub)
  NC.allTimes = matrix(0, nrow = N, ncol = K)
  for (isim in 1:nTimes) {
    NC.tmp = matrix(NA, N, K)
    iter = 0
    for (sub.ix in 1:N) {
      iter = iter + 1
      if(iter <= 4){
        Da = Da_uniform  
        c.prop = round(rdirichlet(n = 1, alpha = Da), digits = 2)
        NC.tmp[iter, ] = round(sweep(c.prop, 1, n.cell.bySub[iter], FUN = "*"))
      }else if(iter > 4 & iter <= 8){
        Da = Da_uniform 
        c.prop = round(rdirichlet(n = 1, alpha = Da), digits = 2)
        NC.tmp[iter, ] = round(sweep(c.prop, 1, n.cell.bySub[iter], FUN = "*"))
      }else if(iter >8 & iter <= 12){
        Da = Da_diff  
        c.prop = round(rdirichlet(n = 1, alpha = Da), digits = 2)
        NC.tmp[iter, ] = round(sweep(c.prop, 1, n.cell.bySub[iter], FUN = "*"))
      }else if(iter > 12){
        Da = Da_uniform
        c.prop = round(rdirichlet(n = 1, alpha = Da), digits = 2)
        NC.tmp[iter, ] = round(sweep(c.prop, 1, n.cell.bySub[iter], FUN = "*"))
      }
      
    }
    #NC.tmp[NC.tmp ==0 ] = 5
    colnames(NC.tmp) = paste0("CellType_", 1:K);
    rownames(NC.tmp) = paste0("Subject_", 1:N)
    
    NC.allTimes = NC.allTimes + NC.tmp
  }
  NC = round(NC.allTimes/nTimes)
  NC[NC <= 1] = 5
  return(NC)
}

