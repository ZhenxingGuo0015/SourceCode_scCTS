Plot_ROC <- function(flag, Pvals, models,ltypes, cols = NULL, leg){
  library(ROCR)
  ROC.DE <- function(DE.gs, pval) {
    pred = prediction(1-pval, DE.gs)
    perf = performance(pred,"tpr","fpr")
    perf
  }
  
  AUC.DE <- function(pval, DE.gs) {
    ii = !is.na(pval)
    library(ROCR)
    pred = prediction(1-pval[ii], DE.gs[ii])
    auc.tmp = performance(pred,"auc")
    auc = as.numeric(auc.tmp@y.values)
    auc
  }
  
  nmodel = length(Pvals)
  if(is.null(cols)){
    cols = 1:nmodel
  }
  # all.flag = rep(flag, ncol(Pvals[[1]]))
  
  all.flag = as.vector(flag)
  auc = rep(NA, nmodel)
  for(i in 1:nmodel){
    this.pval = as.vector(Pvals[[i]])
    idx = is.na(this.pval)
    roc = ROC.DE(all.flag[!idx], this.pval[!idx])
    ###
    auc[i] = round(AUC.DE(this.pval, all.flag), 3)
    ###
    par(cex.axis=1.7)
    if (i == 1){
      plot(roc, xlim = c(0,1), lty = ltypes[i], col = cols[i], ylim = c(0,1),
           lwd = 5)
    }else{
      plot(roc, add = TRUE, col = cols[i], lty = ltypes[i],
           lwd = 5)

    }
  }
  if(leg)
    legend("bottomright", legend = paste0(models, ": ", auc), 
           col = cols, lty=ltypes, lwd = 1.5,cex = 0.8, bty = "n")
}





get_AUC <- function(flag, Pvals, models){
  AUC.DE <- function(pval, DE.gs) {
    ii = !is.na(pval)
    library(ROCR)
    pred = prediction(1-pval[ii], DE.gs[ii])
    auc.tmp = performance(pred,"auc")
    auc = as.numeric(auc.tmp@y.values)
    auc
  }
  
  nmodel = length(Pvals)

  all.flag = as.vector(flag)
  auc = rep(NA, nmodel)
  for(i in 1:nmodel){
    this.pval = as.vector(Pvals[[i]])
    idx = is.na(this.pval)
    auc[i] = round(AUC.DE(this.pval, all.flag), 3)
  }
  names(auc) = models
  return(auc)
}


Plot_pattern <- function(res.sim, dat.info.sum, CT.number = 0,
                         rowInd.fixed = NULL){
  ### plot patterns in Beta, Mu, FC for specified celltype
  K = ncol(res.sim$FullParas$D_gk)
  N = length(unique(res.sim$Counts.wide$subject))
  Esti.FC = 2^dat.info.sum$effect_info
  sim.Beta = res.sim$FullParas$beta
  for (ict in 1:K) {
    if(ict == 1){
      mu_gik = res.sim$FullParas$mu[, , 1]
      FC = Esti.FC[, 1, ]
      Beta = sim.Beta[, , 1]
    }else{
      mu_gik = cbind(mu_gik, res.sim$FullParas$mu[, , ict])
      FC = cbind(FC, Esti.FC[, ict, ])
      Beta = cbind(Beta, sim.Beta[, , ict])
    }
  }
  
  colnames(mu_gik) =
    colnames(FC) =
    colnames(Beta) = paste0(rep(paste0("CT_", 1:K), each = N),
                            rep(paste0("@Subj_", 1:N), K) )
  
  
  idx.Allmarker = which(rowSums(res.sim$FullParas$D_gk) > 0)
  
  idx.Ind.marker = vector("list", length = K)
  for (ik in 1:K) {
    idx.Ind.marker[[ik]] = which((res.sim$FullParas$D_gk[, ik]) > 0)
  }
  
  if(CT.number == 0){
    idx = idx.Allmarker
  }else {
    idx = idx.Ind.marker[[CT.number]]
  }
  library(RColorBrewer)
  
  coul <- colorRampPalette(
    brewer.pal(8, "Spectral"))(10)
  ######
  sub.Beta        = data.matrix(Beta[idx, ])
  distance    = dist(sub.Beta)
  cluster     = hclust(distance, method="ward.D2")
  dendrogram  = as.dendrogram(cluster)
  Rowv        = rowMeans(sub.Beta, na.rm = T)
  dendrogram  = reorder(dendrogram, Rowv)
  reorderfun = function(d,w) { d }
  
  
  if(length(rowInd.fixed) == 0){
    
    # m_gk
    m_gk.beta = m_gk; m_gk.beta[sim.Beta[,1,] == 0] = 0
    rowInd = rev(order.dendrogram(dendrogram))
    sub.m_gk.beta = m_gk.beta[idx, ]
    sub.m_gk.beta_ordered <- sub.m_gk.beta[rowInd, ]
    heatmap(sub.m_gk.beta_ordered, #scale="none",
            Colv=NA, Rowv=NA,  col = coul,
            cexCol  = 0.8,
            #reorderfun=reorderfun,
            main = " ")
    title(main = "m_gk")
    
    ## Beta
    rowInd = rev(order.dendrogram(dendrogram))
    sub.Beta_ordered <- sub.Beta[rowInd, ]
    heatmap(sub.Beta_ordered, #scale="none",
            Colv=NA, Rowv=NA,  col = coul,
            #reorderfun=reorderfun,
            main = " ")
    title(main = "Beta")
    
    
    ## mu
    sub.mu_gik = mu_gik[idx, ]
    sub.mu_gik_ordered = sub.mu_gik[rowInd, ]
    heatmap(sub.mu_gik_ordered,  col = coul,
            Colv=NA,Rowv=NA,
            #reorderfun=reorderfun,
            main = " ")
    title(main = "mu")
    
    
    # # FC
    # sub.FC = FC[idx, ]
    # sub.FC_ordered = sub.FC[rowInd, ]
    # heatmap(sub.FC_ordered,  col = coul,
    #         Colv=NA,Rowv=NA,
    #         #reorderfun=reorderfun,
    #         main = " ")
    # title(main = "Fold Change")
    geneID = rownames(sub.Beta)[rowInd]
    return(geneID)
    
  }else{
    
    ## Beta
    sub.Beta_ordered <- data.matrix(Beta[rowInd.fixed, ])
    heatmap(sub.Beta_ordered, #scale="none",
            Colv=NA,Rowv=NA,  col = coul,
            # reorderfun=reorderfun, 
            main = " ")
    title(main = "Beta")
    
    
    ## mu
    #sub.mu_gik = mu_gik[idx, ]
    sub.mu_gik_ordered = mu_gik[rowInd.fixed, ]
    #  sub.mu_gik[rowInd, ]
    heatmap(sub.mu_gik_ordered,  col = coul,
            Colv=NA,Rowv=NA,
            # reorderfun=reorderfun,
            main = " ")
    title(main = "mu")
    
    
    # FC
    # #sub.FC = FC[idx, ]
    # sub.FC_ordered = FC[rowInd.fixed, ]
    #   #sub.FC[rowInd, ]
    # heatmap(sub.FC_ordered,  col = coul,
    #         Colv=NA,Rowv=NA,
    #         #reorderfun=reorderfun,
    #         main = " ")
    # title(main = "Fold Change")
    # 
  }
  
  
}




CP.Groups <- function(cellProp, AllPossComb){
  ### extract cell type proportion groups
  AllPossComb.CP = matrix(0, nrow = nrow(AllPossComb), ncol = ncol(AllPossComb))
  for (ic in 1:nrow(AllPossComb)) {
    CTs = AllPossComb[ic, ]
    AllPossComb.CP[ic, ] = round(cellProp[CTs]/sum(cellProp[CTs]),2)
  }
  colnames(AllPossComb.CP) = paste0("Celltype_", 1:4)
  rownames(AllPossComb.CP) = paste0("Comb_", 1:nrow(AllPossComb))
  
  distance    = dist(AllPossComb.CP)
  cluster     = hclust(distance, method="ward.D2")
  dendrogram  = as.dendrogram(cluster)
  Rowv        = rowMeans(AllPossComb.CP, na.rm = T)
  dendrogram  = reorder(dendrogram, Rowv)
  reorderfun = function(d,w) { d }
  rowInd = rev(order.dendrogram(dendrogram))
  AllPossComb.CP_ordered <- AllPossComb.CP[rowInd, ]
  
  heatmap(AllPossComb.CP_ordered, #scale="none",
          Colv=NA,Rowv=NA,  
          reorderfun=reorderfun, main = " ",
          cexCol = 1.1)
  return(AllPossComb.CP_ordered)
}




Check_Z <- function(res.sim){
  ## 1. check simulated data: rowmeans(Z_gik) proportional to D_gk*q_gk
  D_gk = res.sim$FullParas$D_gk
  q_gk = res.sim$FullParas$q_gk
  Bernouli.P = D_gk*q_gk
  par(mfrow = c(ncol(D_gk), 2), mai = c(0.4, 0.4, 0.2, 0.2))
  for (ct.ix in 1:ncol(D_gk)) {
    this.P = as.vector(Bernouli.P[, ct.ix]) 
    # N = 20
    plot(this.P, rowMeans(res.sim$FullParas$Z[,, ct.ix]), 
         cex = 0.5, col = "#00000050", xlab = " ", ylab = " ")
    abline(0,1, col = "red")
    title(main = paste0("N = ", N),  line = 0.3)
    title(xlab = "D_gk*q_gk", ylab = "Z",  line = 1.8)
    
    # N = 50
    this.P = as.vector(Bernouli.P[, ct.ix])
    tmp = matrix(rep(this.P, each = 50), nrow = length(this.P), byrow = TRUE)
    this.Z = matrix(rbinom(length(tmp), size = 1, prob = as.vector(tmp)), 
                    nrow = nrow(tmp), byrow = FALSE)
    plot(this.P, rowMeans(this.Z), 
         cex = 0.5, col = "#00000050",
         xlab = " ", ylab = " ")
    abline(0,1, col = "red")
    title( main = "N = 50", line = 0.3)
    title(xlab = "D_gk*q_gk", ylab = "Z",  line = 1.8)
  }
}


Eva_estimates <- function(res.sim, res.marker){
  ##### 
  K = length(res.marker)
  par(mfrow = c(2, K), mai = c(0.4, 0.5, 0.2, 0.2))
  MSE.m_gk = rep(NA, K)
  MSE.q_gk = rep(NA, K)
  for (ct.ix in 1:K) {
    # 1. evaluation of estimate for m_gk
    keep.idx = rownames(res.marker[[ct.ix]]$param) ### 
    marker.idx = paste0("gene_", which(res.sim$FullParas$D_gk[, ct.ix] == 1))
    com.id = intersect(keep.idx, marker.idx)
    
    ### effect size of genes with lower prevalence tend to be over-underestimated
    True_m = m_gk[com.id, ct.ix]; True_q = res.sim$FullParas$q_gk[com.id, ct.ix]
    Esti_m = res.marker[[ct.ix]]$param[com.id, "delta"]
    point.col = ifelse(Esti_m < 0.2, "red", "black")
    
    plot(True_m, Esti_m, cex = 0.6, pch = 16, xlab = " ", 
         ylab = " ", col = "#00000050", xlim = c(0, max( c(Esti_m, True_m))), 
         ylim = c(0, max( c(Esti_m, True_m))))
    points(True_m, Esti_m, col = point.col, pch = 16, cex =0.6)
    abline(0, 1, col = "green3")
    title(xlab = "True LFC", ylab = "Estimated LFC", line = 1.8)
    title(main = paste0("Celltype_", ct.ix), line = 0.3)
    r = round(cor(True_m, Esti_m), 3)
    legend("topleft", legend = paste0("r = ", r), bty = "n")
    
    #### MSE 
    MSE.m_gk[ct.ix] = mean((Esti_m - True_m)^2)
  }
  
  for (ct.ix in 1:K) {
    # 2. q_gk
    keep.idx = rownames(res.marker[[ct.ix]]$param) ### 
    marker.idx = paste0("gene_", which(res.sim$FullParas$D_gk[, ct.ix] == 1))
    com.id = intersect(keep.idx, marker.idx)
    
    sim.q = res.sim$FullParas$q_gk[com.id, ct.ix]

    plot(sim.q, #res.sim$FullParas$q_gk[com.id, ct.ix], 
         res.marker[[ct.ix]]$param[com.id, "q"],
         cex = 0.6, pch = 16, col = "#00000050", 
         xlim = c(0,1), ylim = c(0,1),
         xlab = " ", ylab = " ")
    r = round(cor(sim.q, res.marker[[ct.ix]]$param[com.id, "q"]),3)
    abline(0, 1, col = "red")
    title( xlab = "True prevalence", ylab = "Estimated prevalence", line = 1.8)
    legend("topleft", legend = paste0("r = ", r), bty = "n")
    
    MSE.q_gk[ct.ix] =  mean((res.marker[[ct.ix]]$param[com.id, "q"] - sim.q)^2)
  }
  
  return(res = list(MSE.m_gk = MSE.m_gk, MSE.q_gk = MSE.q_gk))
   
}



