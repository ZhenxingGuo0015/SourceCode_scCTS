library(matrixTests)
library(MASS)
library(MASSExtra)
library(Matrix)
library(doParallel)

### naive wilcoxon test
### the key function row_wilcoxon_twosample comes from package matrixTests
### INPUT: Y, gene expression matrix (Gene * sample)
### INPUT: CT, index of target cell type (TRUE: cell is from target cell type)
### PURPOSE: One side test (CT > nonCT)
### OUTPUT: Wilcoxon test statistic and p-values

wilcoxon_Test <- function(Y, CT) {
  nonCT = !CT ## I put this outside the loop to make it a little faster
  wilcoxon.res <- row_wilcoxon_twosample(x = as.matrix(Y[,CT]) ,
                                         y = as.matrix(Y[,nonCT]), 
                                         alternative = 'greater')
  res.output <- wilcoxon.res[,c('statistic','pvalue')]
  return(res.output)
}

### MarkerFinder
### INPUT: Y, gene expression matrix (normalized to 10k)
### INPUT: celltype: cell type info of cells 
### PURPOSE: calculate effect size and corresponding variance
### OUTPUT: All necessary sample level summary statistics 

### This function is not clean. Some of the calculated statistics are not used
### in current manuscript. 
### Such statistics will be marked with "IGNORE!" ####

markerFinder <- function(Y, celltype){
  #### obtain marker gene within each subject
  gene.num <- nrow(Y)  # gene number
  uniqueCT <- sort(unique(celltype)) # names of cell types
  cell.num <- length(uniqueCT) # number of cell types, NOT number of cells
  feature.names <- rownames(Y) # gene names
  
  #################################
  res <- list() # store result
  matrix.tmp <- matrix(NA, nrow = gene.num, ncol = cell.num) # create a matrix 
  colnames(matrix.tmp) <- uniqueCT # column is cell type
  rownames(matrix.tmp) <- feature.names # row is genes
  
  # create null variables to store different statistics
  res[['expr']] = res[['expr_var']] = res[['expr_se']] = 
    res[['log_expr']] = res[['log_expr_se']] = 
    res[['effect']] = res[['effect_se']] =  
    res[['wilcox.stat']] = res[['wilcox.pval']] = matrix.tmp
  
  # do the loop for each cell type
  for(ct.ix in uniqueCT){ 
    cat(ct.ix, sep = "\n")
    cell.ix <- (celltype == ct.ix) # index for target cell type
    expr.ct.tmp <- Y[,cell.ix]  # gene expression of cells from target cell type
    expr.remain.tmp <- Y[,!cell.ix]  # gene expression of cells from non-targer cell types
    
    res[['expr']][,ct.ix] <- rowMeans(expr.ct.tmp)  # mean expression of target cell type
    res[['expr_var']][,ct.ix] <- apply(expr.ct.tmp,1,var) # variance of expression for target cell type
    res[['expr_se']][,ct.ix] <- sqrt( res[['expr_var']][,ct.ix] / sum(cell.ix)) # sample error for target cell type
    
    res[['log_expr']][,ct.ix] <- log2(res[['expr']][,ct.ix] + 1) # log2 transformed gene expression + psudo count 1
    res[['log_expr_se']][,ct.ix] <- 
      sqrt(res[['expr_var']][,ct.ix]/(sum(cell.ix)*((res[['expr']][,ct.ix]+1)*log(2))^2)) 
    
    test.res <- wilcoxon_Test(Y=Y, CT=cell.ix) # wilcoxon rank sum test one vs. others (a common used/traditional way)
    
    res[['wilcox.stat']][,ct.ix] = test.res[,'statistic']
    res[['wilcox.pval']][,ct.ix] = test.res[,'pvalue']
    
    
  }
  
  ### Loop again, for each cell type calculate effect size (a.k.a log2 fold change defined in Eq 2 in manuscript) 
  ### and its corresponding variance defined in Eq 6 in Supplementary file 
  ### (looks like Eq6 is not correct, it should be summation + not minus - between the two terms) ####
  #par(mfrow=c(2,2))
  for(ct.ix in uniqueCT){
    res[['effect']][,ct.ix] <- res[['log_expr']][, ct.ix] - 
      log2(rowMeans(res[['expr']][,! uniqueCT %in% ct.ix], na.rm = T) +1) 
    
    se.part1 <- res[['log_expr_se']][,ct.ix]^2 
    # se.part2.luxiao <- 1/(log(2)*(rowSums(res[['expr']][,! uniqueCT %in% ct.ix],na.rm = T) +
    #                          (cell.num-1))
    #                       )^2 * rowSums(res[['expr_se']][,! uniqueCT %in% ct.ix]^2)   ##### By Luxiao, questionable

    se.part2 <- 1/(log(2)*(rowSums(res[['expr']][,! uniqueCT %in% ct.ix],na.rm = T) + 0.001)
                   )^2 * rowSums(res[['expr_se']][,! uniqueCT %in% ct.ix]^2)   ##### Modified by Zhenxing


    #### se.part2.luxiao and se.part2 has large difference
    #par(mfrow=c(1,2))
    # plot(se.part2.luxiao, se.part2, xlim = c(0, 0.05), ylim = c(0,0.05),
    #      cex = 0.5, pch = 16, main = ct.ix)
    # abline(0,1, col = "red")
    # 
    # plot(sqrt(se.part1 + se.part2.luxiao),
    #      sqrt(se.part1 + se.part2),
    #      cex = 0.5, pch = 16, xlim = c(0, 0.2), ylim = c(0, 0.2))
    # abline(0,1, col = "red")

    res[['effect_se']][,ct.ix] <- sqrt(se.part1 + se.part2) ### I took a square here, in following analysis we need variance 
    
  }
  

  return(res)
}

### loop_samples
### INPUT: dat.exp, gene expression
### INPUT: subject_info, index of each sample/subject
### INPUT: celltype_info, cell type of each cells
### INPUT: numCores, core number for parallel computation
### INPUT: logcpm.inpu, whether the normalized counts is log
### PURPOSE: estimate variance and effect size for each celltype, each subject
### OUTPUT: a list, summary statistic for each sample/subject
loop_samples <- function(dat.exp, subject_info, celltype_info, 
                         numCores=1, logcpm.input=FALSE){
  
  registerDoParallel(numCores) 
  
  subject_num <- length(unique(subject_info))
  subjects <- sort(unique(subject_info))
  
  res <- list()
  res.store <- foreach (subject.ix = 1:length(subjects) ) %dopar%  {
    
    cat('\n', subjects[subject.ix],':\n')
    target_samples <- which(subject_info==subjects[subject.ix])
    dat.exp_tmp <- dat.exp[,target_samples]
    if(logcpm.input){
      dat.exp_tmp <- exp(dat.exp_tmp) - 1
    }
    celltype_tmp <- celltype_info[target_samples]
    
    res.tmp <- markerFinder(Y = dat.exp_tmp, celltype = celltype_tmp)
    res.tmp
  }
  res <- res.store
  names(res) <- as.character(subjects)
  return(res)
}


### re-construct the data generated by loop_sample
### make the list to an array
marker_summarization <- function(marker.info){
  nGenes <- nrow(marker.info[[1]]$expr)
  nCT <- ncol(marker.info[[1]]$expr)
  nInds <- length(marker.info)
  
  array.tmp <- array(NA, dim=c(nGenes, nCT, nInds))
  dimnames(array.tmp) <- list(rownames(marker.info[[1]]$expr),
                              colnames(marker.info[[1]]$expr),
                              names(marker.info))
  
  expr_info = expr_var_info = 
    log_expr_info = log_expr_se_info =
    effect_info = effect_se_info = 
    wilcox.stat_info = wilcox.pval_info = wilcox.fdr_info = array.tmp

  subjects.tmp <- names(marker.info)
  for(sub.ix in subjects.tmp){
    
    expr_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$expr)
    
    expr_var_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$expr_var)
    
    log_expr_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$log_expr)
    
    log_expr_se_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$log_expr_se)
    
    effect_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$effect)
    
    effect_se_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$effect_se)
    
    wilcox.stat_info[,,sub.ix] <- 
      as.matrix(marker.info[[sub.ix]]$wilcox.stat)
    
    wilcox.pval_info[,,sub.ix] <- as.matrix(marker.info[[sub.ix]]$wilcox.pval)
    
    for(cell.ix in 1:nCT){
      wilcox.fdr_info[,cell.ix,sub.ix] <- 
        p.adjust(c(wilcox.pval_info[,cell.ix,sub.ix]), 'fdr')
    }
  } 
  
  res <- list("expr_info"=expr_info,  "expr_var_info" = expr_var_info, 
              "log_expr_info" = log_expr_info, "log_expr_se_info" = log_expr_se_info,
              "effect_info" = effect_info,  'effect_se_info' = effect_se_info, 
              "wilcox.stat_info" = wilcox.stat_info, 'wilcox.pval_info' = wilcox.pval_info,
              'wilcox.fdr_info' = wilcox.fdr_info
             )
  return(res)
  
}











