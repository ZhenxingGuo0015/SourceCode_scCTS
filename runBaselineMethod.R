###

runBaselineMethod <- function(
  sce,
  use.norm.rep=NULL,
  subject.rep='subject',
  celltype.rep='celltype',
  per.subject=TRUE,
  method=c('wilcox', 'twelch', 'DEseq2', 'NSforest', 'FEAST', 'scGeneFit'),
  celltype.ngenes=NULL,
  numCores=NULL,
  python.path=NULL
  ){
  # check arguments
  stopifnot(all(c(subject.rep, celltype.rep) %in% names(colData(sce))),
            is(per.subject, "logical"))
  # set parallel computation
  numCores.used <- set.parallel.computation(numCores)
  
  # check python path
  if (method %in% baselines.python() & is.null(python.path)){
    cli_abort("{.var {method}} requires a python path.")
  }
  if (!(method %in% baselines.python()) & !is.null(python.path)){
    cli_alert_warning("Invalid argument python.path={.file {python.path}} for {.var {method}}.")
  }
  # get (normalized) count matrix
  if (method %in% baselines.counts()){
    Y <- get.expression.matrix(sce, use.raw=TRUE)
  }else{
    Y <- get.expression.matrix(sce, use.raw=FALSE, use.norm.rep=use.norm.rep)
  }
  # get cell types
  celltypes <- colData(sce)[[celltype.rep]] # convert to char list
  # check the list of numbers of selected markers for each cell type
  if (method %in% baselines.fixnumber()){
    if (is.null(celltype.ngenes)){
      cli_abort("{.var {method}} requires predefined numbers of markers for each cell type.")
    }
    isin.celltypes <- names(celltype.ngenes) %in% celltypes
    if (!all(isin.celltypes)){
      invalid.cts <- names(celltype.ngenes)[!isin.celltypes]
      cli_abort("Invalid number(s) of features to be selected: {.var {invalid.cts}}")
    }
    if (!is.integer(unlist(celltype.ngenes, use.names=FALSE))){
      cli_abort("Invalid number(s) of features to be selected: {.var {celltype.ngenes}}")
    }
  }
  if (!(method %in% baselines.fixnumber()) & !is.null(celltype.ngenes)){
    cli_alert_warning("Invalid argument celltype.ngenes={.var {celltype.ngenes}} for {.var {method}}.")
  }
  
  
  if (per.subject){
    cli_h1("Subject-level {.emph {method}} method")
    
    subjects <- colData(sce)[[subject.rep]] # convert to char list
    unique.subjects <- sort(unique(subjects))
    
    sub.res.list <- list()
    cli_progress_bar("Analyzing each subject", total = length(unique.subjects), type = "tasks")
    for (sub in unique.subjects){
      # cli_h2("Current subject: {sub}")
      sub.Y <- Y[,subjects == sub]
      sub.cts <- celltypes[subjects == sub]
      sub.result = switch(
        method,
        "wilcox" = BaselineMethod.wilcox(sub.Y, sub.cts, numCores.used),
        "twelch" = BaselineMethod.twelch(sub.Y, sub.cts, numCores.used),
        "DEseq2" = BaselineMethod.DEseq2(sub.Y, sub.cts, numCores.used),
        "NSforest" = BaselineMethod.NSforest(sub.Y, sub.cts, celltype.ngenes, python.path),
        "FEAST" = BaselineMethod.FEAST(sub.Y, sub.cts, celltype.ngenes, numCores.used),
        "scGeneFit" = BaselineMethod.scGeneFit(sub.Y, sub.cts, celltype.ngenes, python.path),
        stop(str_glue("No method matched for {method}"))
      )
      if (!(method %in% baselines.fixnumber())){
        # for DE tests, set the name of the last dim of each array as subject
        # name, and then store
        sub.res.list[[sub]] <- lapply(sub.result, function(arr){dimnames(arr)[3] <- sub;arr})
      }else{
        sub.res.list[[sub]] <- sub.result
      }
      cli_progress_update()
    }
    if (!(method %in% baselines.fixnumber())){
      # combine results for multiple subjects
      DE.res <- list()
      for (name in names(sub.res.list[[1]])){
        DE.res[[name]] <- abind(lapply(sub.res.list, function(lst) lst[[name]]))
      }
      return(DE.res)
    }else{
      return(sub.res.list)
    }
  }else{
    cli_h1("Population-level {.emph {method}} method")
    cli_text("{.emph Note: this mode needs batch-effects correction in advance.}")
    subjects <- colData(sce)[[subject.rep]]
    all.result = switch(
      method,
      "wilcox" = BaselineMethod.wilcox(Y, celltypes, numCores.used),
      "twelch" = BaselineMethod.twelch(Y, celltypes, numCores.used),
      "DEseq2" = BaselineMethod.DEseq2(Y, celltypes, numCores.used),
      "NSforest" = BaselineMethod.NSforest(Y, celltypes, celltype.ngenes, python.path),
      "FEAST" = BaselineMethod.FEAST(Y, celltypes, celltype.ngenes, numCores.used),
      "scGeneFit" = BaselineMethod.scGeneFit(Y, celltypes, celltype.ngenes, python.path),
      stop(str_glue("No method matched for {method}"))
    )
    if (!(method %in% baselines.fixnumber())){
      all.result <- lapply(all.result, function(arr){dimnames(arr)[3] <- "all";arr})
    }
    return(all.result)
  }
}


# returns methods implemented in python
baselines.python <- function(){
  return(c("NSforest", "scGeneFit"))
}

# returns methods that require a predefined number of markers
baselines.fixnumber <- function(){
  return(c("NSforest", "FEAST", "scGeneFit"))
}

# returns methods that require raw counts as inputs
baselines.counts <- function(){
  return(c("DEseq2", "FEAST"))
}



#' Naive Wilcoxon test
#'
#' @param expr A gene by cell matrix storing the expression values
#' @param celltypes A vector indicating cell types of each cell
#' @param nCores.used The number of cores actually used
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#'   is genes, the second dim is cell types, and the last dim is subjects.
#' @importFrom matrixTests row_wilcoxon_twosample
#' @importFrom foreach %dopar% foreach
#' @importFrom data.table rbindlist
#' @importFrom iterators iter
#' @importFrom parallel splitIndices
#' @importFrom stats p.adjust
#'
BaselineMethod.wilcox <- function(expr, celltypes, nCores.used){
  # both lognorm and norm data are acceptable for wilcoxon test
  unique.celltypes <- sort(unique(celltypes))
  
  # create 3-dim empty arrays to store results for a single subject (the last dim is 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(unique.celltypes), 1),
                     dimnames = list(rownames(expr), unique.celltypes))
  wilcox.res <- list('wilcox.stat_info' = array.tmp,
                     'wilcox.pval_info' = array.tmp,
                     'wilcox.fdr_info' = array.tmp)
  
  # compare each cell type with other cells
  for(ucelltype in unique.celltypes){
    # cli_text("Comparing {.val {ucelltype}} with other cells...")
    is.current.celltype <- (celltypes == ucelltype)
    not.current.celltype <- !is.current.celltype
    # split genes into chunks
    idx <- NULL  # to prevent "no visible binding for global variable" when checking
    chunks.res.list <- foreach(idx = iter(splitIndices(nrow(expr), nCores.used))) %dopar%{
      chunk.res <- row_wilcoxon_twosample(
        expr[idx,is.current.celltype], expr[idx,not.current.celltype], alternative = 'greater'
      )
      chunk.res[,c('statistic','pvalue')]
    }
    ct.res <- rbindlist(chunks.res.list)
    wilcox.res$wilcox.stat_info[,ucelltype,1] <- ct.res[['statistic']]
    wilcox.res$wilcox.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    wilcox.res$wilcox.fdr_info[,ucelltype,1] <- p.adjust(ct.res[['pvalue']], 'fdr')
  }
  return(wilcox.res)
}

#' t-test
#'
#' @inheritParams BaselineMethod.wilcox
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first dim
#'   is genes, the second dim is cell types, and the last dim is subjects.
#' @importFrom matrixTests row_t_welch
#' @importFrom foreach %dopar% foreach
#' @importFrom data.table rbindlist
#' @importFrom iterators iter
#' @importFrom parallel splitIndices
#' @importFrom stats p.adjust
#'
#'
BaselineMethod.twelch <- function(expr, celltypes, nCores.used){
  
  unique.celltypes <- sort(unique(celltypes))
  # create 3-dim empty arrays to store results for a single subject (i.e, the last dim is 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(unique.celltypes), 1),
                     dimnames = list(rownames(expr), unique.celltypes))
  twelch.res <- list('twelch.stat_info' = array.tmp,
                     'twelch.pval_info' = array.tmp,
                     'twelch.fdr_info' = array.tmp)
  
  # compare each cell type with other cells
  for(ucelltype in unique.celltypes){
    # cli_text("Comparing {.val {ucelltype}} with other cells...")
    is.current.celltype <- (celltypes == ucelltype)
    not.current.celltype <- !is.current.celltype
    # split genes into chunks
    idx <- NULL  # to prevent "no visible binding for global variable" when checking
    chunks.res.list <- foreach(idx = iter(splitIndices(nrow(expr), nCores.used))) %dopar%{
      chunk.res <- row_t_welch(
        expr[idx,is.current.celltype], expr[idx,not.current.celltype], alternative = 'greater'
      )
      chunk.res[,c('statistic','pvalue')]
    }
    ct.res <- rbindlist(chunks.res.list)
    twelch.res$twelch.stat_info[,ucelltype,1] <- ct.res[['statistic']]
    twelch.res$twelch.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    twelch.res$twelch.fdr_info[,ucelltype,1] <- p.adjust(ct.res[['pvalue']], 'fdr')
  }
  return(twelch.res)
}



#' ZINB-WaVE + DESeq2
#'
#' @inheritParams BaselineMethod.wilcox
#'
#' @return A list of 3-dim arrays. Each array corresponds to a stat. The first
#'   dim is genes, the second dim is cell types, and the last dim is a single
#'   subject.
#' @importFrom BiocParallel MulticoreParam
#' @importFrom S4Vectors SimpleList
#' @import SingleCellExperiment
#'
BaselineMethod.DEseq2 <- function(expr, celltypes, nCores.used){
  for (pkg in c("zinbwave", "DESeq2")){
    if(!requireNamespace(pkg)){
      cli_abort("This function requires the {.pkg {pkg}} package.")
    }
  }
  BPPARAM <- MulticoreParam(nCores.used)
  
  ucelltypes <- unique(celltypes)
  cleaned.celltypes <- make.names(celltypes)
  cleaned.ucelltypes <- make.names(ucelltypes)
  # create 3-dim empty arrays to store results for a single subject (i.e, the last dim is 1)
  array.tmp <- array(NA, dim=c(length(rownames(expr)), length(ucelltypes), 1),
                     dimnames = list(rownames(expr), ucelltypes))
  DESeq2.res <- list('DESeq2.stat_info' = array.tmp,
                     'DESeq2.pval_info' = array.tmp,
                     'DESeq2.fdr_info' = array.tmp,
                     'DESeq2.log2FC_info' = array.tmp)
  
  # build a SingleCellExperiment object
  # core <- SingleCellExperiment(assays=list(counts=expr),
  #                              colData=data.frame(celltype = factor(cleaned.celltypes))
  #                              )
  
  ### SingleCellExperiment() always reports error, change it to SummarizedExperiment()
  core <- SummarizedExperiment(assays=SimpleList(counts=expr),
                               colData=DataFrame(celltype = factor(cleaned.celltypes)))
  
  # ZINB-WaVE, specify K = 0 to only compute observational weights
  zinb <- zinbwave::zinbwave(core, 
                             K=0, observationalWeights=TRUE, BPPARAM=BPPARAM, epsilon=1e12)
  mode(assay(zinb)) <- "integer"  # to prevent "converting counts to integer mode"
  # DESeq2
  for (cleaned.uct in cleaned.ucelltypes){
    # build two groups
    is.current.celltype <- (cleaned.celltypes == cleaned.uct)
    colData(zinb)['group'] <- as.factor(ifelse(is.current.celltype, cleaned.uct, "others"))
    
    dds <- DESeq2::DESeqDataSet(zinb, design = ~ group)
    dds <- DESeq2::DESeq(
      dds, sfType="poscounts", useT=TRUE, minmu=1e-6, minRep=Inf, fitType='local',
      parallel=T, BPPARAM=BPPARAM, quiet=TRUE
    )
    ct.res <- DESeq2::results(object = dds, contrast = c("group", cleaned.uct, "others"),
                              alpha = 0.05, pAdjustMethod ='fdr', altHypothesis="greater",
                              parallel=T, BPPARAM=BPPARAM)
    # store results
    ucelltype <- ucelltypes[cleaned.ucelltypes == cleaned.uct]
    DESeq2.res$DESeq2.stat_info[,ucelltype,1] <- ct.res[['stat']]
    DESeq2.res$DESeq2.pval_info[,ucelltype,1] <- ct.res[['pvalue']]
    DESeq2.res$DESeq2.fdr_info[,ucelltype,1] <- ct.res[['padj']]
    DESeq2.res$DESeq2.log2FC_info[,ucelltype,1] <- ct.res[['log2FoldChange']]
  }
  return(DESeq2.res)
}





#' NS-Forest
#'
#' Running this function will create a folder ./NSForest_outputs/ in the current
#' working directory. For each cell type, the genes should be sorted first by
#' binary scores and then by feature importances from the random forest model.
#' The number of actually selected features may be less than the given number of
#' selected features, since \emph{negative markers} defined in the
#' \href{https://doi.org/10.1101/gr.275569.121}{paper} are filtered.
#'
#' @inheritParams BaselineMethod.wilcox
#' @param celltype.ngenes A named list. The names are cell types, and the values
#'   are number of features selected for that cell type.
#' @param python.path The path to the python.
#'
#' @return A named list. Names are unique cell types. Values are selected
#'   features for that cell type.
#'
#' @import dplyr
#' @importFrom readr read_csv
#' @importFrom purrr map2
#'
BaselineMethod.NSforest <- function(expr, celltypes, celltype.ngenes, python.path){
  # import python packages
  pkg <- "reticulate"
  if(!requireNamespace(pkg)){
    cli_abort("This function requires the {.pkg {pkg}} package.")
  }
  reticulate::use_python(python.path)
  nsforest <- reticulate::import("nsforest")
  ad <- reticulate::import("anndata")
  
  ucelltypes <- unique(celltypes)
  X <- reticulate::r_to_py(t(as.data.frame(expr)))
  obs <- reticulate::r_to_py(data.frame(celltype=celltypes, row.names=colnames(expr)))
  # create the anndata object
  adata <- ad$AnnData(X=X, obs=obs)
  adata$var_names <- reticulate::r_to_py(rownames(expr))
  # run NSForest for all cell types
  # the detailed results are saved in ./NSForest_outputs/NSForest_supplementary.csv
  reticulate::py_capture_output(nsforest$NSForest(
    adata, cluster_header='celltype', n_top_genes=adata$n_vars, n_binary_genes=adata$n_vars
  ))
  supp <- read_csv("./NSForest_outputs/NSForest_supplementary.csv", show_col_types=FALSE)
  
  process_cluster <- function(cluster_name, size, data) {
    filtered_data <- data %>%
      filter(get("clusterName") == cluster_name) %>%
      arrange(desc(get("binary_score")), desc(get("rf_feature_importance"))) %>%
      slice_head(n = size)  # use get() to prevent the no visible binding issue
    
    return(filtered_data$binary_genes)
  }
  
  NSforest.res <- map2(.x = names(celltype.ngenes),
                       .y = celltype.ngenes,
                       .f = function(name, size) process_cluster(name, size, supp))
  names(NSforest.res) <- names(celltype.ngenes)
  return(NSforest.res)
}


#' FEAST
#'
#' @inheritParams BaselineMethod.wilcox
#' @param celltype.ngenes A named list. Names are cell types, and the values
#'   are number of features selected for that cell type.
#'
#' @return A list of each cell type's highly variable genes. The number of HVGs
#'   are equal to the given number.
#' @importFrom purrr quietly
#'
BaselineMethod.FEAST <- function(expr, celltypes, celltype.ngenes, nCores.used){
  library(FEAST)
  library(cli)
  library(clintools)
  pkg <- "FEAST"
  if(!requireNamespace(pkg)){
    cli_abort("This function requires the {.pkg {pkg}} package.")
  }
  unique.celltypes <- sort(unique(celltypes))
  FEAST.res <- list()
  for(ucelltype in unique.celltypes){
    ngenes <- celltype.ngenes[[ucelltype]]
    tryCatch({
      idxs <- quietly(FEAST::FEAST_fast)(
        expr[,celltypes == ucelltype], nProc=nCores.used
      )$result  # quietly returns a list of all messages and original return values
      FEAST.res[[ucelltype]] <- rownames(expr)[idxs[1:ngenes]]
    }, error=function(e){print(e);cli_alert_warning("Error occured when processing {.var {ucelltype}}, continue...")})
  }
  return(FEAST.res)
}


#' scGeneFit
#'
#' @inheritParams BaselineMethod.NSforest
#' @return A named list. Names are unique cell types. Values are selected
#'   features for that cell type.
#'
BaselineMethod.scGeneFit <- function(expr, celltypes, celltype.ngenes, python.path){
  # import python packages
  pkg <- "reticulate"
  if(!requireNamespace(pkg)){
    cli_abort("This function requires the {.pkg {pkg}} package.")
  }
  reticulate::use_python(python.path)
  scGeneFit <- reticulate::import("scGeneFit")
  
  unique.celltypes <- sort(unique(celltypes))
  scGeneFit.res <- list()
  for(ucelltype in unique.celltypes){
    ngenes <- celltype.ngenes[[ucelltype]]
    X <- reticulate::r_to_py(t(expr))
    cell.labels <- reticulate::r_to_py(ifelse(celltypes == ucelltype, ucelltype, "others"))
    # + 1 since in python indices start from 0
    reticulate::py_capture_output(
      idxs <- scGeneFit$functions$get_markers(X, cell.labels, num_markers=ngenes) + 1
    )
    scGeneFit.res[[ucelltype]] <- rownames(expr)[idxs]
  }
  return(scGeneFit.res)
}
