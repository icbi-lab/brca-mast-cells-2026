
### Data access and preprocessing functions

load_metabric_data = function(dir){
  mb <- list()
  # expression data: log2 intensity levels (Illumina HT-12 v3 microarray)
  mb$expr <- read.delim(file.path(dir, "brca_metabric", "data_mrna_illumina_microarray.txt"), sep = "\t")
  colnames(mb$expr) <- sub("\\.", "", colnames(mb$expr))
  # expression data: log2 z-scores relative to all samples (Illumina HT-12 v3 microarray)
  mb$expr_zscores <- read.delim(file.path(dir, "brca_metabric", "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt"), sep = "\t")
  colnames(mb$expr_zscores) <- sub("\\.", "", colnames(mb$expr_zscores))
  # somatic mutation data from targeted sequencing
  mb$mutations <- read.delim(file.path(dir, "brca_metabric", "data_mutations.txt"), sep = "\t", skip = 1)
  mb$mutations$Tumor_Sample_Barcode <- sub("-", "", mb$mutations$Tumor_Sample_Barcode)
  # copy-number alterations (-2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification)
  mb$cna <- read.delim(file.path(dir, "brca_metabric", "data_cna.txt"), sep = "\t")
  colnames(mb$cna) <- sub("\\.", "", colnames(mb$cna))
  mb
}

load_metabric_metadata = function(dir){
  # clinical data
  data_patient <- read.delim(file.path(dir, "brca_metabric", "data_clinical_patient.txt"), sep = "\t", skip = 4)
  data_sample <- read.delim(file.path(dir, "brca_metabric", "data_clinical_sample.txt"), sep = "\t", skip = 4)
  metadata <- dplyr::full_join(data_patient, data_sample)
  rownames(metadata) <- sub("-", "", metadata$PATIENT_ID)
  # follow-up data
  rueda_st5 <- read.delim(file.path(dir, "rueda_clinical", "Rueda_ST5.txt"))
  rownames(rueda_st5) <- sub("-", "", rueda_st5$METABRIC.ID)
  colnames(rueda_st5) <- gsub(".", "_", toupper(colnames(rueda_st5)), fixed = TRUE)
  metadata <- data.frame(metadata, rueda_st5[rownames(metadata), setdiff(colnames(rueda_st5), colnames(metadata))])
  metadata
}

merge_rows <- function(data, ids, FUN = sum){
  stopifnot(nrow(data) == length(ids))
  ix <- ids %in% ids[duplicated(ids)]
  if (any(ix)){
    data_unique <- data[!ix,,drop=FALSE]
    data_multi <- data[ix,,drop=FALSE] |> as.data.frame()
    data_multi$ID <- ids[ix]
    data_multi_merged <- data_multi |> group_by(ID) |> dplyr::summarise(across(everything(), FUN), .groups = "drop")
    ids_merged <- data_multi_merged$ID
    data_multi_merged$ID <- NULL
    data <- rbind(data_unique, as(data_multi_merged, class(data_unique)[1]))
  }
  rownames(data) <- c(ids[!ix], ids_merged)
  if (!is.null(rownames(data))){
    data <- data[unique(ids),]
  }
  return(data)
}

se_merge_rows <- function(se, ids, assays = assays){
  merged_assays <- list()
  for (a in names(assays)){
    merged_assays[[a]] <- merge_rows(SummarizedExperiment::assay(se, a), ids = ids, FUN = assays[[a]])
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = merged_assays,
    colData = SummarizedExperiment::colData(se),
  )
}

getMSigDB <- function(collections = c("H", "C2|CP:KEGG"), species = "human", db_species = "HS", id_type = "gene_symbol", format = c("list", "dataframe", "GeneSetCollection"), ...){
  
  stopifnot(requireNamespace("msigdbr"))
  # msigdbr::msigdbr_species()
  # msigdbr::msigdbr_collections()
  
  if (!is.null(collections)){
    
    collections <- strsplit(collections, split = "|", fixed = TRUE)
    categories <- sapply(collections, function(x) x[1] )
    subcategories <- lapply(setNames(collections, categories), function(x) unique(x[-1]) )
    subcategories[sapply(subcategories, length) == 0] <- NA
    df <- stack(subcategories)[,c(2,1)]
    colnames(df) <- c("category", "subcategory")
    
    
    gslist <- lapply(1:nrow(df), function(i){
      categ <- as.character(df[i,1])
      subcateg <- df[i,2]
      if (is.na(subcateg)) subcateg <- NULL
      msigdbr::msigdbr(species = species, collection = categ, db_species = db_species, subcollection = subcateg)
    })
    
    genesets <- Reduce(x = gslist, f = rbind)
    
  } else {
    
    genesets <- msigdbr::msigdbr(species = species)
  }
  
  
  df <- unique(data.frame(genesets[,c(id_type, "gs_name")]))
  colnames(df) <- c("gene", "term")
  GS <- convertGeneSets(df, from = "dataframe", to = format, ...)
  
  GS
}

### Analysis functions

get_association <- function(data, test_var, ref_var, co_vars = NULL, strat_vars = NULL, scale=FALSE, verbose=TRUE, ...){
  
  # get variables
  test_type <- get_var_type(data, test_var)
  ref_type <- get_var_type(data, ref_var)
  
  results <- list(fit=NA, res=NA)
  
  # prepare data
  if (ref_type$type == "survival"){
    testdf <- dplyr::select(data, all_of(unique(setdiff(c(ref_var, test_var, co_vars), ""))))
    testdf <- testdf[rowSums(is.na(testdf)) == 0,,drop=FALSE]
  } else {
    testdf <- dplyr::select(data, all_of(unique(setdiff(c(test_var, ref_var, co_vars), "")))) |> tidyr::drop_na()
  }
  if (nrow(testdf) < 4) return(results)
  if (length(unique(testdf[[test_var]])) < 2) return(results)
  if (length(unique(testdf[[ref_var]])) < 2) return(results)
  if (scale == TRUE){
    if (test_type$type == "numeric"){
      testdf[[test_var]] <- as.numeric(scale(testdf[[test_var]]))
    }
    if (ref_type$type == "numeric"){
      testdf[[ref_var]] <- as.numeric(scale(testdf[[ref_var]]))
    }
  }
  
  # run tests
  if (ref_type$type == "survival"){
    if (verbose) print("Survival: survival ~ test_var + co_vars")
    results <- test_surv(testdf, test_var, ref_var, co_vars, strat_vars)
  } 
  
  if (ref_type$type %in% c("ordinal", "nominal")){
    if (test_type$type == "numeric"){
      if (verbose) print("LM: test_var ~ ref_var + co_vars")
      results <- test_lm(testdf, test_var, ref_var, co_vars)
    } else {
      if (verbose) print("MLR: test_var ~ ref_var + co_vars")
      results <- test_multinom(testdf, test_var, ref_var, co_vars)
    }
  }
  
  if (ref_type$type == "numeric"){
    if (test_type$type == "numeric"){
      if (verbose) print("LM: test_var ~ ref_var + co_vars")
      results <- test_lm(testdf, test_var, ref_var, co_vars)
    } else {
      if (verbose) print("LM: ref_var ~ test_var + co_vars")
      results <- test_lm(testdf, ref_var, test_var, co_vars)
    }
  }
  
  results
}


test_surv <- function(df, test_var, ref_var, co_vars = NULL, strat_vars = NULL, ...){
  strat_vars <- strat_vars[strat_vars %in% co_vars]
  if (length(strat_vars) == 0){
    f <- as.formula(paste0(ref_var, " ~ ", paste(c(test_var, co_vars), collapse = " + ")))
  } else {
    f <- as.formula(paste0(ref_var, " ~ ", paste(c(test_var, setdiff(co_vars, strat_vars), paste0("strata(", strat_vars, ")")), collapse = " + ")))
  }
  print(f)
  fit = survival::coxph(f, data = df)
  tidy_fit <- broom::tidy(fit) |> as.data.frame()
  res <- tidy_fit %>%
    filter(grepl(paste0("^", test_var), term)) %>%
    select(term, estimate, std.error, statistic, p.value)
  res$estimate <- exp(res$estimate)
  results <- list(fit=fit, res=res)
  results
}

test_lm <- function(df, test_var, ref_var, co_vars = NULL, ...){
  f <- as.formula(paste0(test_var, " ~ ", paste(c(ref_var, co_vars), collapse = " + ")))
  fit <- lm(f, data = df)
  tidy_fit <- broom::tidy(fit) |> as.data.frame()
  res <- tidy_fit %>%
    filter(grepl(paste0("^", ref_var), term)) %>%
    select(term, estimate, std.error, statistic, p.value)
  results <- list(fit=fit, res=res)
  results
}

test_multinom <- function(df, test_var, ref_var, co_vars = NULL, ...){
  f <- as.formula(paste0(test_var, " ~ ", paste(c(ref_var, co_vars), collapse = " + ")))
  fit <- nnet::multinom(f, data = df, trace = FALSE)
  tidy_fit <- broom::tidy(fit) |> as.data.frame()
  res <- tidy_fit %>%
    filter(grepl(paste0("^", ref_var), term)) %>%
    select(term, estimate, std.error, statistic, p.value)
  results <- list(fit=fit, res=res)
  results
}

get_n_unique <- function(data, vars){
  sapply(setNames(vars, vars), function(var){
    length(unique(data[[var]]))
  })
}

stratify <- function (x, n = 2, ties = "lower", optimal = FALSE, surv = NULL, ...){
  if (optimal == FALSE) {
    if ("data.frame" %in% class(x)) {
      ix <- sapply(x, is.numeric)
      ix[colnames(x) %in% c("surv")] <- FALSE
      x_strat <- as.data.frame(lapply(x[, ix, drop = FALSE], function(y) {stratify_var(x = y, n = n, ties = ties)}))
      x[, ix] <- x_strat
    }
    else {
      x <- stratify_var(x = x, n = n, ties = ties, ...)
    }
  }
  else {
    if (is.null(dim(x))) 
      x <- data.frame(x = x)
    if (is.null(surv)) 
      surv <- x$surv
    stopifnot(!is.null(surv))
    ix <- sapply(x, is.numeric)
    ix[colnames(x) %in% c("surv")] <- FALSE
    tmp <- data.frame(x[, ix, drop = FALSE], time = surv[,1], status = surv[, 2])
    cutpoints <- survminer::surv_cutpoint(tmp, time = "time", event = "status", variables = setdiff(colnames(tmp), c("time", "status")), progressbar = FALSE, ...)
    cutdata <- select(data.frame(survminer::surv_categorize(x = cutpoints)), -c(time, status))
    tmp <- as.data.frame(lapply(setNames(colnames(cutdata), colnames(cutdata)), function(var) {
      val <- signif(cutpoints$cutpoint[var, "cutpoint"], 2)
      stratvar <- ifelse(cutdata[[var]] == "low", paste0("<", val), paste0(">", val))
      relevel(factor(stratvar), ref = paste0("<", val))
      }))
    x[, ix] <- tmp
    if (ncol(x) == 1) 
      x <- x[[1]]
  }
  x
}

stratify_var <- function (x, n = 2, ties = "lower", ...){
  if (n == 0) 
    return(rep(mean(x), n = length(x)))
  if (n == 1) 
    return(x)
  qs <- quantile(x, seq_along(1:(n - 1))/n, na.rm = TRUE)
  if (n == 2) {
    classes.names <- c("Lo", "Hi")
  }
  else if (n == 3) {
    classes.names <- c("Lo", "Mid", "Hi")
  }
  else {
    names(qs) <- paste0(as.character(round(100 * seq_along(1:(n - 1))/n, 1)), "%")
    tmpnames <- c()
    if (n > 2) {
      for (i in 1:(length(qs) - 1)) {
        tmpnames[i] <- paste0(names(qs)[i], "-", names(qs)[i + 1])
      }
    }
    classes.names <- c(paste0("<", names(qs)[1]), tmpnames, paste0(">", names(qs)[length(qs)]))
  }
  classes <- rep(classes.names[1], length(x))
  if (ties == "higher") {
    for (i in 1:length(qs)) {
      classes[x >= qs[i]] <- classes.names[i + 1]
    }
  }
  else if (ties == "lower") {
    for (i in 1:length(qs)) {
      classes[x > qs[i]] <- classes.names[i + 1]
    }
  }
  if (n == 2) {
    classes <- relevel(factor(classes, levels = classes.names), ref = "Lo")
  }
  else {
    classes <- factor(classes, levels = classes.names, ordered = TRUE)
  }
  names(classes) <- names(x)
  classes[is.na(x)] <- NA
  return(classes)
}

surv_cut <- function(df, time_max = Inf){
  ix <- df$time > time_max
  ix[is.na(ix)] <- FALSE
  df[ix, "status"] <- 0
  df[ix, "time"] <- time_max
  df
}

deconvolute_cibersort <- function(gene_expression_matrix, arrays, absolute = FALSE, abs_method = "sig.score", ...){
  quantile_norm <- arrays
  testit::assert("CIBERSORT.R is provided", exists("cibersort_binary", envir = immunedeconv:::config_env))
  testit::assert("CIBERSORT signature matrix is provided", exists("cibersort_mat", envir = immunedeconv:::config_env))
  source(get("cibersort_binary", envir = immunedeconv:::config_env))
  tmp_mat <- tempfile()
  readr::write_tsv(as_tibble(gene_expression_matrix, rownames = "gene_symbol"), file = tmp_mat)
  arguments <- rlang::dots_list(get("cibersort_mat", envir = immunedeconv:::config_env), 
                                tmp_mat, perm = 0, QN = quantile_norm, absolute = absolute, 
                                abs_method = abs_method, ..., .homonyms = "last")
  call <- rlang::call2(CIBERSORT, !!!arguments)
  res <- eval(call)
  res <- res %>% t()
  return(res)
}

findMarkers <- function(seu, logfc.thresh = 0.5, min.pct = 0.1, test.use = "wilcox", group.by = NULL, group.oi = NULL, group.ref = NULL, assay = NULL, slot = "data", return.thresh = 0.05, fdr = FALSE, n_cores = 1, ...){
  
  if (!is.null(group.by)){
    Idents(seu) <- seu@meta.data[[group.by]]
  }
  
  if (fdr == TRUE){
    logfc.thresh = 0
    min.pct = 0
  }
  
  if (n_cores > 1){
    future::plan("multisession", workers = n_cores)
    on.exit(future::plan("sequential"))
  }
  
  if (is.null(group.oi) & is.null(group.ref)){
    de_markers <- FindAllMarkers(object = seu, assay = assay, logfc.threshold = logfc.thresh, test.use = test.use, slot = slot, min.pct = min.pct, min.diff.pct = -Inf, ...)
    
  } else {
    de_markers <- FindMarkers(object = seu, ident.1 = group.oi, ident.2 = group.ref,
                              assay = assay, logfc.threshold = logfc.thresh, test.use = test.use, slot = slot, min.pct = min.pct, min.diff.pct = -Inf, ...)
    
    de_markers$gene <- rownames(de_markers)
    
    if (is.null(group.ref)){
      de_markers$cluster <- group.oi
    } else {
      de_markers$cluster <- paste0(group.oi, ".vs.", group.ref)
    }
  }

  de_markers$pct.diff <- de_markers$pct.1 - de_markers$pct.2
  de_markers$cluster <- as.character(de_markers$cluster)
  
  if (fdr == TRUE){
    de_markers$fdr <- p.adjust(de_markers$p_val, method = "fdr")
  }
  
  if (all(c("p_val", "avg_log2FC", "pct.diff") %in% colnames(de_markers))){
    de_markers <- de_markers |> dplyr::arrange(p_val, desc(abs(avg_log2FC)), desc(pct.diff))
  }
  de_markers
}

runLIMMA <- function (data, design, formula = ~1, contrasts = NULL, trend = TRUE, 
                      robust = FALSE, p.adj.method = "fdr", do_log = FALSE, normalize = FALSE, shrink=FALSE,
                      norm.method = "vsn", ...){
  stopifnot(requireNamespace("limma", quietly = TRUE))
  if (ncol(data) != nrow(design)) {
    stop("Error: The number of columns in data must match the number of rows in the design data frame.")
  }
  if (!is.null(colnames(data)) & !is.null(rownames(design))) {
    if (!all(colnames(data) == rownames(design))) {
      warning("Warning: Mismatch in 'colnames(data)' and 'rownames(design)'!")
    }
  }

  if (do_log == TRUE) {
    data <- log2(data + 1)
  }
  if (normalize == TRUE & tolower(norm.method) != "vsn") {
    if (do_log == TRUE) {
      stop("Error: Do not combine log-transformation and vsn!")
    }
    data <- limma::normalizeVSN(data)
  }
  else if (normalize == TRUE) {
    data <- limma::normalizeBetweenArrays(data, method = norm.method)
  }
  formula <- update(formula, ~0 + .)
  mm <- model.matrix(formula, design)
  if (qr(mm)$rank < ncol(mm)) {
    warning("Model matrix is rank-deficient.")
  }
  contrasts_named <- lapply(contrasts, function(contr) {
    if (length(contr) == 1) {
      contr
    }
    else {
      paste(c(paste0(contr[1], contr[2]), paste0(contr[1], contr[3])), collapse = " - ")
    }
  })
  args <- c(contrasts_named, list(levels = mm))
  contrasts_limma <- do.call(what = limma::makeContrasts, args = args)
  fit <- limma::lmFit(data, mm)
  cfit <- limma::contrasts.fit(fit = fit, contrasts = contrasts_limma)
  efit <- limma::eBayes(fit = cfit, trend = trend, robust = robust, ...)
  coefnames <- colnames(efit$coefficients)
  results <- lapply(setNames(coefnames, coefnames), function(tmpcoef) limma::topTable(efit, coef = tmpcoef, number = Inf, confint = TRUE, adjust.method = p.adj.method))
  results <- lapply(results, function(tmpres) data.frame(row.names = rownames(tmpres), 
                                                         id = rownames(tmpres), mean = tmpres$AveExpr, stat = tmpres$t, 
                                                         log2FC = tmpres$logFC, CIlower = tmpres$CI.L, CIupper = tmpres$CI.R, pvalue = tmpres$P.Value, padj = tmpres$adj.P.Val))
  
  if (shrink){
    library(ashr)
    results <- lapply(results, function(tmpres){
      tmpres$SE <- (tmpres$CIupper - tmpres$CIlower)/3.92
      ash_out <- ashr::ash(tmpres$log2FC, tmpres$SE, mixcomp="halfuniform")
      tmpres$log2FCshrink <- get_pm(ash_out)
      tmpres
    })
  }
  
  list(data = data, fit = efit, results = results)
}

runDESeq2 <- function(data, design = NULL, formula = ~ 1, contrasts = NULL, lrt_reduced = NULL,
                      prefilter = NULL, postfilter = NULL, min_counts = 5, min_samples = 2,
                      ctrlgenes = NULL, sizefactors = NULL,
                      RUV = list(), SVA = list(),
                      cooksCutoff = TRUE, independentFiltering = TRUE,
                      alpha = 0.05, ordered = TRUE, df = TRUE, ncores = NULL,
                      shrink = TRUE, ihw = TRUE, vst = FALSE, rlog = FALSE,
                      minReplicatesForReplace = 7, fitType = "parametric", ...){
  

  stopifnot(requireNamespace("DESeq2", quietly = TRUE))
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
  if (ifelse(is.null(ncores), TRUE, ncores > 1)) stopifnot(requireNamespace("BiocParallel", quietly = TRUE))
  if (ihw == TRUE) stopifnot(requireNamespace("IHW", quietly = TRUE))
  if (shrink == TRUE) stopifnot(requireNamespace("ashr", quietly = TRUE))
  if (length(RUV) > 0) stopifnot(requireNamespace("RUVSeq", quietly = TRUE))
  if (length(SVA) > 0) stopifnot(requireNamespace("sva", quietly = TRUE))
  
  library(DESeq2)
  
  results <- list()
  
  # Parallel setup ----
  if (is.null(ncores)) ncores <- min(c(10, max(c(1, length(contrasts)))))
  bppar <- NULL
  if (ncores > 1){
    if (tolower(.Platform$OS.type) == "windows"){
      BiocParallel::register(BiocParallel::SnowParam(workers = ncores))
    } else {
      BiocParallel::register(BiocParallel::MulticoreParam(workers = ncores))
    }
    bppar <- BiocParallel::bpparam()
    message(paste0("Parallel on ", ncores, " cores..."))
    
    lapply_fun <- function(X, FUN, BPPARAM = NULL, ...){
      BiocParallel::bplapply(X, FUN = FUN, BPPARAM = BPPARAM, ...)
    }
    
  } else {
    lapply_fun <- function(X, FUN, BPPARAM = NULL, ...){
      lapply(X, FUN = FUN)
    }
  }
  
  
  # Input data ----
  if ("SummarizedExperiment" %in% class(data)){
    cat(crayon::blue("Using 'SummarizedExperiment' as input for 'DESeqDataSet'\n"))
    if (is.null(design)) design <- data.frame(row.names = colnames(assays(data)[[1]]))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSet(data, colData = design, design = formula)
    
  } else if (any(class(data) %in% c("matrix", "data.frame"))){
    cat(crayon::blue("Using raw counts matrix as input for 'DESeqDataSetFromMatrix'\n"))
    if (is.null(design)) design <- data.frame(row.names = colnames(data))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromMatrix(data, colData = design, design = formula)
    
  } else if ("list" %in% class(data)){
    cat(crayon::blue("Using tximport list as input for 'DESeqDataSetFromTximport'\n"))
    if (is.null(design)) design <- data.frame(row.names = colnames(data$abundance))
    for (i in seq_along(data)){
      if (!is.null(ncol(data[[i]]))){
        data[[i]] <- data[[i]][,rownames(design)]
      }
    }
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromTximport(data, colData = design, design = formula)
    
  } else if ("DESeqDataSet" %in% class(data)){
    cat(crayon::blue("Using a 'DESeqDataSet' object as input\n"))
    design <- SummarizedExperiment::colData(data)
    dds <- data
    
  } else {
    stop("Wrong input data format!")
  }
  
  # Pre-filtering ----
  if (!is.null(prefilter)){
    if (length(prefilter) == nrow(dds) & is.logical(prefilter)){
      dds <- dds[naf(prefilter),]
    } else if (is.character(prefilter)){
      dds <- dds[intersect(prefilter, rownames(dds)),]
    } else {
      stop("Index for pre-filtering is of wrong length/format!")
    }
  }
  
  if (!is.null(min_counts) & !is.null(min_samples)){
    dds <- dds[rowSums(DESeq2::counts(dds) > min_counts, na.rm = TRUE) >= min_samples,]
  }
  
  # RUVSeq ----
  if (length(RUV) > 0){
    if (!class(RUV) == "list") stop("Error: 'RUV' must be a named list.")
    if ("n" %in% names(RUV)){
      k <- RUV$n
    } else {
      k <- 1
    }
    
    if ("empirical" %in% names(RUV)){
      ruv <- RUVSeq::RUVg(DESeq2::counts(dds, normalized = FALSE), RUV$empirical, k = k)
      
    } else if ("group" %in% names(RUV)){
      ruv <- RUVSeq::RUVs(DESeq2::counts(dds, normalized = FALSE), scIdx = RUVSeq::makeGroups(as.character(SummarizedExperiment::colData(dds)[[RUV$group]])), k = k)
      
    } else {
      stop("Error: Please provide either 'empirical' or 'group' to run RUVSeq!")
    }
    
    stopifnot(ncol(ruv$W) == k)
    if (!is.null(rownames(ruv$W))) stopifnot(all.equal(rownames(ruv$W), colnames(dds)))
    if (ncol(ruv$W) > 1){
      colnames(ruv$W) <- paste0("RUV", 1:ncol(ruv$W))
    } else {
      colnames(ruv$W) <- "RUV"
    }
    
    design_ruv <- data.frame(SummarizedExperiment::colData(dds) |> as.data.frame(), ruv$W)
    SummarizedExperiment::colData(dds) <- S4Vectors::DataFrame(design_ruv)
    fruv <- as.formula(paste0("~ . + ", paste0(colnames(ruv$W), collapse = " + ")))
    dds <- DESeq2::DESeqDataSet(dds, design = update.formula(old = formula, new = fruv))
  }
  
  # Explicit calculation of sizefactors ----
  if (!is.null(sizefactors)){
    dds$sizeFactor <- sizefactors
  } else if (!is.null(ctrlgenes)){
    if (!all(ctrlgenes %in% rownames(dds))){
      stop("Not all ctrlgenes present in data!")
    }
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = rownames(dds) %in% ctrlgenes)
  }
  
  
  # SVA ----
  if (length(SVA) > 0){
    if (!class(SVA) == "list") stop("Error: 'SVA' must be a named list.")
    if ("n" %in% names(SVA)){
      n <- SVA$n
    } else {
      n <- NULL
    }
    
    if (is.null(DESeq2::sizeFactors(dds))){
      dds <- DESeq2::estimateSizeFactors(dds)
    }
    
    if ("reduced" %in% names(SVA) & "formula" %in% class(SVA$reduced)){
      mm_full <- model.matrix(dds@design, SummarizedExperiment::colData(dds))
      mm_red <- model.matrix(SVA$reduced, SummarizedExperiment::colData(dds))
      svafit <- sva::svaseq(DESeq2::counts(dds, normalized = TRUE), mod = mm_full, mod0 = mm_red)
      
      if (svafit$n.sv > 0){
        if (!is.null(rownames(svafit$sv))) stopifnot(all.equal(rownames(svafit$sv), colnames(dds)))
        if (ncol(svafit$sv) > 1){
          colnames(svafit$sv) <- paste0("SV", 1:ncol(svafit$sv))
        } else {
          colnames(svafit$sv) <- "SV"
        }
        design_sva <- data.frame(SummarizedExperiment::colData(dds) |> as.data.frame(), svafit$sv)
        SummarizedExperiment::colData(dds) <- S4Vectors::DataFrame(design_sva)
        fsva <- as.formula(paste0("~ . + ", paste0(colnames(svafit$sv), collapse = " + ")))
        dds <- DESeq2::DESeqDataSet(dds, design = update.formula(old = formula, new = fsva))
      }
      
    } else {
      stop("Error: Please provide a reduced model formula to run SVA!")
    }
    
  }
  
  # Model fitting ----
  dds <- DESeq2::DESeq(dds, fitType = fitType, minReplicatesForReplace = minReplicatesForReplace, parallel = (ncores > 1), BPPARAM = bppar, ...)
  
  # Post-filtering ----
  if (!is.null(postfilter)){
    if (is.logical(postfilter)) stop("Error: Please provide gene IDs to 'postfilter'!")
    dds <- dds[intersect(postfilter, rownames(dds)),]
  }
  
  # Normalized counts ----
  results$normcounts <- DESeq2::counts(dds, normalized = TRUE)
  if (vst == TRUE) results$vst <- SummarizedExperiment::assays(DESeq2::vst(dds, blind = TRUE))[[1]]
  if (rlog == TRUE) results$rlog <- SummarizedExperiment::assays(DESeq2::rlog(dds, blind = TRUE))[[1]]
  
  # Contrasts ----
  if (!is.null(contrasts)){
    
    if (ihw == TRUE) filterFun <- IHW::ihw else filterFun <- NULL
    
    results$results <- lapply_fun(contrasts, BPPARAM = bppar, FUN = function(tmp){
      
      if (is.null(filterFun)){
        mle <- DESeq2::results(dds, contrast = tmp, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering, alpha = alpha)
      } else {
        mle <- DESeq2::results(dds, filterFun = filterFun, contrast = tmp, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering, alpha = alpha)
      }
      
      if (shrink == TRUE){
        mmse <- DESeq2::lfcShrink(dds, contrast = tmp, res = mle, type = "ashr", svalue = TRUE)
        mmse <- mmse[rownames(mle),]
        mle$svalue <- mmse$svalue
        mle$log2FCshrink <- mmse$log2FoldChange
      }
      
      mle$gene <- rownames(mle)
      mle <- as.data.frame(mle)
      mle <- dplyr::relocate(.data = mle, gene)
      mle <- dplyr::rename(.data = mle, log2FC = log2FoldChange, log2FCse = lfcSE)
      if (ordered == TRUE) mle <- dplyr::arrange(mle, pvalue)
      mle
    })
    
    # Make dataframes
    select_res <- function(contr, res, what){
      df <- res$results[[contr]]
      df$contrast <- contr
      df |> dplyr::select(gene, !!what, contrast)
    }
    
    results$log2FC <- lapply(names(contrasts), select_res, res = results, what = "log2FC") |> Reduce(f = rbind) |>
      tidyr::pivot_wider(names_from = "contrast", values_from = "log2FC") |> as.data.frame() |> col2rownames(gene)
    results$log2FC <-  results$log2FC[rownames(results$normcounts),, drop = FALSE]
    
    results$padj <- lapply(names(contrasts), select_res, res = results, what = "padj") |> Reduce(f = rbind) |>
      tidyr::pivot_wider(names_from = "contrast", values_from = "padj") |> as.data.frame() |> col2rownames(gene)
    results$padj <-  results$padj[rownames(results$normcounts),, drop = FALSE]
  }
  
  if (!is.null(lrt_reduced)){
    dds_rt <- DESeq2::DESeq(dds, test = "LRT", reduced = lrt_reduced,
                            fitType = fitType, minReplicatesForReplace = minReplicatesForReplace,
                            parallel = (ncores > 1), BPPARAM = bppar, ...)
    lrt <- DESeq2::results(dds)
    lrt$log2FoldChange <- NULL
    lrt$lfcSE <- NULL
    lrt$gene <- rownames(lrt)
    lrt <- as.data.frame(lrt)
    lrt <- dplyr::relocate(.data = lrt, gene)
    if (ordered == TRUE) lrt <- dplyr::arrange(lrt, pvalue)
    results$LRT <- lrt
  }
  
  # Results ----
  results$dds <- dds
  results$design <- as.data.frame(SummarizedExperiment::colData(dds))
  results
  
}

runGSEA <- function(data, genesets = NULL, rank_by = "stat", rank_type = "", id_col = NULL, direction_by = NULL, min_size = 3, max_size = 2000, seed = 0, as.df = TRUE, ...){
  
  stopifnot(requireNamespace("clusterProfiler"))
  stopifnot(requireNamespace("fgsea"))
  
  if (is.null(genesets)) {
    genesets <- NA
  }
  if (!"data.frame" %in% class(genesets)) {
    genesets <- convertGeneSets(genesets, to = "data.frame")
  }
  
  ranks <- get_feature_ranks(data, rank_by = rank_by, rank_type = rank_type, id_col = id_col, direction_by = direction_by)
  
  set.seed(seed)
  results <- clusterProfiler::GSEA(ranks, seed = TRUE, eps = 0, minGSSize = min_size, maxGSSize = max_size, TERM2GENE = genesets, pvalueCutoff = 1, pAdjustMethod = "fdr", ...)
  
  if (as.df == TRUE) {
    results <- as.data.frame(results)
    results <- dplyr::select(results, -Description)
    results <- dplyr::mutate(results, core_enrichment = gsub("/", ", ", core_enrichment))
    results <- dplyr::rename(results, term = ID, ES = enrichmentScore, padj = p.adjust, qval = qvalue, gene_set_size = setSize, enriched_genes = core_enrichment)
    results$n_enriched_genes <- sapply(strsplit(results$enriched_genes, split = ", ", fixed = TRUE), FUN = length)
    results <- dplyr::relocate(results, n_enriched_genes, .after = gene_set_size)
  }
  results
}

get_feature_ranks <- function (data, rank_by = "stat", rank_type = "", id_col = NULL, direction_by = NULL, ...){
  data <- as.data.frame(data)
  stopifnot(length(rank_by) == length(rank_type))
  keep <- unique(c(rank_by, id_col, direction_by))
  if (!all(keep %in% colnames(data))) {
    warning(paste0("Column ", setdiff(keep, colnames(data)), " not found in data."))
  }
  keep <- keep[keep %in% colnames(data)]
  df <- data[, keep, drop = FALSE]
  df <- df[!is.na(df[[1]]), , drop = FALSE]
  for (i in seq_along(rank_by)) {
    if (rank_type[i] == "+") {
      if (any(df[[rank_by[i]]] < 0)) {
        stop(paste0("Error: Encountered negative values in column ", 
                    rank_by[i]))
      }
      df[[rank_by[i]]] <- df[[rank_by[i]]] * sign(df[[direction_by]])
    }
    else if (rank_type[i] == "p") {
      if (any(df[[rank_by[i]]] < 0)) {
        stop(paste0("Error: Encountered negative values in column ", 
                    rank_by[i]))
      }
      df[[rank_by[i]]] <- -log10(pmax(df[[rank_by[i]]], .Machine$double.xmin)) * sign(df[[direction_by]])
    }
  }
  if (!is.null(id_col)) {
    df[[id_col]] <- factor(df[[id_col]], levels = sort(unique(df[[id_col]]), decreasing = TRUE), ordered = TRUE)
  }
  df <- dplyr::arrange_all(df)
  df <- df[rev(1:nrow(df)), , drop = FALSE]
  df$rankid <- as.character(df[[1]])
  df$rankid[is.na(df$rankid)] <- "na"
  df$dup <- df$rankid %in% df$rankid[duplicated(df$rankid)]
  df$newstat <- df[[1]]
  if (!is.null(id_col)) {
    df[[id_col]] <- as.character(df[[id_col]])
  }
  for (rankid in unique(subset(df, dup == TRUE)$rankid)) {
    tmp <- df[df$rankid == rankid, , drop = FALSE]
    ix <- range(which(rankid == df$rankid))
    tmp2 <- df[c(ix[1] - 1, ix[2] + 1), , drop = FALSE]
    tmp2[is.na(tmp2)] <- mean(abs(diff(df[[1]])))
    maxdiff <- min(abs(tmp2[[1]] - mean(tmp[[1]]))) * 0.4
    tmpranks <- rev(1:nrow(tmp))
    tmpranks <- (tmpranks - median(tmpranks))/max(abs(tmpranks))
    tmp$newstat <- tmp[[1]] + maxdiff * tmpranks
    df[df$rankid == rankid, ]$newstat <- tmp$newstat
  }
  if (is.null(id_col)) {
    ranks <- setNames(df$newstat, rownames(df))
  }
  else {
    ranks <- setNames(df$newstat, df[[id_col]])
  }
  ranks
}

aggregate_counts <- function(seu, group_by, assay = "RNA", layer = "counts", min_cells = 10, ...){
  
  expr <- PseudobulkExpression(seu, assays = assay, group.by = group_by, layer = layer, method = "aggregate", normalization.method = NULL, ...)
  expr <- expr[[assay]]
  
  metadf <- seu@meta.data %>%
    dplyr::group_by(.data[[group_by]]) %>%
    summarise(
      n_cells = dplyr::n(),
      across(
        everything(),
        ~ if (dplyr::n_distinct(.) == 1) unique(.) else NA,
        .names = "{.col}"
      ), .groups = "drop") |> as.data.frame()
  
  rownames(metadf) <- metadf[[group_by]] |> gsub(pattern = "_", replacement = "-")
  metadf <- metadf[,colSums(!is.na(metadf)) > 0]
  metadf <- metadf |> dplyr::filter(n_cells >= min_cells)
  expr <- expr[,colnames(expr) %in% rownames(metadf)]
  metadf <- metadf[colnames(expr),]
  
  list(expr=expr, meta=metadf)
}

get_contable <- function(df, x, y){
  
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  res <- df |>
    dplyr::group_by(!!x, !!y) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = !!y,
      values_from = n,
      values_fill = 0
    )
  res <- as.data.frame(res)
  rownames(res) <- res[[rlang::as_string(x)]]
  res[[rlang::as_string(x)]] <- NULL
  if (is.factor(df[[rlang::as_string(y)]])) {
    res <- res[, levels(df[[rlang::as_string(y)]]) , drop = FALSE]
  }
  if (is.factor(df[[rlang::as_string(x)]])) {
    res <- res[levels(df[[rlang::as_string(x)]]), , drop = FALSE]
  }
  stopifnot(sum(res) == nrow(tidyr::drop_na(df[, c(rlang::as_string(x), rlang::as_string(y))])))
  return(res)
}


### Plotting functions

ggboxplot <- function(df, x, y, q=0, ymax=NULL, points_color=NULL, point_size=2.5, colors=NULL, fontsize=16, linewidth=NULL, varwidth=FALSE, legend=FALSE, rot_x = 0, y_at_zero=TRUE){
  stopifnot(requireNamespace("ggbeeswarm"))
  gg <- ggplot(df, aes_string(x, y))
  if (!is.null(points_color)){
    gg <- gg + ggbeeswarm::geom_quasirandom(aes_string(colour=points_color), size=point_size, varwidth = varwidth, stroke = NA)
    if (is.null(colors)){
      n <- length(unique(df[[x]][!is.na(df[[x]])]))
      colors <- get_colors(n)
      colors <- colorspace::darken(colors, amount = 0.15)
    }
    gg <- gg + scale_colour_manual(values = colors)
  } else {
    gg <- gg + ggbeeswarm::geom_quasirandom(varwidth = varwidth, size=point_size, stroke = NA, colour="grey30")
  }
  if (is.null(linewidth)){
    linewidth <- fontsize/22
  }
  if (is.null(ymax)){
    ymax <- quantile(df[[y]], q, na.rm = TRUE)
  }
  gg <- gg + geom_boxplot(colour = "black", linewidth = linewidth, fatten = 1.5, varwidth = varwidth, fill = NA, outliers = FALSE, staplewidth = 0.4)
  gg <- gg + scale_y_continuous(limits = c(ifelse(y_at_zero, 0, NA), ymax), expand = ggplot2::expansion(mult = c(0,0)), oob = scales::censor)
  gg <- gg + add_theme(base_size = fontsize, legend = legend, rot_x = rot_x) + coord_cartesian(clip = "off")
  gg
}

get_stat_df <- function(df, x, y, padj_thres = 0.05, spacing = 0.05, ymax=NULL, max_add = NULL, na.omit = FALSE, q = 0.99, TESTFUN = function(x1, x2, ...){ wilcox.test(x1, x2)$p.value }, ...){
  statdf <- get_all_comparisons(df, x, y, as_list = FALSE, na.omit = na.omit)
  colnames(statdf) <- c("x","xend")
  statdf$p <- as.data.frame(t(statdf)) |> sapply(function(groups){
    g1 <- df[[y]][df[[x]] == groups[1]]
    g2 <- df[[y]][df[[x]] == groups[2]]
    TESTFUN(g1, g2)
  })
  statdf$padj <- p.adjust(statdf$p, method = "BH")
  statdf$comparison <- paste0(statdf$xend, "_vs_", statdf$x)
  xall <- unique(union(statdf$x, statdf$xend)) |> sort()
  if (is.null(levels(df[[x]]))){
    xlev <- setNames(as.numeric(factor(xall, ordered = TRUE, levels = xall)), xall)
  } else {
    statdf$x <- factor(statdf$x, levels = levels(df[[x]]), exclude = "", ordered = is.ordered(df[[x]]))
    statdf$xend <- factor(statdf$xend, levels = levels(df[[x]]), exclude = "", ordered = is.ordered(df[[x]]))
    xlev <- setNames(as.numeric(factor(xall, ordered = TRUE, levels = levels(df[[x]]))), xall)
  }
  statdf$xmean <- rowMeans(cbind(xlev[as.character(statdf$x)], xlev[as.character(statdf$xend)]))
  statdf <- subset(statdf, padj <= padj_thres)
  if (!is.null(max_add) & nrow(statdf) > 1){
    statdf <- statdf[1:min(nrow(statdf), max_add),]
  }
  statdf$yoffset <- rank(1 - statdf$p, ties.method = "random") |> scale(center = FALSE) |> as.numeric()
  statdf$yoffset <- (statdf$yoffset - min(statdf$yoffset))
  statdf$yoffset <- statdf$yoffset/max(statdf$yoffset) * (length(statdf$yoffset) - 1)
  if (nrow(statdf) == 1){
    statdf$yoffset <- 0
  }
  ymax2 <- max(quantile(df[[y]], q, na.rm = TRUE), na.rm = TRUE)
  if (is.null(ymax)){
    ymax <- ymax2
  }
  statdf$yoffset <- ymax*spacing*statdf$yoffset
  statdf$y <- ymax2 + statdf$yoffset
  statdf$text <- format.pval(statdf$padj, digits = 3, eps = 10^-9)
  statdf
}

add_comparisons <- function(gg, padj_thres = 0.05, fontsize = 10, max_add = NULL, spacing = NULL, label_offset=0.02, linewidth=NULL, linetype="solid", y_at_zero=TRUE, ymax=NULL, na.omit = FALSE, TESTFUN = function(x, y, ...){ t.test(x, y)$p.value }, ...){
  x <- gg$mapping[["x"]] |> rlang::as_name()
  y <- gg$mapping[["y"]] |> rlang::as_name()
  
  if (is.null(spacing)){
    spacing <- fontsize/200
  }
  statdf <- get_stat_df(gg$data, x, y, padj_thres = padj_thres, spacing = spacing, max_add = max_add, ymax=ymax, na.omit = na.omit, TESTFUN = TESTFUN, ...)
  if (nrow(statdf) < 1) return(gg)
  if (is.null(linewidth)){
    linewidth <- sqrt(fontsize/ggplot2::.pt) / 4
  }
  ar <- grid::arrow(90, length = grid::unit(0.01,"npc"), ends = "both")
  if (is.null(ymax)){
    ymax <- max(statdf$y, na.rm = TRUE) * 1.1
  }
  if (y_at_zero){
    ydist <- ymax
  } else {
    ydist <- ymax - min(gg$data[[y]], na.rm = TRUE)
  }
  gg <- gg +
    geom_segment(statdf, mapping = aes_string(x = "x", xend = "xend", y = "y", yend = "y"),
                 arrow = ar, linewidth = linewidth, lineend = "square", linetype = linetype, inherit.aes = FALSE) +
    geom_text(data = statdf, mapping = aes(x = xmean, y = y, label = text), hjust = 0.5, nudge_y = ydist*label_offset, inherit.aes = FALSE, size = fontsize/ggplot2::.pt)
  gg <- gg + scale_y_continuous(limits = c(ifelse(y_at_zero, 0, NA), ymax), expand = ggplot2::expansion(mult = c(0,0)), oob = scales::censor)
  gg
}

add_theme <- function(axis=TRUE, ticks=TRUE, rot_x=0, rot_y=0, grid=FALSE, base_size = 16, ticks_length = 0.15, ticks_color="black", grid_color = "black", legend_line_color="white", legend=TRUE){
  
  th <- ggplot2::theme_bw(base_size = base_size)
  th <- ggplot2::`%+replace%`(th, ggplot2::theme(panel.border = element_blank(),
                                                 legend.ticks = element_line(colour = legend_line_color),
                                                 panel.grid.major = element_line(colour = grid_color),
                                                 panel.grid.minor = element_blank(),
                                                 axis.line = element_line(),
                                                 axis.ticks = element_line(colour = ticks_color, lineend = "square"),
                                                 axis.ticks.length =unit(ticks_length, "cm"),
                                                 legend.text = element_text(size = base_size),
                                                 legend.justification = c(0,0.5),
                                                 legend.background = element_blank(),
                                                 plot.title = element_text(size = base_size, margin = margin(b = base_size/2)),
                                                 axis.text = element_text(colour = "black", margin = margin(base_size/5, base_size/5, base_size/5, base_size/5))))
  
  if (axis == FALSE){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(axis.line = element_blank()))
  }
  if (ticks == FALSE){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(axis.ticks = element_blank()))
  }
  if (grid == FALSE){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(panel.grid.major = element_blank()))
  }
  if (rot_x < 0){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(axis.text.x = element_text(angle = rot_x, hjust = 0, vjust = 0.5)))
  }
  if (rot_x > 0){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(axis.text.x = element_text(angle = rot_x, hjust = 1, vjust = 0.5)))
  }
  if (rot_y != 0){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(axis.text.y = element_text(angle = rot_y, hjust = 1, vjust = 0.5)))
  }
  if (legend == FALSE){
    th <- ggplot2::`%+replace%`(th, ggplot2::theme(legend.position = "none"))
  }
  
  th
}

PlotDotplot <- function(seu, group.by = NULL, features = NULL, fontsize=15, dot.scale = 5, col.mid = NULL, scale = FALSE, rescale=FALSE, guide_title="mean", col.low = "grey95", col.high = "#0000b5", col.min = NA, col.max = NA, size.min=NA, size.max=NA, lwd = 0.5, ...){
  
  gg <- DotPlot(seu, group.by = group.by, cols = c("blue", "grey97", "red"), features = features, dot.scale = dot.scale, scale = scale, ...) +
    theme(axis.text.x = element_text(size = fontsize, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = fontsize), legend.box = "horizontal")
  
  if (rescale){
    df <- gg$data |> dplyr::rename(cell=id, id=features.plot)
    df <- df |>
      group_by(id) |>
      mutate(color_val = scales::rescale(avg.exp, to = c(0, 1))) %>%
      ungroup()
    
    gg <- ggplot(df, aes(x = id, y = cell, size = pct.exp, color = color_val)) + theme_basic(fontsize = fontsize) +
      geom_point() +
      theme(axis.text.x = element_text(size = fontsize, angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = fontsize),
            legend.box = "horizontal")
    
  }
  
  gg <- gg + xlab("") + ylab("")
  if (is.null(col.mid)){
    gg <- gg + scale_color_gradient(low = col.low, high = col.high, limits = c(col.min, col.max), oob = scales::squish)
  } else {
    gg <- gg + scale_color_gradient2(low = col.low, mid = col.mid, high = col.high, midpoint = 0, oob = scales::squish)
  }
  gg <- gg + scale_size(
    range = c(0, dot.scale),
    limits = c(size.min, size.max),
    guide = "legend"
  )
  
  gg <- gg + guides(color = guide_colorbar(title = guide_title, frame.colour = "black", ticks.colour = "black", frame.linewidth = lwd), size = guide_legend(title = "cells (%)"))
  gg
}

ggforest <- function(df,
                     aes_plot=aes(x = estimate, y = var),
                     aes_point=aes(color = signif, shape = signif),
                     aes_errorbar=aes(xmin = conf.low, xmax = conf.high, color = signif),
                     xlim=c(0.3, 3), point_size=2.5, errorbar_height=0.5, breaks=scales::log_breaks(base = 10), position=position_identity(), ...){
  
  inrange <- df$conf.low >= xlim[1] & df$conf.high <= xlim[2]
  df_errorbar <- df[inrange, ,drop=FALSE]
  df_nobar <- df[!inrange, ,drop=FALSE]
  df_errorbar$var <- factor(df_errorbar$var, levels = levels(df$var))
  df_nobar$var    <- factor(df_nobar$var,    levels = levels(df$var))
  
  gg <- ggplot(df, mapping=aes_plot) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
    geom_point(mapping=aes_point, size = point_size, show.legend = TRUE, position=position) +
    geom_errorbarh(data=df_errorbar, mapping=aes_errorbar, height = errorbar_height, position=position, show.legend = FALSE) +
    geom_errorbarh(data=df_nobar, mapping=aes_errorbar, height = 0, linetype = "dashed", position=position, show.legend = FALSE)
  
  gg <- gg +
    scale_color_manual(name="signif.", values = c("black", "#8E0E0A"), limits = c("not signif.", "signif."), drop=FALSE) +
    scale_shape_manual(name="signif.", values = c(16, 15), limits = c("not signif.", "signif."), drop=FALSE)
  
  gg <- gg +
    scale_x_log10(limits = xlim, oob=scales::squish, breaks = breaks,, labels = scales::label_number()) +
    labs(x = "HR") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t = 15, b = 15),
          axis.text = element_text(colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.title.y = element_text(margin = margin(r = 20)))
  
  gg
}

ggforest_annotated <- function(df, xlim = c(0.6,1.4), breaks = c(0.6, 0.8, 1, 1.2, 1.4), point_size=3.6, errorbar_height=0.55, fontsize=3, title=NULL, cex_title=1.5, ylab=NULL, panel_widths = c(0.5, 0.4, 1)){
  p_varname <- ggplot(df, aes(y = var, x = 1, label = var)) +
    geom_text(hjust = 0, size = fontsize) +
    xlim(0.5, 2) +
    ylab(ylab) + 
    theme_void() + theme(axis.title.y = element_text(angle = 90, size = fontsize*ggplot2:::.pt))
  
  p_HR <- ggplot(df, aes(y = var, x = 1, label = label)) +
    geom_text(hjust = 0, size = fontsize) +
    xlim(0.5, 3) +  
    theme_void()
  
  gg <- ggforest(df, aes_plot = aes(x = estimate, y = var), point_size = point_size, errorbar_height = errorbar_height, xlim = xlim, breaks = breaks)
  p_forest <- gg + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank()
  )
  
  gg <- wrap_plots(p_varname, p_HR, p_forest, ncol = 3, widths = panel_widths) & coord_cartesian(clip = "off")
  gg <- gg & theme(plot.title = element_text(size = fontsize*ggplot2:::.pt * cex_title, hjust = 0.5), plot.margin = margin(t = 5))
  gg + plot_annotation(title = title)
}

convertGeneSets <- function(genesets, from = NULL, to = "list", term = term, gene = gene){
  
  stopifnot(requireNamespace("GSEABase"))
  term <- rlang::enquo(term)
  gene <- rlang::enquo(gene)
  
  if (is.null(from)){
    if ("data.frame" %in% class(genesets)) from <- "dataframe"
    if ("list" %in% class(genesets)) from <- "list"
    if ("GeneSetCollection" %in% class(genesets)) from <- "GeneSetCollection"
  }
  from <- trimws(tolower(from[1]))
  to <- trimws(tolower(to[1]))
  
  gsc <- FALSE
  if (to %in% c("list", "l")) to <- "list"
  if (to %in% c("dataframe", "df", "data.frame", "table")) to <- "dataframe"
  if (to %in% "genesetcollection"){
    to <- "list"
    gsc <- TRUE
  }
  
  if (from == "GeneSetCollection" & to == "GeneSetCollection") return(genesets)
  if (from == "list" & to == "dataframe"){
    genesets <- utils::stack(genesets)[,c(2,1)]
    colnames(genesets) <- c(rlang::as_name(term), rlang::as_name(gene))
  }
  if (from == "dataframe" & to == "list"){
    genesets <- dplyr::select(genesets, !!gene, !!term)
    genesets <- unstack(unique(genesets))
  }
  if (gsc == TRUE){
    genesets <- lapply(setNames(seq_along(genesets), names(genesets)), function(i){
      GSEABase::GeneSet(genesets[[i]], setName = names(genesets[i]))
    })
    genesets <- GSEABase::GeneSetCollection(genesets)
  }
  genesets
}

prettify_gsea_term <- function(x) {
  library(stringr)
  if (is.factor(x)) x <- as.character(x)
  pref_pat <- regex("^(GO(?::)?_?(BP|MF|CC)|GOBP|GOMF|GOCC|HALLMARK|REACTOME|KEGG)[ _]*", ignore_case = TRUE)
  pref_raw <- str_extract(x, pref_pat)
  std_prefix <- dplyr::case_when(
    str_detect(pref_raw, regex("GO(?::)?_?BP|GOBP", ignore_case = TRUE)) ~ "GOBP",
    str_detect(pref_raw, regex("GO(?::)?_?MF|GOMF", ignore_case = TRUE)) ~ "GOMF",
    str_detect(pref_raw, regex("GO(?::)?_?CC|GOCC", ignore_case = TRUE)) ~ "GOCC",
    str_detect(pref_raw, regex("HALLMARK", ignore_case = TRUE))           ~ "HALLMARK",
    str_detect(pref_raw, regex("REACTOME", ignore_case = TRUE))           ~ "REACTOME",
    str_detect(pref_raw, regex("KEGG", ignore_case = TRUE))               ~ "KEGG",
    TRUE ~ NA_character_
  )
  body <- ifelse(is.na(pref_raw), x, str_remove(x, pref_pat))
  body <- body |>
    str_replace_all("[_.]+", " ") |>
    str_squish() |>
    str_to_lower()
  small <- c("a","an","the","and","or","nor","but","for","to","at","by","of","on","in","per","via","with","from","into","over","as")
  tc <- function(s) {
    s_tc <- tools::toTitleCase(s)
    words <- str_split(s_tc, " ", simplify = FALSE)[[1]]
    if (length(words)) {
      for (i in seq_along(words)) {
        wl <- str_to_lower(words[i])
        if (i != 1 && wl %in% small) words[i] <- wl
      }
    }
    paste(words, collapse = " ")
  }
  body <- vapply(body, tc, character(1L))
  acronyms <- c(
    "DNA","RNA","mRNA","tRNA","rRNA","snRNA","lncRNA", "E2F", "ATP", "ADP", "AMP", "G2M", "ECM", "B Cell", "T Cell", "IGA", "CCR", "CCL",
    "MHC","HLA","TCR","BCR","GPCR","MAPK","JAK","STAT", "UV",
    "EGFR","VEGF","PI3K","AKT","MTOR","IFN","TNF", "TNFA", "IL",
    "ERK","PKC","AP-1","CREB","WNT","TGF","SMAD","HIF1A","CXCL","GABA",
    "NFKB","NFkB"
  )
  for (a in acronyms) {
    patt <- regex(paste0("\\b", a, "\\b"), ignore_case = TRUE)
    body <- str_replace_all(body, patt, a)
  }
  roman <- c("I","II","III","IV","V","VI","VII","VIII","IX","X")
  patt_roman <- regex(paste0("\\b(", paste(roman, collapse="|"), ")\\b"), ignore_case = TRUE)
  body <- str_replace_all(body, patt_roman, toupper)
  out <- ifelse(is.na(std_prefix), body, paste0(std_prefix, " ", body))
  out
}

ggsea_split_dotplot <- function(data, pal, title=NULL, top_n_up=15, top_n_down=15, xmin=1.5, xmax=3, fontsize=12, plim=c(0,NA), size_range=c(0.01, 5), maxchar=NULL, margins=ggplot2::margin(t = 8, r = 8, l = 8, b = 8)){
  data <- dplyr::arrange(data, desc(NES))
  data <- rbind(dplyr::slice_head(data, n = top_n_up), dplyr::slice_tail(data, n = top_n_down))
  data$term <- as.character(data$term) |> prettify_gsea_term()
  data$term <- factor(data$term, levels = rev(data$term), ordered = TRUE)
  if (!is.null(maxchar)) data$term <- cutstr(data$term, maxchar = maxchar)
  
  data$category <- factor(data$category, levels = unique(data$category), ordered = TRUE)
  y_scale <- scale_y_discrete(limits = levels(data$term), drop = FALSE, position = "left")
  
  col_scale <- scale_color_manual(values = pal, breaks = names(pal), limits = names(pal), drop = FALSE)
  if (is.na(plim[2])){
    plim[2] <- ceiling(max(-log10(data$padj)))
  }
  size_scale <- scale_size_continuous(limits = plim, range = size_range)
  
  gg1 <- ggplot(subset(data, NES < 0), aes(x = NES, y = term)) + theme_basic(grid.color = "grey90")
  gg1 <- gg1 + geom_point(aes(color = category, size = -log10(padj)), show.legend = TRUE)
  gg1 <- gg1 + y_scale + col_scale + size_scale + xlim(-xmax,-xmin) + theme(plot.margin = margins)
  gg1 <- gg1 + ggtitle(title)
  
  gg2 <- ggplot(subset(data, NES > 0), aes(x = NES, y = term)) + theme_basic(grid.color = "grey90")
  gg2 <- gg2 + geom_point(aes(color = category, size = -log10(padj)), show.legend = TRUE) #+ xlab("up")
  gg2 <- gg2 + y_scale + col_scale + size_scale + xlim(xmin,xmax) + theme(plot.margin = margins) + coord_cartesian(clip = "off")
  
  gg <- (wrap_plots(gg1, gg2, nrow = 1, axes = "collect", axis_titles="collect", guides = "collect") & ylab(""))

  et <- element_text(size = fontsize)
  gg <- gg & theme(text = et, legend.title = et, axis.title = et,  line = element_line(linewidth = 0.5, lineend = "square"), axis.ticks = element_line(linewidth = 0.5, lineend = "butt"))
  gg <- gg & guides(color = guide_legend(override.aes = list(size = 5, label = "")))
  gg
}

get_top_genes <- function(data, x="log2FC", y="padj", n=20){
  data <- subset(data, name != "")
  data$xtmp <- data[[x]]
  data$ytmp <- -log10(data[[y]])
  data$score <- abs(as.numeric(scale(data$xtmp, center = FALSE))) + abs(as.numeric(scale(data$ytmp, center = FALSE)))
  top_up <- subset(data, xtmp > 0) |> dplyr::arrange(desc(score)) |> head(n=n) |> pull(name)
  top_dn <- subset(data, xtmp < 0) |> dplyr::arrange(desc(score)) |> head(n=n) |> pull(name)
  union(top_up, top_dn)
}

roundup <- function(x, digits = 0, ...){
  ceiling(x * 10^digits, ...) / 10^digits
}

rounddown <- function(x, digits = 0, ...){
  floor(x * 10^digits, ...) / 10^digits
}

ggvolcano <- function(data, x = NULL, y = NULL, color = NULL, label = NULL,
                      shape = NULL, stroke = NA, nlabels = NULL, fontsize = 16, 
                      labface = "plain", lab_size = 12, repel = 2, attract = NULL, 
                      box.padding = 0.5, max_overlaps = Inf, seed = 123, ptres = 0.05, 
                      clip = FALSE, symlim = TRUE, expand = c(0, 0), nbreaks_x = 7, 
                      nbreaks_y = 7, xlim = NULL, ylim = NULL, color_up = "#eb9d0e", 
                      color_down = "#146bc7", color_nonsig = "#4d4d4d", label_up = "up", 
                      label_down = "down", label_nonsig = "not signif.", autolabel = NULL, 
                      show_nonsig = TRUE, segment.alpha = 0.6, max.time = 30, max.iter = 10^6, 
                      show_grid = TRUE, title = NULL, title_size = NULL, point_size = 2, raster = FALSE, raster_dpi = 600,
                      scale_size = FALSE, axis_size = NULL, leg_size = NULL, leg_key_size = 4,
                      lwd = 0.8, at_zero = FALSE, clip_frame = "off", ...){
  
  `%<>%` <- magrittr::`%<>%`
  
  data <- as.data.frame(data)
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  if (rlang::quo_is_null(x)) {
    x <- rlang::sym(grep("lfc|log2FoldChange|logFC|log2FC|nes",  names(data), value = TRUE, ignore.case = TRUE)[1])
  }
  if (rlang::quo_is_null(y)) {
    y <- rlang::sym(grep("padj|fdr", names(data), value = TRUE, ignore.case = TRUE)[1])
  }
  data$x <- data[[rlang::as_name(x)]]
  data$y <- -log10(data[[rlang::as_name(y)]])
  data <- data[!is.na(data$x) & !is.na(data$y), ]
  data$xtmp <- data$x
  data$xtmp[is.infinite(data$xtmp)] <- max(abs(data$xtmp[!is.infinite(data$xtmp)])) * 
    sign(data$xtmp[is.infinite(data$xtmp)])
  data$ytmp <- data$y
  data$ytmp[is.infinite(data$ytmp)] <- max(abs(data$ytmp[!is.infinite(data$ytmp)])) * 
    sign(data$ytmp[is.infinite(data$ytmp)])
  data$score <- abs(as.numeric(scale(data$xtmp, center = FALSE))) + 
    abs(as.numeric(scale(data$ytmp, center = FALSE)))
  data$score[is.na(data$score)] <- 0
  if (!is.null(autolabel)) {
    label_up <- sub("_vs_.*", "", autolabel)
    label_down <- sub(".*_vs_", "", autolabel)
  }
  data$class <- label_nonsig
  data$class[data[[rlang::as_name(y)]] <= ptres & data$x > 0] <- label_up
  data$class[data[[rlang::as_name(y)]] <= ptres & data$x < 0] <- label_down
  data$score[data$class == label_nonsig] <- data$score[data$class == label_nonsig] * 0.001
  if (is.null(title_size)) 
    title_size <- fontsize
  if (is.null(axis_size)) 
    axis_size <- fontsize
  if (is.null(leg_size)) 
    leg_size <- fontsize
  if (is.null(lab_size)) 
    lab_size <- fontsize
  shape <- rlang::enquo(shape)
  label <- rlang::enquo(label)
  if (rlang::quo_is_null(label)) {
    data[["label"]] <- rownames(data)
  }
  else {
    data[["label"]] <- data[[rlang::as_name(label)]]
  }
  data <- data[order(data$score, decreasing = TRUE), ]
  if (is.null(nlabels)) {
    nlabels <- min(20, ceiling(nrow(data)/10))
  }
  if (is.infinite(nlabels)) {
    nlabels <- nrow(data)
  }
  data$do_label <- FALSE
  nlabels_left <- nlabels_right <- 0
  if (nrow(subset(data, x < 0)) > 0) 
    nlabels_left <- ceiling(nlabels/2 * max(subset(data, x < 0)$score, na.rm = TRUE)/max(data$score, na.rm = TRUE))
  if (nrow(subset(data, x > 0)) > 0) 
    nlabels_right <- ceiling(nlabels/2 * max(subset(data, x > 0)$score, na.rm = TRUE)/max(data$score, na.rm = TRUE))
  if (is.na(nlabels_left)) 
    nlabels_left <- 0
  if (is.na(nlabels_right)) 
    nlabels_right <- 0
  data$do_label[which(data$x < 0)[1:nlabels_left]] <- TRUE
  data$do_label[which(data$x > 0)[1:nlabels_right]] <- TRUE
  data$do_label[is.na(data$do_label)] <- FALSE
  if (sum(data$do_label) < nlabels) {
    data$do_label[!data$do_label][1:(nlabels - sum(data$do_label))] <- TRUE
  }
  data$label[!data$do_label] <- ""
  data$do_label[data$label == ""] <- FALSE
  data$do_label[data$class == label_nonsig] <- FALSE
  color <- rlang::enquo(color)
  color_user_def <- rlang::quo_is_null(color)
  if (color_user_def) {
    color <- rlang::sym("class")
  }
  sigdata <- subset(data, class != label_nonsig)
  xylimits <- list(xlim = getLimits(sigdata$xtmp, clip = clip, 
                                    expand = expand[1]), ylim = getLimits(sigdata$ytmp, clip = clip, expand = expand[2], negative = FALSE))
  if (symlim == TRUE) {
    xylimits$xlim <- c(min = -max(abs(xylimits$xlim)), max = max(abs(xylimits$xlim)))
  }
  data$xorg <- data$x
  data$yorg <- data$y
  if (!is.null(xlim)) 
    xylimits$xlim <- setNames(xlim, c("min", "max"))
  if (!is.null(ylim)) 
    xylimits$ylim <- setNames(ylim, c("min", "max"))
  xclip_min <- any(naf(data$xorg < xylimits$xlim["min"]))
  xclip_max <- any(naf(data$xorg > xylimits$xlim["max"]))
  xbreaks <- (scales::pretty_breaks(n = nbreaks_x))(xylimits$xlim, n = nbreaks_x)
  if (clip & xclip_min) {
    xylimits$xlim["min"] <- min(xbreaks)
    xclip_min <- any(data$xorg < xylimits$xlim["min"])
  }
  if (clip & xclip_max) {
    xylimits$xlim["max"] <- max(xbreaks)
    xclip_max <- any(data$xorg > xylimits$xlim["max"])
  }
  xylimits$xlim <- xylimits$xlim + c(-diff(xylimits$xlim), diff(xylimits$xlim)) * c(!xclip_min, !xclip_max) * 0.02
  names(xbreaks) <- as.character(xbreaks)
  if (xclip_min) {
    names(xbreaks)[1] <- paste0("<", xbreaks[1])
  }
  if (xclip_max) {
    names(xbreaks)[length(xbreaks)] <- paste0(">", xbreaks[length(xbreaks)])
  }
  yclip_min <- any(naf(data$yorg < xylimits$ylim["min"]))
  yclip_max <- any(naf(data$yorg > xylimits$ylim["max"]))
  ybreaks <- (scales::pretty_breaks(n = nbreaks_y))(xylimits$ylim, n = nbreaks_y)
  if (clip & yclip_min) {
    xylimits$ylim["min"] <- min(ybreaks)
    yclip_min <- any(data$xorg < xylimits$ylim["min"])
  }
  if (clip & yclip_max) {
    xylimits$ylim["max"] <- max(ybreaks)
    yclip_max <- any(data$xorg > xylimits$ylim["max"])
  }
  xylimits$ylim <- xylimits$ylim + c(-diff(xylimits$ylim), diff(xylimits$ylim)) * c(!yclip_min & !at_zero, !yclip_max) * c(0.01, 0.05)
  names(ybreaks) <- as.character(ybreaks)
  if (yclip_min) {
    names(ybreaks)[1] <- paste0("<", ybreaks[1])
  }
  if (yclip_max) {
    names(ybreaks)[length(ybreaks)] <- paste0(">", ybreaks[length(ybreaks)])
  }
  data$x[data$x < xylimits$xlim["min"]] <- xylimits$xlim["min"]
  data$x[data$x > xylimits$xlim["max"]] <- xylimits$xlim["max"]
  data$y[data$y < xylimits$ylim["min"]] <- xylimits$ylim["min"]
  data$y[data$y > xylimits$ylim["max"]] <- xylimits$ylim["max"]
  data <- data[order(data$score, decreasing = FALSE), ]
  
  gg <- data %>% ggplot2::ggplot(ggplot2::aes(x = x, y = y, label = label, color = !!color, shape = !!shape, ...))
  if (show_grid) {
    panel_grid <- ggplot2::element_line(linewidth = lwd, 
                                        color = rgb(0.9, 0.9, 0.9))
  }
  else {
    panel_grid <- ggplot2::element_blank()
  }
  
  gg %<>% +ggplot2::theme_bw(base_size = 20)
  gg %<>% +ggplot2::theme(text = ggplot2::element_text(color = "black", size = fontsize),
                          rect = ggplot2::element_rect(color = "black", linewidth = lwd),
                          line = ggplot2::element_line(linewidth = lwd), 
                          legend.text = ggplot2::element_text(color = "black", size = leg_size),
                          legend.title = ggplot2::element_text(color = "black", size = leg_size),
                          legend.key.spacing.x = ggplot2::unit(leg_size, "pt"),
                          legend.key.spacing.y = ggplot2::unit(leg_size * 0.6, "pt"),
                          panel.grid.minor = ggplot2::element_blank(), 
                          panel.grid.major = panel_grid,
                          panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = ifelse(clip_frame == "on", lwd * 2, lwd)),
                          strip.background = ggplot2::element_blank(), 
                          strip.text = ggplot2::element_text(color = "black", size = title_size), 
                          axis.ticks = ggplot2::element_line(color = "black", linewidth = lwd), 
                          axis.line = ggplot2::element_blank(), plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm"),
                          plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, lineheight = 1.5), axis.title = ggplot2::element_text(size = axis_size, face = "plain"), axis.text = ggplot2::element_text(size = axis_size, color = "black"))

  if (!is.null(ptres)) {
    gg %<>% +ggplot2::geom_hline(yintercept = -log10(ptres), 
                                 linetype = "dashed", color = rgb(0.3, 0.3, 0.3))
  }
  if (scale_size == FALSE) {
    if (!raster){
      gg %<>% + ggplot2::geom_point(size = point_size, alpha = 0.8, stroke = stroke)
    } else {
      gg %<>% + ggrastr::geom_point_rast(size = point_size, alpha = 0.8, stroke = stroke, raster.dpi = raster_dpi)
    }
    
  }
  else {
    gg %<>% +ggplot2::geom_point(aes(size = score), alpha = 0.8, stroke = stroke)
    gg %<>% +ggplot2::scale_size_continuous(range = c(point_size/5, point_size * 2), guide = "none")
  }
  if (color_user_def) {
    color_vals <- setNames(c(color_up, color_down, color_nonsig), c(label_up, label_down, label_nonsig))
    if (show_nonsig) {
      color_breaks <- c(label_up, label_nonsig, label_down)
    }
    else {
      color_breaks <- c(label_up, label_down)
    }
    gg %<>% +ggplot2::scale_colour_manual(values = color_vals, breaks = color_breaks)
    gg %<>% +ggplot2::labs(title = title, y = paste0("-log10 ", rlang::as_name(y)), x = rlang::as_name(x), color = NULL)
  }
  else {
    gg %<>% +ggplot2::labs(title = title, y = paste0("-log10 ", rlang::as_name(y)), x = rlang::as_name(x))
  }
  gg %<>% +ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0)), limits = xylimits$xlim, breaks = xbreaks, labels = names(xbreaks))
  gg %<>% +ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0)), limits = xylimits$ylim, breaks = ybreaks, labels = names(ybreaks))
  
  if (is.null(attract)) 
    attract <- sqrt(repel)
  gg %<>% +ggrepel::geom_text_repel(data = subset(data, do_label == TRUE),
                                    fontface = labface,
                                    size = lab_size/ggplot2:::.pt, 
                                    seed = seed,
                                    xlim = xylimits$xlim - c(-diff(xylimits$xlim), diff(xylimits$xlim)) * 0.18,
                                    ylim = xylimits$ylim - c(-diff(xylimits$ylim) * 0.15, diff(xylimits$ylim) * 0.02),
                                    force = repel, force_pull = attract, max.overlaps = max_overlaps, 
                                    point.padding = 0.35, box.padding = box.padding, max.time = max.time, 
                                    max.iter = max.iter, min.segment.length = 0, vjust = 0, 
                                    color = rgb(0, 0, 0), segment.alpha = segment.alpha)
  gg %<>% +ggplot2::guides(size = "none", label = "none", color = ggplot2::guide_legend(override.aes = list(size = leg_key_size)))
  
  gg %<>% +ggplot2::coord_cartesian(clip = clip_frame)
  return(gg)
}

cxheatmap <- function(data, rowdf = NULL, coldf = NULL, scale = FALSE, cluster_rows = NULL, cluster_cols = NULL,
                      rowdf_side = "left", coldf_side = "top", rowdf_legend = TRUE, coldf_legend = TRUE,
                      legend_border = "black", anno_border = "black",
                      fontsize = 12, rowcex = NULL, colcex = 1,
                      rownames_width = 0.3, colnames_width = 0.3,
                      heatpal = NULL, border = NULL, title = NULL, colors = NULL,
                      inf = F, na = 0, mat = NULL, markoob = FALSE, markshape = 4, marksize = NULL, na_col = "grey", maxchar = 35, ...){
  
  ### Data ----
  datacall <- substitute(data)
  if (is.null(title)){
    if (grepl("row|col", scale, ignore.case = TRUE)){
      title <- "z-score"
    } else {
      title <- deparse1(datacall)
    }
  }
  if (title == FALSE) title <- " "
  
  heatdata <- eval(datacall, envir = parent.frame())
  heatdata <- data.matrix(heatdata)
  heatdata <- matScale(heatdata, rows = grepl("row", scale, ignore.case = TRUE), cols = grepl("col", scale, ignore.case = TRUE))
  
  ### Clustering ----
  clust <- clusterData(heatdata, rows = cluster_rows, cols = cluster_cols, inf = inf, na = na)
  if (is.null(clust$rows)) clust$rows <- FALSE
  if (is.null(clust$cols)) clust$cols <- FALSE
  
  ### Colors ----
  
  # scale colors
  if (is.null(heatpal)){
    heatpal_colors <- getColorScale(heatdata)
    heatpal <- circlize::colorRamp2(breaks = heatpal_colors, colors = names(heatpal_colors))
  }
  
  # annotation colors
  docol <- setdiff(unlist(lapply(list(coldf, rowdf), colnames)), names(colors))
  addcol <- NULL
  if (length(docol) > 0){
    if (!is.null(coldf)) all <- coldf
    if (!is.null(rowdf)) all <- rowdf
    if (!is.null(coldf) & !is.null(rowdf)) all <- dplyr::full_join(coldf, rowdf, by = character())
    addcol <- getColors(all[,docol,drop = FALSE])
  }
  colors <- c(colors, addcol)
  colors <- lapply(colors, function(tmp) tmp[!is.na(tmp) & !is.na(names(tmp))])
  
  # cell border
  if (is.null(border)){
    if (nrow(heatdata) < 100 & ncol(heatdata) < 100){
      border <- grid::gpar(col = rgb(1,1,1), lwd = grid::unit(1, "pt"))
    } else {
      border <- grid::gpar(col = NA)
    }
  } else {
    if (any(border == TRUE)){
      border <- grid::gpar(col = rgb(1,1,1), lwd = grid::unit(1, "pt"))
    } else if (length(border) > 1){
      border <- grid::gpar(col = border[is.na(as.numeric(border))], lwd = grid::unit(as.numeric(border)[!is.na(as.numeric(border))], "pt"))
    } else {
      border <- grid::gpar(col = NA)
    }
  }
  
  ### Annotations ----
  
  # Legends
  legend_params <- list(title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                        legend_height = grid::unit(0.2, "npc"),
                        border = legend_border,
                        labels_gp = grid::gpar(fontsize = fontsize))
  
  # Row annotation
  rowAnn <- NULL
  if (!is.null(rowdf)) rowAnn <- getCXanno(df = rowdf[rownames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           anno_border = anno_border,
                                           side = rowdf_side,
                                           legend = rowdf_legend,
                                           legend_params = legend_params)
  
  # Column annotation
  colAnn <- NULL
  if (!is.null(coldf)) colAnn <- getCXanno(coldf[colnames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           anno_border = anno_border,
                                           side = coldf_side,
                                           legend = coldf_legend,
                                           legend_params = legend_params)
  
  
  # Cell annotation
  if (markoob == TRUE & is.null(mat)){
    mat <- matrix(data = FALSE, nrow = nrow(heatdata), ncol = ncol(heatdata), dimnames = dimnames(heatdata))
    mat[heatdata < min(heatpal_colors)] <- TRUE
    mat[heatdata > max(heatpal_colors)] <- TRUE
  }
  
  if (is.null(marksize)) marksize <- fontsize * 0.6
  cellFUN <- NULL
  if (!is.null(mat)){
    if (!is.logical(mat)) stop("'Mat' must be a logical indicator of whether cells should be marked!")
    cellmat <- mat[rownames(heatdata), colnames(heatdata)]
    cellFUN <- function(j, i, x, y, width, height, fill){
      if (naf(cellmat[i,j] == TRUE)){ grid::grid.points(x, y, pch = markshape, size = unit(marksize, "pt")) }
    }
  }
  
  ### Heatmap ----
  dimnames(heatdata) <- lapply(dimnames(heatdata), function(x) cutstr(x, maxchar = maxchar))
  
  if (is.null(rowcex) & nrow(heatdata) > 10*ncol(heatdata)) rowcex <- 1/log10(nrow(heatdata))
  if (is.null(rowcex)) rowcex <- 1
  if (is.null(colcex)) colcex <- 1
  
  hm <- ComplexHeatmap::Heatmap(name = title,
                                matrix = heatdata,
                                row_names_max_width = grid::unit(rownames_width, "npc"),
                                column_names_max_height = grid::unit(colnames_width, "npc"),
                                column_title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                                rect_gp = border,
                                na_col = na_col,
                                left_annotation = rowAnn,
                                top_annotation = colAnn,
                                row_names_gp = grid::gpar(fontsize = fontsize * rowcex),
                                column_names_gp = grid::gpar(fontsize = fontsize * colcex),
                                col = heatpal,
                                heatmap_legend_param = legend_params,
                                cluster_rows = clust$rows,
                                cluster_columns = clust$cols,
                                cell_fun = cellFUN,
                                ...)
  
  hm
}

getCXanno <- function(df = NULL, side = "top", colors = NULL, fontsize = 12, gap = 0, legend = FALSE, anno_border = NULL, legend_params = NULL, ...){
  
  if (is.null(df)) return(list())
  
  if (is.null(anno_border)){
    gp <- grid::gpar()
  } else {
    gp <- grid::gpar(col = anno_border)
  }
  
  args <- list(df = df,
               name = side,
               annotation_label = colnames(df),
               which = ifelse(side %in% c("top", "bottom"), "column", "row"),
               col = colors[names(colors) %in% colnames(df)],
               show_annotation_name = TRUE,
               gap = grid::unit(gap, "cm"),
               border = FALSE,
               gp = gp,
               annotation_name_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
               simple_anno_size_adjust = TRUE,
               annotation_name_side = ifelse(side %in% c("top", "bottom"), "right", "top"),
               annotation_legend_param = legend_params,
               show_legend = legend,
               ...)
  
  do.call(ComplexHeatmap::HeatmapAnnotation, args)
}

clusterData <- function(data, method = "hclust", rows = NULL, cols = NULL, inf = NULL, na = NULL, ...){
  
  # input arguments
  if (is.null(rows)) rows <- nrow(data) < 1000
  if (is.null(cols)) cols <- ncol(data) < 1000
  if (is.null( na))  na <- ""
  if (is.null(inf)) inf <- ""
  
  # +/-Inf value handling
  if (any(!is.finite(nat(data)))){
    if (naf(inf == FALSE)) data <- subInf(data)
    if (is.na(inf)) data[!is.finite(data)] <- inf
  }
  
  # NA value handling
  if (any(is.na(data))){
    if (is.numeric(na)) data[is.na(data)] <- na
  }
  
  # clustering
  tmp <- list(rows = data, cols = t(data))
  res <- lapply(tmp[c(rows, cols)], function(data){
    # NA value handling
    if (na == "omit") data <- na.omit(data)
    clust <- NULL
    if (method == "hclust") clust <- dendsort::dendsort(stats::hclust(stats::dist(data)))
    clust
  })
  res
}


theme_basic <- function(fontsize = 18, fontfamily = "", lwd = NULL, color = "black", legend.ticks.color = "black", grid.color = NULL,
                        title.face = "bold", title.cex = 1, axis.title.cex = 1, axis.text.cex = 0.8, rot.x = FALSE, rot.y = FALSE, grid.cex = 0.7, ...){
  
  if (is.null(lwd)) lwd <- fontsize/22
  
  th1 <- ggplot2::theme_bw(base_size = fontsize,
                           base_family = fontfamily,
                           base_line_size = lwd,
                           base_rect_size = lwd)
  
  th2 <- ggplot2::theme(line = ggplot2::element_line(colour = color, linewidth = lwd),
                        rect = ggplot2::element_rect(colour = color, fill = NA, linewidth = lwd),
                        text = ggplot2::element_text(colour = color, size = fontsize, family = fontfamily),
                        title = ggplot2::element_text(colour = color, size = fontsize, family = fontfamily),
                        axis.text = ggplot2::element_text(colour = color, size = ggplot2::rel(axis.text.cex)),
                        axis.title = ggplot2::element_text(colour = color, size = ggplot2::rel(axis.title.cex)),
                        axis.ticks = ggplot2::element_line(colour = color),
                        axis.line = ggplot2::element_line(colour = color, lineend = "square"),
                        legend.text = ggplot2::element_text(colour = color),
                        legend.title = ggplot2::element_text(colour = color, hjust = 0),
                        legend.background = ggplot2::element_blank(),
                        legend.ticks = ggplot2::element_line(linewidth = lwd/2, colour = legend.ticks.color),
                        legend.frame = ggplot2::element_rect(linewidth = lwd/2, colour = legend.ticks.color),
                        panel.border = ggplot2::element_blank(),
                        panel.grid = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(colour = color, face = title.face, hjust = 0.5, size = ggplot2::rel(title.cex)),
                        strip.background = ggplot2::element_blank(),
                        strip.text = ggplot2::element_text(colour = color, size = fontsize, family = fontfamily),
                        strip.clip = "off",
                        ...)
  
  if (!is.null(grid.color)){
    th2 <- ggplot2::`%+replace%`(th2, ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = ggplot2::rel(grid.cex), colour = grid.color)))
  }
  
  if (rot.x){
    th2 <- ggplot2::`%+replace%`(th2, ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = 1, vjust = 0.5, angle = 90)))
  }
  
  if (rot.y){
    th2 <- ggplot2::`%+replace%`(th2, ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0.5, vjust = 1, angle = 90)))
  }
  
  th <- ggplot2::`%+replace%`(th1, th2)
  th
}



saveplot <- function(p, file = NULL, dev = "png", width = 3000, height = 2500, dpi = 300, units = "px", bg = "white", ggsave = TRUE, PRINTFUN = print, ...){
  stopifnot(dev %in% c("png", "pdf", "tiff", "svg", "jpeg", "jpg", "cairo_ps", "ps"))
  
  # filename
  if (is.null(file)) file <- rlang::as_name(rlang::enquo(p))
  if (length(file) > 1){ file <- do.call(file.path, as.list(file)) }
  if (nat(baseext(file) != dev)) file <- paste0(file, ".", dev)
  
  get_inches <- function(x, dpi, units){
    if (units == "px"){
      x <- x/dpi
    }
    if (units == "cm"){
      x <- x/dpi / 2.54
    }
    x
  }
  
  # GGSAVE WRAPPER
  # only for single ggplots
  if ("gg" %in% class(p) & ggsave == TRUE & !(dev %in% c("pdf", "ps", "cairo_ps")) ){
    ggplot2::ggsave(filename = file,
                    plot = p,
                    device = dev,
                    width = width,
                    height = height,
                    units = units,
                    dpi = dpi,
                    bg = bg,
                    ...)
    
  } else {
    
    # PNG
    if (dev == "png"){
      if (is.null(bg)) bg <- "white"
      grDevices::png(filename = file,
                     width = width,
                     height = height,
                     units = units,
                     bg = bg,
                     res = dpi,
                     type = "cairo",
                     ...)
      PRINTFUN(p)
      dev.off()
    }
    
    # PDF
    if (dev == "pdf"){
      # multi-page pdf if list of plots (with lenght > 1)
      if (!"list" %in% class(p)) p <- list(p)
      grDevices::pdf(file = file,
                     width = get_inches(width, dpi = dpi, units = units),
                     height = get_inches(height, dpi = dpi, units = units),
                     bg = bg,
                     onefile = TRUE,
                     ...)
      tmp <- lapply(p, PRINTFUN)
      dev.off()
    }
    
    # SVG
    if (dev == "svg"){
      grDevices::svg(filename = file,
                     width = get_inches(width, dpi = dpi, units = units),
                     height = get_inches(height, dpi = dpi, units = units),
                     onefile = TRUE,
                     ...)
      PRINTFUN(p)
      dev.off()
    }
    
    # TIFF
    if (dev == "tiff"){
      grDevices::tiff(filename = file,
                      width = width,
                      height = height,
                      res = dpi,
                      bg = bg,
                      type = "cairo",
                      units = units,
                      ...)
      PRINTFUN(p)
      dev.off()
    }
    
    # JPG
    if (dev %in% c("jpeg", "jpg")){
      grDevices::jpeg(filename = file,
                      width = width,
                      height = height,
                      res = dpi,
                      bg = bg,
                      type = "cairo",
                      units = units,
                      ...)
      PRINTFUN(p)
      dev.off()
    }
    
    # POSTSCRIPT
    if (dev %in% c("ps", "cairo_ps")){
      grDevices::cairo_ps(filename = file,
                          width = get_inches(width, dpi = dpi, units = units),
                          height = get_inches(height, dpi = dpi, units = units),
                          fallback_resolution = dpi,
                          bg = bg,
                          ...)
      PRINTFUN(p)
      dev.off()
    }
  }
  
  invisible(p)
}


### Utility functions

download <- function(file, url, overwrite=FALSE, method = "wget"){
  if (!file.exists(file) | overwrite == TRUE){
    download.file(url = url, destfile = file, method = method)
  }
}

padjust <- function(p, method = "fdr", ...){
  porig <- p
  is.mat <- !is.null(dim(porig))
  if (is.mat) p <- as.vector(data.matrix(p))
  padj <- stats::p.adjust(p, method = method, ...)
  if (is.mat) padj <- matrix(padj, nrow = nrow(porig), dimnames = dimnames(porig))
  if ("data.frame" %in% class(porig)) padj <- as.data.frame(padj)
  padj
}

matScale <- function(data, rows = FALSE, cols = FALSE, FUN = scale, ...){
  
  data.org <- data
  ix <- sapply(as.data.frame(data), is.numeric)
  data <- data[,ix, drop = FALSE]
  names.org <- dimnames(data)
  
  if (rows == TRUE & cols == TRUE){
    d <- dim(data)
    v <- as.vector(data)
    s <- as.numeric(FUN(v, ...))
    data <- matrix(s, nrow = d[1], ncol = d[2])
    dimnames(data) <- names.org
    
  } else if (rows == TRUE){
    data <- t(apply(data, 1, function(tmp) as.numeric(FUN(tmp, ...))))
    
  } else if (cols == TRUE){
    data <- apply(data, 2, function(tmp) as.numeric(FUN(tmp, ...)))
  }
  
  dimnames(data) <- names.org
  data <- cbind(data, data.org[,!ix, drop = FALSE])
  if (!is.null(colnames(data.org))) data <- data[,colnames(data.org)]
  
  stopifnot( all.equal(dim(data.org), dim(data)) )
  return(data)
}

get_var_type <- function(data, var){
  
  x <- data[[var]]
  n_unique <- length(unique(na.omit(x)))
  
  # determine type
  type <- if (inherits(x, "Surv")) {
    "survival"
  } else if (is.numeric(x)) {
    "numeric"
  } else if (is.factor(x)) {
    if (is.ordered(x)) {
      "ordinal"
    } else {
      "nominal"
    }
  } else if (is.character(x)) {
    "nominal"
  } else {
    "other"
  }
  
  return(list(
    type = type,
    n = n_unique
  ))
}

var_as_factor <- function(df, var, sort_n=FALSE, add_n=FALSE, na.omit=FALSE, levels=NULL, na.str="NA"){
  
  # char
  xorig <- df[[var]]
  is_num <- is.numeric(df[[var]]) | mean(!is.na(as.numeric(df[[var]])), na.rm = TRUE) > 0.5
  if (is_num){
    df[[var]] <- as.character(df[[var]])
  }
  
  # NA
  if (na.omit == TRUE){
    df <- df[!is.na(df[[var]]),,drop=FALSE]
  } else {
    df[[var]][is.na(df[[var]])] <- na.str
  }
  
  # factor
  f_levels <- unique(df[[var]])
  tab_var <- table(df[[var]]) |> sort(decreasing = TRUE)
  
  if (is_num){
    f_levels <- sort(unique(xorig), decreasing = FALSE) |> as.character()
    if (any(is.na(xorig)) & !na.omit){
      f_levels <- c(f_levels, na.str)
    }
    
  } else if (sort_n == TRUE){
    f_levels <- tab_var |> names()
  }
  
  if (add_n){
    for (f in f_levels){
      df[[var]][df[[var]] == f] <- paste0(f, " (n=", tab_var[f], ")")
    }
    f_levels <- paste0(f_levels, " (n=", tab_var[f_levels], ")")
  }
  
  if (!is.null(levels)){
    f_levels <- levels
  }
  
  df[[var]] <- factor(df[[var]], ordered = TRUE, levels = f_levels)
  df
}

writeTables <- function(data, file, rowNames = TRUE, adjwidths = TRUE,
                        header = "bold",
                        scale_styles = NULL,
                        condition_styles = NULL,
                        highlight_styles = NULL,
                        check = FALSE, ...){
  
  stopifnot(requireNamespace("openxlsx"))
  
  # helper functions
  .get_scale_style <- function(style){
    if (is.null(names(style))){
      colors <- style
      limits <- NULL
    } else {
      limits <- sort(style)
      colors <- names(limits)
      
    }
    list(style = colors, rule = limits)
  }
  
  .get_condition_style <- function(style){
    style <- as.list(style)
    rule <- style$rule
    style$rule <- NULL
    style <- do.call(openxlsx::createStyle, style)
    list(style = style, rule = rule)
  }
  
  .get_highlight_style <- function(style){
    values <- style[["values"]]
    style["values"] <- NULL
    style <- do.call(openxlsx::createStyle, style)
    list(style = style, values = values)
  }
  
  if (!"list" %in% class(data)){
    data <- setNames(list(data), gsub("\\..*$", "", basename(file)))
  }
  
  # adjust sheet names
  newnames <- cutstr(names(data), maxchar = 29)
  if (any(duplicated(newnames))){
    newnames <- cutstr(names(data), maxchar = 26)
    newnames <- dedupl(newnames)
  }
  names(data) <- newnames
  
  # define styles
  if (!is.null(header)){
    header_style <- openxlsx::createStyle(textDecoration = header)
  }
  
  # create workbook
  wb <- openxlsx::createWorkbook()
  
  # add sheets
  table <- invisible(lapply(names(data), function(tmpname){
    
    tmpdata <- as.data.frame(data[[tmpname]])
    openxlsx::addWorksheet(wb, tmpname)
    
    # write data to sheet
    openxlsx::writeData(wb, sheet = tmpname, x = tmpdata, rowNames = rowNames, headerStyle = header_style, ...)
    
    if (adjwidths == TRUE){
      openxlsx::setColWidths(wb, sheet = tmpname, cols = 1:(ncol(tmpdata) + as.numeric(rowNames)), widths = "auto")
    }
    
    # apply styles: color scales
    if (length(scale_styles) > 0){
      for (col in names(scale_styles)){
        if (!col %in% colnames(tmpdata)) next
        col_style <- .get_scale_style(scale_styles[[col]])
        openxlsx::conditionalFormatting(wb,
                                        sheet = tmpname,
                                        cols = which(colnames(tmpdata) %in% col) + as.numeric(rowNames),
                                        rows = 1 + 1:nrow(tmpdata), # all rows
                                        rule = col_style$rule,
                                        style = col_style$style,
                                        type = "colourScale")
      }
    }
    
    # apply styles: conditional highlighting
    if (length(condition_styles) > 0){
      for (col in names(condition_styles)){
        if (!col %in% colnames(tmpdata)) next
        col_style <- .get_condition_style(condition_styles[[col]])
        openxlsx::conditionalFormatting(wb,
                                        sheet = tmpname,
                                        cols = which(colnames(tmpdata) %in% col) + as.numeric(rowNames),
                                        rows = 1 + 1:nrow(tmpdata), # all rows
                                        rule = col_style$rule,
                                        style = col_style$style)
      }
    }
    
    # apply styles: value highlighting
    if (length(highlight_styles) > 0){
      for (col in names(highlight_styles)){
        if (!col %in% colnames(tmpdata)) next
        col_style <- .get_highlight_style(highlight_styles[[col]])
        openxlsx::addStyle(wb,
                           sheet = tmpname,
                           cols = which(colnames(tmpdata) %in% col) + as.numeric(rowNames),
                           rows = 1 + which(tmpdata[[col]] %in% col_style$values), # select rows
                           style = col_style$style)
      }
    }
    
    tmpdata
  }))
  
  # save file
  if (baseext(file) != "xlsx") file <- paste0(file, ".xlsx")
  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
  
  # check if writing to excel changed any of the values
  if (check == TRUE){
    table_check <- readTables(file, rowNames = rowNames)
    if (is.data.frame(table_check)) table_check <- list(table_check)
    for (i in 1:length(table_check)){
      print(all.equal(table[[i]], table_check[[i]]))
    }
  }
}

getLimits <- function(x, clip = TRUE, expand = 0.1, negative = TRUE){
  
  x <- x[!is.na(x)]
  x <- x + x*expand
  
  if (clip == TRUE){
    h <- hist(x, plot = FALSE, breaks = 30)
    xd <- (h$counts > 3) | (rev(cumsum(rev(h$counts))) > 8)
    xmin <- h$breaks[which(xd)[1]]
    xmax <- rev(h$breaks)[which(rev(xd))[1]]
    
  } else {
    xmin <- NA
    xmax <- NA
  }
  
  if (is.na(xmax)) xmax <- max(x) %>% roundup(., roundup(-log10(abs(.)))) # upper
  if (is.na(xmin)) xmin <- min(x) %>% rounddown(., roundup(-log10(abs(.)))) # lower
  
  if (is.na(xmin)){xmin <- -0.1 * xmax}
  if (is.na(xmax)){xmax <- -0.1 * xmin}
  
  if (is.na(xmax) & is.na(xmin)){
    xmin <- -1
    xmax <- 1
  }
  
  res <- c("min" = xmin, "max" = xmax)
  if (negative == FALSE) res[res < 0] <- 0
  
  res
}

get_all_comparisons <- function(df, x, y, as_list=FALSE, na.omit = FALSE){
  
  groups <- unique(df[[x]]) |> as.character()
  if (na.omit == TRUE){
    groups <- groups[!is.na(groups)]
  }
  groups[is.na(groups)] <- "NA"
  statdf <- expand.grid(groups, groups) |> apply(1, sort) |> t()
  l <- statdf |> apply(1, unique) |> sapply(length)
  statdf <- unique(statdf[l == 2,]) |> data.frame()
  colnames(statdf) <- c("group1","group2")
  
  if (as_list==TRUE){
    statdf <- as.data.frame(t(statdf[,c(2,1)])) |> as.list()
    names(statdf) <- sapply(statdf, paste, collapse="_vs_")
  }
  
  statdf
}

getColorScale <- function(data, ...){
  data <- data.matrix(data)
  data[!is.finite(data)] <- NA
  lims <- range(data, na.rm = TRUE)
  if (min(data, na.rm = TRUE) >= 0){ # positive only
    col <- setNames(c(0, lims[2]), c("white", "red"))
  } else {
    col <- setNames(c(lims[1], 0, lims[2]), c("blue", "white", "red"))
  }
  col
}

get_colors <- function(n, type = "discrete"){
  
  if (type == "discrete"){
    if (n < 3){
      colors <- c("#593eb8", "#ffbb00")
    } else if (n <= length(colorblindr::palette_OkabeIto)){
      colors <- paletteer::paletteer_d("colorblindr::OkabeIto", n) |> as.character()
    } else {
      paldf <- subset(paletteer::palettes_d_names, length <= n & type == "qualitative") |> dplyr::arrange(desc(length))
      colors <- paletteer::paletteer_d(paste0(paldf$package[1], "::", paldf$palette[1]), n) |> as.character()
    }
  }
  
  colors
}

get_mhc_genes <- function(){
  stopifnot(requireNamespace("org.Hs.eg.db"))
  stopifnot(requireNamespace("AnnotationDbi"))
  ids <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
  annot <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ids, column = "GENENAME", keytype = "SYMBOL", multiVals = "first")
  mhc <- annot[grep("histocompatibility", annot, ignore.case = TRUE)]
  mhcI <- names(mhc[grepl("class I,", mhc, ignore.case = TRUE) & !grepl("pseudogene", mhc, ignore.case = TRUE)]) |> grep(pattern = "HLA-", value = TRUE)
  mhcII <- names(mhc[grepl("class II", mhc, ignore.case = TRUE) & !grepl("pseudogene", mhc, ignore.case = TRUE)]) |> grep(pattern = "HLA-", value = TRUE)
  list("classI" = mhcI, "classII" = mhcII)
}

nat <- function (data, ...){
  data[is.na(data)] <- TRUE
  data
}

naf <- function (data, ...){
  data[is.na(data)] <- FALSE
  data
}

cutstr <- function (x, maxchar = 100, add = "...", add_incl = TRUE){
  n <- length(unique(x))
  lev <- levels(x)
  ordered <- FALSE
  if ("ordered" %in% class(x)){
    ordered <- TRUE
  }
  x <- as.character(x)
  naix <- is.na(x)
  x[naix] <- "NA-tmp"
  x <- as.character(x)
  ix <- nchar(x) > maxchar
  x[ix] <- substr(x[ix], 1, ifelse(add_incl, maxchar - nchar(add), maxchar))
  x[ix] <- paste0(x[ix], add)
  x[naix] <- NA
  
  if (!is.null(lev)){
    ix <- nchar(lev) > maxchar
    lev[ix] <- substr(lev[ix], 1, ifelse(add_incl, maxchar - nchar(add), maxchar))
    lev[ix] <- paste0(lev[ix], add)
    x <- factor(x, levels = lev, ordered = ordered)
  }
  
  stopifnot(length(unique(x)) == n)
  x
}

col2rownames <- function(data, col = id, sep = "_"){
  
  col <- rlang::enquo(col)
  if (!rlang::as_name(col) %in% colnames(data)){
    warning(paste0("Warning: Column ", rlang::as_name(col), " not found!"))
    return(data)
  }
  
  names <- dplyr::select(data, !!col)
  names <- apply(names, 1, paste0, collapse = sep)
  rownames(data) <- names
  data <- dplyr::select(data, -!!col)
  data
}

rownames2col <- function(data, col = id, keep = FALSE){
  
  col <- rlang::enquo(col)
  names <- rownames(data)
  if (is.null(names)){
    warning("Warning: No rownames found.")
    return(data)
  }
  if (keep == FALSE) rownames(data) <- NULL
  
  i <- rlang::as_name(col)
  if (i %in% colnames(data)){
    data[,i] <- names
  } else {
    if (class(names) %in% unique(sapply(data, class))){
      data <- cbind(names, data)
    } else {
      data <- data.frame(names, data)
    }
    
    colnames(data)[1] <- i
  }
  
  data
}

fp <- function(...){
  args <- list(...)
  args <- lapply(args, paste0, collapse = "")
  do.call(file.path, args)
}

baseext <- function(path, ...){
  path <- basename(path)
  ext <- gsub(x = path, pattern = ".*\\.", replacement = "")
  ext[!grepl(pattern = ".", x = path, fixed = TRUE)] <- ""
  ext
}

ce <- function(...){
  args <- rlang::enquos(...)
  sapply(args, rlang::as_name)
}

`%L>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe)
  return(res)
}

`%S>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe, simplify = TRUE)
  return(res)
}

listpipe <- function(lhs, rhs, env, simplify = FALSE, qsub = FALSE, ...){
  
  ### LHS (list) ----
  lhs <- eval(lhs, envir = env)
  if (is.null(names(lhs)) & (is.character(lhs) | is.integer(lhs))) names(lhs) <- as.character(lhs)
  
  ### RHS (function) ----
  FUN <- rlang::call_standardise(rhs)
  args <- rlang::call_args(FUN)
  
  anon_f <- any(unlist(sapply(args, grepl, pattern = "function")))
  if (anon_f){
    FUN <- args[[length(args)]]
    FUN <- as.call(parse(text = as.character(FUN)))
    args <- unlist(args[-(length(args)-0:1)], recursive = FALSE)
  }
  
  if (length(args) > 0){
    names(args)[names(args)==""] <- NA
    dix <- which(sapply(args, function(tmp) all(tmp == ".") ))
    
    dot_args <- names(args)[dix]
    dot_args <- dot_args[!is.na(dot_args)]
    if (length(dot_args) == 0) dot_args <- NULL
    
    drm <- dix[is.na(names(args)[dix]) | is.null(names(args)[dix])]
    if (any(drm)) args <- args[-drm]
    if (length(args) == 0) args <- NULL
    
  } else dot_args <- NULL
  
  args <- args[!is.na(names(args))]
  rm_args <- lapply(setNames(names(args), names(args)), function(tmp) rlang::zap() )
  FUN <- rlang::call_modify(FUN, ... = rlang::zap())
  FUN <- rlang::call_modify(FUN, !!!rm_args)
  
  ftxt <- as.character(FUN)
  ftxt <- ftxt[ftxt != "."]
  FUN <- as.call(parse(text = ftxt))
  
  if (length(dot_args) != 0) args[dot_args] <- NULL
  
  ### Iterate function over list ----
  itfun <- function(tmp, dot_args, args, FUN, env){
    tmp_args <- list(tmp)
    if (!is.null(dot_args)) tmp_args <- rep(tmp_args, length(dot_args))
    all_args <- c(setNames(tmp_args, dot_args), args)
    FUN <- rlang::call_modify(FUN, !!!all_args)
    eval(FUN, envir = env)
  }
  
if (simplify == TRUE){
    res <- sapply(setNames(seq_along(lhs), names(lhs)), function(i){
      assign(x = ".name", value = names(lhs)[i], envir = env)
      itfun(lhs[[i]], dot_args, args, FUN, env)
    }, ...)
  } else {
    res <- lapply(setNames(seq_along(lhs), names(lhs)), function(i){
      assign(x = ".name", value = names(lhs)[i], envir = env)
      itfun(lhs[[i]], dot_args, args, FUN, env)
    }, ...)
    
  }
  
  suppressWarnings(rm(".name", envir = env))
  return(res)
}

