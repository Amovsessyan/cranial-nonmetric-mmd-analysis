#####################################################################
## Cranial non-metric traits, MMD, Mantel test,
## subset screening, Honest permutation procedures,
## sequential expansion of the Honest core, and leave-one-out analysis.
##
## Repository version:
## - input data are read from CSV files in data/
## - no embedded data tables
## - no figures
## - numerical results are saved as CSV files in results/
##
## Required input files:
##   data/ural_trait_frequencies.csv
##   data/siberia_trait_frequencies.csv
##   data/ural_fst_matrix.csv
##   data/siberia_fst_matrix.csv
##
## Trait-frequency CSV files must contain:
##   population, N, and 30 cranial non-metric trait columns.
##
## FST CSV files must contain:
##   population as the first column, followed by a square FST matrix
##   with matching population names.
#####################################################################

suppressPackageStartupMessages({
  library(Matrix)
})

set.seed(42)

# =========================
# SETTINGS
# =========================

angular_transform   <- "Freeman"
truncate_negative   <- FALSE
k_max               <- 12

nr_perm_mantel_all  <- 4999
nr_perm_mantel_best <- 999
nr_perm_pipeline    <- 499
nr_perm_loo         <- 1999

max_comb_per_k      <- 200000
max_comb_pipeline   <- 10000

chunk_cols          <- 20000

# Population order used for aligning morphology and FST matrices.
# This prevents accidental row-order mismatches in the input CSV files.
ural_order <- c(
  "Finns",
  "Karelians",
  "Mordovians",
  "Maris",
  "Komis",
  "Udmurts",
  "Mansis",
  "Khanty",
  "Selkups"
)

siberia_order <- c(
  "Mongolian",
  "Aleutian",
  "Buryat",
  "Chukchi",
  "Evenki",
  "Khanty",
  "Naukan",
  "Selkup",
  "Tuvinians",
  "Yakut"
)

# =========================
# PATHS
# =========================

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

data_dir <- file.path(project_root, "data")
out_dir  <- file.path(project_root, "results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ural_traits_file    <- file.path(data_dir, "ural_trait_frequencies.csv")
siberia_traits_file <- file.path(data_dir, "siberia_trait_frequencies.csv")
ural_fst_file       <- file.path(data_dir, "ural_fst_matrix.csv")
siberia_fst_file    <- file.path(data_dir, "siberia_fst_matrix.csv")

required_files <- c(
  ural_traits_file,
  siberia_traits_file,
  ural_fst_file,
  siberia_fst_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input files are missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Current project folder: ", project_root)
message("Data dir:               ", data_dir)
message("Output dir:             ", out_dir)

# =========================
# HELPERS: READING AND PARSING
# =========================

parse_num <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("NA", "NaN", "", "NULL")] <- NA_character_
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

trim_all <- function(x) {
  trimws(as.character(x))
}

read_morph_csv <- function(file) {
  df <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  
  colnames(df) <- trim_all(colnames(df))
  
  if (!all(c("population", "N") %in% colnames(df))) {
    stop("Morphology file must contain columns named 'population' and 'N': ", file)
  }
  
  pops <- trim_all(df$population)
  
  n_vec <- parse_num(df$N)
  names(n_vec) <- pops
  
  trait_names <- setdiff(colnames(df), c("population", "N"))
  
  P <- as.matrix(df[, trait_names, drop = FALSE])
  P <- apply(P, 2, parse_num)
  
  rownames(P) <- pops
  colnames(P) <- trait_names
  
  P <- as.matrix(P)
  
  P[P < 0.005] <- 0.005
  P[P > 0.995] <- 0.995
  
  list(P = P, n = n_vec)
}

read_fst_csv <- function(file) {
  df <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  
  colnames(df) <- trim_all(colnames(df))
  
  if (!("population" %in% colnames(df))) {
    stop("FST file must contain a first column named 'population': ", file)
  }
  
  pops <- trim_all(df$population)
  mat_cols <- setdiff(colnames(df), "population")
  
  fst <- as.matrix(df[, mat_cols, drop = FALSE])
  fst <- apply(fst, 2, parse_num)
  
  rownames(fst) <- pops
  colnames(fst) <- mat_cols
  
  fst <- as.matrix(fst)
  diag(fst) <- 0
  
  fst
}

match_region_data <- function(morph_list,
                              fst_mat,
                              expected_order,
                              expected_traits,
                              region_name) {
  missing_morph <- setdiff(expected_order, rownames(morph_list$P))
  missing_n     <- setdiff(expected_order, names(morph_list$n))
  missing_fst_r <- setdiff(expected_order, rownames(fst_mat))
  missing_fst_c <- setdiff(expected_order, colnames(fst_mat))
  
  if (length(missing_morph) > 0) {
    stop(region_name, ": missing morphology rows: ", paste(missing_morph, collapse = ", "))
  }
  
  if (length(missing_n) > 0) {
    stop(region_name, ": missing sample sizes for: ", paste(missing_n, collapse = ", "))
  }
  
  if (length(missing_fst_r) > 0) {
    stop(region_name, ": missing FST rows: ", paste(missing_fst_r, collapse = ", "))
  }
  
  if (length(missing_fst_c) > 0) {
    stop(region_name, ": missing FST columns: ", paste(missing_fst_c, collapse = ", "))
  }
  
  P <- morph_list$P[expected_order, , drop = FALSE]
  n_vec <- morph_list$n[expected_order]
  fst <- fst_mat[expected_order, expected_order, drop = FALSE]
  
  if (ncol(P) != expected_traits) {
    stop(region_name, ": expected ", expected_traits, " traits, got ", ncol(P))
  }
  
  if (!all(rownames(P) == rownames(fst))) {
    stop(region_name, ": population order mismatch between morphology and FST rows.")
  }
  
  if (!all(rownames(P) == colnames(fst))) {
    stop(region_name, ": population order mismatch between morphology and FST columns.")
  }
  
  list(P = P, n = n_vec, fst = fst)
}

# =========================
# CORE FUNCTIONS
# =========================

upper_vec <- function(M) {
  M[upper.tri(M)]
}

theta_matrix <- function(P, n_vec, angular = "Freeman") {
  if (angular == "Anscombe") {
    return(asin(1 - 2 * P))
  }
  
  K <- round(sweep(P, 1, n_vec, `*`))
  denom <- n_vec + 1
  
  A <- asin(sqrt(sweep(K, 1, denom, `/`)))
  B <- asin(sqrt(sweep(K + 1, 1, denom, `/`)))
  
  0.5 * (A + B)
}

mmd_matrix_fast <- function(P,
                            n_vec,
                            cols = NULL,
                            angular = "Freeman",
                            truncate = FALSE) {
  if (!is.null(cols)) {
    P <- P[, cols, drop = FALSE]
  }
  
  Tn <- ncol(P)
  
  if (Tn < 1) {
    stop("Need at least one trait.")
  }
  
  th <- theta_matrix(P, n_vec, angular = angular)
  
  norms <- rowSums(th^2)
  G <- th %*% t(th)
  
  sqdist_mean <- (outer(norms, norms, `+`) - 2 * G) / Tn
  
  v <- 1 / (n_vec + 0.5)
  var_mat <- outer(v, v, `+`)
  
  M <- sqdist_mean - var_mat
  diag(M) <- 0
  
  if (truncate) {
    M[M < 0] <- 0
  }
  
  M
}

mmd_dist_vec <- function(P,
                         n_vec,
                         cols = NULL,
                         angular = "Freeman",
                         truncate = FALSE) {
  upper_vec(mmd_matrix_fast(P, n_vec, cols, angular, truncate))
}

corr_mmd_fst <- function(P,
                         n_vec,
                         fst_upper,
                         cols = NULL,
                         angular = "Freeman",
                         truncate = FALSE) {
  d <- mmd_dist_vec(P, n_vec, cols, angular, truncate)
  
  keep <- is.finite(d) & is.finite(fst_upper)
  
  if (sum(keep) < 3) return(NA_real_)
  if (sd(d[keep]) == 0 || sd(fst_upper[keep]) == 0) return(NA_real_)
  
  cor(d[keep], fst_upper[keep], method = "pearson")
}

mantel_perm_mmd <- function(P,
                            n_vec,
                            fst_mat,
                            cols = NULL,
                            angular = "Freeman",
                            truncate = FALSE,
                            perms = 999,
                            show_progress = TRUE) {
  n <- nrow(P)
  fst_upper <- upper_vec(fst_mat)
  
  r_obs <- corr_mmd_fst(P, n_vec, fst_upper, cols, angular, truncate)
  
  if (!is.finite(r_obs)) {
    return(list(r = NA_real_, p = NA_real_))
  }
  
  more <- 0L
  
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = perms, style = 3)
  }
  
  for (b in seq_len(perms)) {
    perm <- sample.int(n)
    
    Pp <- P[perm, , drop = FALSE]
    rownames(Pp) <- rownames(P)
    
    np <- n_vec[perm]
    names(np) <- names(n_vec)
    
    r_perm <- corr_mmd_fst(Pp, np, fst_upper, cols, angular, truncate)
    
    if (is.finite(r_perm) && abs(r_perm) >= abs(r_obs)) {
      more <- more + 1L
    }
    
    if (show_progress) {
      setTxtProgressBar(pb, b)
    }
  }
  
  if (show_progress) {
    close(pb)
  }
  
  pval <- (more + 1) / (perms + 1)
  
  list(r = r_obs, p = pval)
}

# =========================
# LEAVE-ONE-OUT ANALYSIS
# =========================

loo_mantel <- function(P,
                       n_vec,
                       fst,
                       cols = NULL,
                       angular = "Freeman",
                       truncate = FALSE,
                       perms = 1999) {
  out <- data.frame(
    excluded = character(0),
    r = numeric(0),
    p = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (pop in rownames(P)) {
    keep <- setdiff(rownames(P), pop)
    
    P_sub <- P[keep, , drop = FALSE]
    n_sub <- n_vec[keep]
    fst_sub <- fst[keep, keep, drop = FALSE]
    
    mb <- mantel_perm_mmd(
      P = P_sub,
      n_vec = n_sub,
      fst_mat = fst_sub,
      cols = cols,
      angular = angular,
      truncate = truncate,
      perms = perms,
      show_progress = FALSE
    )
    
    out <- rbind(
      out,
      data.frame(
        excluded = pop,
        r = mb$r,
        p = mb$p,
        stringsAsFactors = FALSE
      )
    )
  }
  
  out
}

sample_combos_idx <- function(p, k, max_comb, seed = 42) {
  n_total <- suppressWarnings(choose(p, k))
  
  if (is.finite(n_total) && n_total <= max_comb) {
    return(utils::combn(seq_len(p), k))
  }
  
  set.seed(seed)
  
  mat <- replicate(
    max_comb,
    sort(sample.int(p, k)),
    simplify = "matrix"
  )
  
  uniq <- unique(t(mat))
  t(uniq)
}

build_sparse_Bk <- function(Ck, p) {
  k <- nrow(Ck)
  nsubs <- ncol(Ck)
  
  i <- as.vector(Ck)
  j <- rep(seq_len(nsubs), each = k)
  
  Matrix::sparseMatrix(
    i = i,
    j = j,
    x = 1,
    dims = c(p, nsubs),
    giveCsparse = TRUE
  )
}

best_stat_fixed_k <- function(P,
                              n_vec,
                              fst_upper,
                              Ck,
                              Bk,
                              angular = "Freeman",
                              truncate = FALSE,
                              chunk = 20000) {
  n <- nrow(P)
  k <- nrow(Ck)
  nsubs <- ncol(Ck)
  
  idx_i <- rep(seq_len(n), times = n)
  idx_j <- rep(seq_len(n), each = n)
  
  mask <- idx_i < idx_j
  
  ii <- idx_i[mask]
  jj <- idx_j[mask]
  
  th <- theta_matrix(P, n_vec, angular = angular)
  
  MD <- (th[ii, , drop = FALSE] - th[jj, , drop = FALSE])^2
  
  var_pair <- 1 / (n_vec[ii] + 0.5) + 1 / (n_vec[jj] + 0.5)
  
  y <- fst_upper
  yc <- y - mean(y)
  yss <- sum(yc^2)
  
  bestT <- -Inf
  best_cols <- NULL
  
  for (start in seq(1, nsubs, by = chunk)) {
    end <- min(nsubs, start + chunk - 1)
    cols_idx <- start:end
    
    SUM <- MD %*% Bk[, cols_idx, drop = FALSE]
    
    X <- SUM / k
    X <- sweep(X, 1, var_pair, "-")
    
    if (truncate) {
      X <- pmax(0, X)
    }
    
    xm <- colMeans(X)
    Xc <- sweep(X, 2, xm, "-")
    
    num <- as.numeric(crossprod(yc, Xc))
    den <- sqrt(yss * colSums(Xc^2))
    
    r <- num / den
    r[!is.finite(r)] <- NA_real_
    
    loc <- which.max(r)
    Tloc <- r[loc]
    
    if (is.finite(Tloc) && Tloc > bestT) {
      bestT <- Tloc
      best_cols <- Ck[, cols_idx[loc]]
    }
  }
  
  list(T = bestT, best_cols = best_cols)
}

build_sparse_B_list <- function(combos_by_k, p) {
  B_list <- list()
  
  for (k_name in names(combos_by_k)) {
    Ck <- combos_by_k[[k_name]]
    
    if (is.null(Ck) || ncol(Ck) == 0) next
    
    k <- as.integer(k_name)
    nsubs <- ncol(Ck)
    
    i <- as.vector(Ck)
    j <- rep(seq_len(nsubs), each = k)
    
    B_list[[k_name]] <- Matrix::sparseMatrix(
      i = i,
      j = j,
      x = 1,
      dims = c(p, nsubs),
      giveCsparse = TRUE
    )
  }
  
  B_list
}

pipeline_best_stat <- function(P,
                               n_vec,
                               fst_upper,
                               angular,
                               truncate,
                               B_list,
                               combos_by_k,
                               chunk = 20000) {
  n <- nrow(P)
  
  idx_i <- rep(seq_len(n), times = n)
  idx_j <- rep(seq_len(n), each = n)
  
  mask <- idx_i < idx_j
  
  ii <- idx_i[mask]
  jj <- idx_j[mask]
  
  th <- theta_matrix(P, n_vec, angular = angular)
  
  MD <- (th[ii, , drop = FALSE] - th[jj, , drop = FALSE])^2
  
  var_pair <- 1 / (n_vec[ii] + 0.5) + 1 / (n_vec[jj] + 0.5)
  
  y <- fst_upper
  yc <- y - mean(y)
  yss <- sum(yc^2)
  
  bestT <- -Inf
  best_k <- NA_integer_
  best_cols <- NULL
  
  for (k_name in names(B_list)) {
    Bk <- B_list[[k_name]]
    Ck <- combos_by_k[[k_name]]
    
    k <- as.integer(k_name)
    nsubs <- ncol(Bk)
    
    for (start in seq(1, nsubs, by = chunk)) {
      end <- min(nsubs, start + chunk - 1)
      cols_idx <- start:end
      
      SUM <- MD %*% Bk[, cols_idx, drop = FALSE]
      
      X <- SUM / k
      X <- sweep(X, 1, var_pair, "-")
      
      if (truncate) {
        X <- pmax(0, X)
      }
      
      xm <- colMeans(X)
      Xc <- sweep(X, 2, xm, "-")
      
      num <- as.numeric(crossprod(yc, Xc))
      den <- sqrt(yss * colSums(Xc^2))
      
      r <- num / den
      r[!is.finite(r)] <- NA_real_
      
      loc <- which.max(r)
      Tloc <- r[loc]
      
      if (is.finite(Tloc) && Tloc > bestT) {
        bestT <- Tloc
        best_k <- k
        best_cols <- Ck[, cols_idx[loc]]
      }
    }
  }
  
  list(T = bestT, best_k = best_k, best_cols = best_cols)
}

honest_pipeline_p <- function(P,
                              n_vec,
                              fst_mat,
                              combos_by_k,
                              angular,
                              truncate,
                              perms = 499,
                              seed = 42,
                              chunk = 20000) {
  set.seed(seed)
  
  fst_upper <- upper_vec(fst_mat)
  p <- ncol(P)
  
  B_list <- build_sparse_B_list(combos_by_k, p)
  
  obs <- pipeline_best_stat(
    P,
    n_vec,
    fst_upper,
    angular = angular,
    truncate = truncate,
    B_list = B_list,
    combos_by_k = combos_by_k,
    chunk = chunk
  )
  
  T_obs <- obs$T
  
  more <- 0L
  bestk_perm <- integer(perms)
  
  pb <- txtProgressBar(min = 0, max = perms, style = 3)
  
  n <- nrow(P)
  
  for (b in seq_len(perms)) {
    perm <- sample.int(n)
    
    Pp <- P[perm, , drop = FALSE]
    rownames(Pp) <- rownames(P)
    
    np <- n_vec[perm]
    names(np) <- names(n_vec)
    
    st <- pipeline_best_stat(
      Pp,
      np,
      fst_upper,
      angular = angular,
      truncate = truncate,
      B_list = B_list,
      combos_by_k = combos_by_k,
      chunk = chunk
    )
    
    if (is.finite(st$T) && st$T >= T_obs) {
      more <- more + 1L
    }
    
    bestk_perm[b] <- st$best_k
    
    setTxtProgressBar(pb, b)
  }
  
  close(pb)
  
  p_honest <- (more + 1) / (perms + 1)
  
  list(
    T_obs = T_obs,
    p_honest = p_honest,
    best_k_obs = obs$best_k,
    best_cols_obs = obs$best_cols,
    bestk_perm = bestk_perm
  )
}

# =========================
# SEQUENTIAL EXPANSION OF HONEST CORE
# =========================

run_sequential_expansion <- function(P,
                                     n_vec,
                                     fst,
                                     core_traits,
                                     k_target,
                                     traits,
                                     angular_transform,
                                     truncate_negative,
                                     max_comb_pipeline,
                                     nr_perm_pipeline,
                                     chunk_cols) {
  out <- list()
  
  base_cols <- match(core_traits, traits)
  
  mb_core <- mantel_perm_mmd(
    P,
    n_vec,
    fst,
    cols = base_cols,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_mantel_best,
    show_progress = FALSE
  )
  
  out[[1]] <- data.frame(
    k = length(core_traits),
    statistic = mb_core$r,
    p_value = mb_core$p,
    traits = paste(core_traits, collapse = "; "),
    mode = "core only (fixed set)",
    stringsAsFactors = FALSE
  )
  
  if (length(core_traits) < k_target) {
    for (m in (length(core_traits) + 1):k_target) {
      extras_needed <- m - length(core_traits)
      
      remaining <- setdiff(traits, core_traits)
      
      core_idx <- match(core_traits, traits)
      rem_idx <- match(remaining, traits)
      
      C_ex <- sample_combos_idx(
        length(rem_idx),
        extras_needed,
        max_comb_pipeline,
        seed = 20000 + m
      )
      
      C_ex_global <- matrix(
        rem_idx[C_ex],
        nrow = nrow(C_ex),
        ncol = ncol(C_ex)
      )
      
      C_full <- rbind(
        matrix(core_idx, nrow = length(core_idx), ncol = ncol(C_ex_global)),
        C_ex_global
      )
      
      combos_this <- list()
      combos_this[[as.character(m)]] <- C_full
      
      st <- honest_pipeline_p(
        P = P,
        n_vec = n_vec,
        fst_mat = fst,
        combos_by_k = combos_this,
        angular = angular_transform,
        truncate = truncate_negative,
        perms = nr_perm_pipeline,
        seed = 42,
        chunk = chunk_cols
      )
      
      out[[length(out) + 1]] <- data.frame(
        k = m,
        statistic = st$T_obs,
        p_value = st$p_honest,
        traits = paste(traits[st$best_cols_obs], collapse = "; "),
        mode = "honest core + extras",
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, out)
}

# =========================
# PER-REGION ANALYSIS
# =========================

run_region_analysis <- function(region_name, P, n_vec, fst) {
  traits <- colnames(P)
  fst_upper <- upper_vec(fst)
  p <- ncol(P)
  k_max2 <- min(k_max, p)
  
  cat("\n============================================================\n")
  cat("REGION:", region_name, "\n")
  cat("Populations:", nrow(P), "\n")
  cat("Traits:", p, "\n")
  cat("============================================================\n")
  
  cat("=== ALL-TRAITS MANTEL (MMD) ===\n")
  cat("Permutations:", nr_perm_mantel_all, "\n")
  
  mantel_all <- mantel_perm_mmd(
    P,
    n_vec,
    fst,
    cols = NULL,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_mantel_all,
    show_progress = TRUE
  )
  
  cat("\nALL traits: r =", round(mantel_all$r, 7), ", p =", mantel_all$p, "\n\n")
  
  # Screening
  cat("=== SCREENING BEST SUBSET BY k (FAST, SAMPLED) ===\n")
  
  summary_table <- data.frame(
    k = integer(0),
    n_subsets = integer(0),
    r_max = numeric(0),
    p_best = numeric(0),
    traits = character(0),
    stringsAsFactors = FALSE
  )
  
  for (k in seq_len(k_max2)) {
    Ck <- sample_combos_idx(
      p,
      k,
      max_comb_per_k,
      seed = 100 + k
    )
    
    Bk <- build_sparse_Bk(Ck, p)
    
    cat("k =", k, "/", k_max2, " | subsets =", ncol(Ck), " ... ")
    
    bestk <- best_stat_fixed_k(
      P,
      n_vec,
      fst_upper,
      Ck,
      Bk,
      angular = angular_transform,
      truncate = truncate_negative,
      chunk = chunk_cols
    )
    
    mb <- mantel_perm_mmd(
      P,
      n_vec,
      fst,
      cols = bestk$best_cols,
      angular = angular_transform,
      truncate = truncate_negative,
      perms = nr_perm_mantel_best,
      show_progress = FALSE
    )
    
    best_traits <- traits[bestk$best_cols]
    
    summary_table <- rbind(
      summary_table,
      data.frame(
        k = k,
        n_subsets = ncol(Ck),
        r_max = bestk$T,
        p_best = mb$p,
        traits = paste(best_traits, collapse = "; "),
        stringsAsFactors = FALSE
      )
    )
    
    cat("r_max =", round(bestk$T, 4), " | p_best =", mb$p, "\n")
  }
  
  best_k_screen <- summary_table$k[which.max(summary_table$r_max)]
  
  best_traits_screen <- strsplit(
    summary_table$traits[which.max(summary_table$r_max)],
    ";\\s*"
  )[[1]]
  
  # Honest max-over-k
  cat("\n=== HONEST max-over-k p-value (permute + reselect over ALL k) ===\n")
  cat("max_comb_pipeline =", max_comb_pipeline, " per k | perms =", nr_perm_pipeline, "\n\n")
  
  combos_by_k_pipe <- list()
  
  for (k in seq_len(k_max2)) {
    combos_by_k_pipe[[as.character(k)]] <- sample_combos_idx(
      p,
      k,
      max_comb_pipeline,
      seed = 5000 + k
    )
    
    cat(
      "PIPE k=", k,
      " subsets=", ncol(combos_by_k_pipe[[as.character(k)]]),
      "\n",
      sep = ""
    )
  }
  
  cat("\n")
  
  honest_allk <- honest_pipeline_p(
    P = P,
    n_vec = n_vec,
    fst_mat = fst,
    combos_by_k = combos_by_k_pipe,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_pipeline,
    seed = 42,
    chunk = chunk_cols
  )
  
  honest_k <- honest_allk$best_k_obs
  honest_cols <- honest_allk$best_cols_obs
  honest_traits <- traits[honest_cols]
  
  cat("\nHONEST max-over-k:\n")
  cat("  best k   =", honest_k, "\n")
  cat("  T_obs    =", round(honest_allk$T_obs, 4), "\n")
  cat("  p_honest =", honest_allk$p_honest, "\n")
  cat("  traits:\n")
  print(honest_traits)
  
  # Honest fixed k=12
  k_fixed <- min(12, k_max2)
  
  cat("\n=== HONEST FIXED k=", k_fixed, " (permute + reselect within k fixed) ===\n", sep = "")
  
  combos_fixed <- list()
  
  combos_fixed[[as.character(k_fixed)]] <- sample_combos_idx(
    p,
    k_fixed,
    max_comb_pipeline,
    seed = 9000 + k_fixed
  )
  
  honest_kfixed <- honest_pipeline_p(
    P = P,
    n_vec = n_vec,
    fst_mat = fst,
    combos_by_k = combos_fixed,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_pipeline,
    seed = 42,
    chunk = chunk_cols
  )
  
  traits_kfixed <- traits[honest_kfixed$best_cols_obs]
  
  cat("\nHONEST fixed k:\n")
  cat("  k        =", k_fixed, "\n")
  cat("  T_obs    =", round(honest_kfixed$T_obs, 4), "\n")
  cat("  p_honest =", honest_kfixed$p_honest, "\n")
  cat("  traits:\n")
  print(traits_kfixed)
  
  # Leave-one-out analysis
  cat("\n=== LEAVE-ONE-OUT ANALYSIS ===\n")
  
  loo_all <- loo_mantel(
    P = P,
    n_vec = n_vec,
    fst = fst,
    cols = NULL,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_loo
  )
  
  loo_honest_max <- loo_mantel(
    P = P,
    n_vec = n_vec,
    fst = fst,
    cols = honest_allk$best_cols_obs,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_loo
  )
  
  loo_honest_fixed <- loo_mantel(
    P = P,
    n_vec = n_vec,
    fst = fst,
    cols = honest_kfixed$best_cols_obs,
    angular = angular_transform,
    truncate = truncate_negative,
    perms = nr_perm_loo
  )
  
  loo_all$panel <- "All traits"
  loo_honest_max$panel <- paste0("Honest max-over-k, k=", honest_k)
  loo_honest_fixed$panel <- "Honest fixed k=12"
  
  loo_table <- rbind(
    loo_all,
    loo_honest_max,
    loo_honest_fixed
  )
  
  cat("\nLeave-one-out summary:\n")
  print(loo_table)
  
  # Sequential expansion of Honest core
  k_target <- k_fixed
  core_traits <- honest_traits
  
  cat("\n=== HONEST core+extras to k=", k_target, " ===\n", sep = "")
  
  honest_seq <- run_sequential_expansion(
    P = P,
    n_vec = n_vec,
    fst = fst,
    core_traits = core_traits,
    k_target = k_target,
    traits = traits,
    angular_transform = angular_transform,
    truncate_negative = truncate_negative,
    max_comb_pipeline = max_comb_pipeline,
    nr_perm_pipeline = nr_perm_pipeline,
    chunk_cols = chunk_cols
  )
  
  traits_core12 <- character(0)
  
  row12 <- honest_seq[honest_seq$k == k_fixed, , drop = FALSE]
  
  if (nrow(row12) > 0 && row12$mode[1] == "honest core + extras") {
    traits_core12 <- strsplit(row12$traits[1], ";\\s*")[[1]]
  }
  
  if (nrow(row12) > 0) {
    cat("\nSEQUENTIAL EXPANSION TO k=12:\n")
    print(honest_seq[, c("k", "statistic", "p_value", "mode")], row.names = FALSE)
  }
  
  # Save tables
  prefix <- gsub("\\s+", "_", region_name)
  
  write.csv(
    summary_table,
    file.path(out_dir, paste0(prefix, "_trait_congruence_k_table.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  write.csv(
    honest_seq,
    file.path(out_dir, paste0(prefix, "_sequential_expansion.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  write.csv(
    loo_table,
    file.path(out_dir, paste0(prefix, "_leave_one_out.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  selected_traits <- data.frame(
    panel = c(
      "Best screening subset",
      "Honest max-over-k",
      "Honest fixed k=12",
      "Sequential expansion final k=12"
    ),
    traits = c(
      paste(best_traits_screen, collapse = "; "),
      paste(honest_traits, collapse = "; "),
      paste(traits_kfixed, collapse = "; "),
      if (length(traits_core12)) paste(traits_core12, collapse = "; ") else ""
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    selected_traits,
    file.path(out_dir, paste0(prefix, "_selected_traits.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  article_summary <- data.frame(
    analysis = c(
      "All traits Mantel",
      "Best screening subset",
      "Honest max-over-k",
      "Honest fixed k=12",
      "Sequential expansion final k=12"
    ),
    k = c(
      p,
      best_k_screen,
      honest_k,
      k_fixed,
      k_fixed
    ),
    statistic = c(
      mantel_all$r,
      max(summary_table$r_max),
      honest_allk$T_obs,
      honest_kfixed$T_obs,
      if (nrow(row12)) row12$statistic[1] else NA_real_
    ),
    p = c(
      mantel_all$p,
      summary_table$p_best[which.max(summary_table$r_max)],
      honest_allk$p_honest,
      honest_kfixed$p_honest,
      if (nrow(row12)) row12$p_value[1] else NA_real_
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    article_summary,
    file.path(out_dir, paste0(prefix, "_article_summary.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
  
  list(
    region = region_name,
    mantel_all = mantel_all,
    summary_table = summary_table,
    best_k_screen = best_k_screen,
    best_traits_screen = best_traits_screen,
    honest_allk = honest_allk,
    honest_traits = honest_traits,
    honest_kfixed = honest_kfixed,
    traits_kfixed = traits_kfixed,
    honest_seq = honest_seq,
    loo_table = loo_table
  )
}

# =========================
# LOAD DATA
# =========================

ural_m <- read_morph_csv(ural_traits_file)
ural_f <- read_fst_csv(ural_fst_file)

siberia_m <- read_morph_csv(siberia_traits_file)
siberia_f <- read_fst_csv(siberia_fst_file)

ural <- match_region_data(
  morph_list = ural_m,
  fst_mat = ural_f,
  expected_order = ural_order,
  expected_traits = 30,
  region_name = "Ural"
)

siberia <- match_region_data(
  morph_list = siberia_m,
  fst_mat = siberia_f,
  expected_order = siberia_order,
  expected_traits = 30,
  region_name = "Siberia"
)

# =========================
# RUN BOTH REGIONS
# =========================

res_ural <- run_region_analysis(
  "Ural",
  ural$P,
  ural$n,
  ural$fst
)

res_siberia <- run_region_analysis(
  "Siberia",
  siberia$P,
  siberia$n,
  siberia$fst
)

# =========================
# SAVE COMBINED SUMMARIES
# =========================

article_rows <- rbind(
  data.frame(
    region = "Ural",
    analysis = c(
      "All traits Mantel",
      "Best screening subset",
      "Honest max-over-k",
      "Honest fixed k=12",
      "Sequential expansion final k=12"
    ),
    k = c(
      ncol(ural$P),
      res_ural$best_k_screen,
      res_ural$honest_allk$best_k_obs,
      12,
      12
    ),
    statistic = c(
      res_ural$mantel_all$r,
      max(res_ural$summary_table$r_max),
      res_ural$honest_allk$T_obs,
      res_ural$honest_kfixed$T_obs,
      res_ural$honest_seq$statistic[res_ural$honest_seq$k == 12][1]
    ),
    p = c(
      res_ural$mantel_all$p,
      res_ural$summary_table$p_best[which.max(res_ural$summary_table$r_max)],
      res_ural$honest_allk$p_honest,
      res_ural$honest_kfixed$p_honest,
      res_ural$honest_seq$p_value[res_ural$honest_seq$k == 12][1]
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    region = "Siberia",
    analysis = c(
      "All traits Mantel",
      "Best screening subset",
      "Honest max-over-k",
      "Honest fixed k=12",
      "Sequential expansion final k=12"
    ),
    k = c(
      ncol(siberia$P),
      res_siberia$best_k_screen,
      res_siberia$honest_allk$best_k_obs,
      12,
      12
    ),
    statistic = c(
      res_siberia$mantel_all$r,
      max(res_siberia$summary_table$r_max),
      res_siberia$honest_allk$T_obs,
      res_siberia$honest_kfixed$T_obs,
      res_siberia$honest_seq$statistic[res_siberia$honest_seq$k == 12][1]
    ),
    p = c(
      res_siberia$mantel_all$p,
      res_siberia$summary_table$p_best[which.max(res_siberia$summary_table$r_max)],
      res_siberia$honest_allk$p_honest,
      res_siberia$honest_kfixed$p_honest,
      res_siberia$honest_seq$p_value[res_siberia$honest_seq$k == 12][1]
    ),
    stringsAsFactors = FALSE
  )
)

write.csv(
  article_rows,
  file.path(out_dir, "article_summary_both_regions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

loo_both_regions <- rbind(
  cbind(region = "Ural", res_ural$loo_table),
  cbind(region = "Siberia", res_siberia$loo_table)
)

write.csv(
  loo_both_regions,
  file.path(out_dir, "leave_one_out_both_regions.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("Outputs saved to:\n")
cat("  ", out_dir, "\n", sep = "")
cat("============================================================\n")