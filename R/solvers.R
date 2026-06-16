# R Solvers for sadMCDM Package
# Requires: lpSolve

#' @importFrom lpSolve lp
NULL

# --- 1. AHP Helper Functions ---

calculate_geometric_mean_weights <- function(mat) {
  n <- nrow(mat)
  if (n <= 2) {
    return(list(weights = rep(1/n, n), cr = 0.0))
  }
  
  # Row geometric means
  geom_means <- exp(rowMeans(log(mat)))
  
  # Normalize weights
  weights <- geom_means / sum(geom_means)
  
  # Principal eigenvalue lambda_max
  weighted_sum <- mat %*% weights
  eigenvalues <- weighted_sum / weights
  lambda_max <- mean(eigenvalues)
  
  # Consistency Index (CI)
  ci <- (lambda_max - n) / (n - 1)
  
  # Consistency Ratio (CR)
  random_index <- c(0.0, 0.0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49)
  ri <- if (n <= length(random_index)) random_index[n] else 1.49
  cr <- if (ri > 0) ci / ri else 0.0
  
  list(weights = as.vector(weights), cr = as.double(cr))
}

map_ratio_to_score <- function(ratio) {
  if (ratio <= 1.05) return(1)
  if (ratio <= 1.2) return(2)
  if (ratio >= 3.0) return(9)
  score <- 2 + (ratio - 1.2) * (7.0 / 1.8)
  return(as.integer(round(score)))
}

map_ratio_to_saaty <- function(val_i, val_j, is_benefit=TRUE) {
  if (val_i == val_j) return(1.0)
  if (!is_benefit) {
    temp <- val_i
    val_i <- val_j
    val_j <- temp
  }
  if (val_j == 0) {
    ratio <- if (val_i > 0) 9.0 else 1.0
  } else {
    ratio <- val_i / val_j
  }
  if (ratio < 1.0) {
    inv_ratio <- 1.0 / ratio
    score <- map_ratio_to_score(inv_ratio)
    return(1.0 / score)
  } else {
    score <- map_ratio_to_score(ratio)
    return(as.double(score))
  }
}

# --- 2. solve_ahp ---

solve_ahp <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  # Elicit Criteria Weights
  crit_matrix_raw <- preference_data$criteria_matrix
  if (is.null(crit_matrix_raw)) {
    crit_matrix <- matrix(1.0, nrow=n, ncol=n)
  } else {
    crit_matrix <- matrix(as.numeric(unlist(crit_matrix_raw)), nrow=n, ncol=n, byrow=TRUE)
  }
  
  # Enforce reciprocal
  for (i in 1:n) {
    crit_matrix[i, i] <- 1.0
    if (i < n) {
      for (j in (i+1):n) {
        if (crit_matrix[i, j] == 0) {
          crit_matrix[i, j] <- 1.0
        }
        crit_matrix[j, i] <- 1.0 / crit_matrix[i, j]
      }
    }
  }
  
  crit_res <- calculate_geometric_mean_weights(crit_matrix)
  crit_weights <- crit_res$weights
  crit_cr <- crit_res$cr
  
  # Elicit Alternative Priorities
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  alt_crs <- list()
  
  alt_matrices_pref <- preference_data$alternatives_matrices
  if (is.null(alt_matrices_pref)) alt_matrices_pref <- list()
  
  for (j in 1:n) {
    cid <- criteria_ids[j]
    cid_str <- as.character(cid)
    
    if (cid_str %in% names(alt_matrices_pref)) {
      alt_matrix <- matrix(as.numeric(unlist(alt_matrices_pref[[cid_str]])), nrow=m, ncol=m, byrow=TRUE)
    } else {
      # Pre-fill based on consequence ratio
      alt_matrix <- matrix(1.0, nrow=m, ncol=m)
      is_benefit <- (criteria_types[j] == "benefit")
      for (r in 1:m) {
        val_r <- matrix_data[[r]][[j]]
        for (c in 1:m) {
          val_c <- matrix_data[[c]][[j]]
          alt_matrix[r, c] <- map_ratio_to_saaty(val_r, val_c, is_benefit)
        }
      }
    }
    
    # Enforce reciprocal
    for (r in 1:m) {
      alt_matrix[r, r] <- 1.0
      if (r < m) {
        for (c in (r+1):m) {
          if (alt_matrix[r, c] == 0) {
            alt_matrix[r, c] <- 1.0
          }
          alt_matrix[c, r] <- 1.0 / alt_matrix[r, c]
        }
      }
    }
    
    alt_res <- calculate_geometric_mean_weights(alt_matrix)
    norm_matrix[, j] <- alt_res$weights
    alt_crs[[cid_str]] <- alt_res$cr
  }
  
  # Compute global score
  global_scores <- as.vector(norm_matrix %*% crit_weights)
  
  # Ranks (descending order, ties resolved by order)
  ranks <- rank(-global_scores, ties.method = "first")
  
  list(
    weights = crit_weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks),
    criteria_cr = crit_cr,
    alternatives_cr = alt_crs
  )
}

# --- 3. solve_bwm ---

solve_bwm_weights_lp <- function(n, best_idx, worst_idx, best_to_others, others_to_worst) {
  # Variables: w_1, ..., w_n, xi (n + 1 variables)
  obj <- c(rep(0, n), 1)
  
  # Constraints: Sum w_i = 1, plus inequalities
  # w_i in [0,1], xi >= 0 (handled by lpSolve non-negative bounds except upper bound)
  # We will define explicitly the upper bounds of w_i <= 1
  # and constraints matrix
  max_constraints <- 1 + n + 4 * n
  const_mat <- matrix(0, nrow = max_constraints, ncol = n + 1)
  const_dir <- rep("", max_constraints)
  const_rhs <- rep(0, max_constraints)
  
  # 1. Sum w_i = 1
  const_mat[1, 1:n] <- 1.0
  const_dir[1] <- "="
  const_rhs[1] <- 1.0
  
  row_idx <- 2
  
  # 2. Upper bounds w_i <= 1
  for (i in 1:n) {
    const_mat[row_idx, i] <- 1.0
    const_dir[row_idx] <- "<="
    const_rhs[row_idx] <- 1.0
    row_idx <- row_idx + 1
  }
  
  # 3. Inequalities
  for (j in 1:n) {
    a_Bj <- best_to_others[j]
    a_jW <- others_to_worst[j]
    
    if (j != best_idx) {
      # w_B - a_Bj * w_j - xi <= 0
      const_mat[row_idx, best_idx] <- 1.0
      const_mat[row_idx, j] <- -a_Bj
      const_mat[row_idx, n + 1] <- -1.0
      const_dir[row_idx] <- "<="
      const_rhs[row_idx] <- 0.0
      row_idx <- row_idx + 1
      
      # -w_B + a_Bj * w_j - xi <= 0
      const_mat[row_idx, best_idx] <- -1.0
      const_mat[row_idx, j] <- a_Bj
      const_mat[row_idx, n + 1] <- -1.0
      const_dir[row_idx] <- "<="
      const_rhs[row_idx] <- 0.0
      row_idx <- row_idx + 1
    }
    
    if (j != worst_idx) {
      # w_j - a_jW * w_W - xi <= 0
      const_mat[row_idx, j] <- 1.0
      const_mat[row_idx, worst_idx] <- -a_jW
      const_mat[row_idx, n + 1] <- -1.0
      const_dir[row_idx] <- "<="
      const_rhs[row_idx] <- 0.0
      row_idx <- row_idx + 1
      
      # -w_j + a_jW * w_W - xi <= 0
      const_mat[row_idx, j] <- -1.0
      const_mat[row_idx, worst_idx] <- a_jW
      const_mat[row_idx, n + 1] <- -1.0
      const_dir[row_idx] <- "<="
      const_rhs[row_idx] <- 0.0
      row_idx <- row_idx + 1
    }
  }
  
  # Trim matrices
  const_mat <- const_mat[1:(row_idx-1), , drop=FALSE]
  const_dir <- const_dir[1:(row_idx-1)]
  const_rhs <- const_rhs[1:(row_idx-1)]
  
  res <- lpSolve::lp("min", obj, const_mat, const_dir, const_rhs)
  
  if (res$status == 0) {
    weights <- res$solution[1:n]
    xi <- res$solution[n + 1]
    weights[weights < 0] <- 0.0
    weights[weights > 1] <- 1.0
    weights <- weights / sum(weights)
    return(list(weights = weights, xi = xi, success = TRUE))
  } else {
    return(list(weights = rep(1/n, n), xi = 0.0, success = FALSE))
  }
}

solve_bwm <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  # Map indices (accepting both 0-based and 1-based, converting to 1-based)
  best_idx <- as.integer(preference_data$best_idx)
  if (!is.na(best_idx) && best_idx < n) best_idx <- best_idx + 1
  if (is.na(best_idx) || best_idx <= 0 || best_idx > n) best_idx <- 1
  
  worst_idx <- as.integer(preference_data$worst_idx)
  if (!is.na(worst_idx) && worst_idx < n) worst_idx <- worst_idx + 1
  if (is.na(worst_idx) || worst_idx <= 0 || worst_idx > n) worst_idx <- n
  
  best_to_others <- as.numeric(preference_data$best_to_others)
  if (length(best_to_others) != n) best_to_others <- rep(1.0, n)
  
  others_to_worst <- as.numeric(preference_data$others_to_worst)
  if (length(others_to_worst) != n) others_to_worst <- rep(1.0, n)
  
  crit_res <- solve_bwm_weights_lp(n, best_idx, worst_idx, best_to_others, others_to_worst)
  crit_weights <- crit_res$weights
  xi <- crit_res$xi
  crit_success <- crit_res$success
  
  # Elicit Alternative Priorities (identical to AHP)
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  alt_crs <- list()
  
  alt_matrices_pref <- preference_data$alternatives_matrices
  if (is.null(alt_matrices_pref)) alt_matrices_pref <- list()
  
  for (j in 1:n) {
    cid <- criteria_ids[j]
    cid_str <- as.character(cid)
    
    if (cid_str %in% names(alt_matrices_pref)) {
      alt_matrix <- matrix(as.numeric(unlist(alt_matrices_pref[[cid_str]])), nrow=m, ncol=m, byrow=TRUE)
    } else {
      alt_matrix <- matrix(1.0, nrow=m, ncol=m)
      is_benefit <- (criteria_types[j] == "benefit")
      for (r in 1:m) {
        val_r <- matrix_data[[r]][[j]]
        for (c in 1:m) {
          val_c <- matrix_data[[c]][[j]]
          alt_matrix[r, c] <- map_ratio_to_saaty(val_r, val_c, is_benefit)
        }
      }
    }
    
    # Enforce reciprocal
    for (r in 1:m) {
      alt_matrix[r, r] <- 1.0
      if (r < m) {
        for (c in (r+1):m) {
          if (alt_matrix[r, c] == 0) {
            alt_matrix[r, c] <- 1.0
          }
          alt_matrix[c, r] <- 1.0 / alt_matrix[r, c]
        }
      }
    }
    
    alt_res <- calculate_geometric_mean_weights(alt_matrix)
    norm_matrix[, j] <- alt_res$weights
    alt_crs[[cid_str]] <- alt_res$cr
  }
  
  # Compute global score
  global_scores <- as.vector(norm_matrix %*% crit_weights)
  ranks <- rank(-global_scores, ties.method = "first")
  
  list(
    weights = crit_weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks),
    criteria_success = crit_success,
    consistency_xi = xi,
    alternatives_cr = alt_crs
  )
}

# --- 4. solve_bwt ---

interpolate_bisection <- function(val, col_min, col_max, x_05, is_benefit) {
  if (col_max == col_min) return(1.0)
  
  # Safeguard bounds
  x_05 <- max(col_min + 1e-9, min(col_max - 1e-9, x_05))
  val <- max(col_min, min(col_max, val))
  
  if (is_benefit) {
    if (val <= x_05) {
      denom <- x_05 - col_min
      return(0.5 * (val - col_min) / denom)
    } else {
      denom <- col_max - x_05
      return(0.5 + 0.5 * (val - x_05) / denom)
    }
  } else {
    if (val <= x_05) {
      denom <- x_05 - col_min
      return(1.0 - 0.5 * (val - col_min) / denom)
    } else {
      denom <- col_max - x_05
      return(0.5 * (col_max - val) / denom)
    }
  }
}

solve_bwt_weights_lp <- function(n, best_idx, worst_idx, best_to_others, others_to_worst) {
  # BWT minimax LP has the identical structure as BWM weights LP
  solve_bwm_weights_lp(n, best_idx, worst_idx, best_to_others, others_to_worst)
}

solve_bwt <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  best_idx <- as.integer(preference_data$best_idx)
  if (!is.na(best_idx) && best_idx < n) best_idx <- best_idx + 1
  if (is.na(best_idx) || best_idx <= 0 || best_idx > n) best_idx <- 1
  
  worst_idx <- as.integer(preference_data$worst_idx)
  if (!is.na(worst_idx) && worst_idx < n) worst_idx <- worst_idx + 1
  if (is.na(worst_idx) || worst_idx <= 0 || worst_idx > n) worst_idx <- n
  
  best_to_others <- as.numeric(preference_data$best_to_others)
  if (length(best_to_others) != n) best_to_others <- rep(1.0, n)
  
  others_to_worst <- as.numeric(preference_data$others_to_worst)
  if (length(others_to_worst) != n) others_to_worst <- rep(1.0, n)
  
  crit_res <- solve_bwt_weights_lp(n, best_idx, worst_idx, best_to_others, others_to_worst)
  crit_weights <- crit_res$weights
  xi <- crit_res$xi
  crit_success <- crit_res$success
  
  # Alternative Priorities
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  bisection_midpoints <- preference_data$bisection_midpoints
  if (is.null(bisection_midpoints)) bisection_midpoints <- list()
  
  for (j in 1:n) {
    cid <- criteria_ids[j]
    cid_str <- as.character(cid)
    
    col <- sapply(matrix_data, function(row) row[[j]])
    col_min <- min(col)
    col_max <- max(col)
    
    default_mid <- (col_min + col_max) / 2.0
    x_05 <- if (cid_str %in% names(bisection_midpoints)) as.numeric(bisection_midpoints[[cid_str]]) else default_mid
    if (is.na(x_05)) x_05 <- default_mid
    
    is_benefit <- (criteria_types[j] == "benefit")
    for (r in 1:m) {
      val <- matrix_data[[r]][[j]]
      norm_matrix[r, j] <- interpolate_bisection(val, col_min, col_max, x_05, is_benefit)
    }
  }
  
  global_scores <- as.vector(norm_matrix %*% crit_weights)
  ranks <- rank(-global_scores, ties.method = "first")
  
  list(
    weights = crit_weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks),
    criteria_success = crit_success,
    consistency_xi = xi
  )
}

# --- 5. solve_macbeth ---

solve_macbeth_lp <- function(comparison_matrix) {
  n <- nrow(comparison_matrix)
  if (n <= 1) {
    return(list(weights = rep(1.0, n), success = TRUE))
  }
  
  # Setup LP: variables s_1, ..., s_n, delta (n + 1 variables)
  # Maximize delta => Minimize -delta
  obj <- c(rep(0, n), -1.0)
  
  # Find anchors
  net_pref <- rep(0, n)
  for (i in 1:n) {
    for (j in 1:n) {
      val <- comparison_matrix[i, j]
      if (!is.null(val) && !is.na(val)) {
        if (val > 0) {
          net_pref[i] <- net_pref[i] + 1
          net_pref[j] <- net_pref[j] - 1
        } else if (val < 0) {
          net_pref[i] <- net_pref[i] - 1
          net_pref[j] <- net_pref[j] + 1
        }
      }
    }
  }
  
  best_idx <- which.max(net_pref)
  worst_idx <- which.min(net_pref)
  if (best_idx == worst_idx) {
    worst_idx <- if (best_idx == n) 1 else best_idx + 1
  }
  
  # Count constraints
  # s[best] = 1 (eq)
  # s[worst] = 0 (eq)
  # bounds (s_i <= 1) (ub)
  # comparison conditions
  max_constr <- 2 + n + n*n
  const_mat <- matrix(0, nrow=max_constr, ncol=n+1)
  const_dir <- rep("", max_constr)
  const_rhs <- rep(0, max_constr)
  
  # s[best] = 1
  const_mat[1, best_idx] <- 1.0
  const_dir[1] <- "="
  const_rhs[1] <- 1.0
  
  # s[worst] = 0
  const_mat[2, worst_idx] <- 1.0
  const_dir[2] <- "="
  const_rhs[2] <- 0.0
  
  row_idx <- 3
  
  # bounds s_i <= 1
  for (i in 1:n) {
    const_mat[row_idx, i] <- 1.0
    const_dir[row_idx] <- "<="
    const_rhs[row_idx] <- 1.0
    row_idx <- row_idx + 1
  }
  
  # comparisons
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      val <- comparison_matrix[i, j]
      if (is.null(val) || is.na(val)) next
      
      if (val == 0) {
        # s_i - s_j = 0
        const_mat[row_idx, i] <- 1.0
        const_mat[row_idx, j] <- -1.0
        const_dir[row_idx] <- "="
        const_rhs[row_idx] <- 0.0
        row_idx <- row_idx + 1
      } else if (val > 0) {
        # s_j - s_i + val*delta <= 0
        step <- as.integer(val)
        const_mat[row_idx, i] <- -1.0
        const_mat[row_idx, j] <- 1.0
        const_mat[row_idx, n+1] <- step
        const_dir[row_idx] <- "<="
        const_rhs[row_idx] <- 0.0
        row_idx <- row_idx + 1
      }
    }
  }
  
  const_mat <- const_mat[1:(row_idx-1), , drop=FALSE]
  const_dir <- const_dir[1:(row_idx-1)]
  const_rhs <- const_rhs[1:(row_idx-1)]
  
  # Set delta bounds >= 0.001 (we add as constraint delta >= 0.001)
  # delta >= 0.001 => -delta <= -0.001
  delta_mat <- c(rep(0, n), -1.0)
  const_mat <- rbind(const_mat, delta_mat)
  const_dir <- c(const_dir, "<=")
  const_rhs <- c(const_rhs, -0.001)
  
  res <- lpSolve::lp("min", obj, const_mat, const_dir, const_rhs)
  
  if (res$status == 0) {
    scores <- res$solution[1:n]
    scores[scores < 0] <- 0.0
    scores[scores > 1] <- 1.0
    total <- sum(scores)
    weights <- if (total > 0) scores / total else rep(1/n, n)
    return(list(weights = weights, success = TRUE))
  } else {
    scores <- net_pref - min(net_pref)
    total <- sum(scores)
    weights <- if (total > 0) scores / total else rep(1/n, n)
    return(list(weights = weights, success = FALSE))
  }
}

get_macbeth_scores <- function(mat, n, is_benefit) {
  if (is.null(mat)) {
    return(if (is_benefit) seq(0, 1, length.out=n) else seq(1, 0, length.out=n))
  }
  
  # Setup LP for score values s in [0,1]
  obj <- c(rep(0, n), -1.0) # maximize delta
  
  worst_idx <- if (is_benefit) 1 else n
  best_idx <- if (is_benefit) n else 1
  
  max_constr <- 2 + n + n*n + 1
  const_mat <- matrix(0, nrow=max_constr, ncol=n+1)
  const_dir <- rep("", max_constr)
  const_rhs <- rep(0, max_constr)
  
  # s[best] = 1
  const_mat[1, best_idx] <- 1.0
  const_dir[1] <- "="
  const_rhs[1] <- 1.0
  
  # s[worst] = 0
  const_mat[2, worst_idx] <- 1.0
  const_dir[2] <- "="
  const_rhs[2] <- 0.0
  
  row_idx <- 3
  for (i in 1:n) {
    const_mat[row_idx, i] <- 1.0
    const_dir[row_idx] <- "<="
    const_rhs[row_idx] <- 1.0
    row_idx <- row_idx + 1
  }
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      val <- mat[i, j]
      if (is.null(val) || is.na(val)) next
      
      if (val == 0) {
        const_mat[row_idx, i] <- 1.0
        const_mat[row_idx, j] <- -1.0
        const_dir[row_idx] <- "="
        const_rhs[row_idx] <- 0.0
        row_idx <- row_idx + 1
      } else if (val > 0) {
        step <- as.integer(val)
        const_mat[row_idx, i] <- -1.0
        const_mat[row_idx, j] <- 1.0
        const_mat[row_idx, n+1] <- step
        const_dir[row_idx] <- "<="
        const_rhs[row_idx] <- 0.0
        row_idx <- row_idx + 1
      }
    }
  }
  
  const_mat <- const_mat[1:(row_idx-1), , drop=FALSE]
  const_dir <- const_dir[1:(row_idx-1)]
  const_rhs <- const_rhs[1:(row_idx-1)]
  
  # delta >= 0.001
  delta_mat <- c(rep(0, n), -1.0)
  const_mat <- rbind(const_mat, delta_mat)
  const_dir <- c(const_dir, "<=")
  const_rhs <- c(const_rhs, -0.001)
  
  res <- lpSolve::lp("min", obj, const_mat, const_dir, const_rhs)
  if (res$status == 0) {
    scores <- res$solution[1:n]
    scores[scores < 0] <- 0.0
    scores[scores > 1] <- 1.0
    return(scores)
  } else {
    return(if (is_benefit) seq(0, 1, length.out=n) else seq(1, 0, length.out=n))
  }
}

solve_macbeth <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  crit_matrix_raw <- preference_data$criteria_matrix
  if (is.null(crit_matrix_raw)) {
    crit_matrix <- matrix(1.0, nrow=n, ncol=n)
  } else {
    crit_matrix <- matrix(as.numeric(unlist(crit_matrix_raw)), nrow=n, ncol=n, byrow=TRUE)
  }
  
  crit_res <- solve_macbeth_lp(crit_matrix)
  crit_weights <- crit_res$weights
  crit_success <- crit_res$success
  
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  levels_success <- list()
  levels_matrices_pref <- preference_data$levels_matrices
  if (is.null(levels_matrices_pref)) levels_matrices_pref <- list()
  
  for (j in 1:n) {
    cid <- criteria_ids[j]
    cid_str <- as.character(cid)
    
    col <- sapply(matrix_data, function(row) row[[j]])
    col_min <- min(col)
    col_max <- max(col)
    
    levels <- seq(col_min, col_max, length.out=5)
    
    if (cid_str %in% names(levels_matrices_pref)) {
      levels_matrix <- matrix(as.numeric(unlist(levels_matrices_pref[[cid_str]])), nrow=5, ncol=5, byrow=TRUE)
      levels_success[[cid_str]] <- TRUE # will check solver
    } else {
      levels_matrix <- NULL
      levels_success[[cid_str]] <- TRUE
    }
    
    scores_levels <- get_macbeth_scores(levels_matrix, 5, criteria_types[j] == "benefit")
    
    for (r in 1:m) {
      val <- matrix_data[[r]][[j]]
      if (col_max == col_min) {
        norm_matrix[r, j] <- 1.0
      } else {
        norm_matrix[r, j] <- approx(levels, scores_levels, xout=val, rule=2)$y
      }
    }
  }
  
  global_scores <- as.vector(norm_matrix %*% crit_weights)
  ranks <- rank(-global_scores, ties.method = "first")
  
  list(
    weights = crit_weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks),
    criteria_success = crit_success,
    levels_success = levels_success
  )
}

# --- 6. solve_electre ---

solve_electre <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  matrix <- matrix(0.0, nrow=m, ncol=n)
  for (i in 1:m) {
    for (j in 1:n) {
      matrix[i, j] <- matrix_data[[i]][[j]]
    }
  }
  
  version <- preference_data$electre_version
  if (is.null(version)) version <- "I"
  
  # weights
  weights <- rep(0.0, n)
  if (version != "IV") {
    scores <- preference_data$weights
    if (is.null(scores)) scores <- list()
    total_score <- 0.0
    for (i in 1:n) {
      cid <- criteria_ids[i]
      score <- if (as.character(cid) %in% names(scores)) as.numeric(scores[[as.character(cid)]]) else 10.0
      weights[i] <- score
      total_score <- total_score + score
    }
    weights <- if (total_score > 0) weights / total_score else rep(1/n, n)
  } else {
    weights <- rep(1/n, n)
  }
  
  # ranges
  ranges <- rep(0.0, n)
  for (j in 1:n) {
    col <- matrix[, j]
    r <- max(col) - min(col)
    ranges[j] <- if (r > 0) r else 1.0
  }
  
  # thresholds
  q <- rep(0.0, n)
  p <- rep(1e-9, n)
  v <- rep(Inf, n)
  thresholds_pref <- preference_data$thresholds
  if (is.null(thresholds_pref)) thresholds_pref <- list()
  
  for (j in 1:n) {
    cid_str <- as.character(criteria_ids[j])
    if (cid_str %in% names(thresholds_pref)) {
      t <- thresholds_pref[[cid_str]]
      if ("q" %in% names(t)) q[j] <- as.numeric(t$q)
      if ("p" %in% names(t)) p[j] <- as.numeric(t$p)
      if ("v" %in% names(t)) v[j] <- as.numeric(t$v)
    }
  }
  
  # Normalize matrix
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  for (j in 1:n) {
    col <- matrix[, j]
    c_min <- min(col)
    c_max <- max(col)
    denom <- c_max - c_min
    if (denom == 0) denom <- 1.0
    if (criteria_types[j] == "benefit") {
      norm_matrix[, j] <- (col - c_min) / denom
    } else {
      norm_matrix[, j] <- (c_max - col) / denom
    }
  }
  
  global_scores <- rep(0.0, m)
  ranks <- rep(1, m)
  extra_data <- list()
  
  if (version %in% c("I", "IS")) {
    concordance <- matrix(0.0, nrow=m, ncol=m)
    discordance <- matrix(0.0, nrow=m, ncol=m)
    
    c_threshold <- if ("concordance_threshold" %in% names(preference_data)) as.numeric(preference_data$concordance_threshold) else 0.6
    d_threshold <- if ("discordance_threshold" %in% names(preference_data)) as.numeric(preference_data$discordance_threshold) else 0.4
    
    for (a in 1:m) {
      for (b in 1:m) {
        if (a == b) next
        c_sum <- 0.0
        max_disc <- 0.0
        vetoed <- FALSE
        
        for (j in 1:n) {
          diff <- norm_matrix[a, j] - norm_matrix[b, j]
          
          if (version == "I") {
            if (diff >= 0) {
              c_sum <- c_sum + weights[j]
            } else {
              val <- -diff
              if (val > max_disc) max_disc <- val
            }
          } else { # IS
            norm_q <- q[j] / ranges[j]
            norm_p <- p[j] / ranges[j]
            norm_v <- v[j] / ranges[j]
            
            if (diff >= -norm_q) {
              c_sum <- c_sum + weights[j]
            } else if (diff < -norm_p) {
              c_sum <- c_sum + 0.0
            } else {
              c_sum <- c_sum + weights[j] * (norm_p + diff) / (norm_p - norm_q)
            }
            
            if (diff < -norm_v) {
              vetoed <- TRUE
            }
          }
        }
        concordance[a, b] <- c_sum
        if (version == "I") {
          discordance[a, b] <- max_disc
        } else {
          discordance[a, b] <- if (vetoed) 1.0 else 0.0
        }
      }
    }
    
    # Outranking graph
    outranks <- matrix(FALSE, nrow=m, ncol=m)
    for (a in 1:m) {
      for (b in 1:m) {
        if (a != b) {
          if (version == "I") {
            outranks[a, b] <- (concordance[a, b] >= c_threshold) && (discordance[a, b] <= d_threshold)
          } else {
            outranks[a, b] <- (concordance[a, b] >= c_threshold) && (discordance[a, b] == 0.0)
          }
        }
      }
    }
    
    # Kernel
    kernel <- c()
    for (i in 1:m) {
      outranked <- FALSE
      for (j in 1:m) {
        if (j != i && outranks[j, i]) {
          outranked <- TRUE
          break
        }
      }
      if (!outranked) {
        kernel <- c(kernel, i)
      }
    }
    if (length(kernel) == 0) kernel <- seq_len(m)
    
    global_scores <- as.vector(rowSums(outranks))
    ranks <- rank(-global_scores, ties.method = "first")
    
    extra_data <- list(
      concordance = concordance,
      discordance = discordance,
      kernel = alternatives_ids[kernel]
    )
    
  } else if (version == "II") {
    # Distillation
    c_strong <- if ("concordance_strong" %in% names(preference_data)) as.numeric(preference_data$concordance_strong) else 0.7
    c_weak <- if ("concordance_weak" %in% names(preference_data)) as.numeric(preference_data$concordance_weak) else 0.5
    d_strong <- if ("discordance_strong" %in% names(preference_data)) as.numeric(preference_data$discordance_strong) else 0.3
    d_weak <- if ("discordance_weak" %in% names(preference_data)) as.numeric(preference_data$discordance_weak) else 0.5
    
    concordance <- matrix(0.0, nrow=m, ncol=m)
    discordance <- matrix(0.0, nrow=m, ncol=m)
    for (a in 1:m) {
      for (b in 1:m) {
        if (a == b) next
        c_sum <- 0.0
        max_disc <- 0.0
        for (j in 1:n) {
          diff <- norm_matrix[a, j] - norm_matrix[b, j]
          if (diff >= 0) {
            c_sum <- c_sum + weights[j]
          } else {
            val <- -diff
            if (val > max_disc) max_disc <- val
          }
        }
        concordance[a, b] <- c_sum
        discordance[a, b] <- max_disc
      }
    }
    
    S_strong <- matrix(FALSE, nrow=m, ncol=m)
    S_weak <- matrix(FALSE, nrow=m, ncol=m)
    for (a in 1:m) {
      for (b in 1:m) {
        if (a == b) next
        S_strong[a, b] <- (concordance[a, b] >= c_strong) && (discordance[a, b] <= d_strong)
        S_weak[a, b] <- (concordance[a, b] >= c_weak) && (discordance[a, b] <= d_weak)
      }
    }
    
    # Helper recursive distillation
    distill <- function(alternatives, is_forward=TRUE) {
      if (length(alternatives) == 0) return(c())
      sub_m <- length(alternatives)
      deg <- rep(0.0, sub_m)
      for (i in 1:sub_m) {
        a <- alternatives[i]
        for (j in 1:sub_m) {
          b <- alternatives[j]
          if (i != j) {
            if (S_strong[a, b]) deg[i] <- deg[i] + 1
            if (is_forward && S_weak[a, b]) deg[i] <- deg[i] + 0.5
          }
        }
      }
      
      target_val <- if (is_forward) max(deg) else min(deg)
      best_subset <- alternatives[which(deg == target_val)]
      
      if (length(best_subset) == length(alternatives)) {
        return(best_subset)
      } else {
        return(distill(best_subset, is_forward))
      }
    }
    
    alts <- seq_len(m)
    ranking_order <- c()
    while (length(alts) > 0) {
      best_alts <- distill(alts, is_forward=TRUE)
      for (ba in best_alts) {
        ranking_order <- c(ranking_order, ba)
        alts <- alts[alts != ba]
      }
    }
    
    for (rank_pos in 1:m) {
      idx <- ranking_order[rank_pos]
      ranks[idx] <- rank_pos
      global_scores[idx] <- as.double(m - rank_pos + 1)
    }
    
    extra_data <- list(
      concordance = concordance,
      discordance = discordance
    )
    
  } else if (version %in% c("III", "IV")) {
    credibility <- matrix(0.0, nrow=m, ncol=m)
    for (a in 1:m) {
      for (b in 1:m) {
        if (a == b) next
        c_sum <- 0.0
        for (j in 1:n) {
          diff <- if (criteria_types[j] == "benefit") matrix[a, j] - matrix[b, j] else matrix[b, j] - matrix[a, j]
          if (diff >= -q[j]) {
            c_sum <- c_sum + weights[j]
          } else if (diff < -p[j]) {
            c_sum <- c_sum + 0.0
          } else {
            c_sum <- c_sum + weights[j] * (p[j] + diff) / (p[j] - q[j])
          }
        }
        
        d_j <- rep(0.0, n)
        for (j in 1:n) {
          diff <- if (criteria_types[j] == "benefit") matrix[a, j] - matrix[b, j] else matrix[b, j] - matrix[a, j]
          if (diff >= -p[j]) {
            d_j[j] <- 0.0
          } else if (diff < -v[j]) {
            d_j[j] <- 1.0
          } else {
            d_j[j] <- (-diff - p[j]) / (v[j] - p[j])
          }
        }
        
        rho <- c_sum
        for (j in 1:n) {
          if (d_j[j] > c_sum) {
            rho <- rho * (1.0 - d_j[j]) / (1.0 - c_sum)
          }
        }
        credibility[a, b] <- rho
      }
    }
    
    # Distill based on credibility
    # Sort by Net flows as a proxy for distillation ranking to keep it simple and clean
    net_flows <- rowSums(credibility) - colSums(credibility)
    ranks <- rank(-net_flows, ties.method = "first")
    global_scores <- as.vector(net_flows)
    
    extra_data <- list(
      credibility = credibility
    )
    
  } else if (version == "TRI") {
    profiles_list <- preference_data$profiles
    if (is.null(profiles_list) || length(profiles_list) == 0) {
      # fallback
      k_categories <- if ("num_categories" %in% names(preference_data)) as.integer(preference_data$num_categories) else 3
      profiles <- matrix(0.0, nrow=k_categories-1, ncol=n)
      for (k in 1:(k_categories-1)) {
        for (j in 1:n) {
          col <- matrix[, j]
          profiles[k, j] <- min(col) + (max(col) - min(col)) * (k / k_categories)
        }
      }
    } else {
      profiles <- matrix(as.numeric(unlist(profiles_list)), ncol=n, byrow=TRUE)
    }
    
    num_profiles <- nrow(profiles)
    rho_a_b <- matrix(0.0, nrow=m, ncol=num_profiles)
    lambda_t <- if ("lambda_threshold" %in% names(preference_data)) as.numeric(preference_data$lambda_threshold) else 0.6
    
    for (i in 1:m) {
      for (h in 1:num_profiles) {
        c_a_b <- 0.0
        for (j in 1:n) {
          val_a <- matrix[i, j]
          val_b <- profiles[h, j]
          diff_a_b <- if (criteria_types[j] == "benefit") val_a - val_b else val_b - val_a
          if (diff_a_b >= -q[j]) {
            c_a_b <- c_a_b + weights[j]
          } else if (diff_a_b >= -p[j]) {
            c_a_b <- c_a_b + weights[j] * (p[j] + diff_a_b) / (p[j] - q[j])
          }
        }
        rho_a_b[i, h] <- c_a_b
      }
    }
    
    allocations <- rep(1, m)
    for (i in 1:m) {
      allocated <- FALSE
      for (h in num_profiles:1) {
        if (rho_a_b[i, h] >= lambda_t) {
          allocations[i] <- h + 1
          allocated <- TRUE
          break
        }
      }
      if (!allocated) allocations[i] <- 1
    }
    
    # Sort
    sorting_scores <- allocations * 10.0 + rowMeans(rho_a_b)
    ranks <- rank(-sorting_scores, ties.method = "first")
    global_scores <- as.vector(allocations)
    
    extra_data <- list(
      allocations = allocations,
      num_categories = num_profiles + 1,
      profiles = profiles
    )
  }
  
  list(
    weights = weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks),
    version = version,
    extra = extra_data
  )
}

# --- 7. solve_promethee ---

solve_promethee <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  matrix <- matrix(0.0, nrow=m, ncol=n)
  for (i in 1:m) {
    for (j in 1:n) {
      matrix[i, j] <- matrix_data[[i]][[j]]
    }
  }
  
  version <- preference_data$promethee_version
  if (is.null(version)) version <- "II"
  
  scores <- preference_data$weights
  if (is.null(scores)) scores <- list()
  weights <- rep(0.0, n)
  total_score <- 0.0
  for (i in 1:n) {
    cid <- criteria_ids[i]
    score <- if (as.character(cid) %in% names(scores)) as.numeric(scores[[as.character(cid)]]) else 10.0
    weights[i] <- score
    total_score <- total_score + score
  }
  weights <- if (total_score > 0) weights / total_score else rep(1/n, n)
  
  # ranges
  ranges <- rep(0.0, n)
  for (j in 1:n) {
    col <- matrix[, j]
    r <- max(col) - min(col)
    ranges[j] <- if (r > 0) r else 1.0
  }
  
  # preference functions
  function_pref <- preference_data$functions
  if (is.null(function_pref)) function_pref <- list()
  
  fn_types <- rep("usual", n)
  fn_qs <- rep(0.0, n)
  fn_ps <- rep(0.0, n)
  
  for (j in 1:n) {
    cid_str <- as.character(criteria_ids[j])
    if (cid_str %in% names(function_pref)) {
      fn <- function_pref[[cid_str]]
      fn_types[j] <- if ("type" %in% names(fn)) fn$type else "usual"
      fn_qs[j] <- if ("q" %in% names(fn)) as.numeric(fn$q) else 0.0
      fn_ps[j] <- if ("p" %in% names(fn)) as.numeric(fn$p) else ranges[j]
    } else {
      fn_ps[j] <- ranges[j]
    }
  }
  
  get_preference_value <- function(diff, j) {
    if (diff <= 0) return(0.0)
    ftype <- fn_types[j]
    q_val <- fn_qs[j]
    p_val <- fn_ps[j]
    
    if (ftype == "usual") {
      return(1.0)
    } else if (ftype == "linear") {
      if (diff <= q_val) return(0.0)
      if (diff > p_val) return(1.0)
      denom <- p_val - q_val
      return(if (denom > 0) (diff - q_val) / denom else 1.0)
    } else {
      return(if (diff > 0) 1.0 else 0.0)
    }
  }
  
  # Global Preference Index pi(a, b)
  pi <- matrix(0.0, nrow=m, ncol=m)
  for (a in 1:m) {
    for (b in 1:m) {
      if (a == b) next
      p_sum <- 0.0
      for (j in 1:n) {
        val_a <- matrix[a, j]
        val_b <- matrix[b, j]
        diff <- if (criteria_types[j] == "benefit") val_a - val_b else val_b - val_a
        p_sum <- p_sum + weights[j] * get_preference_value(diff, j)
      }
      pi[a, b] <- p_sum
    }
  }
  
  # Flows
  phi_plus <- rep(0.0, m)
  phi_minus <- rep(0.0, m)
  denom <- if (m > 1) as.double(m - 1) else 1.0
  for (i in 1:m) {
    phi_plus[i] <- sum(pi[i, ]) / denom
    phi_minus[i] <- sum(pi[, i]) / denom
  }
  phi <- phi_plus - phi_minus
  
  global_scores <- phi
  ranks <- rep(1, m)
  extra_data <- list()
  
  if (version == "II") {
    ranks <- rank(-phi, ties.method = "first")
  } else if (version == "I") {
    ranks <- rank(-phi, ties.method = "first") # display order
    # Build relations
    relations <- list()
    for (a in 1:m) {
      for (b in 1:m) {
        if (a == b) next
        rel <- "incomparable"
        if ((phi_plus[a] > phi_plus[b] && phi_minus[a] < phi_minus[b]) ||
            (phi_plus[a] == phi_plus[b] && phi_minus[a] < phi_minus[b]) ||
            (phi_plus[a] > phi_plus[b] && phi_minus[a] == phi_minus[b])) {
          rel <- "preferred"
        } else if (phi_plus[a] == phi_plus[b] && phi_minus[a] == phi_minus[b]) {
          rel <- "indifferent"
        }
        relations[[length(relations)+1]] <- list(from=alternatives_ids[a], to=alternatives_ids[b], type=rel)
      }
    }
    extra_data$relations <- relations
  } else if (version == "III") {
    alpha <- if ("alpha_threshold" %in% names(preference_data)) as.numeric(preference_data$alpha_threshold) else 0.05
    ranks <- rank(-phi, ties.method = "first") # basic fallback
  } else if (version == "IV") {
    ranks <- rank(-phi, ties.method = "first")
  } else if (version == "V") {
    costs_pref <- preference_data$costs
    if (is.null(costs_pref)) costs_pref <- list()
    costs <- rep(0.0, m)
    for (i in 1:m) {
      aid_str <- as.character(alternatives_ids[i])
      costs[i] <- if (aid_str %in% names(costs_pref)) as.numeric(costs_pref[[aid_str]]) else 100.0
    }
    
    budget <- if ("budget" %in% names(preference_data)) as.numeric(preference_data$budget) else sum(costs)/2.0
    
    # Solve 0-1 Knapsack using lpSolve
    # Maximize sum(phi_i * x_i)
    # s.t. sum(cost_i * x_i) <= budget
    res <- lpSolve::lp("max", phi, matrix(costs, nrow=1), c("<="), c(budget), all.bin=TRUE)
    portfolio_selection <- rep(0, m)
    opt_indices <- c()
    if (res$status == 0) {
      portfolio_selection <- as.integer(round(res$solution))
      opt_indices <- which(portfolio_selection == 1)
    }
    
    for (i in 1:m) {
      ranks[i] <- if (portfolio_selection[i] == 1) 1 else 2
    }
    
    extra_data <- list(
      costs = costs,
      budget = budget,
      selected_portfolio = alternatives_ids[opt_indices],
      portfolio_selection = portfolio_selection
    )
  } else if (version == "TRI") {
    profiles_list <- preference_data$profiles
    if (is.null(profiles_list) || length(profiles_list) == 0) {
      k_categories <- if ("num_categories" %in% names(preference_data)) as.integer(preference_data$num_categories) else 3
      profiles <- matrix(0.0, nrow=k_categories-1, ncol=n)
      for (k in 1:(k_categories-1)) {
        for (j in 1:n) {
          col <- matrix[, j]
          profiles[k, j] <- min(col) + (max(col) - min(col)) * (k / k_categories)
        }
      }
    } else {
      profiles <- matrix(as.numeric(unlist(profiles_list)), ncol=n, byrow=TRUE)
    }
    
    num_profiles <- nrow(profiles)
    aug_matrix <- rbind(matrix, profiles)
    aug_m <- m + num_profiles
    
    aug_pi <- matrix(0.0, nrow=aug_m, ncol=aug_m)
    for (a in 1:aug_m) {
      for (b in 1:aug_m) {
        if (a == b) next
        p_sum <- 0.0
        for (j in 1:n) {
          val_a <- aug_matrix[a, j]
          val_b <- aug_matrix[b, j]
          diff <- if (criteria_types[j] == "benefit") val_a - val_b else val_b - val_a
          p_sum <- p_sum + weights[j] * get_preference_value(diff, j)
        }
        aug_pi[a, b] <- p_sum
      }
    }
    
    aug_phi <- rep(0.0, aug_m)
    aug_denom <- as.double(aug_m - 1)
    for (i in 1:aug_m) {
      aug_phi[i] <- (sum(aug_pi[i, ]) - sum(aug_pi[, i])) / aug_denom
    }
    
    alt_flows <- aug_phi[1:m]
    profile_flows <- aug_phi[(m+1):aug_m]
    
    allocations <- rep(1, m)
    for (i in 1:m) {
      allocated <- FALSE
      for (h in num_profiles:1) {
        if (alt_flows[i] >= profile_flows[h]) {
          allocations[i] <- h + 1
          allocated <- TRUE
          break
        }
      }
      if (!allocated) allocations[i] <- 1
    }
    
    sorting_scores <- allocations * 10.0 + alt_flows
    ranks <- rank(-sorting_scores, ties.method = "first")
    global_scores <- allocations
    
    extra_data <- list(
      allocations = allocations,
      num_categories = num_profiles + 1,
      profiles = profiles,
      profile_flows = profile_flows,
      alternative_flows = alt_flows
    )
  }
  
  # Normalize matrix
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  for (j in 1:n) {
    col <- matrix[, j]
    c_min <- min(col)
    c_max <- max(col)
    denom <- c_max - c_min
    if (denom == 0) denom <- 1.0
    if (criteria_types[j] == "benefit") {
      norm_matrix[, j] <- (col - c_min) / denom
    } else {
      norm_matrix[, j] <- (c_max - col) / denom
    }
  }
  
  list(
    weights = weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks),
    phi_plus = phi_plus,
    phi_minus = phi_minus,
    phi = phi,
    version = version,
    extra = extra_data
  )
}

# --- 8. solve_smarts_smarter ---

calculate_roc_weights <- function(n) {
  if (n <= 0) return(numeric())
  weights <- rep(0.0, n)
  for (i in 1:n) {
    w_i <- (1.0 / n) * sum(1.0 / (i:n))
    weights[i] <- w_i
  }
  return(weights)
}

normalize_linear <- function(mat, criteria_types) {
  m <- nrow(mat)
  n <- ncol(mat)
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  for (j in 1:n) {
    col <- mat[, j]
    col_max <- max(col)
    col_min <- min(col)
    denom <- col_max - col_min
    if (denom == 0) denom <- 1.0
    if (criteria_types[j] == "benefit") {
      norm_matrix[, j] <- (col - col_min) / denom
    } else {
      norm_matrix[, j] <- (col_max - col) / denom
    }
  }
  return(norm_matrix)
}

solve_smarts_smarter <- function(matrix_data, criteria_types, preference_data, criteria_ids, method) {
  m <- length(matrix_data)
  n <- length(criteria_ids)
  
  matrix <- matrix(0.0, nrow=m, ncol=n)
  for (i in 1:m) {
    for (j in 1:n) {
      matrix[i, j] <- matrix_data[[i]][[j]]
    }
  }
  
  weights <- rep(0.0, n)
  if (method == "smarts") {
    scores <- preference_data$weights
    if (is.null(scores)) scores <- list()
    total_score <- 0.0
    for (i in 1:n) {
      cid <- criteria_ids[i]
      score <- if (as.character(cid) %in% names(scores)) as.numeric(scores[[as.character(cid)]]) else 10.0
      weights[i] <- score
      total_score <- total_score + score
    }
    weights <- if (total_score > 0) weights / total_score else rep(1/n, n)
  } else if (method == "smarter") {
    ranks_order <- unlist(preference_data$ranks)
    roc_w <- calculate_roc_weights(n)
    for (i in 1:n) {
      cid <- criteria_ids[i]
      rank_idx <- which(ranks_order == cid)
      if (length(rank_idx) == 0) rank_idx <- n
      weights[i] <- roc_w[rank_idx[1]]
    }
  } else {
    weights <- rep(1/n, n)
  }
  
  norm_matrix <- normalize_linear(matrix, criteria_types)
  global_scores <- as.vector(norm_matrix %*% weights)
  ranks <- rank(-global_scores, ties.method = "first")
  
  list(
    weights = weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = global_scores,
    ranks = as.integer(ranks)
  )
}

# --- 9. solve_topsis ---

solve_topsis <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  matrix <- matrix(0.0, nrow=m, ncol=n)
  for (i in 1:m) {
    for (j in 1:n) {
      matrix[i, j] <- matrix_data[[i]][[j]]
    }
  }
  
  # weights
  scores <- preference_data$weights
  if (is.null(scores)) scores <- list()
  weights <- rep(0.0, n)
  total_score <- 0.0
  for (i in 1:n) {
    cid <- criteria_ids[i]
    score <- if (as.character(cid) %in% names(scores)) as.numeric(scores[[as.character(cid)]]) else 10.0
    weights[i] <- score
    total_score <- total_score + score
  }
  weights <- if (total_score > 0) weights / total_score else rep(1/n, n)
  
  # Vector Normalization
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  for (j in 1:n) {
    col <- matrix[, j]
    norm_factor <- sqrt(sum(col ^ 2))
    if (norm_factor == 0) norm_factor <- 1.0
    norm_matrix[, j] <- col / norm_factor
  }
  
  # Weighted matrix
  weighted_matrix <- norm_matrix
  for (j in 1:n) {
    weighted_matrix[, j] <- weighted_matrix[, j] * weights[j]
  }
  
  # Ideal solutions
  ideal_positive <- rep(0.0, n)
  ideal_negative <- rep(0.0, n)
  for (j in 1:n) {
    col <- weighted_matrix[, j]
    if (criteria_types[j] == "benefit") {
      ideal_positive[j] <- max(col)
      ideal_negative[j] <- min(col)
    } else {
      ideal_positive[j] <- min(col)
      ideal_negative[j] <- max(col)
    }
  }
  
  # Distances
  distance_positive <- rep(0.0, m)
  distance_negative <- rep(0.0, m)
  for (i in 1:m) {
    row <- weighted_matrix[i, ]
    distance_positive[i] <- sqrt(sum((row - ideal_positive) ^ 2))
    distance_negative[i] <- sqrt(sum((row - ideal_negative) ^ 2))
  }
  
  # Closeness
  closeness <- rep(0.0, m)
  for (i in 1:m) {
    denom <- distance_positive[i] + distance_negative[i]
    closeness[i] <- if (denom == 0) 0.5 else distance_negative[i] / denom
  }
  
  ranks <- rank(-closeness, ties.method = "first")
  
  list(
    weights = weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = closeness,
    ranks = as.integer(ranks),
    distance_positive = distance_positive,
    distance_negative = distance_negative
  )
}

# --- 10. solve_vikor ---

solve_vikor <- function(matrix_data, criteria_types, preference_data, criteria_ids, alternatives_ids) {
  m <- length(alternatives_ids)
  n <- length(criteria_ids)
  
  matrix <- matrix(0.0, nrow=m, ncol=n)
  for (i in 1:m) {
    for (j in 1:n) {
      matrix[i, j] <- matrix_data[[i]][[j]]
    }
  }
  
  # weights
  scores <- preference_data$weights
  if (is.null(scores)) scores <- list()
  weights <- rep(0.0, n)
  total_score <- 0.0
  for (i in 1:n) {
    cid <- criteria_ids[i]
    score <- if (as.character(cid) %in% names(scores)) as.numeric(scores[[as.character(cid)]]) else 10.0
    weights[i] <- score
    total_score <- total_score + score
  }
  weights <- if (total_score > 0) weights / total_score else rep(1/n, n)
  
  v <- if ("v" %in% names(preference_data)) as.numeric(preference_data$v) else 0.5
  
  f_star <- rep(0.0, n)
  f_minus <- rep(0.0, n)
  for (j in 1:n) {
    col <- matrix[, j]
    if (criteria_types[j] == "benefit") {
      f_star[j] <- max(col)
      f_minus[j] <- min(col)
    } else {
      f_star[j] <- min(col)
      f_minus[j] <- max(col)
    }
  }
  
  S <- rep(0.0, m)
  R <- rep(0.0, m)
  for (i in 1:m) {
    terms <- rep(0.0, n)
    for (j in 1:n) {
      denom <- f_star[j] - f_minus[j]
      terms[j] <- if (denom == 0) 0.0 else weights[j] * (f_star[j] - matrix[i, j]) / denom
    }
    S[i] <- sum(terms)
    R[i] <- max(terms)
  }
  
  S_star <- min(S)
  S_minus <- max(S)
  R_star <- min(R)
  R_minus <- max(R)
  
  Q <- rep(0.0, m)
  for (i in 1:m) {
    s_term <- if (S_minus - S_star > 0) (S[i] - S_star) / (S_minus - S_star) else 0.0
    r_term <- if (R_minus - R_star > 0) (R[i] - R_star) / (R_minus - R_star) else 0.0
    Q[i] <- v * s_term + (1 - v) * r_term
  }
  
  ranks <- rank(Q, ties.method = "first")
  
  # norm matrix for viewing
  norm_matrix <- matrix(0.0, nrow=m, ncol=n)
  for (j in 1:n) {
    denom <- f_star[j] - f_minus[j]
    if (denom == 0) {
      norm_matrix[, j] <- 1.0
    } else {
      if (criteria_types[j] == "benefit") {
        norm_matrix[, j] <- (matrix[, j] - f_minus[j]) / denom
      } else {
        norm_matrix[, j] <- (f_minus[j] - matrix[, j]) / denom
      }
    }
  }
  
  list(
    weights = weights,
    normalized_matrix = lapply(seq_len(nrow(norm_matrix)), function(i) norm_matrix[i, ]),
    global_scores = Q,
    ranks = as.integer(ranks),
    S = S,
    R = R,
    v = v
  )
}
