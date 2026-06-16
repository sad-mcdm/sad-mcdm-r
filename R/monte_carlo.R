# Monte Carlo simulation for sadMCDM package

run_monte_carlo <- function(matrix_data, criteria_types, weights, variations_pct, num_simulations, method, preference_data, criteria_ids) {
  # Convert matrix_data list of lists to R matrix
  m <- length(matrix_data)
  n <- length(criteria_ids)
  mat <- matrix(0.0, nrow=m, ncol=n)
  for (i in 1:m) {
    for (j in 1:n) {
      mat[i, j] <- matrix_data[[i]][[j]]
    }
  }
  
  col_maxs <- apply(mat, 2, max)
  col_mins <- apply(mat, 2, min)
  ranges <- col_maxs - col_mins
  deltas <- (variations_pct / 100.0) * ranges
  
  # Pre-extract helper structures for bisection/macbeth
  bisection_mids <- rep(0.0, n)
  macbeth_level_scores <- list()
  
  if (method == "bwt") {
    mids <- preference_data$bisection_midpoints
    if (is.null(mids)) mids <- list()
    for (j in 1:n) {
      cid_str <- as.character(criteria_ids[j])
      default_mid <- (col_mins[j] + col_maxs[j]) / 2.0
      val <- if (cid_str %in% names(mids)) as.numeric(mids[[cid_str]]) else default_mid
      bisection_mids[j] <- if (is.na(val)) default_mid else val
    }
  } else if (method == "macbeth") {
    levels_pref <- preference_data$levels_matrices
    if (is.null(levels_pref)) levels_pref <- list()
    for (j in 1:n) {
      cid_str <- as.character(criteria_ids[j])
      is_benefit <- (criteria_types[j] == "benefit")
      if (cid_str %in% names(levels_pref)) {
        levels_matrix <- matrix(as.numeric(unlist(levels_pref[[cid_str]])), nrow=5, ncol=5, byrow=TRUE)
        macbeth_level_scores[[j]] <- get_macbeth_scores(levels_matrix, 5, is_benefit)
      } else {
        macbeth_level_scores[[j]] <- if (is_benefit) seq(0, 1, length.out=5) else seq(1, 0, length.out=5)
      }
    }
  }
  
  winner_counts <- rep(0, m)
  rank_counts <- matrix(0, nrow=m, ncol=m)
  global_scores_sum <- rep(0.0, m)
  
  for (s in 1:num_simulations) {
    # Generate uniform noise and perturb matrix
    perturbed <- mat
    for (j in 1:n) {
      noise <- runif(m, min = -deltas[j], max = deltas[j])
      perturbed[, j] <- mat[, j] + noise
      # Clip
      perturbed[, j] <- pmax(col_mins[j], pmin(col_maxs[j], perturbed[, j]))
    }
    
    # Solve based on method
    if (method %in% c("topsis", "vikor", "electre", "promethee")) {
      # Reconstruct matrix_data as list of lists
      perturbed_list <- lapply(1:m, function(r) as.list(perturbed[r, ]))
      
      # Reconstruct weights dict
      weights_dict <- list()
      for (j in 1:n) {
        weights_dict[[as.character(criteria_ids[j])]] <- weights[j] * 100.0
      }
      pref_copy <- preference_data
      pref_copy$weights <- weights_dict
      
      if (method == "topsis") {
        res <- solve_topsis(perturbed_list, criteria_types, pref_copy, criteria_ids, seq_len(m))
      } else if (method == "vikor") {
        res <- solve_vikor(perturbed_list, criteria_types, pref_copy, criteria_ids, seq_len(m))
      } else if (method == "electre") {
        res <- solve_electre(perturbed_list, criteria_types, pref_copy, criteria_ids, seq_len(m))
      } else { # promethee
        res <- solve_promethee(perturbed_list, criteria_types, pref_copy, criteria_ids, seq_len(m))
      }
      
      scores <- res$global_scores
      ranks_iter <- res$ranks
      
    } else {
      norm_perturbed <- matrix(0.0, nrow=m, ncol=n)
      
      if (method == "bwt") {
        for (j in 1:n) {
          is_benefit <- (criteria_types[j] == "benefit")
          x_05 <- bisection_mids[j]
          for (i in 1:m) {
            norm_perturbed[i, j] <- interpolate_bisection(perturbed[i, j], col_mins[j], col_maxs[j], x_05, is_benefit)
          }
        }
      } else if (method == "macbeth") {
        for (j in 1:n) {
          levels <- seq(col_mins[j], col_maxs[j], length.out=5)
          scores_levels <- macbeth_level_scores[[j]]
          for (i in 1:m) {
            if (col_maxs[j] == col_mins[j]) {
              norm_perturbed[i, j] <- 1.0
            } else {
              norm_perturbed[i, j] <- approx(levels, scores_levels, xout=perturbed[i, j], rule=2)$y
            }
          }
        }
      } else {
        norm_perturbed <- normalize_linear(perturbed, criteria_types)
      }
      
      scores <- as.vector(norm_perturbed %*% weights)
      ranks_iter <- rank(-scores, ties.method = "first")
    }
    
    global_scores_sum <- global_scores_sum + scores
    for (i in 1:m) {
      r_pos <- ranks_iter[i]
      r_pos <- max(1, min(m, r_pos))
      rank_counts[i, r_pos] <- rank_counts[i, r_pos] + 1
      if (r_pos == 1) {
        winner_counts[i] <- winner_counts[i] + 1
      }
    }
  }
  
  probabilities_first <- winner_counts / num_simulations
  rank_probabilities <- rank_counts / num_simulations
  average_scores <- global_scores_sum / num_simulations
  
  average_ranks <- rep(0.0, m)
  for (i in 1:m) {
    avg_r <- sum((1:m) * rank_counts[i, ]) / num_simulations
    average_ranks[i] <- avg_r
  }
  
  list(
    first_place_probabilities = as.vector(probabilities_first),
    rank_probabilities = lapply(1:m, function(i) rank_probabilities[i, ]),
    average_scores = as.vector(average_scores),
    average_ranks = as.vector(average_ranks),
    deltas = as.vector(deltas)
  )
}
