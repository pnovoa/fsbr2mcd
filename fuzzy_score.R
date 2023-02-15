rm(list = ls())
library(tidyverse)
library(patchwork)
source("aw_methods.R")
source("fuzzy_analysis.R")

# Creating case study data

set.seed(1321)

nsol <- 897
ncrit <- 4
max_eval <- 10.

sol_names <- paste0("S", 1:nsol)
crit_names <- paste0("C", 1:ncrit)

eval_matrix = matrix(runif(nsol*ncrit)*10., byrow = TRUE, nrow = nsol)

rownames(eval_matrix) <- sol_names
colnames(eval_matrix) <- crit_names

# Compute the score intervals

generate_polyhedron_vertices <- function(ncrit) {
  m <- matrix(
    rep(1 / 1:ncrit, ncrit),
    nrow = ncrit,
    byrow = T
  )
  
  m[lower.tri(m)] <- 0.
  
  return(m)
}

vert_matrix <- generate_polyhedron_vertices(ncrit = ncrit)

vert_eval_matrix <- eval_matrix %*% vert_matrix

score_interval_matrix <- t(apply(vert_eval_matrix, MARGIN = 1, function(sol_evals) c(min(sol_evals), max(sol_evals))))

colnames(score_interval_matrix) <- c("LB", "UB")

# Applying weight approximation methods

## Computing the weights

weights_settings_list <- list(
  "EW" = ew_weights(ncrit),
  "ROC" = roc_weights(ncrit),
  "RR" = rr_weights(ncrit),
  "RS" = rs_weights(ncrit)
)

score_by_weights_setting_matrix <- sapply(weights_settings_list, function(w_setting) eval_matrix %*% w_setting)


# Mean weights score 

comb_weight_score_matrix <- cbind(vert_eval_matrix[,-ncrit], score_by_weights_setting_matrix)

MWS <- apply(comb_weight_score_matrix, MARGIN = 1, function(sol_w_s) mean(sol_w_s))

MWS2 <- apply(score_by_weights_setting_matrix, MARGIN = 1, function(sol_w_s) mean(sol_w_s))

core_matrix <- cbind(score_by_weights_setting_matrix, MWS, MWS2)


crisp_ranked_solution_matrix <- apply(core_matrix, MARGIN = 2, function(CORE){
  
  names(sort(rank(-CORE)))
  
})



fuzzy_ranked_solution_matrix <- apply(core_matrix, MARGIN = 2, function(CORE){
  fn_def_matrix <- cbind(score_interval_matrix, CORE)
  fuzzy_solution_sort(fn_def_matrix)
})

colnames(fuzzy_ranked_solution_matrix) <- paste0("F_", colnames(core_matrix))

rank_label_matrix <- cbind(crisp_ranked_solution_matrix, fuzzy_ranked_solution_matrix)

rank_matrix <- apply(rank_label_matrix, MARGIN = 2, function(x){
  tb_rnk <- tibble(
    Rank = 1:nsol,
    Label = factor(x, levels = sol_names)
  ) %>% arrange(Label)
  
  tb_rnk$Rank
  
})

df_rank_matrix <- as_tibble(rank_matrix) %>%
  mutate(Solution = factor(sol_names, levels=sol_names))

corr_matrix <- sapply(df_rank_matrix %>% select(-Solution), function(m1){
  sapply(df_rank_matrix %>% select(-Solution), function(m2){
    cor(m1, m2, method = "kendall")
  })
})


df_label_rank_matrix <- as_tibble(rank_label_matrix) %>%
  mutate(across(everything(), ~ factor(.x, levels=sol_names))) %>%
  mutate(Rank = 1:nsol)

topn <- 50

match_matrix <- sapply(df_label_rank_matrix %>% select(-Rank), function(m1){
  sapply(df_label_rank_matrix %>% select(-Rank), function(m2){
    length(intersect(m1[1:topn], m2[1:topn]))/topn
  })
})


# Plots

library(corrplot)

p1 <- corrplot(corr_matrix)

p2 <- corrplot(match_matrix)

print(p1)

print(p2)
