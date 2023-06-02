rm(list = ls())
library(tidyverse)
library(patchwork)
source("aw_methods.R")
source("fuzzy_analysis.R")
source("plots.R")
library(patchwork)
library(GGally)

eval_matrix <- read.table(file = "matrix_aircrafts.csv") %>% as.matrix()
nsol <- nrow(eval_matrix)
ncrit <- ncol(eval_matrix)

sol_names <- paste0("S", 1:nsol)
crit_names <- paste0("C", 1:ncrit)

rownames(eval_matrix) <- sol_names
colnames(eval_matrix) <- crit_names




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

# Interval center score

IC <- apply(score_interval_matrix, MARGIN=1, FUN=function(x) mean(x))

# Mean weights score 

comb_weight_score_matrix <- cbind(vert_eval_matrix[,-ncrit], score_by_weights_setting_matrix)

#MWS <- apply(comb_weight_score_matrix, MARGIN = 1, function(sol_w_s) mean(sol_w_s))

MWS <- apply(cbind(score_by_weights_setting_matrix, IC), MARGIN = 1, function(sol_w_s) mean(sol_w_s))

core_matrix <- cbind(IC, score_by_weights_setting_matrix, MWS)


crisp_ranked_solution_matrix <- apply(core_matrix, MARGIN = 2, function(CORE){
  
  names(sort(rank(-CORE)))
  
})



fuzzy_ranked_solution_matrix <- apply(core_matrix, MARGIN = 2, function(CORE){
  fn_def_matrix <- cbind(score_interval_matrix, CORE)
  fuzzy_solution_sort(fn_def_matrix)
})

colnames(fuzzy_ranked_solution_matrix) <- paste0("F", colnames(core_matrix))

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

corr_matrix <- sapply(df_rank_matrix %>% select(starts_with("F")), function(m1){
  sapply(df_rank_matrix %>% select(starts_with("F")), function(m2){
    cor(m1, m2, method = "kendall")
  })
})

newnames <- str_replace(colnames(corr_matrix), "F", "")
colnames(corr_matrix) <- newnames
rownames(corr_matrix) <- newnames


df_label_rank_matrix <- as_tibble(rank_label_matrix) %>%
  mutate(across(everything(), ~ factor(.x, levels=sol_names))) %>%
  mutate(Rank = 1:nsol)

topn <- 5

match_matrix <- sapply(df_label_rank_matrix %>% select(-Rank) %>% select(starts_with("F")), function(m1){
  sapply(df_label_rank_matrix %>% select(-Rank) %>% select(starts_with("F")), function(m2){
    length(intersect(m1[1:topn], m2[1:topn]))/topn
  })
})
colnames(match_matrix) <- newnames
rownames(match_matrix) <- newnames

# Plots

library(corrplot)

p1 <- ggcorr(df_rank_matrix %>% 
               select(starts_with("F")) %>% 
               rename_with(~str_replace(.x, "F", ""), everything()), 
             method = c("pairwise", "kendall"), 
             low = "#DBDBDB", mid = "#5C5E60", high = "#DBDBDB", 
             label = TRUE, 
             label_round = 3) + 
  ggtitle("a) Correlation") + 
  theme(title = element_text(face = "bold", size = 12))

#p1 <- p1 + scale_fill_viridis_c(option = "C", )

# ggsave(filename = "corr_plot.pdf", p1, width = 5, height = 5)

p2 <- ggcorr(data = NULL, cor_matrix = match_matrix, limits = FALSE,
             label = TRUE, 
             low = "#212016", mid = "#5C5E60", high = "#DBDBDB", 
             label_round = 3) + 
  ggtitle("b) Matching rate for the top 5 solutions") + 
  theme(title = element_text(face = "bold", size = 12))

ggsave(filename = "corr_match_plots.pdf", p1 + p2, width = 9, height = 4)

f_names <- colnames(core_matrix)
letter_vec <- letters[1:length(f_names)]
letter_vec <- paste0(letter_vec, ") ", f_names)
names(letter_vec) <- f_names

# fuzzy_plots <- lapply(f_names, function(core_name){
#   sort_sols <- df_label_rank_matrix[, paste0("F", core_name)] %>% as_vector()
#   CORE <- core_matrix[,core_name]
#   m <- cbind(score_interval_matrix, CORE)
#   plot_fuzzy_scores(m, sort_sols, letter_vec[core_name])
# })

fuzzy_plots <- lapply(f_names, function(core_name){
  sort_sols <- sol_names #df_label_rank_matrix[, paste0("F", core_name)] %>% as_vector()
  CORE <- core_matrix[,core_name]
  m <- cbind(score_interval_matrix, CORE)
  plot_fuzzy_scores(m, sort_sols, letter_vec[core_name])
})


p_all <- wrap_plots(fuzzy_plots, ncol = 3)
ggsave(filename = "plot_fuzzy.pdf", plot = p_all, width = 7, height = 9)

fuzzy_plots <- lapply(f_names, function(core_name){
  sort_sols <- sol_names #df_label_rank_matrix[, paste0("F", core_name)] %>% as_vector()
  CORE <- core_matrix[,core_name]
  Rank <- df_label_rank_matrix[, paste0("F", core_name)] %>% as_vector()
  m <- cbind(score_interval_matrix, CORE, Rank)
  plot_fuzzy_scores_label(m, sort_sols, letter_vec[core_name])
})


p_all <- wrap_plots(fuzzy_plots, ncol = 3)
ggsave(filename = "plot_fuzzy_labels.pdf", plot = p_all, width = 7, height = 9)



weights_labels <- apply(vert_matrix, MARGIN = 2, function(c) paste0("(", paste(round(c,2), collapse = ", "), ")", collapse = ""))
colnames(vert_eval_matrix) <- paste0("VE", 1:ncrit)

p_int <- plot_intervals_and_weights(cbind(vert_eval_matrix, score_interval_matrix), weights_labels) + 
  ggtitle("b) Distribution of extreme weights") + 
  theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 9))  + 
  guides(fill=guide_legend(nrow=3, byrow=FALSE))


weights_labels <- sapply(weights_settings_list, function(c) paste0("(", paste(round(c,2), collapse = ", "), ")", collapse = ""))
weights_labels <- paste0(names(weights_settings_list), "=", weights_labels)

approx_weights_matrix <- cbind(score_by_weights_setting_matrix, score_interval_matrix)

colnames(approx_weights_matrix) <- c(paste0("VE", 1:length(weights_labels)) , "LB", "UB")

p_aw <- plot_intervals_and_weights(approx_weights_matrix, weights_labels) + 
  ggtitle("c) Distribution of approximated weights")  +
  theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 9))  + 
  guides(fill=guide_legend(nrow=4, byrow=FALSE))



p_int_only <- plot_intervals_only(score_interval_matrix) + 
  ggtitle("a) Score intervals")

#p_int_all <- p_int_only / p_int / p_aw

p_int_all <- p_int_only / (p_int | p_aw)

ggsave(filename = "plot_intervals_all.pdf", plot = p_int_all, width = 9, height = 9)