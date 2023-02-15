require(tidyverse)

plot_fuzzy_scores <- function(interval_matrix, sorted_sol, fn_name){
  
  pfns <- interval_matrix %>%
    as_tibble(rownames="Solution") %>%
    select(Solution, LB, CORE, UB) %>%
    mutate(Solution = factor(Solution, levels=rev(sorted_sol))) %>%
    pivot_longer(cols = -Solution, names_to = "Point", values_to = "Score") %>%
    mutate(
      MF = ifelse(Point == "CORE", 1, 0)
    ) %>% mutate(x = Score, y=as.numeric(Solution) + ifelse(Point == "CORE", 1, 0)) %>% 
    ggplot(aes(x = Score, y = Solution)) +
    geom_polygon(aes(y=y, fill = Solution, group = Solution), color="black") + 
    scale_y_continuous(breaks = seq(1,nsol), labels = rev(sorted_sol)) + 
    theme(legend.position = "none") +
    ggtitle(fn_name) +
    scale_fill_viridis_d(option = "C")
  
  return(pfns)
  
}


plot_intervals <- function(assess_matrix, weights_labels){
  
  cc <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=100))
  
  porig <- assess_matrix %>%
    as_tibble(rownames="Solution") %>%
    mutate(Solution = factor(Solution, levels=rev(rownames(assess_matrix)))) %>%
    rename_with(starts_with("VE"), .fn = ~ weights_labels) %>%
    select(Solution, contains("("), LB, UB) %>%
    pivot_longer(cols = contains("("), names_to = "Weight", values_to = "Score") %>%
    mutate(Weight = factor(Weight, levels=weights_labels)) %>%
    ggplot(aes(x=LB, y=Solution)) +
    geom_segment(aes(xend=UB, yend=Solution), linewidth=0.5, linetype="solid") +
    geom_point(aes(x=Score, fill=Weight), shape=21, stroke=0.9, size=2) +
    xlab("Score") +
    scale_fill_viridis_d(option = "C")
    #scale_fill_brewer(palette = "YlOrRd")
  
  print(porig)
}
