tmp <- readRDS("./results/data/alpha_2_3_15/scenario1_misspecified_alpha_2_3_15.rds")

p1 <- ggplot(tmp, aes(x=n_counts, y=es_correct, color=method)) + 
  stat_summary() + 
  scale_x_log10()

p2 <- ggplot(tmp, aes(x=n_counts, y=es_incorrect, color=method)) + 
  stat_summary() + 
  scale_x_log10()

grid.arrange(p1, p2)



tmp2 <- readRDS("./results/data/alpha_2_3_15/scenario2_narrow_calibrated_alpha_2_3_15.rds")

p3 <- ggplot(tmp2, aes(x=n_counts, y=es_correct, color=method)) + 
  stat_summary() + 
  scale_x_log10()

p4 <- ggplot(tmp2, aes(x=n_counts, y=es_incorrect, color=method)) + 
  stat_summary() + 
  scale_x_log10()

grid.arrange(p3, p4)


tmp3 <- readRDS("./results/data/alpha_2_3_15/scenario3_dispersed_calibrated_alpha_2_3_15.rds")

p5 <- ggplot(tmp3, aes(x=n_counts, y=es_correct, color=method)) + 
  stat_summary() + 
  scale_x_log10()

p6 <- ggplot(tmp3, aes(x=n_counts, y=es_incorrect, color=method)) + 
  stat_summary() + 
  scale_x_log10()

grid.arrange(p5, p6)


tmp4 <- readRDS("./results/data/alpha_2_3_15/scenario4_less_misspecified_alpha_2_3_15.rds")

p7 <- ggplot(tmp4, aes(x=n_counts, y=es_correct, color=method)) + 
  stat_summary() + 
  scale_x_log10()

p8 <- ggplot(tmp4, aes(x=n_counts, y=es_incorrect, color=method)) + 
  stat_summary() + 
  scale_x_log10()

grid.arrange(p7, p8)