source("code/genus_process.R")


# significance
# for the comparison of 3 or more groups use the kruskall-wallis test, for 2 groups the wilcoxon test
  
sig_genera <-
  composite %>%  
  nest(data = -taxonomy) %>% 
  mutate(test = map(.x = data, 
                    ~ wilcox.test(rel_abund ~ srn, data = .x) %>% tidy())) %>% 
  unnest(test) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH")) %>% 
  filter(p.adjust < 0.05) %>% 
  select(taxonomy, p.adjust)
  
# iN CASE NEEDED, THE LIMIT OF DETECTION
# calculate the seqs per sample and then the limit of detection (LOD)
nseqs_per_sample <-
  shared %>% 
  group_by(group) %>% 
  summarize(N = sum(count), .groups = "drop") %>% 
  count(N) %>% # all the samples should have the same sampling depth
  pull(N)

stopifnot(length(nseqs_per_sample) == 1)

lod <- 100*1/nseqs_per_sample
# the obtained number is the lod for this dataset. 
# One use you could do with it, is to replace the zeros with  2/3*lod, then, you can log transform the data without loosing shape of the data as it's below the lod.


composite %>% 
  inner_join(sig_genera, by = "taxonomy") %>% 
  mutate(rel_abund = 100 * (rel_abund + lod/200),
         taxonomy = str_replace(taxonomy, "(.*)", "*\\1*"),
         taxonomy = str_replace(taxonomy, 
                                "\\*(.*)_unclassified\\*", 
                                "Unclassified<br>*\\1*"),
         srn = factor(srn, levels = c(T, F))) %>% 
  ggplot(aes(x = rel_abund, y = taxonomy, color = srn, fill = srn)) +
  geom_vline(xintercept = lod, linetype = "dashed", color = "gray 80") +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8,
                                              jitter.width = 0.3), 
              alpha = 0.75,
              shape = 21) +
  stat_summary(fun.data = median_hilow,
               fun.args = list(conf.int = 0.5),
               geom = "pointrange",
               position = position_dodge(width = 0.8),
               color = "black",
               show.legend = FALSE) +
  scale_x_log10() +
  scale_color_manual(NULL, 
                     breaks = c(F, T),
                     values = c("gray", "steelblue"),
                     labels = c("Healthy", "SRN")) +
  scale_fill_manual(NULL, 
                     breaks = c(F, T),
                     values = c("gray", "steelblue"),
                     labels = c("Healthy", "SRN")) +
  labs(x = "Relative abundance (%)",
       y = NULL) +
  theme_classic() +
  theme(axis.text.y = element_markdown())

ggsave("figures/significant_genera.png", width = 6, height = 5)




test <-
  composite %>% 
  inner_join(sig_genera, by = "taxonomy") %>% 
  select(group, taxonomy, rel_abund, fit_result, srn) %>% 
  pivot_wider(names_from = taxonomy, values_from = rel_abund) %>% 
  pivot_longer(cols = -c(group, srn), names_to =  "metric", values_to = "score") %>% 
  filter(metric == "fit_result")


get_sens_spec <-
  function(treshold, score, actual, direction) {
    
    # treshold <- 100
    # score <- test$score
    # actual <- test$srn
    # direction <- "greaterthan"
    
    predicted <- if(direction == "greaterthan") {score > treshold} else{score < treshold}
    
    tp <- sum(predicted & actual)
    tn <- sum(!predicted & !actual)
    fp <- sum(predicted & !actual)
    fn <- sum(!predicted & actual)
    
    specificity <- tn / (tn + fp)
    sensitivity <- tp / (tp + fn)
    
    tibble("specificity" = specificity,
           "sensitivity" = sensitivity)
  }
    


get_roc_data <- 
  function(x, direction) {
  
  # x <- test
  # direction <- "greaterthan"
  
  tresholds = unique(x$score) %>%  sort()
  
  map_dfr(.x = tresholds, ~get_sens_spec(.x, x$score, x$srn, direction))# %>% 
    #rbind(c(specificity = 0, sensitivity = 1))
  
}


# get_sens_spec(100, test$score, test$srn, "greaterthan")
# get_roc_data(test, "greaterthan")


roc_data <-
  composite %>% 
  inner_join(sig_genera, by = "taxonomy") %>% 
  select(group, taxonomy, rel_abund, fit_result, srn) %>% 
  pivot_wider(names_from = taxonomy, values_from = rel_abund) %>% 
  pivot_longer(cols = -c(group, srn), names_to =  "metric", values_to = "score") %>% 
  # filter(metric == "fit_result") %>% 
  nest(data = -metric) %>% 
  mutate(direction = ifelse(metric == "Lachnospiraceae_unclassified",
                            "lessthan", "greaterthan")) %>% 
  mutate(roc_data = map2(.x = data, .y = direction, ~get_roc_data(.x, .y))) %>% 
  unnest(roc_data) %>% 
  select(metric, specificity, sensitivity)

roc_data %>% 
  ggplot(aes(x = 1-specificity, y= sensitivity, color = metric)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, color = "gray75", linetype = "dashed") +
  theme_classic() +
  coord_fixed()

ggsave("figures/roc_figure.png", height = 4, width = 6)




