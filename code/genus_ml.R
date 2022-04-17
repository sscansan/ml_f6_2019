source("code/genus_process.R")

library(mikropml)


srn_genus_data <-
  composite %>% 
  select(group, taxonomy, rel_abund, srn) %>% 
  pivot_wider(names_from = taxonomy, values_from = rel_abund) %>% 
  select(-group) %>% 
  mutate(srn = if_else(srn, "srn", "healthy")) %>% 
  select(srn, everything())

srn_genus_results <-
  run_ml(srn_genus_data,
       method = "glmnet", 
       outcome_colname = "srn",
       kfold = 5, #k splits
       cv_folds = 100, #ctoss validations
       training_frac = 0.8, # 80/20 Split
       seed = 14091989)

srn_genus_results
