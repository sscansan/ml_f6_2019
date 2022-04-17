source("code/genus_process.R")

library(mikropml)


srn_genus_data <-
  composite %>% 
  select(group, taxonomy, rel_abund, srn) %>% 
  pivot_wider(names_from = taxonomy, values_from = rel_abund) %>% 
  select(-group) %>% 
  mutate(srn = if_else(srn, "srn", "healthy")) %>% 
  select(srn, everything())

srn_genus_nopp <-
  run_ml(srn_genus_data,
         method = "glmnet", 
         outcome_colname = "srn",
         kfold = 5, #k splits
         cv_times = 100, #ctoss validations
         training_frac = 0.8, # 80/20 Split
         seed = 14091989)

# Preprocessing and re run
srn_genus_preprocess <- 
  preprocess_data(srn_genus_data, outcome_colname = "srn")$dat_transformed

srn_genus_results <-
  run_ml(srn_genus_preprocess,
       method = "glmnet", 
       outcome_colname = "srn",
       kfold = 5, #k splits
       cv_times = 100, #ctoss validations
       training_frac = 0.8, # 80/20 Split
       seed = 14091989)

# Hyperparametrers adjustment



