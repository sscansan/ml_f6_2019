source("code/genus_process.R")

library(mikropml)
library(tictoc)
library(furrr)

# plan("sequential") #not parallel
# plan("multicore") # does not work with windows or RStudio
plan("multisession")

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

test_hp <- list(alpha = c(0),
                lambda = c(0.1, 0.5, 1, 2, 3, 4, 5, 10)) 

get_srn_genus_results <- function(seed) {
  run_ml(srn_genus_preprocess,
       method = "glmnet", 
       outcome_colname = "srn",
       kfold = 5, #k splits
       cv_times = 100, #ctoss validations
       training_frac = 0.8, # 80/20 Split
       hyperparameters = test_hp,
       seed = seed)
}

# Hyperparametrers adjustment
# glmnet has two hyperparameters, lambda and alpha.
# for different values of lambda (weighing and regularisation)
# alpha = 0 - ridge regression (L2)
# 0 < alpha < 1 - elastic net 
# alpha = 1 - lasso regression (L1)
 

tic()
iterative_run_ml_result <- future_map(1:100, get_srn_genus_results,
                                      .options = furrr_options(seed = TRUE))
toc()

performance <- 
  iterative_run_ml_result %>% 
  map(pluck, "trained_model") %>% 
  combine_hp_performance()

plot_hp_performance(performance$dat, lambda, AUC)

performance$dat %>% 
  group_by(alpha, lambda) %>% 
  summarize(mean_AUC = mean(AUC),
            lquartile = quantile(AUC, prob = 0.25),
            uquartile = quantile(AUC, prob = 0.75),
            .groups = "drop") %>%
  # top_n(n=3, mean_AUC)
  ggplot(aes(x = lambda, y = mean_AUC, color = as.character(alpha))) +
  geom_line()


plan("sequential")
