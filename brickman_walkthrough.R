source("brickman_walkthrough_help.R")

set.seed(607)

### DATA PREPARATION

# INITIAL SF OBJECT - SUBSTITUTE YOUR PRESENCE/ABSENCE DATA HERE
#   PRESENCE: factor (0 or 1)
#   MONTH: numeric (1 through 12)

library(ecomon) # for example presence data - delete if not using ecomon dataset

pa_data <- ecomon::read_staged(species = "calfin", form = "sf") |>
  transmute(PRESENCE = (total_m2 > 10000) |> as.numeric() |> as.factor(),
            MONTH = lubridate::month(date))

# COVARIATES
VARS <- c("Bathy_depth", "Xbtm", "MLD", "Sbtm", "SSS", "SST", "Tbtm", "U", "V")

# pairing brickman present data to presence data - this is the input dataset
model_data <- brickman::extract_points(brickman::compose_filename("PRESENT"), 
                                       vars = VARS, 
                                       pts = pa_data, 
                                       complete = TRUE,
                                       simplify_names = TRUE) |>
  bind_cols(pa_data) |>
  select(lat, lon, PRESENCE, MONTH, all_of(VARS)) |>
  rowwise() |>
  mutate_at(VARS[!VARS == "Bathy_depth"], ~.x[[MONTH]]) |>
  ungroup() |>
  mutate(MONTH = as.factor(MONTH))

### PERFORMING MODELLING 

# performing initial split 
data_split <- initial_split(model_data, prop = 3/4, strata = PRESENCE)
training_data <- training(data_split)
testing_data <- testing(data_split)

# example recipe 
recipe <- recipe(PRESENCE ~ ., data = training_data) |>
  update_role(lat, lon, U, V, new_role = "ID") |>
  step_mutate(Vel = sqrt(U^2 + V^2), role = "predictor") |>
  step_corr(all_numeric_predictors(), threshold = .95) |>
  step_zv(all_predictors()) |>
  step_normalize(all_numeric_predictors())

# example model: random forest
model <- rand_forest(trees = 15) |>
  set_engine("ranger") |>
  set_mode("classification")

# workflow 
workflow <- workflow() |>
  add_recipe(recipe) |>
  add_model(model) |>
  fit(training_data)

# results
test_results <- augment(workflow, testing_data)

### COLLECTING METRICS 

# overall metrics
pa_metrics <- metric_set(roc_auc, sens, spec, accuracy)
pa_metrics(test_results, 
           truth = PRESENCE, 
           estimate = .pred_class, 
           .pred_1,
           event_level = "second")

# collecting AUC by month
auc_monthly <- count(test_results, MONTH) |>
  bind_cols(AUC = split(test_results, test_results$MONTH) |> 
              lapply(function(x) roc_auc_vec(x$PRESENCE,
                                             x$.pred_1,
                                             event_level = "second")) |>
      unlist()) |>
  mutate(MONTH = as.numeric(MONTH))

ggplot(data = auc_monthly, 
       mapping = aes(x = MONTH, y = AUC)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name = "Month", 
                     breaks = 1:12, 
                     labels = c("Jan", "Feb", "Mar", "Apr", 
                                "May", "Jun", "Jul", "Aug",
                                "Sep", "Oct", "Nov", "Dec")) +
  scale_y_continuous(name = "AUC", limits = c(.5, 1)) +
  ggtitle("AUC by Month") +
  theme_classic() + 
  theme(panel.grid.major.y = element_line())

### GENERATING PREDICTIONS

downsample_val <- 3

# predictions for the most extreme climate situation: RCP85 2075. 
rcp85_2075 <- get_predictions(wkf = workflow, 
                              brickman_vars = VARS, 
                              year = 2075, 
                              scenario = "RCP85", 
                              augment_preds = FALSE, 
                              verbose = FALSE, 
                              downsample = downsample_val) 

# predictions for the least extreme future climate situation: RCP45 2055
rcp45_2055 <- get_predictions(wkf = workflow, 
                              brickman_vars = VARS,
                              year = 2055, 
                              scenario = "RCP45", 
                              augment_preds = FALSE, 
                              verbose = FALSE, 
                              downsample = downsample_val)

# present day predictions
present_preds <- get_predictions(wkf = workflow, 
                                 brickman_vars = VARS, 
                                 year = NA, 
                                 scenario = "PRESENT",
                                 augment_preds = FALSE,
                                 verbose = FALSE,
                                 downsample = downsample_val)

### RETRIEVING PLOTS

# raw plots for RCP85 2075
raw_plots_rcp85_2075 <- get_value_plots(preds_list = rcp85_2075, 
                                        title = "RCP85 2075 Predicted Presence Probability",
                                        pt_size = .3, 
                                        xlim = NULL,
                                        ylim = NULL)

# another example of raw plots, now for RCP45 2055 
# This plot is cropped to the Gulf of Maine/Gulf of Saint Lawrence
raw_plots_rcp45_2055 <- get_value_plots(preds_list = rcp45_2055, 
                                        title = "RCP45 2055 Predicted Presence Probability",
                                        pt_size = .3, 
                                        xlim = NULL,
                                        ylim = NULL)

# difference plots for RCP85 2075 relative to present day
difference_plots <- get_value_plots(preds_list = rcp85_2075,
                                    title = "RCP85 2075 Change in Presence Probability",
                                    pt_size = .3,
                                    xlim = NULL, 
                                    ylim = NULL,
                                    comparison_list = present_preds)

# threshold plots for RCP85 2075, relative to present day
threshold_plots <- get_threshold_plots(preds_list = rcp85_2075, 
                                       comparison_list = present_preds, 
                                       threshold = .5,
                                       title = "Shifts in presence between Present Day and RCP85 2075",
                                       pt_size = .3,
                                       xlim = NULL,
                                       ylim = NULL)




