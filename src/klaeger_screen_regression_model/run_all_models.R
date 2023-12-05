# This script builds out the required directory structure for the model results
# and then runs all the single model scripts.

library(here)

dir.create(here('results/rand_forest_regression_models/'), 
					 showWarnings = FALSE, 
					 recursive = TRUE)

source(here('src/klaeger_screen_regression_model/single_model_scripts/rand_forest_model_CAF.R'))
source(here('src/klaeger_screen_regression_model/single_model_scripts/rand_forest_model_P0422.R'))
source(here('src/klaeger_screen_regression_model/single_model_scripts/rand_forest_model_P0411.R'))

dir.create(here('results/xgboost_regression_models/'), 
					 showWarnings = FALSE, 
					 recursive = TRUE)

source(here('src/klaeger_screen_regression_model/single_model_scripts/xgboost_model_CAF.R'))
source(here('src/klaeger_screen_regression_model/single_model_scripts/xgboost_model_P0422.R'))
source(here('src/klaeger_screen_regression_model/single_model_scripts/xgboost_model_P0411.R'))

dir.create(here('results/glm_regression_models/'), 
					 showWarnings = FALSE, 
					 recursive = TRUE)

source(here('src/klaeger_screen_regression_model/single_model_scripts/glm_model_CAF.R'))
source(here('src/klaeger_screen_regression_model/single_model_scripts/glm_model_P0422.R'))
source(here('src/klaeger_screen_regression_model/single_model_scripts/glm_model_P0411.R'))
