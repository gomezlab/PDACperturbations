# This script simply runs all of the model building files in the single model scripts folder.

source('src/klaeger_synergy_binary_model/single_model_scripts/rand_forest_model_CAF.R')
source('src/klaeger_synergy_binary_model/single_model_scripts/rand_forest_model_P1004.R')
source('src/klaeger_synergy_binary_model/single_model_scripts/rand_forest_model_P1304.R')

source('src/klaeger_synergy_binary_model/single_model_scripts/xgboost_model_CAF.R')
source('src/klaeger_synergy_binary_model/single_model_scripts/xgboost_model_P1004.R')
source('src/klaeger_synergy_binary_model/single_model_scripts/xgboost_model_P1304.R')

source('src/klaeger_synergy_binary_model/single_model_scripts/svm_model_CAF.R')
source('src/klaeger_synergy_binary_model/single_model_scripts/svm_model_P1004.R')
source('src/klaeger_synergy_binary_model/single_model_scripts/svm_model_P1304.R')