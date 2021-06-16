#!/usr/bin/env bash

sbatch -c 64 --time=24:00:00 --wrap="./rand_forest_model_CAF.R"
sbatch -c 64 --time=24:00:00 --wrap="./rand_forest_model_P1004.R"
sbatch -c 64 --time=24:00:00 --wrap="./rand_forest_model_P1304.R"

sbatch -c 64 --time=24:00:00 --wrap="./svm_model_CAF.R"
sbatch -c 64 --time=24:00:00 --wrap="./svm_model_P1004.R"
sbatch -c 64 --time=24:00:00 --wrap="./svm_model_P1304.R"

sbatch -c 64 --time=24:00:00 --wrap="./xgboost_model_CAF.R"
sbatch -c 64 --time=24:00:00 --wrap="./xgboost_model_P1004.R"
sbatch -c 64 --time=24:00:00 --wrap="./xgboost_model_P1304.R"
