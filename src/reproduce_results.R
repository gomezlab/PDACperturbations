#!/usr/bin/env Rscript

library(rmarkdown)
library(here)

render(here('src/process_screen_data/process_screen_YehLabHTS.Rmd'))
render(here('src/process_klaeger_data/klaeger_data_processing.Rmd'))
render(here('src/process_screen_klaeger_for_ML/process_for_ML.Rmd'))

source(here('src/klaeger_screen_regression_model/run_all_models.R'))
render(here('src/klaeger_screen_regression_model/assess_regression_models.Rmd'))

source(here('src/klaeger_screen_binary_model/run_all_models.R'))
render(here('src/klaeger_screen_binary_model/assess_binary_models.Rmd'))

render(here('src/klaeger_screen_binary_predictions/build_klaeger_screen_binary_predictions.Rmd'))

render(here('src/process_screen_data/build_screen_EDA_figures.Rmd'))