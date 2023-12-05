#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

tic()

doParallel::registerDoParallel(cores=detectCores() - 2)

cell_line_compound_splits = read_rds(here('results/klaeger_screen_regression_CV_split.rds'))

xgb_spec <- boost_tree(
	trees = tune(), 
	tree_depth = tune(),
	min_n = tune(), 
	loss_reduction = tune(),                     
	sample_size = tune(), 
	mtry = tune(),         
	learn_rate = tune(),
) %>% 
	set_engine("xgboost") %>% 
	set_mode("regression")

xgb_grid <- grid_latin_hypercube(
	trees(),
	tree_depth(),
	min_n(),
	loss_reduction(),
	sample_size = sample_prop(),
	#mtry is affected by the size of the data sets, so we need to provide a sample
	#data set for the mtry ranges to be set
	finalize(mtry(), cell_line_compound_splits[["P0119-T1 CAF"]]),
	learn_rate(),
	size = 100
)

xgb_wf <- workflow() %>%
	add_formula(viability ~ .) %>%
	add_model(xgb_spec)


xgb_res_P0411 <- tune_grid(
	xgb_wf,
	resamples = cell_line_compound_splits[["P0411-T1"]],
	grid = xgb_grid,
	control = control_grid(save_pred = TRUE)
) %>% write_rds(here('results/xgboost_regression_models/P0411.rds'), compress = 'gz')

toc()