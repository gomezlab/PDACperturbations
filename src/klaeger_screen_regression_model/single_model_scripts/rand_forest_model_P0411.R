#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

tic()

doParallel::registerDoParallel(cores=detectCores() - 2)

cell_line_compound_splits = read_rds(here('results/klaeger_screen_regression_CV_split.rds'))

rand_forest_spec <- rand_forest(
	trees = tune(), 
	mtry = tune(),
	min_n = tune()
) %>% set_engine("ranger") %>% 
	set_mode("regression")

rand_forest_grid <- grid_latin_hypercube(
	trees(c(1000,5000)),
	min_n(),
	finalize(mtry(),cell_line_compound_splits[["P0119-T1 CAF"]]),
	size = 100
)

rand_forest_wf <- workflow() %>%
	add_formula(viability ~ .) %>%
	add_model(rand_forest_spec)


rand_forest_res_P0411 <- tune_grid(
	rand_forest_wf,
	resamples = cell_line_compound_splits[["P0411-T1"]],
	grid = rand_forest_grid,
	control = control_grid(save_pred = TRUE)
) %>% write_rds(here('results/rand_forest_regression_models/P0411.rds'), compress = 'gz')

toc()