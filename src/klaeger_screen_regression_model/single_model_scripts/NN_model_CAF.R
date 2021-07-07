#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

tic()
doParallel::registerDoParallel(cores=detectCores())

cell_line_compound_splits = read_rds(here('results/klaeger_screen_regression_CV_split.rds'))

nn_spec = mlp(
	epochs = tune(),
	hidden_units = tune(),
	dropout = tune()
) %>% set_engine("keras", verbose = 0) %>%
	set_mode("regression")

nn_grid = grid_latin_hypercube(
	epochs(),
	hidden_units(),
	dropout(),
	size = 100
)

nnet_wf <- workflow() %>%
	add_formula(viability ~ .) %>%
	add_model(nn_spec)

tune_grid(
	nnet_wf,
	resamples = cell_line_compound_splits$CAF,
	grid = nn_grid,
	control = control_grid(save_pred = TRUE)) %>% 
	write_rds(here('results/NN_regression_models/CAF.rds'), compress = 'gz')

toc()
