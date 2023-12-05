#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

tic()

cell_line_compound_splits = read_rds(here('results/klaeger_screen_classification_90_CV_split.rds'))

nn_spec = mlp(
	epochs = tune(),
	hidden_units = tune(),
	dropout = tune()
) %>% set_engine("keras", verbose = 0) %>%
	set_mode("classification")

nn_grid = grid_latin_hypercube(
	epochs(),
	hidden_units(),
	dropout(c(0,0.25)),
	size = 100
)

nnet_wf <- workflow() %>%
	add_formula(viability_90 ~ .) %>%
	add_model(nn_spec)

tune_grid(
	nnet_wf,
	resamples = cell_line_compound_splits[["P0422-T1"]],
	grid = nn_grid,
	control = control_grid(save_pred = TRUE)) %>% 
	write_rds(here('results/NN_below90_models/P0422.rds'), compress = 'gz')

toc()