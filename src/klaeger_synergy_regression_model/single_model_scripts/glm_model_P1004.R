#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

tic()
doParallel::registerDoParallel(cores=detectCores() - 2)

cell_line_compound_splits = read_rds(here('results/klaeger_synergy_regression_CV_split.rds'))

norm_recipe <- 
	recipe(
		viability ~ ., 
		data = cell_line_compound_splits$P1004$splits[[1]]$data
	) %>%
	step_center(all_predictors()) %>%
	step_scale(all_predictors()) %>%
	# estimate the means and standard deviations
	prep(training = cell_line_compound_splits$P1004$splits[[1]]$data, retain = TRUE)

glm_spec <- linear_reg(
	penalty = tune(), 
	mixture = tune()
) %>% set_engine("glmnet") %>% 
	set_mode("regression")

glm_grid <- grid_latin_hypercube(
	penalty(),
	mixture(),
	size = 100
)

glm_wf <- workflow() %>%
	add_recipe(norm_recipe) %>%
	add_model(glm_spec)

tune_grid(
	glm_wf,
	resamples = cell_line_compound_splits$P1004,
	grid = glm_grid,
	control = control_grid(save_pred = TRUE)
) %>% write_rds(here('results/glm_regression_models/P1004.rds'), compress = 'gz')
toc()