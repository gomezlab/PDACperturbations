#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

knitr::opts_knit$set(root.dir = here())

tic()
doParallel::registerDoParallel(cores=detectCores() - 2)

cell_line_compound_splits = read_rds(here('results/klaeger_synergy_regression_CV_split.rds'))

svm_spec <- svm_poly(
	cost = tune(), 
	degree = tune(),
	scale_factor = tune()
) %>% set_engine("kernlab") %>% 
	set_mode("regression")

svm_grid <- grid_latin_hypercube(
	cost(),
	degree(),
	scale_factor(),
	size = 100
)

svm_wf <- workflow() %>%
	add_formula(viability ~ .) %>%
	add_model(svm_spec)


svm_res_P1004 <- tune_grid(
	svm_wf,
	resamples = cell_line_compound_splits$P1004,
	grid = svm_grid,
	control = control_grid(save_pred = TRUE)
) %>% write_rds(here('results/svm_regression_models/P1004.rds'), compress = 'gz')
toc()