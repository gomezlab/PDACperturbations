#!/usr/bin/env Rscript

library(tidyverse)
library(here)
library(tidymodels)
library(tictoc)
library(doParallel)

knitr::opts_knit$set(root.dir = here())

doParallel::registerDoParallel(cores=detectCores() - 2)

cell_line_compound_splits = read_rds(here('results/klaeger_synergy_classification_90_CV_split.rds'))

svm_spec <- svm_poly(
	cost = tune(), 
	degree = tune(),
	scale_factor = tune()
) %>% set_engine("kernlab") %>% 
	set_mode("classification")

svm_grid <- grid_latin_hypercube(
	cost(),
	degree(),
	scale_factor(),
	size = 100
)

svm_wf <- workflow() %>%
	add_formula(viability_90 ~ .) %>%
	add_model(svm_spec)

tic()
svm_res_CAF <- tune_grid(
	svm_wf,
	resamples = cell_line_compound_splits$P1304,
	grid = svm_grid,
	control = control_grid(save_pred = TRUE)
) %>% write_rds(here('results/svm_below90_models/P1304.rds'), compress = 'gz')
toc()
