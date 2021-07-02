# Predicting Pancreatic Cancer Cell Line Response to Kinase Inhibitors

Data and code associated with kinase activation state and cell viabiltiy modeling. Most of the code consists of Rmarkdown documents, with the model testing code saved as R script files. The code requires several packages and is organized into sequential steps:

## Required Packages

I've written a [script](src/package_check.R) that checks for and installs all of the packages required in the repository. I use the pacman package for this purpose and installing pacman if missing is covered in the script. There are also two github based packages:

  * [BerginskiRmisc](https://github.com/mbergins/BerginskiRMisc) for my custom theme and helper scripts
    * The helper script I use calls the "convert" command from imagemagick to trim whitespace around images, so imagemagick will need to be installed as well
  * [DarkKinaseTools](https://github.com/IDG-Kinase/DarkKinaseTools) for kinase lists

The BerginskiRmisc package uses the convert command from Imagemagick to trim whitespace around images, if you don't have Imagemagick installed it should work with a failure message.

## Data Cleaning and Organization

There are two primary data sets in the repository, results of a screening assay and the data downloaded from the supplement of Klaeger et al. These first scripts take each of these data sets and produce R data files that are then used in downstream processing.

* Screen Data: [here](src/process_synergy_data/process_synergy_YehLabHTS.md)
  * Reads in and organizes the screen data
* Klaeger Data: [here](src/process_klaeger_data/klaeger_data_processing.md)
  * Organizes the Klaeger data for downstream processing

### Compound Matching

Most of the compound names in the synergy/Klaeger collections don't match up exactly, so we had to go through and manually match most of the shared compounds. This report has the code used to simplify this search:

* Report [here](src/find_synergy_klaeger_matches/find_synergy_klaeger_matches.md)

### Preparation for Machine Learning

With the compound/drug names matched, I preped the data sets for machine learning (both regression and above/below 90 viability classification). This code also removes any genes which don't vary in the Klaeger collection after the compounds have been filtered to only the matching compounds. The primary output here are machine learning ready data sets (deposited in the [results](results)) and cross-validation split data sets for both regression and binary classification.

* Report [here](src/process_synergy_klaeger_for_ML/process_for_ML.md)

## Model Testing

The model testing code is saved single self contained scripts which fully implement and run the hyperparameter scanning and model testing. The code is organized this way to make it simplier to run the modelling code on the UNC computing cluster, but should also be compatible with any computing environment. This code takes a long time to run. All of the regression testing models are available [here](src/klaeger_synergy_regression_model/single_model_scripts) and the binary above/below 90 are available [here](src/klaeger_synergy_binary_model/single_model_scripts). There are also scripts (search for run_all_models.R) that build out directory infrastructure and run all the models sequentially. 

## Prediction Results

Using the random forest model and associated code, predicting the rest of the Klaeger compounds is [here](src/klaeger_synergy_binary_predictions/build_klaeger_synergy_binary_predictions.md).

## Reproducibility Script

I've attempted to write a single [script](src/reproduce_results.R) that runs all the Rmarkdown files and scripts to completely reproduce the models and figures from the paper. I've only tested the code on Linux, but I see no reason why it wouldn't work on other platforms. Let me know if you attempt to run this script and it fails on your platform.

This script takes a long time to run (XX hours on 64 core cluster node, XX hours on a Ryzen 7 5800x). This is mostly due to the hyperparameter scanning for the regression and binary models. 
