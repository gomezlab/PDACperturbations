# PDACperturbations

Data and code associated with pancreatic cancer drug screening and modeling

## Data Cleaning and Organization

* Synergy Data: [here](src/process_synergy_data/process_synergy.md)
  * Reads in and organizes the synergy screen data and tries to search synergistic relationships using linear models
* Klaeger Data: [here](src/process_klaeger_data/klaeger_data_processing.md)
  * Organizes the Klaeger data for downstream processing

### Compound Matching

Most of the compound names in the synergy/Klaeger collections don't match up exactly, so we had to go through and manually match most of the shared compounds.

* Report [here](src/find_synergy_klaeger_matches/find_synergy_klaeger_matches.md)

### Preparation for Machine Learning

With the compound/drug names matched, I preped the data set for machine learning (both regression and above/below 90 viability classification). This code also removes any genes which don't vary in the Klaeger collection after the compounds have been filtered to only the matching compounds.

* Report [here](src/process_synergy_klaeger_for_ML/process_for_ML.md)

## Modeling Results

* Synergy Screen + Klaeger to Predict Raw Viability Values: [here](src/klaeger_synergy_model/build_klaeger_synergy_model.md)
  * Random forest modeling code to try to predict the cell viability values in the single drug treatments in the synergy data from the Klaeger kinase activation values
* Synergy Screen + Klaeger to Predict Binary Viability Values: [here](src/klaeger_synergy_binary_model/build_klaeger_synergy_binary_model.md)
  * Random forest modeling code to try to predict binarized (below 90 and below 40) cell viability values in the single drug treatments in the synergy data from the Klaeger kinase activation values

## Prediction Results

* Predicting the binary outcome for the rest of the Klaeger compounds: [here](src/klaeger_synergy_binary_predictions/build_klaeger_synergy_binary_predictions.md)
