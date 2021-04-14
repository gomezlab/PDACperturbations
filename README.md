# PDACperturbations

Data and code associated with pancreatic cancer drug screening and modeling

## Data Cleaning

* Synergy Data: [here](src/process_synergy_data/process_synergy.md)
  * Reads in and organizes the synergy screen data and tries to search synergistic relationships using linear models
* Klaeger Data: [here](src/process_klaeger_data/klaeger_data_processing.md)
  * Organizes the Klaeger data for downstream processing

## Modeling Results

* Synergy Screen + Klaeger to Predict Raw Viability Values: [here](src/klaeger_synergy_model/build_klaeger_synergy_model.md)
  * Random forest modeling code to try to predict the cell viability values in the single drug treatments in the synergy data from the Klaeger kinase activation values
* Synergy Screen + Klaeger to Predict Binary Viability Values: [here](src/klaeger_synergy_binary_model/build_klaeger_synergy_binary_model.md)
  * Random forest modeling code to try to predict binarized (below 90 and below 40) cell viability values in the single drug treatments in the synergy data from the Klaeger kinase activation values

## Prediction Results

* Predicting the binary outcome for the rest of the Klaeger compounds: [here](src/klaeger_synergy_binary_predictions/build_klaeger_synergy_binary_predictions.md)
