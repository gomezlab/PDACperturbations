Klaeger Synergy Viability Predictions
================
Matthew Berginski
2021-05-11

# Read In and Combine Klaeger/Synergy Data

The first step for the prediction pipeline is the same as for the model
testing and production pipeline. Namely, read in the pre-processed
Klaeger data and then combine it with the appropriate synergy screen
results. Then, using this combined data set, build models for each of
the cell lines.

First a few facts about the compounds/concentrations used for the model
construction and predictions:

  - \# of Compounds Used for Model Construction: 40
  - \# of Compound/Concentration Combinations in Model Construction: 200
  - \# of Compounds in Prediction Set: 182
  - \# of Compound/Concentration Combinations in Prediction Set: 1456

Assuming we’d like to test all of the drugs (222 total) in combination
(24420 combinations) at 6 concentrations for each compound. The total
number of assays to run would be 879120, which would need 2290 384 well
plates per cell line.

# Modeling

Build two sets of random forest models, one to predict whether a
compound at a concentration will cause cell viability to dip below 90%
and another to predict cell viability below 40%. Unique models are built
for each of the cell lines in the data set.

# Make Predictions Using Models

Now to use the models to make predictions on the remainder of the
compounds in the Klaeger data set and build a few plots showing the
distribution of viability predictions.

## Below 90% Cell Viability

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Below 40% Cell Viability

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Looking for Compounds to Test - Below 90 Predictions

With predictions of the likelihood that the compounds will push cell
viability below 90% or below 40%, the next step is to try to use these
predictions to pick out compound that are likely to have interesting
viability results. As for what interesting means, I’ve collected four
types of interesting results:

  - Compounds with Low Predicted Effect on Viability: Pick out the
    compounds that are predicted to have no/minimal effect on viability.
    This might seem to be the least interesting type of prediction at
    first, but a model is only as good as it’s ability to predict
    negative as well as positive results.
  - Compounds with High Predicted Effect on Viability: Pick out the
    compounds that are predicted to strongly decrease cell viability.
  - Compounds with Large Differences across Concentrations: Pick out
    compounds where the likelihood of affecting cell viability changes
    across the predicted concentrations. These should be compounds where
    low concentrations don’t affect cell viability, while high
    concentrations start to inhibit cell growth.
  - Compounds with Large Differences across Cell Lines: Pick out the
    compounds with large differences between the cell lines. This one
    seems kind of obvious as well, but if the model predicts large
    differences between the cell lines, this is probably worth testing

As for calculating all these, I’ve taken an approach that looks across
all the cell lines simultaneously. Searching for across the board low
effect and high effects is done by looking for the lowest and highest
average probability on a per compound basis. Searching for large
differences across concentrations is done by looking for the compounds
with the highest range between probabilities. The final criteria
(differences between cell lines) is a bit tougher to quantify, but what
I’ve done is:

  - matched up all the predictions across compound and concentration for
    each cell line
  - calculated the absolute value of all the combination of differences
    in probability value
      - We have three cell lines here (CAF, P1004 and P1304), so three
        combinations (CAF/P1004), (CAF/P1304) and (P1004/P1304)
  - find the compounds with the largest mean value across all
    differences

## Compounds with Low Predicted Effect

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/low%20effect%20below%2090-1.png)<!-- -->

## Compounds with High Predicted Effect

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/high%20effect%20below%2090-1.png)<!-- -->

## Compounds with High Predicted Differences Across Concentrations

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/concentration%20diff%20below%2090-1.png)<!-- -->

## Compounds with High Predicted Differences Between Cell Lines

![](build_klaeger_synergy_binary_predictions_files/figure-gfm/cell%20line%20diff%20below%2090-1.png)<!-- -->
