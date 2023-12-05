Klaeger Synergy Viability Binary Predictions
================
Matthew Berginski
2023-11-14

# Read In and Combine Klaeger/Synergy Data

``` r
klaeger_screen_model_data = read_rds(here('results/klaeger_screen_for_classification_90.rds')) %>%
    mutate(viability_90 = as.factor(viability_90))

klaeger_for_prediction = read_rds(here('results/klaeger_wide_for_prediction.rds'))
```

``` r
best_model_configs = read_rds(here('results/best_model_configs/rand_forest.rds'))
names(best_model_configs) <- c("P0119-T1 CAF", "P0422-T1", "P0411-T1")
```

# Modeling

# Make Predictions Using Models

Now to use the models to make predictions on the remainder of the
compounds in the Klaeger data set and build a few plots showing the
distribution of viability predictions.

## Below 90% Cell Viability

![](build_klaeger_screen_binary_predictions_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

![](build_klaeger_screen_binary_predictions_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
CAF_importance = binary_90_models[["P0119-T1 CAF"]] %>% 
    vip(num_features = 30)

CAF_DK = CAF_importance$data$Variable[CAF_importance$data$Variable %in% dark_kinases$symbol]

P0422_importance = binary_90_models[["P0422-T1"]] %>% 
    vip(num_features = 30)

P0422_DK = P0422_importance$data$Variable[P0422_importance$data$Variable %in% dark_kinases$symbol]

P0411_importance = binary_90_models[["P0411-T1"]] %>% 
    vip(num_features = 30)

P0411_DK = P0411_importance$data$Variable[P0411_importance$data$Variable %in% dark_kinases$symbol]

CAF_plot_data = as.data.frame(CAF_importance$data) %>% 
    mutate(Variable = fct_relevel(as.factor(Variable), rev(Variable)),
                 cell_line = "P0119-T1 CAF") #%>% 
    # mutate(cell_line = case_when(
    #   cell_line == "P1004" ~ "P0422-T1",
    #   cell_line == "P1304" ~ "P0411-T1",
    #   cell_line == "CAF" ~ "P0119-T1 CAF"
    # ))
P0422_plot_data = as.data.frame(P0422_importance$data) %>% 
    mutate(Variable = ifelse(Variable == "CSNK2A1;CSNK2A3","CSNK2A(1|3)",Variable)) %>%
    mutate(Variable = fct_relevel(as.factor(Variable), rev(Variable)),
                 cell_line = "P0422-T1") #%>% 
    # mutate(cell_line = case_when(
    #   cell_line == "P1004" ~ "P0422-T1",
    #   cell_line == "P1304" ~ "P0411-T1",
    #   cell_line == "CAF" ~ "P0119-T1 CAF"
    # ))

P0411_plot_data = as.data.frame(P0411_importance$data) %>% 
    mutate(Variable = fct_relevel(as.factor(Variable), rev(Variable)),
                 cell_line = "P0411-T1") #%>% 
    # mutate(cell_line = case_when(
    #   cell_line == "P1004" ~ "P0422-T1",
    #   cell_line == "P1304" ~ "P0411-T1",
    #   cell_line == "CAF" ~ "P0119-T1 CAF"
    # ))
```

``` r
CAF_plot = ggplot(CAF_plot_data, aes(y=Variable,x=Importance)) + 
    geom_col() +
    labs(x="Importance",y="Kinase Target",title = 'P0119-T1 CAF Model') +
    BerginskiRMisc::theme_berginski()#+
 # theme(axis.text=element_text(size=6))

P0422_plot = ggplot(P0422_plot_data, aes(y=Variable,x=Importance)) + 
    geom_col() +
    labs(x="Importance",y="Kinase Target",title = 'P0422-T1 Model') +
    BerginskiRMisc::theme_berginski()

P0411_plot = ggplot(P0411_plot_data, aes(y=Variable,x=Importance)) + 
    geom_col() +
    labs(x="Importance",y="",title = 'P0411-T1 Model') +
    BerginskiRMisc::theme_berginski()

VIP_upset = bind_rows(CAF_plot_data, P0422_plot_data, P0411_plot_data) %>%
    group_by(Variable) %>%
    summarise(cell_lines = list(cell_line)) %>%
  ggplot(aes(x=rev(cell_lines)))+
    geom_bar() + 
    scale_x_upset() +
    scale_y_continuous(breaks=seq(0,14,by=2)) +
    BerginskiRMisc::theme_berginski() +
    labs(x="Cell Lines",y="Number of Kinases")

full_importance_plot = CAF_plot + P0422_plot + P0411_plot
ggsave(here('figures/prediction_results/VIP_plot.png'),width=9,height=4)
BerginskiRMisc::trimImage(here('figures/prediction_results/VIP_plot.png'))

P0411_plot = ggplot(P0411_plot_data, aes(y=Variable,x=Importance)) +
    geom_col() +
    labs(x="Importance",y="Kinase Target",title = 'P0411-T1 Model') +
    BerginskiRMisc::theme_berginski()

VIP_with_upset = (P0411_plot + P0422_plot) / (CAF_plot + VIP_upset)
ggsave(here('figures/prediction_results/VIP_with_upset.png'),height=8,width=9)


VIP_with_upset = ggarrange(P0411_plot, P0422_plot,CAF_plot, VIP_upset, ncol= 2, nrow=2)
ggsave(here('figures/prediction_results/VIP_with_upset.png'),height=8,width=9)
BerginskiRMisc::trimImage(here('figures/prediction_results/VIP_with_upset.png'))
VIP_with_upset 
```

![](build_klaeger_screen_binary_predictions_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Looking for Compounds to Test - Below 90 Predictions

With predictions of the likelihood that the compounds will push cell
viability below 90% or below 40%, the next step is to try to use these
predictions to pick out compound that are likely to have interesting
viability results. As for what interesting means, I’ve collected four
types of interesting results:

- Compounds with Low Predicted Effect on Viability: Pick out the
  compounds that are predicted to have no/minimal effect on viability.
  This might seem to be the least interesting type of prediction at
  first, but a model is only as good as it’s ability to predict negative
  as well as positive results.
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
- calculated the absolute value of all the combination of differences in
  probability value
- We have three cell lines here (CAF, P0422 and P0411), so three
  combinations (CAF/P0422), (CAF/P0411) and (P0422/P0411)
- find the compounds with the largest mean value across all differences

## Compounds with Low Predicted Effect

![](build_klaeger_screen_binary_predictions_files/figure-gfm/low%20effect%20below%2090-1.png)<!-- -->

## CAF Survival Max

## Compounds with High Predicted Effect

![](build_klaeger_screen_binary_predictions_files/figure-gfm/high%20effect%20below%2090-1.png)<!-- -->

## Compounds with High Predicted Differences Across Concentrations

![](build_klaeger_screen_binary_predictions_files/figure-gfm/concentration%20diff%20below%2090-1.png)<!-- -->

## Compounds with High Predicted Differences Between Cell Lines

![](build_klaeger_screen_binary_predictions_files/figure-gfm/cell%20line%20diff%20below%2090-1.png)<!-- -->
