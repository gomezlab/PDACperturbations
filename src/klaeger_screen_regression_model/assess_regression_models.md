Assess Klaeger/Synergy Models
================
Matthew Berginski
2023-11-10

# Convenience Functions

``` r
get_best_RMSE_model <- function(model_metrics) {
    minimum_model_val = max(model_metrics$n)
    model_metrics %>%
        filter(.metric == "rmse",n == minimum_model_val) %>% 
        arrange(mean) %>%
        dplyr::slice(1) %>%
        return()
}

get_best_predictions <- function(model_set,best_config,cell_line,model_type) {
    model_set %>% 
        collect_predictions() %>% 
        filter(.config == best_config) %>%
        select(id,.pred,viability,.config) %>%
        mutate(cell_line = cell_line, model_type = model_type)  %>%
        return()
}

best_predictions = data.frame()
best_RMSE_models = list()
```

# XGBoost Model Assessment

``` r
xgb_res_CAF = read_rds(here('results/xgboost_regression_models/CAF.rds'))
xgb_res_P0422 = read_rds(here('results/xgboost_regression_models/P0422.rds'))
xgb_res_P0411 = read_rds(here('results/xgboost_regression_models/P0411.rds'))

best_RMSE_models$XGBoost$CAF = get_best_RMSE_model(xgb_res_CAF %>% collect_metrics())
best_RMSE_models$XGBoost$P0422 = get_best_RMSE_model(xgb_res_P0422 %>% collect_metrics())
best_RMSE_models$XGBoost$P0411 = get_best_RMSE_model(xgb_res_P0411 %>% collect_metrics())

best_predictions = bind_rows(
    best_predictions,
    get_best_predictions(xgb_res_CAF,best_RMSE_models$XGBoost$CAF$.config[1],"CAF","XGBoost"),
    get_best_predictions(xgb_res_P0422,best_RMSE_models$XGBoost$P0422$.config[1],"P0422","XGBoost"),
    get_best_predictions(xgb_res_P0411,best_RMSE_models$XGBoost$P0411$.config[1],"P0411","XGBoost"),
)
```

# Random Forest Model Assessment

``` r
rand_forest_res_CAF = read_rds(here('results/rand_forest_regression_models/CAF.rds'))
rand_forest_res_P0422 = read_rds(here('results/rand_forest_regression_models/P0422.rds'))
rand_forest_res_P0411 = read_rds(here('results/rand_forest_regression_models/P0411.rds'))

best_RMSE_models$rand_forest$CAF = get_best_RMSE_model(rand_forest_res_CAF %>% collect_metrics())
best_RMSE_models$rand_forest$P0422 = get_best_RMSE_model(rand_forest_res_P0422 %>% collect_metrics())
best_RMSE_models$rand_forest$P0411 = get_best_RMSE_model(rand_forest_res_P0411 %>% collect_metrics())

best_predictions = bind_rows(
    best_predictions,
    get_best_predictions(rand_forest_res_CAF,best_RMSE_models$rand_forest$CAF$.config[1],"CAF","Random Forest"),
    get_best_predictions(rand_forest_res_P0422,best_RMSE_models$rand_forest$P0422$.config[1],"P0422","Random Forest"),
    get_best_predictions(rand_forest_res_P0411,best_RMSE_models$rand_forest$P0411$.config[1],"P0411","Random Forest"),
)
```

# SVM Model Assessment

Blocking SVM model assessment for regression, since I don’t really think
SVM’s were made for regression…

``` r
# svm_res_CAF = read_rds(here('results/svm_regression_models/CAF.rds'))
# svm_res_P0422 = read_rds(here('results/svm_regression_models/P0422.rds'))
# svm_res_P0411 = read_rds(here('results/svm_regression_models/P0411.rds'))
# 
# best_RMSE_models$svm$CAF = get_best_RMSE_model(svm_res_CAF %>% collect_metrics())
# best_RMSE_models$svm$P0422 = get_best_RMSE_model(svm_res_P0422 %>% collect_metrics())
# best_RMSE_models$svm$P0411 = get_best_RMSE_model(svm_res_P0411 %>% collect_metrics())
# 
# best_predictions = bind_rows(
#   best_predictions,
#   get_best_predictions(svm_res_CAF,best_RMSE_models$svm$CAF$.config[1],"CAF","SVM"),
#   get_best_predictions(svm_res_P0422,best_RMSE_models$svm$P0422$.config[1],"P0422","SVM"),
#   get_best_predictions(svm_res_P0411,best_RMSE_models$svm$P0411$.config[1],"P0411","SVM"),
# )
```

# GLM Model Assessment

``` r
glm_res_CAF = read_rds(here('results/glm_regression_models/CAF.rds'))
glm_res_P0422 = read_rds(here('results/glm_regression_models/P0422.rds'))
glm_res_P0411 = read_rds(here('results/glm_regression_models/P0411.rds'))

best_RMSE_models$glm$CAF = get_best_RMSE_model(glm_res_CAF %>% collect_metrics())
best_RMSE_models$glm$P0422 = get_best_RMSE_model(glm_res_P0422 %>% collect_metrics())
best_RMSE_models$glm$P0411 = get_best_RMSE_model(glm_res_P0411 %>% collect_metrics())

best_predictions = bind_rows(
    best_predictions,
    get_best_predictions(glm_res_CAF,best_RMSE_models$glm$CAF$.config[1],"CAF","GLMnet"),
    get_best_predictions(glm_res_P0422,best_RMSE_models$glm$P0422$.config[1],"P0422","GLMnet"),
    get_best_predictions(glm_res_P0411,best_RMSE_models$glm$P0411$.config[1],"P0411","GLMnet"),
)
```

``` r
rmse_text = tribble(
    ~cell_line, ~model_type, ~RMSE,
    "CAF", "GLMnet", best_RMSE_models$glm$CAF$mean,
    "P0422", "GLMnet", best_RMSE_models$glm$P0422$mean,
    "P0411", "GLMnet", best_RMSE_models$glm$P0411$mean,
    "CAF", "Random Forest", best_RMSE_models$rand_forest$CAF$mean,
    "P0422", "Random Forest", best_RMSE_models$rand_forest$P0422$mean,
    "P0411", "Random Forest", best_RMSE_models$rand_forest$P0411$mean,
    "CAF", "XGBoost", best_RMSE_models$XGBoost$CAF$mean,
    "P0422", "XGBoost", best_RMSE_models$XGBoost$P0422$mean,
    "P0411", "XGBoost", best_RMSE_models$XGBoost$P0411$mean,
) %>% mutate(x = range(best_predictions$viability)[2],
                         y = range(best_predictions$viability)[1],
                         RMSE_text = paste0("RMSE: ",round(RMSE,1))) %>%
mutate(cell_line = case_when(
        cell_line == "P0422" ~ "P0422-T1",
        cell_line == "P0411" ~ "P0411-T1",
        cell_line == "CAF" ~ "P0119-T1 CAF"))

 best_predictions = best_predictions %>%
    mutate(cell_line = case_when(
        cell_line == "P0422" ~ "P0422-T1",
        cell_line == "P0411" ~ "P0411-T1",
        cell_line == "CAF" ~ "P0119-T1 CAF"))

ggplot(best_predictions,aes(x=viability, y=.pred)) + 
    ylim(range(best_predictions$viability)) + 
    xlim(range(best_predictions$viability)) + 
    geom_hex() + 
    geom_smooth(method = lm, se = FALSE, alpha = 0.5, color = "green") +
    geom_abline(intercept = 0,slope = 1, alpha = 0.5, linetype = 4) + 
    geom_text(data = rmse_text, 
                        mapping = aes(x = x,y = y,label = RMSE_text),  
                        hjust = "inward", vjust = "inward") +
    labs(x = "Actual Cell Viability", y = "Predicted Viability", fill = "Density of\nPoints") +
    BerginskiRMisc::theme_berginski() +
    facet_grid(vars(model_type), vars(fct_relevel(cell_line,"P0411-T1","P0422-T1","P0119-T1 CAF")))
```

    ## Warning: Removed 56 rows containing non-finite values (`stat_binhex()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 56 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 2 rows containing missing values (`geom_hex()`).

![](assess_regression_models_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
dir.create(here('figures/model_assessment'), recursive = T, showWarnings = F)
ggsave(here('figures/model_assessment/regression_models.png'))
```

    ## Saving 7 x 5 in image

    ## Warning: Removed 56 rows containing non-finite values (`stat_binhex()`).

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 56 rows containing non-finite values (`stat_smooth()`).
    ## Removed 2 rows containing missing values (`geom_hex()`).

``` r
BerginskiRMisc::trimImage(here('figures/model_assessment/regression_models.png'))
```
