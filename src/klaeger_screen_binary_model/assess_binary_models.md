Assess Klaeger/Synergy Models
================
Matthew Berginski
2023-11-10

# Convenience Functions

``` r
process_full_roc_data <- function(all_results, cell_line, required_data_vals = NA) {
    
    #Some of the hyperparameter settings throw errors on specific data sets
    #resulting in fewer than expected results being returned. This variable makes
    #sure that all the tested hyperparameter config have the expected number of
    #predictions returned.
    if (is.na(required_data_vals)) {
        required_data_vals = max(all_results %>% count(.config) %>% pull(n))
    }
    
    roc_vals = data.frame()
    for (model_id in unique(all_results$.config)) {
        these_results = all_results %>% filter(.config == model_id)
        
        if (dim(these_results)[1] < required_data_vals) {
            next;
        }
        pred <- prediction(these_results$.pred_TRUE,these_results$viability_90)
        perf_roc <- performance(pred, measure = "auc")
        prc_roc <- performance(pred, measure = "aucpr")
        
        roc_vals = bind_rows(
            roc_vals,
            data.frame(model_id = model_id,roc = perf_roc@y.values[[1]],prc = prc_roc@y.values[[1]] )
        )
    }
    
    roc_vals = roc_vals %>%
        mutate(roc_rank = dense_rank(desc(roc)),
                     prc_rank = dense_rank(desc(prc))) %>%
        arrange(roc_rank)
}

get_ROC_curve_values <- function(results, cell_line) {
    pred <- prediction(results$.pred_TRUE,results$viability_90)
    perf <- performance(pred,measure = "tpr",x.measure = "fpr")
    return(data.frame(fpr = perf@x.values[[1]],
                                        tpr = perf@y.values[[1]],
                                        cell_line = cell_line))
}

get_PRC_curve_values <- function(results, cell_line) {
    pred <- prediction(results$.pred_TRUE,results$viability_90)
    perf <- performance(pred,measure = "prec",x.measure = "rec")
    return(data.frame(precision = perf@y.values[[1]],
                                        recall = perf@x.values[[1]],
                                        cell_line = cell_line))
}
```

# XGBoost Model Assessment

``` r
xgb_res_CAF = read_rds(here('results/xgboost_below90_models/CAF.rds'))
xgb_res_P0422 = read_rds(here('results/xgboost_below90_models/P0422.rds'))
xgb_res_P0411 = read_rds(here('results/xgboost_below90_models/P0411.rds'))

xgboost_results = list()
xgboost_results[["CAF"]] = process_full_roc_data(xgb_res_CAF %>% collect_predictions(), "CAF")
xgboost_results[["P0422"]] = process_full_roc_data(xgb_res_P0422 %>% collect_predictions(), "P0422")
xgboost_results[["P0411"]] = process_full_roc_data(xgb_res_P0411 %>% collect_predictions(), "P0411")

best_ROC_vals_xgboost = bind_rows(
    xgb_res_CAF %>%
        collect_predictions() %>%
        filter(.config == xgboost_results[["CAF"]]$model_id[1]) %>%
        get_ROC_curve_values("CAF"),
    xgb_res_P0422 %>%
        collect_predictions() %>%
        filter(.config == xgboost_results[["P0422"]]$model_id[1]) %>%
        get_ROC_curve_values("P0422"),
    xgb_res_P0411 %>%
        collect_predictions() %>%
        #There was a problem with the first 4 XGBoost models where the predictions
        #didn't span a wide enough range to get a decent ROC curve estimate. Model
        #rank 5 seemed to work correctly.
        filter(.config == xgboost_results[["P0411"]]$model_id[5]) %>%
        get_ROC_curve_values("P0411"),
) %>% mutate(model_type = "XGBoost")

best_PRC_vals_xgboost = bind_rows(
    xgb_res_CAF %>%
        collect_predictions() %>%
        filter(.config == xgboost_results[["CAF"]]$model_id[1]) %>%
        get_PRC_curve_values("CAF"),
    xgb_res_P0422 %>%
        collect_predictions() %>%
        filter(.config == xgboost_results[["P0422"]]$model_id[1]) %>%
        get_PRC_curve_values("P0422"),
    xgb_res_P0411 %>%
        collect_predictions() %>%
        filter(.config == xgboost_results[["P0411"]]$model_id[5]) %>%
        get_PRC_curve_values("P0411"),
) %>% mutate(model_type = "XGBoost")
```

# Rand Forest Model Assessment

``` r
rand_forest_res_CAF = read_rds(here('results/rand_forest_below90_models/CAF.rds'))
rand_forest_res_P0422 = read_rds(here('results/rand_forest_below90_models/P0422.rds'))
rand_forest_res_P0411 = read_rds(here('results/rand_forest_below90_models/P0411.rds'))

rand_forest_results = list()
rand_forest_results[["CAF"]] = process_full_roc_data(rand_forest_res_CAF %>% collect_predictions(), "CAF")
rand_forest_results[["P0422"]] = process_full_roc_data(rand_forest_res_P0422 %>% collect_predictions(), "P0422")
rand_forest_results[["P0411"]] = process_full_roc_data(rand_forest_res_P0411 %>% collect_predictions(), "P0411")

rand_forest_best_model_parameters = list()
rand_forest_best_model_parameters[["CAF"]] = rand_forest_res_CAF %>% 
    collect_metrics() %>% 
    filter(.config == rand_forest_results[["CAF"]]$model_id[1], .metric == "roc_auc")
rand_forest_best_model_parameters[["P0422"]] = rand_forest_res_P0422 %>% 
    collect_metrics() %>% 
    filter(.config == rand_forest_results[["P0422"]]$model_id[1], .metric == "roc_auc")
rand_forest_best_model_parameters[["P0411"]] = rand_forest_res_P0411 %>% 
    collect_metrics() %>% 
    filter(.config == rand_forest_results[["P0411"]]$model_id[1], .metric == "roc_auc")
dir.create(here('results/best_model_configs/'), recursive = T, showWarnings = F)
write_rds(rand_forest_best_model_parameters,here('results/best_model_configs/rand_forest.rds'))

best_ROC_vals_rand_forest = bind_rows(
    rand_forest_res_CAF %>%
        collect_predictions() %>%
        filter(.config == rand_forest_results[["CAF"]]$model_id[1]) %>%
        get_ROC_curve_values("CAF"),
    rand_forest_res_P0422 %>%
        collect_predictions() %>%
        filter(.config == rand_forest_results[["P0422"]]$model_id[1]) %>%
        get_ROC_curve_values("P0422"),
    rand_forest_res_P0411 %>%
        collect_predictions() %>%
        filter(.config == rand_forest_results[["P0411"]]$model_id[1]) %>%
        get_ROC_curve_values("P0411"),
) %>% mutate(model_type = "Random Forest")

best_PRC_vals_rand_forest = bind_rows(
    rand_forest_res_CAF %>%
        collect_predictions() %>%
        filter(.config == rand_forest_results[["CAF"]]$model_id[1]) %>%
        get_PRC_curve_values("CAF"),
    rand_forest_res_P0422 %>%
        collect_predictions() %>%
        filter(.config == rand_forest_results[["P0422"]]$model_id[1]) %>%
        get_PRC_curve_values("P0422"),
    rand_forest_res_P0411 %>%
        collect_predictions() %>%
        filter(.config == rand_forest_results[["P0411"]]$model_id[1]) %>%
        get_PRC_curve_values("P0411"),
) %>% mutate(model_type = "Random Forest")
```

# SVM Model Assessment

``` r
svm_res_CAF = read_rds(here('results/svm_below90_models/CAF.rds'))
svm_res_P0422 = read_rds(here('results/svm_below90_models/P0422.rds'))
svm_res_P0411 = read_rds(here('results/svm_below90_models/P0411.rds'))

svm_results = list()
svm_results[["CAF"]] = process_full_roc_data(svm_res_CAF %>% collect_predictions(), "CAF")
svm_results[["P0422"]] = process_full_roc_data(svm_res_P0422 %>% collect_predictions(), "P0422")
svm_results[["P0411"]] = process_full_roc_data(svm_res_P0411 %>% collect_predictions(), "P0411")

best_ROC_vals_svm = bind_rows(
    svm_res_CAF %>%
        collect_predictions() %>%
        filter(.config == svm_results[["CAF"]]$model_id[1]) %>%
        get_ROC_curve_values("CAF"),
    svm_res_P0422 %>%
        collect_predictions() %>%
        filter(.config == svm_results[["P0422"]]$model_id[1]) %>%
        get_ROC_curve_values("P0422"),
    svm_res_P0411 %>%
        collect_predictions() %>%
        filter(.config == svm_results[["P0411"]]$model_id[1]) %>%
        get_ROC_curve_values("P0411"),
) %>% mutate(model_type = "SVM")

best_PRC_vals_svm = bind_rows(
    svm_res_CAF %>%
        collect_predictions() %>%
        filter(.config == svm_results[["CAF"]]$model_id[1]) %>%
        get_PRC_curve_values("CAF"),
    svm_res_P0422 %>%
        collect_predictions() %>%
        filter(.config == svm_results[["P0422"]]$model_id[1]) %>%
        get_PRC_curve_values("P0422"),
    svm_res_P0411 %>%
        collect_predictions() %>%
        filter(.config == svm_results[["P0411"]]$model_id[1]) %>%
        get_PRC_curve_values("P0411"),
) %>% mutate(model_type = "SVM")
```

# Model Assessment/Visualization

``` r
all_ROC_vals = bind_rows(
    best_ROC_vals_xgboost,
    best_ROC_vals_rand_forest,
    best_ROC_vals_svm
 ) %>% mutate(cell_line = case_when(
    cell_line == "P0422" ~ "P0422-T1",
    cell_line == "P0411" ~ "P0411-T1",
    cell_line == "CAF" ~ "P0119-T1 CAF"
 ))

all_PRC_vals = bind_rows(
    best_PRC_vals_xgboost,
    best_PRC_vals_rand_forest,
    best_PRC_vals_svm
 )%>% mutate(cell_line = case_when(
    cell_line == "P0422" ~ "P0422-T1",
    cell_line == "P0411" ~ "P0411-T1",
    cell_line == "CAF" ~ "P0119-T1 CAF"
 ))

# ROC_data = tribble(
#   ~"Cell Line",~"Model Type",~"AUROC",
#   "CAF","XGBoost",xgboost_results$roc_CAF[1],
#   "CAF","Rand Forest",rand_forest_results$roc_CAF[1],
#   "CAF","SVM",svm_results$roc_CAF[1],
#   "P0422","XGBoost",xgboost_results$roc_P0422[1],
#   "P0422","Rand Forest",rand_forest_results$roc_P0422[1],
#   "P0422","SVM",svm_results$roc_P0422[1],
#   "P1304","XGBoost",xgboost_results$roc_P1304[1],
#   "P1304","Rand Forest",rand_forest_results$roc_P1304[1],
#   "P1304","SVM",svm_results$roc_P1304[1], 
# ) %>% mutate(AUROC = signif(AUROC,3))
# # gtsave(here('figures/model_assessment/ROC_table.png'))
# 
# PRC_data = tribble(
#   ~"Cell Line",~"Model Type",~"AUPRC",
#   "CAF","XGBoost",xgboost_results$prc_CAF[1],
#   "CAF","Rand Forest",rand_forest_results$prc_CAF[1],
#   "CAF","SVM",svm_results$prc_CAF[1],
#   "P0422","XGBoost",xgboost_results$prc_P0422[1],
#   "P0422","Rand Forest",rand_forest_results$prc_P0422[1],
#   "P0422","SVM",svm_results$prc_P0422[1],
#   "P1304","XGBoost",xgboost_results$prc_P1304[1],
#   "P1304","Rand Forest",rand_forest_results$prc_P1304[1],
#   "P1304","SVM",svm_results$prc_P1304[1], 
# ) %>% mutate(AUPRC = signif(AUPRC,3))
# 
# ROC_data %>% left_join(PRC_data) %>%
#   gt() %>%
#   gtsave(here('figures/model_assessment/results_table.png'))
# BerginskiRMisc::trimImage(here('figures/model_assessment/results_table.png'))

ROC_text = tribble(
    ~cell_line,~text,
    "CAF",paste0("Rand Forest: ",sprintf("%.3f",rand_forest_results$CAF$roc[1]), "\n", 
                             "XGBoost: ",sprintf("%.3f",xgboost_results$CAF$roc[1]), "\n", 
                             "SVM: ",sprintf("%.3f",svm_results$CAF$roc[1])),
    "P0422",paste0("Rand Forest: ",sprintf("%.3f",rand_forest_results$P0422$roc[1]), "\n", 
                                 "XGBoost: ",sprintf("%.3f",xgboost_results$P0422$roc[1]), "\n", 
                                 "SVM: ",sprintf("%.3f",svm_results$P0422$roc[1])),
    "P0411",paste0("Rand Forest: ",sprintf("%.3f",rand_forest_results$P0411$roc[1]), "\n", 
                                 "XGBoost: ",sprintf("%.3f",xgboost_results$P0411$roc[5]), "\n", 
                                 "SVM: ",sprintf("%.3f",svm_results$P0411$roc[1])),
)%>% mutate(cell_line = case_when(
    cell_line == "P0422" ~ "P0422-T1",
    cell_line == "P0411" ~ "P0411-T1",
    cell_line == "CAF" ~ "P0119-T1 CAF"
 ))


ROC_plots = ggplot(all_ROC_vals, aes(x=fpr, y=tpr,color=model_type)) + 
    geom_abline(intercept = 0,slope = 1, alpha=0.5,linetype=2) +
    geom_line(alpha=0.75) + 
    geom_text(mapping=aes(x=1,y=0,label=text),data=ROC_text,color='black',vjust="inward",hjust="inward",size=3.5) +
    # geom_segment(x=0,y=0,xend=1,yend=1) +
    xlim(c(0,1)) + ylim(c(0,1)) +
    labs(x="False Positive Rate",y="True Positive Rate",color="Model Type") +
    BerginskiRMisc::theme_berginski() +
    facet_wrap(~factor(cell_line, levels=c( 'P0411-T1', 'P0422-T1','P0119-T1 CAF')))

PRC_text = tribble(
    ~cell_line,~text,
    "CAF",paste0("Rand Forest: ",sprintf("%.3f",rand_forest_results$CAF$prc[1]), "\n", 
                             "XGBoost: ",sprintf("%.3f",xgboost_results$CAF$prc[1]), "\n", 
                             "SVM: ",sprintf("%.3f",svm_results$CAF$prc[1])),
    "P0422",paste0("Rand Forest: ",sprintf("%.3f",rand_forest_results$P0422$prc[1]), "\n", 
                                 "XGBoost: ",sprintf("%.3f",xgboost_results$P0422$prc[1]), "\n", 
                                 "SVM: ",sprintf("%.3f",svm_results$P0422$prc[1])),
    "P0411",paste0("Rand Forest: ",sprintf("%.3f",rand_forest_results$P0411$prc[1]), "\n", 
                                 "XGBoost: ",sprintf("%.3f",xgboost_results$P0411$prc[5]), "\n", 
                                 "SVM: ",sprintf("%.3f",svm_results$P0411$prc[1])),
)%>% mutate(cell_line = case_when(
    cell_line == "P0422" ~ "P0422-T1",
    cell_line == "P0411" ~ "P0411-T1",
    cell_line == "CAF" ~ "P0119-T1 CAF"
 ))

PRC_guessing_lines = tribble(
    ~cell_line,~yintercept,
    "CAF",mean(as.logical(rand_forest_res_CAF$splits[[1]]$data$viability_90)),
    "P0422",mean(as.logical(rand_forest_res_P0422$splits[[1]]$data$viability_90)),
    "P0411",mean(as.logical(rand_forest_res_P0411$splits[[1]]$data$viability_90)),
)%>% mutate(cell_line = case_when(
    cell_line == "P0422" ~ "P0422-T1",
    cell_line == "P0411" ~ "P0411-T1",
    cell_line == "CAF" ~ "P0119-T1 CAF"
 ))

PRC_plots = ggplot(all_PRC_vals, aes(x=recall, y=precision,color=model_type)) + 
    geom_line(alpha = 0.75) + 
    geom_text(mapping=aes(x=1,y=1,label=text),data=PRC_text,color='black',vjust="inward",hjust="inward",size=3.5) +
    geom_hline(data=PRC_guessing_lines, mapping = aes(yintercept = yintercept), alpha=0.5, linetype=2) +
    labs(x="Recall",y="Precision",color="Model Type") +
    BerginskiRMisc::theme_berginski() +
  #facet_wrap(~cell_line)
    facet_wrap(~factor(cell_line, levels=c( 'P0411-T1', 'P0422-T1','P0119-T1 CAF')))

model_assess_plots = ROC_plots / PRC_plots
dir.create(here('figures/model_assessment'), recursive = T, showWarnings = F)
ggsave(here('figures/model_assessment/roc_prc.png'), width=6*1.4,height=4*1.4)
```

    ## Warning: Removed 3 rows containing missing values (`geom_line()`).

``` r
BerginskiRMisc::trimImage(here('figures/model_assessment/roc_prc.png'))
```