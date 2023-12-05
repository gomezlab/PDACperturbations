Build EDA Figures
================
Matthew Berginski
2023-11-14

# Read In Combined Klaeger/Synergy Data and Organize

This chunk of code reads in the pre-processed Klaeger/Synergy data.

``` r
klaeger_data_matches_full = read_rds(here('results/klaeger_screen_for_regression.rds')) 
```

``` r
drug_viability_rank = klaeger_data_matches_full %>%
    group_by(drug) %>%
    summarise(mean_viability = mean(viability)) %>%
    arrange(desc(mean_viability))

cell_line_viability_rank = klaeger_data_matches_full %>%
    group_by(cell_line) %>%
    summarise(mean_viability = mean(viability)) %>%
    arrange(desc(mean_viability))

klaeger_data_matches_full = klaeger_data_matches_full %>%
    mutate(drug = fct_relevel(as.factor(drug), drug_viability_rank$drug),
                 cell_line = fct_relevel(as.factor(cell_line), cell_line_viability_rank$cell_line))
```

# Exploritory Data Analysis

## Cell Viability Visualizations

``` r
ggplot(klaeger_data_matches_full, aes(x = viability)) + 
    geom_histogram(breaks=c(seq(0,130,by=5))) +
    geom_vline(aes(xintercept = 90), color='green',linewidth=2) +
    labs(x="Cell Viability", y="Number of Drug/Cell Line Combos") +
    theme_berginski()
```

![](build_screen_EDA_figures_files/figure-gfm/distribution%20of%20viability-1.png)<!-- -->

``` r
ggsave(here('figures/EDA_plots/all_cell_via_with_thresh.png'),width=4,height=2.75)
BerginskiRMisc::trimImage(here('figures/EDA_plots/all_cell_via_with_thresh.png'))
```

``` r
library(ggridges)
klaeger_data_matches_full %>%
  mutate(cell_line = fct_relevel(cell_line, c("P0119-T1 CAF","P0422-T1","P0411-T1")))%>%
  ggplot(aes(x=viability,y=cell_line)) +
    geom_density_ridges() +
    # scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = expand_scale(mult = c(0.01, .7))) +
    labs(x="Cell Viability",y='Cell Line') +
    # theme_ridges() +
    theme_berginski()+
  theme(axis.text=element_text(size=9),axis.title=element_text(size=14))
```

    ## Warning: `expand_scale()` was deprecated in ggplot2 3.3.0.
    ## ℹ Please use `expansion()` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Picking joint bandwidth of 3.1

![](build_screen_EDA_figures_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
dir.create(here('figures/EDA_plots'), recursive = T, showWarnings = F)
ggsave(here('figures/EDA_plots/cell_line_viability_plots.png'),height = 6, width = 4)
```

    ## Picking joint bandwidth of 3.1

``` r
trimImage(here('figures/EDA_plots/cell_line_viability_plots.png'))

cell_line_binary_freq<- klaeger_data_matches_full %>% group_by(cell_line) %>% select(viability, cell_line) %>% mutate(total_rows = n())%>%
  filter(viability <90) %>% mutate(perc_below = n())
```

``` r
library(ggridges)
ggplot(klaeger_data_matches_full, aes(x=viability,y=drug)) +
    geom_density_ridges() +
    labs(x="Cell Viability",y='Compound') +
    # scale_y_discrete(expand = expand_scale(mult = c(0.01, .7))) +
    theme_berginski()+
  theme(axis.text=element_text(size=6), axis.title=element_text(size=10))
```

    ## Picking joint bandwidth of 4.77

![](build_screen_EDA_figures_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave(here('figures/EDA_plots/compound_viability_plots.png'),height = 4.6, width=4.5)
```

    ## Picking joint bandwidth of 4.77

``` r
trimImage(here('figures/EDA_plots/compound_viability_plots.png'))
```

``` r
treatment_variability = klaeger_data_matches_full %>% 
    group_by(cell_line,drug,concentration_M) %>% 
    summarise(viability_sd = sd(viability),
                        viability_mean = mean(viability)) %>% 
    arrange(desc(viability_sd))
```

    ## `summarise()` has grouped output by 'cell_line', 'drug'. You can override using
    ## the `.groups` argument.

``` r
ggplot(treatment_variability, aes(x=viability_sd,y=fct_rev(as_factor(cell_line)))) +
    geom_density_ridges() +
    labs(x="Standard Deviation in Cell Viability By Treatment",y='') +
    theme_berginski()
```

    ## Picking joint bandwidth of 0.713

![](build_screen_EDA_figures_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
