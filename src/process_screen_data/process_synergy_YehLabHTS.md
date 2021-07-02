Process YehLabHTS Data for ML
================
Matthew Berginski
2021-06-16

``` r
#loading and then reoutputing the synergy screen results
load(here('data/normalized_data.RData'))

write_rds(drug_results, here('data/normalized_data.rds'), compress = 'gz')
```

``` r
synergy = read_rds(here('data/normalized_data.rds')) %>%
    clean_names() %>%
    filter(cell_line != "P140710N1") %>%
    rename(dose_anchor_m = anchor_dose,
                 dose_compound_m = dose, 
                 viability = normalized) %>% 
    mutate(cell_line = case_when(
        cell_line == "P100422" ~ "P1004",
        cell_line == "P130411" ~ "P1304",
        cell_line == "P170119" ~ "CAF"
    )) %>%
    write_rds(here('results/synergy_combined.rds'), compress = 'gz')
```