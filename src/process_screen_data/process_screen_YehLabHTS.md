Process YehLabHTS Data for ML
================
Matthew Berginski
2023-11-14

``` r
#loading and then reoutputing the screen results
load(here('data/normalized_data.RData'))

write_rds(drug_results, here('data/normalized_data.rds'), compress = 'gz')
```

``` r
screen = read_rds(here('data/normalized_data.rds')) %>%
    clean_names() %>%
    # Remove the NAF line due to problems with line verification
    filter(cell_line != "P140710N1") %>%
    rename(dose_anchor_m = anchor_dose,
                 dose_compound_m = dose, 
                 viability = normalized) %>% 
    mutate(cell_line = case_when(
        cell_line == "P100422" ~ "P0422-T1",
        cell_line == "P130411" ~ "P0411-T1",
        cell_line == "P170119" ~ "P0119-T1 CAF")) %>%
    write_rds(here('results/screen_combined.rds'), compress = 'gz')
```
