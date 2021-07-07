Process Data for ML
================
Matthew Berginski
2021-06-02

``` r
drug_matches = read_csv(here('src/find_synergy_klaeger_matches/klaeger_synergy_drug_matches.csv'))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   klaeger_drugs = col_character(),
    ##   synergy_drugs = col_character()
    ## )

``` r
klaeger = read_rds(here('results/klaeger_full_tidy.rds')) %>%
    mutate(drug = trimws(drug))

klaeger_compound_matches = klaeger %>%
    #filter to just the drugs that are in the match list and also toss out the zero
    #concentration values
    filter(drug %in% drug_matches$klaeger_drugs, concentration_M != 0)

#Make sure all the expected matches are found in the filtered Klaeger set
assertthat::assert_that(all(klaeger_compound_matches$drug %in% drug_matches$klaeger_drugs))
```

    ## [1] TRUE

``` r
klaeger_wide = klaeger_compound_matches %>%
    pivot_wider(names_from = "gene_name", values_from = "relative_intensity")
```

``` r
synergy = read_rds(here('results/synergy_combined.rds')) %>%
    mutate(anchor = trimws(anchor), compound = trimws(compound))

synergy_compound_matches = bind_rows(
    synergy %>%
        filter(compound %in% drug_matches$synergy_drugs, dose_anchor_m == 0, dose_compound_m != 0),
    synergy %>%
        filter(anchor %in% drug_matches$synergy_drugs, dose_compound_m == 0, dose_anchor_m != 0)
) %>% mutate(drug = case_when(
    dose_anchor_m == 0 ~ compound,
    dose_compound_m == 0 ~ anchor
)) %>% mutate(dose = case_when(
    dose_anchor_m == 0 ~ dose_compound_m,
    dose_compound_m == 0 ~ dose_anchor_m
))

#Make sure all the expected matches are found in the filtered Klaeger set
assertthat::assert_that(all(synergy_compound_matches$drug %in% drug_matches$synergy_drugs))
```

    ## [1] TRUE

``` r
synergy_compound_to_join = synergy_compound_matches %>%
    select(drug, dose, cell_line, viability) %>%
    left_join(drug_matches, by=c('drug' = 'synergy_drugs')) %>%
    select(-drug)
```

``` r
klaeger_synergy = klaeger_wide %>%
    left_join(synergy_compound_to_join, by=c('drug'='klaeger_drugs','concentration_M'='dose')) %>%
    select(drug,concentration_M,viability,cell_line,everything()) %>%
    filter(!is.na(viability))

klaeger_synergy_tidy = klaeger_synergy %>%
    pivot_longer(-c('drug','concentration_M','viability','cell_line'), names_to = "gene_name", values_to = 'relative_intensity')

no_gene_variation = klaeger_synergy_tidy %>%
    group_by(gene_name) %>%
    summarise(gene_sd = sd(relative_intensity)) %>%
    filter(gene_sd == 0)

klaeger_synergy_gene_variation = klaeger_synergy %>%
    #toss out the genes with no variation
    select(-one_of(no_gene_variation$gene_name))

write_rds(klaeger_synergy_gene_variation,here('results/klaeger_synergy_for_regression.rds'))

klaeger_synergy_gene_variation %>%
    mutate(viability_90 = viability <= 90) %>%
    select(-viability) %>%
    select(drug,concentration_M,viability_90,cell_line,everything()) %>% 
    write_rds(here('results/klaeger_synergy_for_classification_90.rds'))
```
