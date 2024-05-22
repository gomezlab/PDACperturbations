Process Data for ML
================
Matthew Berginski
2024-03-02

The following section reads in the klaeger data and filters to the
collection of compounds that are also in the screen screen. The data is
also rearranged to match up with the expectation that the relative
intensity values for every gene will be in a seperate column.

``` r
drug_matches = read_csv(here('src/find_screen_klaeger_matches/klaeger_screen_drug_matches.csv'))
```

    ## Rows: 62 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): klaeger_drugs, screen_drugs
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

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

The screen data is processed in much the same way, toss out all the
compounds not in the match list and prep the data for combining with
Klaeger.

``` r
screen = read_rds(here('results/screen_combined.rds')) %>%
    mutate(anchor = trimws(anchor), compound = trimws(compound))

screen_compound_matches = bind_rows(
    screen %>%
        filter(compound %in% drug_matches$screen_drugs, dose_anchor_m == 0, dose_compound_m != 0),
    screen %>%
        filter(anchor %in% drug_matches$screen_drugs, dose_compound_m == 0, dose_anchor_m != 0)
) %>% mutate(drug = case_when(
    dose_anchor_m == 0 ~ compound,
    dose_compound_m == 0 ~ anchor
)) %>% mutate(dose = case_when(
    dose_anchor_m == 0 ~ dose_compound_m,
    dose_compound_m == 0 ~ dose_anchor_m
))

#Make sure all the expected matches are found in the filtered Klaeger set
assertthat::assert_that(all(screen_compound_matches$drug %in% drug_matches$screen_drugs))
```

    ## [1] TRUE

``` r
screen_compound_to_join = screen_compound_matches %>%
    select(drug, dose, cell_line, viability) %>%
    left_join(drug_matches, by=c('drug' = 'screen_drugs')) %>%
    select(-drug)
```

Now we join up the screen and Klaeger, followed by a filter to toss out
any gene that doesn’t vary at all across the matched compounds. Then
save out the data frame with the raw viability values and produce a
second data set with the viability converted into a binary above/below
90 threshold.

``` r
klaeger_screen = klaeger_wide %>%
    left_join(screen_compound_to_join, by=c('drug'='klaeger_drugs','concentration_M'='dose')) %>%
    select(drug,concentration_M,viability,cell_line,everything()) %>%
    filter(!is.na(viability)) %>% 
    arrange(cell_line,drug,concentration_M)

klaeger_screen_tidy = klaeger_screen %>%
    pivot_longer(-c('drug','concentration_M','viability','cell_line'), names_to = "gene_name", values_to = 'relative_intensity')

no_gene_variation = klaeger_screen_tidy %>%
    group_by(gene_name) %>%
    summarise(gene_sd = sd(relative_intensity)) %>%
    filter(gene_sd == 0)

klaeger_screen_gene_variation = klaeger_screen %>%
    #toss out the genes with no variation
    select(-one_of(no_gene_variation$gene_name))

write_rds(klaeger_screen_gene_variation,here('results/klaeger_screen_for_regression.rds'), compress='gz')

klaeger_screen_below_90 = klaeger_screen_gene_variation %>%
    mutate(viability_90 = viability <= 90) %>%
    select(-viability) %>%
    select(drug,concentration_M,viability_90,cell_line,everything()) %>% 
    write_rds(here('results/klaeger_screen_for_classification_90.rds'), compress='gz')

klaeger_for_prediction = klaeger %>% 
    filter(!drug %in% drug_matches$klaeger_drugs, 
                 concentration_M != 0,
                 ! gene_name %in% no_gene_variation$gene_name)

#After filtering out all the genes without variance in the klaeger set, there
#are probably a few drugs that don't hit any of the genes included in the model.
#So I'll find those and filter them in the next step.
klaeger_prediction_no_variation = klaeger_for_prediction %>% 
    group_by(drug) %>% 
    summarise(relative_sd = sd(relative_intensity)) %>% 
    filter(relative_sd == 0)

klaeger_wide_for_prediction = klaeger_for_prediction %>%
    filter(! drug %in% klaeger_prediction_no_variation$drug) %>%
    pivot_wider(names_from = "gene_name", values_from = "relative_intensity") %>%
    write_rds(here('results/klaeger_wide_for_prediction.rds'), compress='gz')
```

Finally, in prep for running cross validation across
leave-one-compound-out:

``` r
cell_line_compound_splits = list()

for (this_cell_line in unique(klaeger_screen_below_90$cell_line)) {
    
    klaeger_data_cell_line = klaeger_screen_below_90 %>%
        filter(cell_line == this_cell_line)
    
    klaeger_data_cell_line = klaeger_data_cell_line %>%
        select(-cell_line,-concentration_M) %>%
        mutate(viability_90 = as.factor(viability_90))
    
    splits = list()
    index = 1
    id = c()
    for (exclude_compound in unique(klaeger_data_cell_line$drug)) {
        assessment_ids = which(klaeger_data_cell_line$drug == exclude_compound)
        analysis_ids = which(klaeger_data_cell_line$drug != exclude_compound)
        
        splits[[index]] = make_splits(list("analysis" = analysis_ids,"assessment" = assessment_ids),
                                                                    klaeger_data_cell_line %>% select(-drug))
        index = index + 1
        
        id = c(id,exclude_compound)
    }
    
    cell_line_compound_splits[[this_cell_line]] = new_rset(
        splits = splits,
        ids = id,
        attrib = paste0("Per compound cv splits for ", this_cell_line),
        subclass = c("vfold_cv", "rset")
    )
}

write_rds(cell_line_compound_splits, here('results/klaeger_screen_classification_90_CV_split.rds'), compress = 'gz')
```

``` r
cell_line_compound_splits = list()

for (this_cell_line in unique(klaeger_screen_gene_variation$cell_line)) {
    
    klaeger_data_cell_line = klaeger_screen_gene_variation %>%
        filter(cell_line == this_cell_line) %>%
        select(-cell_line,-concentration_M)
    
    splits = list()
    index = 1
    id = c()
    for (exclude_compound in unique(klaeger_data_cell_line$drug)) {
        assessment_ids = which(klaeger_data_cell_line$drug == exclude_compound)
        analysis_ids = which(klaeger_data_cell_line$drug != exclude_compound)
        
        splits[[index]] = make_splits(list("analysis" = analysis_ids,"assessment" = assessment_ids),
                                                                    klaeger_data_cell_line %>% select(-drug))
        index = index + 1
        
        id = c(id,exclude_compound)
    }
    
    cell_line_compound_splits[[this_cell_line]] = new_rset(
        splits = splits,
        ids = id,
        attrib = paste0("Per compound cv splits for ", this_cell_line),
        subclass = c("vfold_cv", "rset")
    )
    print(names(cell_line_compound_splits))
}
```

    ## [1] "P0119-T1 CAF"
    ## [1] "P0119-T1 CAF" "P0411-T1"    
    ## [1] "P0119-T1 CAF" "P0411-T1"     "P0422-T1"

``` r
write_rds(cell_line_compound_splits, here('results/klaeger_screen_regression_CV_split.rds'), compress = 'gz')
```
