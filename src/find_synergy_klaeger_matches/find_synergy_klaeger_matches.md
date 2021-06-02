Find Synergy-Klaeger Drug Matches
================
Matthew Berginski
2021-06-02

``` r
klaeger = read_rds(here('results/klaeger_full_tidy.rds')) %>%
    mutate(drug = trimws(drug))

synergy = read_rds(here('results/synergy_combined.rds')) %>%
    mutate(anchor = trimws(anchor), compound = trimws(compound))

klaeger_drugs = unique(klaeger$drug) %>% sort()
synergy_drugs = unique(c(synergy$anchor,synergy$compound)) %>% sort()

perfect_matches = klaeger_drugs[klaeger_drugs %in% synergy_drugs]

#remove the direct perfect matches from the above lists and then dump each list
#out to CSV for manual matching

klaeger_drugs_no_match = klaeger_drugs[! klaeger_drugs %in% perfect_matches]
synergy_drugs_no_match = synergy_drugs[! synergy_drugs %in% perfect_matches]

write_csv(as.data.frame(klaeger_drugs_no_match), here('src/find_synergy_klaeger_matches/no_clear_match_klaeger_drug_list.csv'))
write_csv(as.data.frame(synergy_drugs_no_match), here('src/find_synergy_klaeger_matches/no_clear_match_synergy_drug_list.csv'))
```

From here I opened each individual drug list file and manually matched
the drug names, producing “manual\_match\_list.csv”. Now I’ll open this
file and run a few checks on the manual results and then combine them
with the perfect match list from above.

``` r
manual_matches = read_csv(here('src/find_synergy_klaeger_matches/manual_match_list.csv'))
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   klaeger_drugs = col_character(),
    ##   synergy_drugs = col_character()
    ## )

``` r
assertthat::assert_that(all(manual_matches$synergy_drugs %in% synergy_drugs_no_match))
```

    ## [1] TRUE

``` r
#Dinaciclib is annotated in two ways in the synergy screen, as "Dinaciclib",
#which has a perfect match in the Klaeger set and as "Dinaciclib (SCH727965)"
#which doesn't have a perfect match in the Klaeger set. These are clearly the
#same drug though, so I added "Dinaciclib" to the manual match list, but since
#this is a perfect match in the Klaeger set it was removed from the no match
#list. This causes the assert to fail though, so I've filtered the Dinaciclib to
#"Dinaciclib (SCH727965)" match out for this check.
assertthat::assert_that(all(manual_matches %>% 
                                                            filter(klaeger_drugs != "Dinaciclib") %>% 
                                                            pull(klaeger_drugs) %in% klaeger_drugs_no_match))
```

    ## [1] TRUE

``` r
all_matches = bind_rows(
    manual_matches,
    data.frame(klaeger_drugs = perfect_matches,
                         synergy_drugs = perfect_matches)
) %>% arrange(synergy_drugs) %>%
    write_csv(here('src/find_synergy_klaeger_matches/klaeger_synergy_drug_matches.csv'))

unmatched_klaeger = data.frame(
    klaeger_drugs = klaeger_drugs[!klaeger_drugs %in% all_matches$klaeger_drugs])
unmatched_synergy = data.frame(
    synergy_drugs = synergy_drugs[!synergy_drugs %in% all_matches$synergy_drugs])

write_csv(unmatched_klaeger,here('src/find_synergy_klaeger_matches/unmatched_klaeger.csv'))
write_csv(unmatched_synergy,here('src/find_synergy_klaeger_matches/unmatched_synergy.csv'))
```
