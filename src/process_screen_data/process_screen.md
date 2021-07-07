Process Synergy
================
Matthew Berginski
2021-06-09

``` r
synergy_files = list.files(here('data','SynergyData'), full.names = T)

synergy_data = data.frame()

col_types = cols(
    Ind = col_character(),
    Anchor = col_character(),
    Compound = col_character(),
    `DoseAnchor(M)` = col_double(),
    `DoseCompound(M)` = col_double(),
    plate = col_double(),
    row = col_double(),
    col = col_double(),
    Viability = col_double()
)

for (this_file in synergy_files) {
    cell_line = case_when(
        str_detect(this_file,"CAF") ~ "CAF",
        str_detect(this_file,"NAF") ~ "NAF",
        str_detect(this_file,"P1004") ~ "P1004",
        str_detect(this_file,"P1304") ~ "P1304",
        TRUE ~ "No_Line"
    )
    
    if (cell_line == "No_Line") {
        next;
    }
    
    synergy_data = rbind(synergy_data,
                                             read_delim(this_file, delim=";", col_types = col_types) %>%
                                                mutate(cell_line = cell_line))
}

synergy_data = synergy_data %>% 
    clean_names() %>%
    #there were downstream problems with the NAF line used in the screen, removing
    #from data set
    filter(cell_line != "NAF") %>%
    # select(-ind,-plate,-row,-col) %>%
    write_rds(here('results/synergy_combined.rds'))
```