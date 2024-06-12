Prep Validation Screen Data
================
Matthew Berginski
2023-11-10

# Read in plate data

All the data from the plate reader is collected into a set of excel
spreadsheets. The luminescence values show up at the bottom of the
spreadsheet. I’ve manually matched up the names of the excel files to
the appropriate cell line names and compounds in each plate. There is a
text file with the same information in the same data folder I’ve read
these files from.

The plate is organized with the concentrations decreasing over each
column of the plate, following the concentration used in Klaeger. The
last two columns in each plate are DMSO controls.

``` r
read_xlsx_plate_data <- function(file_name,these_compounds,this_cell_line,skip_lines = 16) {
    
    temp_plate = read_excel(file_name, skip = skip_lines) %>%
        select(-A, -`...2`) %>%
        dplyr::slice(1:6) %>%
        mutate(row = c(seq(6))) %>% 
        mutate(compound = c(rep(these_compounds[1],3),rep(these_compounds[2],3))) %>%
        pivot_longer(-c(compound,row), names_to = "concentration_M",values_to = "lum") %>%
        mutate(cell_line = this_cell_line, plate_name = basename(file_name)) %>%
        identity()
    
    DMSO_vals = temp_plate %>%
        filter(concentration_M == "...11" | concentration_M == "...12") %>%
        mutate(concentration_M = 0) %>%
        mutate(compound = "DMSO")
    
    temp_plate = temp_plate %>%
        filter(concentration_M != "...11", concentration_M != "...12") %>%
        mutate(concentration_M = case_when(
            concentration_M == "...3" ~ 3e-5,
            concentration_M == "...4" ~ 3e-6,
            concentration_M == "...5" ~ 1e-6,
            concentration_M == "...6" ~ 3e-7,
            concentration_M == "...7" ~ 1e-7,
            concentration_M == "...8" ~ 3e-8,
            concentration_M == "...9" ~ 1e-8,
            concentration_M == "...10" ~ 3e-9,
            T ~ NA_real_
        )) %>%
        bind_rows(DMSO_vals)
    
    return(temp_plate)
}

read_csv_plate_data <- function(file_name,these_compounds,this_cell_line) {
    
    temp_plate = read_csv(file_name, col_names = FALSE, show_col_types = F) %>%
        select(-X1, -X12) %>%
        dplyr::slice(2:7) %>%
        mutate(row = 1:6) %>% 
        mutate(compound = c(rep(these_compounds[1],3),rep(these_compounds[2],3))) %>%
        pivot_longer(-c(compound,row), names_to = "concentration_M",values_to = "lum") %>%
        mutate(cell_line = this_cell_line, plate_name = basename(file_name)) %>%
        mutate(lum = as.numeric(lum)) %>% 
        identity()
    
    DMSO_vals = temp_plate %>%
        filter(concentration_M == "X10" | concentration_M == "X11") %>%
        mutate(concentration_M = 0) %>%
        mutate(compound = "DMSO")
    
    temp_plate = temp_plate %>%
        filter(concentration_M != "X10", concentration_M != "X11") %>%
        mutate(concentration_M = case_when(
            concentration_M == "X2" ~ 3e-5,
            concentration_M == "X3" ~ 3e-6,
            concentration_M == "X4" ~ 1e-6,
            concentration_M == "X5" ~ 3e-7,
            concentration_M == "X6" ~ 1e-7,
            concentration_M == "X7" ~ 3e-8,
            concentration_M == "X8" ~ 1e-8,
            concentration_M == "X9" ~ 3e-9,
            T ~ NA_real_
        )) %>%
        bind_rows(DMSO_vals)
    
    return(temp_plate)
}

plate_vals_part3 = bind_rows(
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/119_1.csv'),c("Masitinib","Ripasudil"),"CAF"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/119_2.csv'),c("AT-13148","RGB-286638"),"CAF"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/119_3.csv'),c("PHA-793887","AT-9283"),"CAF"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/119_4.csv'),c("Lestaurtinib","KW-2449"),"CAF"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/119_5.csv'),c("K-252a","PF-03814735"),"CAF"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/119_6.csv'),c("XL-228","JANEX-1"),"CAF"),
    
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/411_1.csv'),c("Masitinib","Ripasudil"),"P0411"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/411_2.csv'),c("AT-13148","RGB-286638"),"P0411"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/411_3.csv'),c("PHA-793887","AT-9283"),"P0411"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/411_4.csv'),c("Lestaurtinib","KW-2449"),"P0411"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/411_5.csv'),c("K-252a","PF-03814735"),"P0411"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/411_6.csv'),c("XL-228","JANEX-1"),"P0411"),
    
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/422_1.csv'),c("Masitinib","Ripasudil"),"P0422"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/422_2.csv'),c("AT-13148","RGB-286638"),"P0422"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/422_3.csv'),c("PHA-793887","AT-9283"),"P0422"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/422_4.csv'),c("Lestaurtinib","KW-2449"),"P0422"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/422_5.csv'),c("K-252a","PF-03814735"),"P0422"),
    read_csv_plate_data(here('data/PDAC_validation_screen_part3/422_6.csv'),c("XL-228","JANEX-1"),"P0422")
)

plate_vals_part2 = bind_rows(
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/119_1.xlsx'),c("Masitinib","Ripasudil"),"CAF"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/119_2.xlsx'),c("AT-13148","RGB-286638"),"CAF"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/119_3.xlsx'),c("PHA-793887","AT-9283"),"CAF"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/119_4.xlsx'),c("Lestaurtinib","KW-2449"),"CAF"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/119_5.xlsx'),c("K-252a","PF-03814735"),"CAF"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/119_6.xlsx'),c("XL-228","JANEX-1"),"CAF"),
    
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/411_1.xlsx'),c("Masitinib","Ripasudil"),"P0411"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/411_2.xlsx'),c("AT-13148","RGB-286638"),"P0411"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/411_3.xlsx'),c("PHA-793887","AT-9283"),"P0411"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/411_4.xlsx'),c("Lestaurtinib","KW-2449"),"P0411"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/411_5.xlsx'),c("K-252a","PF-03814735"),"P0411"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/411_6.xlsx'),c("XL-228","JANEX-1"),"P0411"),
    
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/422_1.xlsx'),c("Masitinib","Ripasudil"),"P0422"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/422_2.xlsx'),c("AT-13148","RGB-286638"),"P0422"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/422_3.xlsx'),c("PHA-793887","AT-9283"),"P0422"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/422_4.xlsx'),c("Lestaurtinib","KW-2449"),"P0422"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/422_5.xlsx'),c("K-252a","PF-03814735"),"P0422"),
    read_xlsx_plate_data(here('data/PDAC_validation_screen_part2/422_6.xlsx'),c("XL-228","JANEX-1"),"P0422")
)
```

    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## New names:
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`

``` r
plate_order = c("Masitinib","Ripasudil","AT-13148","RGB-286638","PHA-793887",
                                "AT-9283","Lestaurtinib","KW-2449","K-252a","PF-03814735",
                                "XL-228","JANEX-1")

plate_vals_part2 = plate_vals_part2 %>%
    mutate(compound = fct_relevel(as.factor(compound), plate_order))

plate_vals_part3 = plate_vals_part3 %>%
    mutate(compound = fct_relevel(as.factor(compound), plate_order))

if (any(is.na(plate_vals_part2$concentration_M)) | any(is.na(plate_vals_part3$concentration_M))) {
    stop("Found missing concentration value, this shouldn't happen.")
}
```

# Per Plate DMSO Boxplots

``` r
all_DMSO_vals = plate_vals_part3 %>%
    filter(compound == "DMSO") %>%
    mutate(plate_name = paste0("p3_", plate_name)) %>%
    bind_rows(
        plate_vals_part2 %>%
            filter(compound == "DMSO") %>%
            mutate(plate_name = paste0("p2_", plate_name))
    )


ggplot(all_DMSO_vals, aes(x=plate_name,y=lum)) + 
    geom_boxplot() +
    # ylim(c(min(plate_vals$lum), NA)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
```

![](prep_validation_data_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

# Per Plate Diagnostic Plotting/Analysis

``` r
DMSO_checks_part2 = plate_vals_part2 %>% 
    bind_rows(plate_vals_part3) %>%
    filter(compound == "DMSO") %>%
    group_by(plate_name) %>%
    nest() %>%
    mutate(DMSO_1 = map_dbl(data, ~ tidy(t.test(lum ~ target_row, 
                                                                                            data = .x %>% 
                                                                                                mutate(target_row = ifelse(row == 1, 1, 0))))$p.value)) %>%
    mutate(DMSO_2 = map_dbl(data, ~ tidy(t.test(lum ~ target_row, 
                                                                                            data = .x %>% 
                                                                                                mutate(target_row = ifelse(row == 2, 1, 0))))$p.value)) %>%
    mutate(DMSO_3 = map_dbl(data, ~ tidy(t.test(lum ~ target_row, 
                                                                                            data = .x %>% 
                                                                                                mutate(target_row = ifelse(row == 3, 1, 0))))$p.value)) %>%
    mutate(DMSO_4 = map_dbl(data, ~ tidy(t.test(lum ~ target_row, 
                                                                                            data = .x %>% 
                                                                                                mutate(target_row = ifelse(row == 4, 1, 0))))$p.value)) %>%
    mutate(DMSO_5 = map_dbl(data, ~ tidy(t.test(lum ~ target_row, 
                                                                                            data = .x %>% 
                                                                                                mutate(target_row = ifelse(row == 5, 1, 0))))$p.value)) %>%
    mutate(DMSO_6 = map_dbl(data, ~ tidy(t.test(lum ~ target_row, 
                                                                                            data = .x %>% 
                                                                                                mutate(target_row = ifelse(row == 6, 1, 0))))$p.value)) %>%
    select(-data) %>%
    pivot_longer(-plate_name, values_to = "p.value", names_to = "DMSO_row") %>%
    mutate(p.value.adj = p.adjust(p.value, method="fdr"))
```

``` r
treatment_lum_summaries_part2 = plate_vals_part2 %>%
    filter(compound != "DMSO") %>% 
    group_by(cell_line, compound, concentration_M) %>%
    dplyr::mutate(sd = sd(lum), 
                        num_samples = n(),
                        mean_lum = mean(lum),
                        CV = sd/mean_lum)

treatment_lum_summaries_part3 = plate_vals_part3 %>% 
    group_by(cell_line, compound, concentration_M) %>%
    dplyr::mutate(sd = sd(lum), 
                        num_samples = n(),
                        mean_lum = mean(lum),
                        CV = sd/mean_lum)

ggplot(treatment_lum_summaries_part2, aes(x=CV)) + geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](prep_validation_data_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggplot(treatment_lum_summaries_part3, aes(x=CV)) + geom_histogram()
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](prep_validation_data_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

# Per Plate Normalization

``` r
plate_DMSO_means_part2 = plate_vals_part2 %>%
    filter(compound == "DMSO") %>%
    group_by(plate_name) %>%
    summarise(DMSO_mean = mean(lum),
                        num_DMSO = n())

plate_DMSO_means_part3 = plate_vals_part3 %>%
    filter(compound == "DMSO") %>%
    group_by(plate_name) %>%
    summarise(DMSO_mean = mean(lum),
                        num_DMSO = n())

if (! all(plate_DMSO_means_part2$num_DMSO == max(plate_DMSO_means_part2$num_DMSO))) {
    stop("Number of DMSO values in a plate isn't correct.")
} else {
    plate_DMSO_means_part2 = plate_DMSO_means_part2 %>% select(-num_DMSO)
}

if (! all(plate_DMSO_means_part3$num_DMSO == max(plate_DMSO_means_part3$num_DMSO))) {
    stop("Number of DMSO values in a plate isn't correct.")
} else {
    plate_DMSO_means_part3 = plate_DMSO_means_part3 %>% select(-num_DMSO)
}
```

``` r
qc_filter <- function(data_set, viability_threshold = 120) {
    qc_fails = data_set %>% 
        filter(viability > viability_threshold) %>% 
        select(row,compound,cell_line,plate_name) %>% 
        unique() %>%
        mutate(qc_fail = T)
    
    data_set = data_set %>%
        left_join(qc_fails) %>%
        mutate(qc_fail = ifelse(is.na(qc_fail), F, T))
    
    return(data_set)
}

dir.create(here('results/validation_results/'), showWarnings = F)

plate_vals_plate_norm_part2 = plate_vals_part2 %>%
    left_join(plate_DMSO_means_part2) %>% 
    mutate(viability = lum/DMSO_mean*100) %>%
    qc_filter() %>%
    write_rds(here('results/validation_results/validation_screen_part2.rds'))
```

    ## Joining with `by = join_by(plate_name)`
    ## Joining with `by = join_by(row, compound, cell_line, plate_name)`

``` r
plate_vals_plate_norm_part3 = plate_vals_part3 %>%
    left_join(plate_DMSO_means_part3) %>% 
    mutate(viability = lum/DMSO_mean*100) %>%
    qc_filter() %>%
    write_rds(here('results/validation_results/validation_screen_part3.rds'))
```

    ## Joining with `by = join_by(plate_name)`
    ## Joining with `by = join_by(row, compound, cell_line, plate_name)`

``` r
treatment_via_summaries_part2 = plate_vals_plate_norm_part2 %>% 
    group_by(compound, concentration_M, cell_line) %>%
    summarise(sd = sd(viability),
                        num_samples = n(),
                        mean_via = mean(viability),
                        CV = sd/mean_via)
```

    ## `summarise()` has grouped output by 'compound', 'concentration_M'. You can override using the `.groups`
    ## argument.

``` r
treatment_via_summaries_part3 = plate_vals_plate_norm_part3 %>% 
    group_by(compound, concentration_M, cell_line) %>%
    summarise(sd = sd(viability),
                        num_samples = n(),
                        mean_via = mean(viability),
                        CV = sd/mean_via)
```

    ## `summarise()` has grouped output by 'compound', 'concentration_M'. You can override using the `.groups`
    ## argument.

``` r
b = ggplot(plate_vals_plate_norm_part2 %>% filter(compound != "DMSO"), aes(x=concentration_M, y = viability, color=qc_fail)) +
    geom_hline(aes(yintercept = 90),alpha=0.5) +
    geom_hline(aes(yintercept = 100),alpha=0.5, linetype = 2) +
    geom_hline(aes(yintercept = 120),alpha=0.5, linetype = 2) +
    scale_x_log10() +
    geom_jitter() +
    geom_smooth() +
    BerginskiRMisc::theme_berginski() +
    facet_grid(rows = vars(cell_line), cols = vars(compound)) +
    ggtitle("Replicate 2")

c = ggplot(plate_vals_plate_norm_part3 %>% filter(compound != "DMSO"), aes(x=concentration_M, y = viability, color=qc_fail)) +
    geom_hline(aes(yintercept = 90),alpha=0.5) +
    geom_hline(aes(yintercept = 100),alpha=0.5, linetype = 2) +
    geom_hline(aes(yintercept = 120),alpha=0.5, linetype = 2) +
    scale_x_log10() +
    geom_jitter() +
    geom_smooth() +
    BerginskiRMisc::theme_berginski() +
    facet_grid(rows = vars(cell_line), cols = vars(compound)) +
    ggtitle("Replicate 3")

d = b / c 

ggsave(here('figures/validation_testing/validation_viability_part2_part3_plate_norm.png'),width=20*0.75,height=15*0.75)
```

    ## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'
    ## `geom_smooth()` using method = 'loess' and formula = 'y ~ x'

``` r
BerginskiRMisc::trimImage(here('figures/validation_testing/validation_viability_part2_part3_plate_norm.png'))
```