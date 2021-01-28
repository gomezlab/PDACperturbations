# Data files

## KlaegerLong.csv
Data processed from the "Kinobeads" worksheet from "Klaeger_Science _ 2017 Supplementary Table 2 Target Lists.xlsx"

Processing and notes:
- 0.0 values replaced by a pseudocount
- NA values replaced by a pseudocount
- pseudocount is calculated as 1/4 of the minimum observed non-zero value
- values were log2 transformed after pseudocounts were applied
- some columns contain groups of kinases that were treated as indistinguishable in terms of kinobead binding
  - e.g., CSNK2A1;CSNK2A3
- kinase names are original and not mapped to any identifier (e.g., HGNC)
- compound names are similarly left in the original form

## KlaegerWide.csv
This file is the wideform version of the above KlaegerLong.csv file, with the new "missing" values generated during the transformation converted to the log2 transformed pseudocount.
