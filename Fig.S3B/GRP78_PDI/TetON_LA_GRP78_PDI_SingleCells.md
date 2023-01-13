Figure S3B: TetON-Lamin A / GRP78 and PDI
================
Sandra Vidak/Gianluca Pegoraro
January 13th 2023

### Introduction

Columbus screen names:

`171006-TetON-LA-GRP78-PDI_20171006_164850`

`180130-40x-TetON-LA-GRP78-PDI_20180130_123005`

### Analysis Setup

Load required packages.

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.3.6      ✔ purrr   0.3.5 
    ## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.1      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.3      ✔ forcats 0.5.2 
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(fs)
library(Hmisc)
```

    ## Loading required package: lattice
    ## Loading required package: survival
    ## Loading required package: Formula
    ## 
    ## Attaching package: 'Hmisc'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     src, summarize
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     format.pval, units

``` r
library(ggthemes)
library(DescTools) # for Dunnett's Test
```

    ## 
    ## Attaching package: 'DescTools'
    ## 
    ## The following objects are masked from 'package:Hmisc':
    ## 
    ##     %nin%, Label, Mean, Quantile

``` r
library(curl)
```

    ## Using libcurl 7.79.1 with LibreSSL/3.3.6
    ## 
    ## Attaching package: 'curl'
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     parse_date

``` r
source("R/Plotters.R") #Functions needed for plotting
```

Set the palette and the running theme for ggplot2.

### Experimental Metadata

Read plate layouts.

``` r
cell_levs <- c("Uninduced", "Induced")

plate_layouts <- read_tsv("metadata/plate_layout.txt") %>%
  filter(!is.na(cell_line)) %>%
  separate(col = cell_line, 
           into = c("cell_line"), 
           remove = T) %>%
  mutate(cell_line = factor(cell_line, levels = cell_levs))

glimpse(plate_layouts)
```

    ## Rows: 36
    ## Columns: 5
    ## $ screen    <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "171006-TetON-L…
    ## $ row       <dbl> 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, …
    ## $ column    <dbl> 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 5, 5, …
    ## $ marker    <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "PDI",…
    ## $ cell_line <fct> Uninduced, Uninduced, Uninduced, Induced, Induced, Induced, …

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

### Download the data if needed

Download and unzip the Columbus results of the experiments from Figshare
if they have not been already downloaded.

``` r
if(!dir.exists("input")) {
  URL <- "https://figshare.com/ndownloader/files/38669273"
  curl_download(URL, "input.zip")
  unzip("input.zip")
}
```

### Read and Process Columbus data

Recursively search the `data` directory and its subdirectories for files
whose name includes the Glob patterns defined in the chunk above, and
read the cell-level Columbus data from the results text files.

``` r
read_columbus_results <- function(path, glob) {
  dir_ls(path = path,
         recurse = T,
         glob = glob)  %>%
    read_tsv(
      id = "file_name"
    ) %>%
    select(
      screen = ScreenName,
      plate = PlateName,
      well = WellName,
      row = Row,
      column = Column,
      nuc_area = `Nuclei Selected - Nucleus Area [px²]`,
      cyto_area = `Nuclei Selected - Cytoplasm Area [px²]`,
      cell_area = `Nuclei Selected - Cell Area [px²]`,
      nuc_marker_int = `Nuclei Selected - Intensity Nucleus BP600/37 Mean`,
      cyto_marker_int = `Nuclei Selected - Intensity Cytoplasm BP600/37 Mean`,
      ratio_marker_int = `Nuclei Selected - Nuc_Cyto_BP600_Ratio`
    )
}

glob_path <- "*- Nuclei Selected[0].txt"
col_tbl <- read_columbus_results("input", glob_path)

glimpse(col_tbl)
```

    ## Rows: 121,105
    ## Columns: 11
    ## $ screen           <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "171006-…
    ## $ plate            <chr> "Plate1", "Plate1", "Plate1", "Plate1", "Plate1", "Pl…
    ## $ well             <chr> "B14", "B14", "B14", "B14", "B14", "B14", "B14", "B14…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 1…
    ## $ nuc_area         <dbl> 2459, 1802, 1688, 1956, 1720, 1766, 1712, 2058, 2168,…
    ## $ cyto_area        <dbl> 5569, 1822, 1769, 13614, 9869, 11284, 4401, 6526, 559…
    ## $ cell_area        <dbl> 8028, 3624, 3457, 15570, 11589, 13050, 6113, 8584, 77…
    ## $ nuc_marker_int   <dbl> 2591.38, 1223.84, 1313.43, 1407.22, 1776.15, 2386.63,…
    ## $ cyto_marker_int  <dbl> 3412.72, 2498.72, 3226.09, 2910.42, 3375.52, 3737.21,…
    ## $ ratio_marker_int <dbl> 0.759330, 0.489789, 0.407128, 0.483511, 0.526184, 0.6…

Join Columbus data with the plate layout information.

``` r
cell_tbl <- col_tbl %>%
  mutate(sum_marker_int = nuc_marker_int + cyto_marker_int) %>%
  inner_join(plate_layouts,
             by = c("row", "column", "screen")) %>%
  select(screen,
         plate,
         well,
         row,
         column,
         cell_line,
         marker,
         nuc_area:sum_marker_int)

glimpse(cell_tbl)
```

    ## Rows: 29,718
    ## Columns: 14
    ## $ screen           <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "171006-…
    ## $ plate            <chr> "Plate1", "Plate1", "Plate1", "Plate1", "Plate1", "Pl…
    ## $ well             <chr> "B3", "B3", "B3", "B3", "B3", "B3", "B3", "B3", "B3",…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,…
    ## $ cell_line        <fct> Uninduced, Uninduced, Uninduced, Uninduced, Uninduced…
    ## $ marker           <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78",…
    ## $ nuc_area         <dbl> 1910, 1719, 1763, 1951, 1947, 1779, 953, 3494, 1633, …
    ## $ cyto_area        <dbl> 4432, 8506, 4300, 6495, 13152, 3002, 5827, 12589, 752…
    ## $ cell_area        <dbl> 6342, 10225, 6063, 8446, 15099, 4781, 6780, 16083, 91…
    ## $ nuc_marker_int   <dbl> 1297.03, 1515.97, 2715.81, 2518.69, 2433.18, 1473.20,…
    ## $ cyto_marker_int  <dbl> 1759.62, 2454.98, 3408.47, 3626.51, 3129.34, 3295.05,…
    ## $ ratio_marker_int <dbl> 0.737107, 0.617509, 0.796783, 0.694521, 0.777538, 0.4…
    ## $ sum_marker_int   <dbl> 3056.65, 3970.95, 6124.28, 6145.20, 5562.52, 4768.25,…

Calculate number of cells and mean per well for all properties. 4
technical replicates for each of the 3 biological replicates -\> 12
wells for marker/cell line combination.

``` r
well_tbl <- cell_tbl %>%
  group_by(screen,
           well,
           row,
           column,
           cell_line,
           marker) %>%
  summarise(cell_n = n(),
            across(nuc_area:sum_marker_int,
                   list(mean = ~ mean(.x, na.rm = T))))

glimpse(well_tbl)
```

    ## Rows: 36
    ## Columns: 14
    ## Groups: screen, well, row, column, cell_line [36]
    ## $ screen                <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "17…
    ## $ well                  <chr> "B3", "B5", "C3", "C5", "D3", "D5", "H3", "H5", …
    ## $ row                   <dbl> 2, 2, 3, 3, 4, 4, 8, 8, 9, 9, 10, 10, 2, 2, 3, 3…
    ## $ column                <dbl> 3, 5, 3, 5, 3, 5, 3, 5, 3, 5, 3, 5, 3, 5, 3, 5, …
    ## $ cell_line             <fct> Uninduced, Uninduced, Uninduced, Uninduced, Unin…
    ## $ marker                <chr> "GRP78", "PDI", "GRP78", "PDI", "GRP78", "PDI", …
    ## $ cell_n                <int> 912, 937, 966, 950, 978, 915, 1012, 1036, 944, 1…
    ## $ nuc_area_mean         <dbl> 1961.951, 1966.489, 1978.992, 1940.893, 2037.984…
    ## $ cyto_area_mean        <dbl> 6709.349, 6424.189, 6476.345, 6522.003, 6361.741…
    ## $ cell_area_mean        <dbl> 8671.299, 8390.678, 8455.336, 8462.896, 8399.725…
    ## $ nuc_marker_int_mean   <dbl> 2070.3435, 1827.7416, 1996.5161, 1815.9125, 1779…
    ## $ cyto_marker_int_mean  <dbl> 3104.7562, 3521.1667, 2966.3184, 3528.3377, 2702…
    ## $ ratio_marker_int_mean <dbl> 0.6742553, 0.5260615, 0.6805383, 0.5259330, 0.66…
    ## $ sum_marker_int_mean   <dbl> 5175.0996, 5350.3126, 4962.8345, 5344.2502, 4481…

Calculate the mean of the technical replicates for each biological
replicate. Now every marker/cell line combination has an n = 3
biological replicates.

``` r
bioreps_tbl <- well_tbl %>%
  group_by(screen,
           cell_line,
           marker) %>%
  summarise(across(cell_n:sum_marker_int_mean,
                    ~ mean(.x, na.rm = T)))

glimpse(bioreps_tbl)
```

    ## Rows: 12
    ## Columns: 11
    ## Groups: screen, cell_line [6]
    ## $ screen                <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "17…
    ## $ cell_line             <fct> Uninduced, Uninduced, Induced, Induced, Uninduce…
    ## $ marker                <chr> "GRP78", "PDI", "GRP78", "PDI", "GRP78", "PDI", …
    ## $ cell_n                <dbl> 952.0000, 934.0000, 986.0000, 1063.6667, 930.666…
    ## $ nuc_area_mean         <dbl> 1992.975, 1943.537, 1691.595, 1643.935, 1837.211…
    ## $ cyto_area_mean        <dbl> 6515.812, 6589.530, 6449.827, 6062.655, 6153.666…
    ## $ cell_area_mean        <dbl> 8508.787, 8533.067, 8141.422, 7706.590, 7990.877…
    ## $ nuc_marker_int_mean   <dbl> 1948.6810, 1829.8501, 1554.7619, 1604.7220, 451.…
    ## $ cyto_marker_int_mean  <dbl> 2924.4471, 3645.2673, 2641.0669, 3838.8207, 591.…
    ## $ ratio_marker_int_mean <dbl> 0.6740483, 0.5110368, 0.5999345, 0.4248802, 0.79…
    ## $ sum_marker_int_mean   <dbl> 4873.1280, 5475.5855, 4195.8288, 5444.2757, 1042…

### Biological Replicates Level plot for Fig.S3B

![](output/FigS3_B-1.png)<!-- -->

### Calculate Dunnett’s test for the continuous variables.

Define a custom function to run a Dunnett post-hoc test only on the Mean
marker intensity sum (Cyto + Nucleus), using the cell line as the
predictor variable, and fixing WT1 as the negative control. The output
of the Dunnett’s test is then rearranged to a tidy table to make it work
with `dplyr`.

``` r
calc_dunnett <- function(df){
  as.data.frame(as.table(DunnettTest(ratio_marker_int_mean ~ cell_line,
                          control = "Uninduced",
                          data = df)$Uninduced)) %>%
    pivot_wider(names_from = Var2, values_from = Freq) %>%
    rename(comparison = Var1)
}
```

Run the custom function on all the data grouped based on the IF marker
and save the data to a .csv file.

``` r
dunnett_test <- bioreps_tbl %>%
  group_by(marker) %>%
  group_modify(~ calc_dunnett(.x))

write_csv(dunnett_test, "output/dunnett_results.csv")

knitr::kable(dunnett_test, digits = 3)
```

| marker | comparison        |   diff | lwr.ci | upr.ci |  pval |
|:-------|:------------------|-------:|-------:|-------:|------:|
| GRP78  | Induced-Uninduced | -0.102 | -0.494 |   0.29 | 0.510 |
| PDI    | Induced-Uninduced | -0.136 | -0.713 |   0.44 | 0.548 |

Document the information about the analysis session

``` r
sessionInfo()
```

    ## R version 4.2.2 (2022-10-31)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] curl_4.3.3        DescTools_0.99.47 ggthemes_4.2.4    Hmisc_4.7-1      
    ##  [5] Formula_1.2-4     survival_3.4-0    lattice_0.20-45   fs_1.5.2         
    ##  [9] forcats_0.5.2     stringr_1.4.1     dplyr_1.0.10      purrr_0.3.5      
    ## [13] readr_2.1.3       tidyr_1.2.1       tibble_3.1.8      ggplot2_3.3.6    
    ## [17] tidyverse_1.3.2  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_4.0.5         lubridate_1.8.0     RColorBrewer_1.1-3 
    ##  [4] httr_1.4.4          tools_4.2.2         backports_1.4.1    
    ##  [7] utf8_1.2.2          R6_2.5.1            rpart_4.1.19       
    ## [10] DBI_1.1.3           colorspace_2.0-3    nnet_7.3-18        
    ## [13] withr_2.5.0         Exact_3.2           tidyselect_1.2.0   
    ## [16] gridExtra_2.3       bit_4.0.4           compiler_4.2.2     
    ## [19] cli_3.4.1           rvest_1.0.3         htmlTable_2.4.1    
    ## [22] expm_0.999-6        xml2_1.3.3          labeling_0.4.2     
    ## [25] scales_1.2.1        checkmate_2.1.0     mvtnorm_1.1-3      
    ## [28] proxy_0.4-27        digest_0.6.30       foreign_0.8-84     
    ## [31] rmarkdown_2.17      base64enc_0.1-3     jpeg_0.1-9         
    ## [34] pkgconfig_2.0.3     htmltools_0.5.3     highr_0.9          
    ## [37] dbplyr_2.2.1        fastmap_1.1.0       htmlwidgets_1.5.4  
    ## [40] rlang_1.0.6         readxl_1.4.1        rstudioapi_0.14    
    ## [43] farver_2.1.1        generics_0.1.3      jsonlite_1.8.3     
    ## [46] vroom_1.6.0         googlesheets4_1.0.1 magrittr_2.0.3     
    ## [49] interp_1.1-3        Matrix_1.5-3        Rcpp_1.0.9         
    ## [52] munsell_0.5.0       fansi_1.0.3         lifecycle_1.0.3    
    ## [55] stringi_1.7.8       yaml_2.3.6          rootSolve_1.8.2.3  
    ## [58] MASS_7.3-58.1       grid_4.2.2          parallel_4.2.2     
    ## [61] crayon_1.5.2        lmom_2.9            deldir_1.0-6       
    ## [64] haven_2.5.1         splines_4.2.2       hms_1.1.2          
    ## [67] knitr_1.40          pillar_1.8.1        boot_1.3-28.1      
    ## [70] gld_2.6.6           reprex_2.0.2        glue_1.6.2         
    ## [73] evaluate_0.17       latticeExtra_0.6-30 data.table_1.14.4  
    ## [76] modelr_0.1.9        png_0.1-7           vctrs_0.5.0        
    ## [79] tzdb_0.3.0          cellranger_1.1.0    gtable_0.3.1       
    ## [82] assertthat_0.2.1    xfun_0.34           broom_1.0.1        
    ## [85] e1071_1.7-12        class_7.3-20        googledrive_2.0.0  
    ## [88] gargle_1.2.1        cluster_2.1.4       ellipsis_0.3.2
