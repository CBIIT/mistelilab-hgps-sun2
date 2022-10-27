Figure S1D: TetON Cells / Hsp110
================
Sandra Vidak/Gianluca Pegoraro
October 27th 2022

### Introduction

Columbus screen names:

`200707-40x-hTERTfibro-Hsp90-Hsp110-SUN1_20200707_115858`

`200712-40x-hTERTfibro-Hsp110_20200714_115633`

`200724-40x-hTERTfibro-Calnexin-Hsp90-Hsp110-SUN`\_20200724_122959\`

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
source("R/Plotters.R") #Functions needed for plotting
```

Set the palette and the running theme for ggplot2.

### Experimental Metadata

Read plate layouts.

``` r
cell_levs <- c("WT", "HGPS1", "HGPS2")

plate_layouts <- read_tsv("metadata/plate_layout.txt") %>%
  filter(!is.na(cell_line)) %>%
  separate(col = cell_line, 
           into = c("cell_line"), 
           remove = T) %>%
  mutate(cell_line = factor(cell_line, levels = cell_levs))

glimpse(plate_layouts)
```

    ## Rows: 27
    ## Columns: 5
    ## $ screen    <chr> "200707-40x-hTERTfibro-Hsp90-Hsp110-SUN1_20200707_115858", "…
    ## $ row       <dbl> 11, 12, 13, 11, 12, 13, 11, 12, 13, 11, 12, 13, 11, 12, 13, …
    ## $ column    <dbl> 3, 3, 3, 10, 10, 10, 17, 17, 17, 3, 3, 3, 9, 9, 9, 15, 15, 1…
    ## $ marker    <chr> "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110", …
    ## $ cell_line <fct> WT, WT, WT, HGPS1, HGPS1, HGPS1, HGPS2, HGPS2, HGPS2, WT, WT…

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

### Read and Process Columbus data

Recursively search the `input` directory and its subdirectories for
files whose name includes the Glob patterns defined in the chunk above,
and read the cell-level Columbus data from the results text files.

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

    ## Rows: 30,633
    ## Columns: 11
    ## $ screen           <chr> "200712-40x-hTERTfibro-Hsp110_20200714_115633", "2007…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "K15", "K15", "K15", "K15", "K15", "K15", "K15", "K15…
    ## $ row              <dbl> 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 1…
    ## $ column           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1…
    ## $ nuc_area         <dbl> 1537, 2401, 1023, 1666, 1624, 1660, 1294, 1305, 1901,…
    ## $ cyto_area        <dbl> 2309, 9142, 5552, 19244, 9325, 12161, 8248, 14982, 10…
    ## $ cell_area        <dbl> 3846, 11543, 6575, 20910, 10949, 13821, 9542, 16287, …
    ## $ nuc_marker_int   <dbl> 490.918, 503.028, 402.938, 495.666, 637.902, 556.551,…
    ## $ cyto_marker_int  <dbl> 850.653, 590.609, 538.651, 616.081, 626.944, 616.197,…
    ## $ ratio_marker_int <dbl> 0.577107, 0.851712, 0.748050, 0.804546, 1.017480, 0.9…

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

    ## Rows: 11,326
    ## Columns: 14
    ## $ screen           <chr> "200712-40x-hTERTfibro-Hsp110_20200714_115633", "2007…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "K15", "K15", "K15", "K15", "K15", "K15", "K15", "K15…
    ## $ row              <dbl> 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 1…
    ## $ column           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1…
    ## $ cell_line        <fct> HGPS2, HGPS2, HGPS2, HGPS2, HGPS2, HGPS2, HGPS2, HGPS…
    ## $ marker           <chr> "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hs…
    ## $ nuc_area         <dbl> 1537, 2401, 1023, 1666, 1624, 1660, 1294, 1305, 1901,…
    ## $ cyto_area        <dbl> 2309, 9142, 5552, 19244, 9325, 12161, 8248, 14982, 10…
    ## $ cell_area        <dbl> 3846, 11543, 6575, 20910, 10949, 13821, 9542, 16287, …
    ## $ nuc_marker_int   <dbl> 490.918, 503.028, 402.938, 495.666, 637.902, 556.551,…
    ## $ cyto_marker_int  <dbl> 850.653, 590.609, 538.651, 616.081, 626.944, 616.197,…
    ## $ ratio_marker_int <dbl> 0.577107, 0.851712, 0.748050, 0.804546, 1.017480, 0.9…
    ## $ sum_marker_int   <dbl> 1341.571, 1093.637, 941.589, 1111.747, 1264.846, 1172…

Calculate number of cells and mean per well for all properties.

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

    ## Rows: 27
    ## Columns: 14
    ## Groups: screen, well, row, column, cell_line [27]
    ## $ screen                <chr> "200707-40x-hTERTfibro-Hsp90-Hsp110-SUN1_2020070…
    ## $ well                  <chr> "K10", "K17", "K3", "L10", "L17", "L3", "M10", "…
    ## $ row                   <dbl> 11, 11, 11, 12, 12, 12, 13, 13, 13, 11, 11, 11, …
    ## $ column                <dbl> 10, 17, 3, 10, 17, 3, 10, 17, 3, 15, 3, 9, 15, 3…
    ## $ cell_line             <fct> HGPS1, HGPS2, WT, HGPS1, HGPS2, WT, HGPS1, HGPS2…
    ## $ marker                <chr> "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110"…
    ## $ cell_n                <int> 279, 438, 323, 268, 513, 322, 277, 391, 311, 419…
    ## $ nuc_area_mean         <dbl> 2058.018, 1941.605, 2160.526, 2062.567, 1952.140…
    ## $ cyto_area_mean        <dbl> 13971.760, 10441.635, 13259.929, 14659.575, 9818…
    ## $ cell_area_mean        <dbl> 16029.778, 12383.240, 15420.455, 16722.142, 1177…
    ## $ nuc_marker_int_mean   <dbl> 375.25204, 510.99237, 290.33364, 535.33806, 893.…
    ## $ cyto_marker_int_mean  <dbl> 314.27688, 454.62393, 268.47418, 448.07094, 841.…
    ## $ ratio_marker_int_mean <dbl> 1.2341478, 1.1522838, 1.1117972, 1.2207303, 1.09…
    ## $ sum_marker_int_mean   <dbl> 689.5289, 965.6163, 558.8078, 983.4090, 1734.793…

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

    ## Rows: 9
    ## Columns: 11
    ## Groups: screen, cell_line [9]
    ## $ screen                <chr> "200707-40x-hTERTfibro-Hsp90-Hsp110-SUN1_2020070…
    ## $ cell_line             <fct> WT, HGPS1, HGPS2, WT, HGPS1, HGPS2, WT, HGPS1, H…
    ## $ marker                <chr> "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110"…
    ## $ cell_n                <dbl> 318.6667, 274.6667, 447.3333, 383.0000, 491.3333…
    ## $ nuc_area_mean         <dbl> 2143.960, 2050.804, 1938.905, 1980.900, 1813.995…
    ## $ cyto_area_mean        <dbl> 13561.596, 14119.334, 10773.492, 11640.079, 9913…
    ## $ cell_area_mean        <dbl> 15705.56, 16170.14, 12712.40, 13620.98, 11727.57…
    ## $ nuc_marker_int_mean   <dbl> 379.5258, 549.1569, 705.3901, 866.8587, 773.2634…
    ## $ cyto_marker_int_mean  <dbl> 340.5177, 435.7027, 616.2005, 1114.2206, 942.164…
    ## $ ratio_marker_int_mean <dbl> 1.1412693, 1.2839610, 1.1873121, 0.7945742, 0.83…
    ## $ sum_marker_int_mean   <dbl> 720.0434, 984.8596, 1321.5906, 1981.8467, 1715.4…

### Biological Replicates Level plots For Figure S1D

![](output/Fig.S1_D-1.png)<!-- -->

### Calculate Dunnett’s test for the continuous variables.

Define a custom function to run a Dunnett post-hoc test only on the Mean
marker intensity sum (Cyto + Nucleus), using the cell line as the
predictor variable, and fixing WT1 as the negative control. The output
of the Dunnett’s test is then rearranged to a tidy table to make it work
with `dplyr`.

``` r
calc_dunnett <- function(df){
  as.data.frame(as.table(DunnettTest(ratio_marker_int_mean ~ cell_line,
                          control = "WT",
                          data = df)$WT)) %>%
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

| marker | comparison |   diff | lwr.ci | upr.ci |  pval |
|:-------|:-----------|-------:|-------:|-------:|------:|
| Hsp110 | HGPS1-WT   |  0.023 | -0.445 |  0.491 | 0.986 |
| Hsp110 | HGPS2-WT   | -0.011 | -0.479 |  0.457 | 0.997 |

Document the information about the analysis session

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
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
    ##  [1] DescTools_0.99.47 ggthemes_4.2.4    Hmisc_4.7-1       Formula_1.2-4    
    ##  [5] survival_3.4-0    lattice_0.20-45   fs_1.5.2          forcats_0.5.2    
    ##  [9] stringr_1.4.1     dplyr_1.0.10      purrr_0.3.5       readr_2.1.3      
    ## [13] tidyr_1.2.1       tibble_3.1.8      ggplot2_3.3.6     tidyverse_1.3.2  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bit64_4.0.5         lubridate_1.8.0     RColorBrewer_1.1-3 
    ##  [4] httr_1.4.4          tools_4.2.1         backports_1.4.1    
    ##  [7] utf8_1.2.2          R6_2.5.1            rpart_4.1.19       
    ## [10] DBI_1.1.3           colorspace_2.0-3    nnet_7.3-18        
    ## [13] withr_2.5.0         Exact_3.2           tidyselect_1.2.0   
    ## [16] gridExtra_2.3       bit_4.0.4           compiler_4.2.1     
    ## [19] cli_3.4.1           rvest_1.0.3         htmlTable_2.4.1    
    ## [22] expm_0.999-6        xml2_1.3.3          labeling_0.4.2     
    ## [25] scales_1.2.1        checkmate_2.1.0     mvtnorm_1.1-3      
    ## [28] proxy_0.4-27        digest_0.6.30       foreign_0.8-83     
    ## [31] rmarkdown_2.17      base64enc_0.1-3     jpeg_0.1-9         
    ## [34] pkgconfig_2.0.3     htmltools_0.5.3     highr_0.9          
    ## [37] dbplyr_2.2.1        fastmap_1.1.0       htmlwidgets_1.5.4  
    ## [40] rlang_1.0.6         readxl_1.4.1        rstudioapi_0.14    
    ## [43] farver_2.1.1        generics_0.1.3      jsonlite_1.8.3     
    ## [46] vroom_1.6.0         googlesheets4_1.0.1 magrittr_2.0.3     
    ## [49] interp_1.1-3        Matrix_1.5-1        Rcpp_1.0.9         
    ## [52] munsell_0.5.0       fansi_1.0.3         lifecycle_1.0.3    
    ## [55] stringi_1.7.8       yaml_2.3.6          rootSolve_1.8.2.3  
    ## [58] MASS_7.3-58.1       grid_4.2.1          parallel_4.2.1     
    ## [61] crayon_1.5.2        lmom_2.9            deldir_1.0-6       
    ## [64] haven_2.5.1         splines_4.2.1       hms_1.1.2          
    ## [67] knitr_1.40          pillar_1.8.1        boot_1.3-28        
    ## [70] gld_2.6.6           reprex_2.0.2        glue_1.6.2         
    ## [73] evaluate_0.17       latticeExtra_0.6-30 data.table_1.14.4  
    ## [76] modelr_0.1.9        png_0.1-7           vctrs_0.5.0        
    ## [79] tzdb_0.3.0          cellranger_1.1.0    gtable_0.3.1       
    ## [82] assertthat_0.2.1    xfun_0.34           broom_1.0.1        
    ## [85] e1071_1.7-12        class_7.3-20        googledrive_2.0.0  
    ## [88] gargle_1.2.1        cluster_2.1.4       ellipsis_0.3.2
