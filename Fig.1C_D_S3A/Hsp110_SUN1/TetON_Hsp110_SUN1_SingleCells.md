Figure 1C and S3A: TetON Cells / Hsp110 and SUN1
================
Sandra Vidak/Gianluca Pegoraro
March 13th 2022

### Introduction

Columbus screen names:

`20191008-40x-TetON-GRP94-Hsp70-Hsp40-Hsp90_20191008_131309`

`20191009-40x-TetON-GRP94-Hsp70-Hsp90-Hsp40_20191009_131506`

`20191010-40x-TetON-Hsp110-GRP94-Calnexin-Hsp70-SUN1-Hsp40-Hsp90_20191010_144805`

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

    ## Using libcurl 7.86.0 with LibreSSL/3.3.6
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

    ## Rows: 12
    ## Columns: 4
    ## $ row       <dbl> 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10
    ## $ column    <dbl> 2, 2, 2, 2, 2, 2, 9, 9, 9, 9, 9, 9
    ## $ marker    <chr> "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110", …
    ## $ cell_line <fct> Uninduced, Uninduced, Uninduced, Induced, Induced, Induced, …

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

### Download the data if needed

Download and unzip the Columbus results of the experiments from Figshare
if they have not been already downloaded.

``` r
if(!dir.exists("input")) {
  URL <- "https://figshare.com/ndownloader/files/38158173"
  curl_download(URL, "input.zip")
  unzip("input.zip")
}
```

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

    ## Rows: 327,191
    ## Columns: 11
    ## $ screen           <chr> "20191009-40x-TetON-Hsp110-Calnexin-Hsc70-SUN1_201910…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "B2", "B2", "B2", "B2", "B2", "B2", "B2", "B2", "B2",…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ nuc_area         <dbl> 1890, 1491, 990, 1136, 3845, 1289, 4611, 1848, 2540, …
    ## $ cyto_area        <dbl> 1075, 2031, 704, 1096, 6995, 1539, 1960, 1438, 2295, …
    ## $ cell_area        <dbl> 2965, 3522, 1694, 2232, 10840, 2828, 6571, 3286, 4835…
    ## $ nuc_marker_int   <dbl> 693.037, 677.667, 572.840, 712.676, 693.425, 629.361,…
    ## $ cyto_marker_int  <dbl> 943.717, 689.404, 782.678, 886.878, 791.174, 934.969,…
    ## $ ratio_marker_int <dbl> 0.734369, 0.982975, 0.731898, 0.803579, 0.876451, 0.6…

Join Columbus data with the plate layout information.

``` r
cell_tbl <- col_tbl %>%
  mutate(sum_marker_int = nuc_marker_int + cyto_marker_int) %>%
  inner_join(plate_layouts,
             by = c("row", "column")) %>%
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

    ## Rows: 55,174
    ## Columns: 14
    ## $ screen           <chr> "20191009-40x-TetON-Hsp110-Calnexin-Hsc70-SUN1_201910…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "B2", "B2", "B2", "B2", "B2", "B2", "B2", "B2", "B2",…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ cell_line        <fct> Uninduced, Uninduced, Uninduced, Uninduced, Uninduced…
    ## $ marker           <chr> "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hsp110", "Hs…
    ## $ nuc_area         <dbl> 1890, 1491, 990, 1136, 3845, 1289, 4611, 1848, 2540, …
    ## $ cyto_area        <dbl> 1075, 2031, 704, 1096, 6995, 1539, 1960, 1438, 2295, …
    ## $ cell_area        <dbl> 2965, 3522, 1694, 2232, 10840, 2828, 6571, 3286, 4835…
    ## $ nuc_marker_int   <dbl> 693.037, 677.667, 572.840, 712.676, 693.425, 629.361,…
    ## $ cyto_marker_int  <dbl> 943.717, 689.404, 782.678, 886.878, 791.174, 934.969,…
    ## $ ratio_marker_int <dbl> 0.734369, 0.982975, 0.731898, 0.803579, 0.876451, 0.6…
    ## $ sum_marker_int   <dbl> 1636.754, 1367.071, 1355.518, 1599.554, 1484.599, 156…

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

    ## Rows: 36
    ## Columns: 14
    ## Groups: screen, well, row, column, cell_line [36]
    ## $ screen                <chr> "20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp4…
    ## $ well                  <chr> "B2", "B9", "C2", "C9", "D2", "D9", "H2", "H9", …
    ## $ row                   <dbl> 2, 2, 3, 3, 4, 4, 8, 8, 9, 9, 10, 10, 2, 2, 3, 3…
    ## $ column                <dbl> 2, 9, 2, 9, 2, 9, 2, 9, 2, 9, 2, 9, 2, 9, 2, 9, …
    ## $ cell_line             <fct> Uninduced, Uninduced, Uninduced, Uninduced, Unin…
    ## $ marker                <chr> "Hsp110", "SUN1", "Hsp110", "SUN1", "Hsp110", "S…
    ## $ cell_n                <int> 1536, 1555, 1566, 1412, 1610, 1450, 1543, 1608, …
    ## $ nuc_area_mean         <dbl> 2659.503, 2635.758, 2641.430, 2604.302, 2436.370…
    ## $ cyto_area_mean        <dbl> 3066.577, 3044.633, 3041.213, 3471.748, 3094.278…
    ## $ cell_area_mean        <dbl> 5726.080, 5680.390, 5682.643, 6076.050, 5530.648…
    ## $ nuc_marker_int_mean   <dbl> 839.8879, 825.5209, 827.7724, 1168.8189, 871.535…
    ## $ cyto_marker_int_mean  <dbl> 1026.4642, 345.2158, 1025.6999, 368.7408, 1068.1…
    ## $ ratio_marker_int_mean <dbl> 0.8226794, 2.4776955, 0.8112201, 3.4418960, 0.82…
    ## $ sum_marker_int_mean   <dbl> 1866.3521, 1170.6654, 1853.4723, 1538.2249, 1939…

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
    ## $ screen                <chr> "20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp4…
    ## $ cell_line             <fct> Uninduced, Uninduced, Induced, Induced, Uninduce…
    ## $ marker                <chr> "Hsp110", "SUN1", "Hsp110", "SUN1", "Hsp110", "S…
    ## $ cell_n                <dbl> 1570.667, 1472.333, 1547.000, 1621.000, 1561.000…
    ## $ nuc_area_mean         <dbl> 2579.101, 2642.569, 2027.667, 2041.796, 2566.296…
    ## $ cyto_area_mean        <dbl> 3067.356, 3278.907, 3599.789, 3346.337, 3022.299…
    ## $ cell_area_mean        <dbl> 5646.457, 5921.477, 5627.456, 5388.133, 5588.595…
    ## $ nuc_marker_int_mean   <dbl> 846.3985, 1114.5744, 819.0576, 1934.2940, 672.04…
    ## $ cyto_marker_int_mean  <dbl> 1040.0937, 370.6190, 854.3445, 552.5658, 842.327…
    ## $ ratio_marker_int_mean <dbl> 0.8182083, 3.1851500, 0.9839436, 4.0455666, 0.80…
    ## $ sum_marker_int_mean   <dbl> 1886.4922, 1485.3912, 1673.4021, 2486.8598, 1514…

### Biological Replicates Level plots for Fig.1C and S3A

![](output/Fig1_C_S3A-1.png)<!-- -->

### Calculate Dunnett’s test for the continuous variables.

Define a custom function to run a Dunnett post-hoc test only on the Mean
marker intensity sum (Cyto + Nucleus), using the cell line as the
predictor variable, and fixing Uninduced as the negative control. The
output of the Dunnett’s test is then rearranged to a tidy table to make
it work with `dplyr`.

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

| marker | comparison        |  diff | lwr.ci | upr.ci |  pval |
|:-------|:------------------|------:|-------:|-------:|------:|
| Hsp110 | Induced-Uninduced | 0.115 |  0.033 |  0.198 | 0.018 |
| SUN1   | Induced-Uninduced | 0.662 | -0.233 |  1.556 | 0.109 |

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
    ##  [5] Formula_1.2-4     survival_3.5-3    lattice_0.20-45   fs_1.5.2         
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
    ## [58] MASS_7.3-58.3       grid_4.2.2          parallel_4.2.2     
    ## [61] crayon_1.5.2        lmom_2.9            deldir_1.0-6       
    ## [64] haven_2.5.1         splines_4.2.2       hms_1.1.2          
    ## [67] knitr_1.40          pillar_1.8.1        boot_1.3-28.1      
    ## [70] gld_2.6.6           reprex_2.0.2        glue_1.6.2         
    ## [73] evaluate_0.17       latticeExtra_0.6-30 data.table_1.14.4  
    ## [76] modelr_0.1.9        png_0.1-7           vctrs_0.5.0        
    ## [79] tzdb_0.3.0          cellranger_1.1.0    gtable_0.3.1       
    ## [82] assertthat_0.2.1    xfun_0.34           broom_1.0.1        
    ## [85] e1071_1.7-12        class_7.3-21        googledrive_2.0.0  
    ## [88] gargle_1.2.1        cluster_2.1.4       ellipsis_0.3.2
