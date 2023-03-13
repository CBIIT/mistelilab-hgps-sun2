Figure 1C: TetON Cells / Hsc70
================
Sandra Vidak/Gianluca Pegoraro
November 16th 2022

### Introduction

Columbus screen names:

`20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp40_20191008_123603`

`20191009-40x-TetON-Hsp110-Calnexin-Hsc70-SUN1_20191009_125328`

`20191010-40x-TetON-Hsc70-Hsp40_20191010_142649`

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

    ## Rows: 6
    ## Columns: 4
    ## $ row       <dbl> 2, 3, 4, 8, 9, 10
    ## $ column    <dbl> 7, 7, 7, 7, 7, 7
    ## $ marker    <chr> "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hsc70"
    ## $ cell_line <fct> Uninduced, Uninduced, Uninduced, Induced, Induced, Induced

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

### Download the data if needed

Download and unzip the Columbus results of the experiments from Figshare
if they have not been already downloaded.

``` r
if(!dir.exists("input")) {
  URL <- "https://figshare.com/ndownloader/files/38158128"
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

    ## Rows: 218,828
    ## Columns: 11
    ## $ screen           <chr> "20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp40_201…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "B16", "B16", "B16", "B16", "B16", "B16", "B16", "B16…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 1…
    ## $ nuc_area         <dbl> 1518, 3332, 2068, 5254, 1348, 3234, 2767, 1619, 1701,…
    ## $ cyto_area        <dbl> 2114, 2523, 2033, 4064, 1970, 2387, 2685, 1039, 6559,…
    ## $ cell_area        <dbl> 3632, 5855, 4101, 9318, 3318, 5621, 5452, 2658, 8260,…
    ## $ nuc_marker_int   <dbl> 439.182, 525.465, 471.506, 488.599, 499.123, 517.741,…
    ## $ cyto_marker_int  <dbl> 497.314, 541.748, 536.463, 467.419, 497.477, 555.837,…
    ## $ ratio_marker_int <dbl> 0.883110, 0.969945, 0.878916, 1.045310, 1.003310, 0.9…

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

    ## Rows: 27,962
    ## Columns: 14
    ## $ screen           <chr> "20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp40_201…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "B7", "B7", "B7", "B7", "B7", "B7", "B7", "B7", "B7",…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,…
    ## $ cell_line        <fct> Uninduced, Uninduced, Uninduced, Uninduced, Uninduced…
    ## $ marker           <chr> "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hsc70",…
    ## $ nuc_area         <dbl> 1290, 1793, 1682, 6450, 2231, 2065, 1620, 1452, 1997,…
    ## $ cyto_area        <dbl> 325, 1529, 2216, 8940, 3403, 825, 2352, 178, 2941, 14…
    ## $ cell_area        <dbl> 1615, 3322, 3898, 15390, 5634, 2890, 3972, 1630, 4938…
    ## $ nuc_marker_int   <dbl> 562.119, 721.327, 681.144, 959.265, 1306.670, 607.091…
    ## $ cyto_marker_int  <dbl> 721.760, 919.993, 674.949, 1077.140, 1437.130, 773.68…
    ## $ ratio_marker_int <dbl> 0.778817, 0.784056, 1.009180, 0.890570, 0.909223, 0.7…
    ## $ sum_marker_int   <dbl> 1283.879, 1641.320, 1356.093, 2036.405, 2743.800, 138…

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

    ## Rows: 18
    ## Columns: 14
    ## Groups: screen, well, row, column, cell_line [18]
    ## $ screen                <chr> "20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp4…
    ## $ well                  <chr> "B7", "C7", "D7", "H7", "I7", "J7", "B7", "C7", …
    ## $ row                   <dbl> 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, 3, 4, 8…
    ## $ column                <dbl> 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, …
    ## $ cell_line             <fct> Uninduced, Uninduced, Uninduced, Induced, Induce…
    ## $ marker                <chr> "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hs…
    ## $ cell_n                <int> 1510, 1511, 1542, 1634, 1626, 1566, 1613, 1614, …
    ## $ nuc_area_mean         <dbl> 2709.498, 2581.174, 2586.439, 1988.895, 2007.636…
    ## $ cyto_area_mean        <dbl> 3135.740, 3246.265, 3071.267, 3373.010, 3382.977…
    ## $ cell_area_mean        <dbl> 5845.238, 5827.439, 5657.706, 5361.905, 5390.613…
    ## $ nuc_marker_int_mean   <dbl> 828.384, 1054.481, 899.332, 1188.169, 1211.778, …
    ## $ cyto_marker_int_mean  <dbl> 962.0593, 1299.0195, 1088.6013, 1703.8028, 1684.…
    ## $ ratio_marker_int_mean <dbl> 0.8655093, 0.8321555, 0.8422823, 0.7075251, 0.73…
    ## $ sum_marker_int_mean   <dbl> 1790.443, 2353.500, 1988.151, 2891.828, 2896.411…

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

    ## Rows: 6
    ## Columns: 11
    ## Groups: screen, cell_line [6]
    ## $ screen                <chr> "20191008-40x-TetON-Hsp110-Hsc70-SUN1-Hsp90-Hsp4…
    ## $ cell_line             <fct> Uninduced, Induced, Uninduced, Induced, Uninduce…
    ## $ marker                <chr> "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hsc70", "Hs…
    ## $ cell_n                <dbl> 1521.000, 1608.667, 1592.000, 1532.333, 1493.667…
    ## $ nuc_area_mean         <dbl> 2625.704, 2024.454, 2561.764, 2110.345, 2569.536…
    ## $ cyto_area_mean        <dbl> 3151.091, 3429.462, 2897.574, 3563.268, 3292.877…
    ## $ cell_area_mean        <dbl> 5776.794, 5453.917, 5459.338, 5673.613, 5862.414…
    ## $ nuc_marker_int_mean   <dbl> 927.3988, 1196.2333, 1256.2622, 2666.8230, 1316.…
    ## $ cyto_marker_int_mean  <dbl> 1116.560, 1683.634, 1431.065, 3308.610, 1461.492…
    ## $ ratio_marker_int_mean <dbl> 0.8466490, 0.7212851, 0.8860138, 0.8196992, 0.91…
    ## $ sum_marker_int_mean   <dbl> 2044.031, 2879.998, 2687.327, 5975.792, 2778.065…

### Biological Replicates Level plot for Fig.1C

![](output/Fig1_C-1.png)<!-- -->

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
| Hsc70  | Induced-Uninduced | -0.029 | -0.277 |  0.219 | 0.763 |

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
    ##  [1] curl_4.3.3        DescTools_0.99.47 ggthemes_4.2.4    Hmisc_4.7-1      
    ##  [5] Formula_1.2-4     survival_3.4-0    lattice_0.20-45   fs_1.5.2         
    ##  [9] forcats_0.5.2     stringr_1.4.1     dplyr_1.0.10      purrr_0.3.5      
    ## [13] readr_2.1.3       tidyr_1.2.1       tibble_3.1.8      ggplot2_3.3.6    
    ## [17] tidyverse_1.3.2  
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
