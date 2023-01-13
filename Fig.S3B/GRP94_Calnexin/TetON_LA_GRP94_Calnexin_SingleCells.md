Single Cell Analysis for TetON-LA Cells
================
Sandra Vidak/Gianluca Pegoraro
January 13th 2023

### Introduction

Columbus screen names:

`171006-TetON-LA-GRP94-Calnexin_20171006_171548`

`180130-40x-TetON-LA-GRP94,Calnexin-SUN1_20180130_114147`

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
    ## $ screen    <chr> "171006-TetON-LA-GRP94-Calnexin_20171006_171548", "171006-Te…
    ## $ row       <dbl> 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, …
    ## $ column    <dbl> 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2, 2, 4, 4, …
    ## $ marker    <chr> "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "Calne…
    ## $ cell_line <fct> Uninduced, Uninduced, Uninduced, Induced, Induced, Induced, …

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

Download and unzip the Columbus results of the experiments from Figshare
if they have not been already downloaded.

``` r
if(!dir.exists("input")) {
  URL <- "https://figshare.com/ndownloader/files/38669300"
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

    ## Rows: 153,645
    ## Columns: 11
    ## $ screen           <chr> "171006-TetON-LA-GRP94-Calnexin_20171006_171548", "17…
    ## $ plate            <chr> "Plate1", "Plate1", "Plate1", "Plate1", "Plate1", "Pl…
    ## $ well             <chr> "B13", "B13", "B13", "B13", "B13", "B13", "B13", "B13…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 1…
    ## $ nuc_area         <dbl> 1686, 1641, 1755, 1689, 3086, 1549, 1787, 2638, 2054,…
    ## $ cyto_area        <dbl> 1493, 6760, 5968, 3782, 7969, 11990, 3861, 8689, 1609…
    ## $ cell_area        <dbl> 3179, 8401, 7723, 5471, 11055, 13539, 5648, 11327, 18…
    ## $ nuc_marker_int   <dbl> 227.922, 338.604, 285.409, 304.580, 340.570, 262.719,…
    ## $ cyto_marker_int  <dbl> 278.468, 324.927, 325.882, 375.914, 343.985, 266.381,…
    ## $ ratio_marker_int <dbl> 0.818486, 1.042090, 0.875805, 0.810239, 0.990071, 0.9…

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

    ## Rows: 31,239
    ## Columns: 14
    ## $ screen           <chr> "171006-TetON-LA-GRP94-Calnexin_20171006_171548", "17…
    ## $ plate            <chr> "Plate1", "Plate1", "Plate1", "Plate1", "Plate1", "Pl…
    ## $ well             <chr> "B2", "B2", "B2", "B2", "B2", "B2", "B2", "B2", "B2",…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ cell_line        <fct> Uninduced, Uninduced, Uninduced, Uninduced, Uninduced…
    ## $ marker           <chr> "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GRP94",…
    ## $ nuc_area         <dbl> 1293, 1621, 1704, 1611, 1954, 1884, 1875, 1683, 1954,…
    ## $ cyto_area        <dbl> 1744, 1846, 5496, 9889, 2665, 8378, 5833, 7026, 6036,…
    ## $ cell_area        <dbl> 3037, 3467, 7200, 11500, 4619, 10262, 7708, 8709, 799…
    ## $ nuc_marker_int   <dbl> 318.493, 301.035, 315.989, 321.253, 239.463, 281.372,…
    ## $ cyto_marker_int  <dbl> 336.944, 353.451, 397.297, 383.383, 302.250, 338.292,…
    ## $ ratio_marker_int <dbl> 0.945240, 0.851702, 0.795348, 0.837942, 0.792269, 0.8…
    ## $ sum_marker_int   <dbl> 655.437, 654.486, 713.286, 704.636, 541.713, 619.664,…

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
    ## $ screen                <chr> "171006-TetON-LA-GRP94-Calnexin_20171006_171548"…
    ## $ well                  <chr> "B2", "B4", "C2", "C4", "D2", "D4", "H2", "H4", …
    ## $ row                   <dbl> 2, 2, 3, 3, 4, 4, 8, 8, 9, 9, 10, 10, 2, 2, 3, 3…
    ## $ column                <dbl> 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, …
    ## $ cell_line             <fct> Uninduced, Uninduced, Uninduced, Uninduced, Unin…
    ## $ marker                <chr> "GRP94", "Calnexin", "GRP94", "Calnexin", "GRP94…
    ## $ cell_n                <int> 1005, 1136, 1012, 953, 964, 1020, 1037, 1015, 94…
    ## $ nuc_area_mean         <dbl> 1968.333, 2026.074, 1966.806, 1931.409, 1891.002…
    ## $ cyto_area_mean        <dbl> 6060.169, 5344.236, 6216.121, 6517.486, 6416.634…
    ## $ cell_area_mean        <dbl> 8028.502, 7370.310, 8182.927, 8448.895, 8307.636…
    ## $ nuc_marker_int_mean   <dbl> 284.8481, 264.8836, 279.8431, 249.0214, 289.7483…
    ## $ cyto_marker_int_mean  <dbl> 308.2083, 248.4560, 306.6213, 225.9160, 317.3632…
    ## $ ratio_marker_int_mean <dbl> 0.9280041, 1.0723685, 0.9188292, 1.1077599, 0.91…
    ## $ sum_marker_int_mean   <dbl> 593.0564, 513.3395, 586.4644, 474.9374, 607.1115…

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
    ## Groups: screen, cell_line [8]
    ## $ screen                <chr> "171006-TetON-LA-GRP94-Calnexin_20171006_171548"…
    ## $ cell_line             <fct> Uninduced, Uninduced, Induced, Induced, Uninduce…
    ## $ marker                <chr> "Calnexin", "GRP94", "Calnexin", "GRP94", "Calne…
    ## $ cell_n                <dbl> 1036.3333, 993.6667, 1002.3333, 1029.6667, 984.0…
    ## $ nuc_area_mean         <dbl> 1959.848, 1942.047, 1681.717, 1675.370, 1899.869…
    ## $ cyto_area_mean        <dbl> 5970.531, 6230.975, 6424.571, 6308.467, 5810.713…
    ## $ cell_area_mean        <dbl> 7930.379, 8173.022, 8106.288, 7983.837, 7710.582…
    ## $ nuc_marker_int_mean   <dbl> 253.5057, 284.8132, 260.3150, 278.1049, 234.4912…
    ## $ cyto_marker_int_mean  <dbl> 235.0675, 310.7309, 244.8077, 320.3673, 192.3320…
    ## $ ratio_marker_int_mean <dbl> 1.0848114, 0.9220026, 1.0728799, 0.8761103, 1.29…
    ## $ sum_marker_int_mean   <dbl> 488.5733, 595.5441, 505.1227, 598.4721, 426.6857…

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

| marker   | comparison        |   diff | lwr.ci | upr.ci |  pval |
|:---------|:------------------|-------:|-------:|-------:|------:|
| Calnexin | Induced-Uninduced | -0.053 | -0.879 |  0.772 | 0.866 |
| GRP94    | Induced-Uninduced | -0.157 | -0.477 |  0.162 | 0.243 |

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
