Figure S1B: TetON Cells / GRP94
================
Sandra Vidak/Gianluca Pegoraro
March 13th 2023

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
cell_levs <- c("Uninduced", "2d", "4d", "6d")

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
    ## $ row       <dbl> 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
    ## $ column    <dbl> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3
    ## $ marker    <chr> "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GRP94…
    ## $ cell_line <fct> Uninduced, Uninduced, Uninduced, 2d, 2d, 2d, 4d, 4d, 4d, 6d,…

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

### Download the data if needed

Download and unzip the Columbus results of the experiments from Figshare
if they have not been already downloaded.

``` r
if(!dir.exists("input")) {
  URL <- "https://figshare.com/ndownloader/files/39550555"
  curl_download(URL, "input.zip")
  unzip("input.zip")
  dir_copy("data", "input")
  dir_delete("data")
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

    ## Rows: 309,346
    ## Columns: 11
    ## $ screen           <chr> "20191008-40x-TetON-GRP94-Hsp70-Hsp40-Hsp90_20191008_…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "B17", "B17", "B17", "B17", "B17", "B17", "B17", "B17…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 1…
    ## $ nuc_area         <dbl> 1658, 1635, 1626, 1158, 1644, 1133, 1401, 1715, 2145,…
    ## $ cyto_area        <dbl> 1121, 1086, 2939, 2697, 3181, 2205, 1100, 2856, 3065,…
    ## $ cell_area        <dbl> 2779, 2721, 4565, 3855, 4825, 3338, 2501, 4571, 5210,…
    ## $ nuc_marker_int   <dbl> 430.258, 409.428, 412.504, 433.524, 597.744, 394.823,…
    ## $ cyto_marker_int  <dbl> 359.518, 361.772, 368.287, 399.813, 443.879, 415.877,…
    ## $ ratio_marker_int <dbl> 1.196760, 1.131730, 1.120060, 1.084320, 1.346640, 0.9…

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

    ## Rows: 55,177
    ## Columns: 14
    ## $ screen           <chr> "20191008-40x-TetON-GRP94-Hsp70-Hsp40-Hsp90_20191008_…
    ## $ plate            <chr> "AssayPlate_PerkinElmer_CellCarrier-384 Ultra", "Assa…
    ## $ well             <chr> "B3", "B3", "B3", "B3", "B3", "B3", "B3", "B3", "B3",…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,…
    ## $ cell_line        <fct> Uninduced, Uninduced, Uninduced, Uninduced, Uninduced…
    ## $ marker           <chr> "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GRP94",…
    ## $ nuc_area         <dbl> 4267, 1509, 1929, 2340, 4192, 1361, 1664, 1436, 1631,…
    ## $ cyto_area        <dbl> 3795, 1463, 923, 2053, 2139, 1660, 1072, 1524, 1790, …
    ## $ cell_area        <dbl> 8062, 2972, 2852, 4393, 6331, 3021, 2736, 2960, 3421,…
    ## $ nuc_marker_int   <dbl> 203.710, 168.824, 253.594, 187.466, 208.701, 166.359,…
    ## $ cyto_marker_int  <dbl> 215.170, 184.591, 200.460, 199.874, 215.050, 213.169,…
    ## $ ratio_marker_int <dbl> 0.946743, 0.914585, 1.265060, 0.937921, 0.970475, 0.7…
    ## $ sum_marker_int   <dbl> 418.880, 353.415, 454.054, 387.340, 423.751, 379.528,…

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
    ## $ screen                <chr> "20191008-40x-TetON-GRP94-Hsp70-Hsp40-Hsp90_2019…
    ## $ well                  <chr> "B3", "C3", "D3", "E3", "F3", "G3", "H3", "I3", …
    ## $ row                   <dbl> 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 3, 4,…
    ## $ column                <dbl> 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, …
    ## $ cell_line             <fct> Uninduced, Uninduced, Uninduced, 2d, 2d, 2d, 4d,…
    ## $ marker                <chr> "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GR…
    ## $ cell_n                <int> 1582, 1525, 1536, 1580, 1714, 1527, 1571, 1496, …
    ## $ nuc_area_mean         <dbl> 2451.439, 2558.076, 2582.745, 2472.744, 2311.251…
    ## $ cyto_area_mean        <dbl> 3081.379, 3213.163, 3133.562, 3059.963, 2891.039…
    ## $ cell_area_mean        <dbl> 5532.819, 5771.239, 5716.307, 5532.706, 5202.289…
    ## $ nuc_marker_int_mean   <dbl> 226.7791, 240.8891, 273.3726, 326.8548, 339.0343…
    ## $ cyto_marker_int_mean  <dbl> 233.3086, 240.8811, 269.6929, 285.0406, 284.1472…
    ## $ ratio_marker_int_mean <dbl> 0.9822478, 1.0102559, 1.0227741, 1.1544212, 1.20…
    ## $ sum_marker_int_mean   <dbl> 460.0878, 481.7702, 543.0655, 611.8955, 623.1815…

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
    ## Groups: screen, cell_line [12]
    ## $ screen                <chr> "20191008-40x-TetON-GRP94-Hsp70-Hsp40-Hsp90_2019…
    ## $ cell_line             <fct> Uninduced, 2d, 4d, 6d, Uninduced, 2d, 4d, 6d, Un…
    ## $ marker                <chr> "GRP94", "GRP94", "GRP94", "GRP94", "GRP94", "GR…
    ## $ cell_n                <dbl> 1547.667, 1607.000, 1553.333, 1438.000, 1624.000…
    ## $ nuc_area_mean         <dbl> 2530.753, 2432.514, 2047.399, 2187.023, 2526.407…
    ## $ cyto_area_mean        <dbl> 3142.701, 3045.919, 3554.789, 3759.386, 2820.448…
    ## $ cell_area_mean        <dbl> 5673.455, 5478.433, 5602.188, 5946.410, 5346.855…
    ## $ nuc_marker_int_mean   <dbl> 247.0136, 332.8849, 402.8297, 500.9772, 212.6786…
    ## $ cyto_marker_int_mean  <dbl> 247.9609, 287.2453, 316.1956, 345.9218, 238.1536…
    ## $ ratio_marker_int_mean <dbl> 1.0050926, 1.1673414, 1.3014019, 1.4768760, 0.90…
    ## $ sum_marker_int_mean   <dbl> 494.9745, 620.1302, 719.0252, 846.8990, 450.8322…

### Biological Replicates Level plots for Fig.S1B

![](output/FigS1_B-1.png)<!-- -->

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

| marker | comparison   |  diff | lwr.ci | upr.ci |  pval |
|:-------|:-------------|------:|-------:|-------:|------:|
| GRP94  | 2d-Uninduced | 0.121 | -0.198 |  0.441 | 0.585 |
| GRP94  | 4d-Uninduced | 0.193 | -0.126 |  0.512 | 0.263 |
| GRP94  | 6d-Uninduced | 0.293 | -0.027 |  0.612 | 0.071 |

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
