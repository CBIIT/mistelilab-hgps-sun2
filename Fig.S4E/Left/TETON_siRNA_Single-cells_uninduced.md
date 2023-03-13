Figure S4E Left: TETON-Progerin Cells
================
Sandra Vidak/Gianluca Pegoraro
March 13th 2023

### Introduction

Columbus screen names:

`191126-60x-TetON-siRNA-GRP78-PDI-plate1_20191126_122500`

`191126-60x-TetON-siRNA-GRP78-PDI-plate2_20191126_133003`

`200219-40x-TetON-siRNA-SUN2-GRP78_20200219_121533`

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
cell_levs <- c("Uninduced_siCTRL", "Uninduced_siSUN1",
               "Uninduced_siSUN2")

plate_layouts <- read_tsv("metadata/plate_layout.txt") %>%
  filter(!is.na(cell_line)) %>%
  mutate(cell_line = factor(cell_line, levels = cell_levs))

glimpse(plate_layouts)
```

    ## Rows: 27
    ## Columns: 5
    ## $ screen    <chr> "191126-60x-TetON-siRNA-GRP78-PDI-plate1_20191126_122500", "…
    ## $ row       <dbl> 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, …
    ## $ column    <dbl> 3, 3, 3, 8, 8, 8, 13, 13, 13, 3, 3, 3, 8, 8, 8, 13, 13, 13, …
    ## $ marker    <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78…
    ## $ cell_line <fct> Uninduced_siCTRL, Uninduced_siCTRL, Uninduced_siCTRL, Unindu…

Plot plate layouts.

![](output/cell-line-layout-1.png)<!-- -->

![](output/ab-layout-1.png)<!-- -->

### Download the data if needed

Download and unzip the Columbus results of the experiments from Figshare
if they have not been already downloaded.

``` r
if(!dir.exists("input")) {
  URL <- "https://figshare.com/ndownloader/files/38843265"
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

    ## Rows: 32,481
    ## Columns: 11
    ## $ screen           <chr> "191126-60x-TetON-siRNA-GRP78-PDI-plate1_20191126_122…
    ## $ plate            <chr> "Plate 1", "Plate 1", "Plate 1", "Plate 1", "Plate 1"…
    ## $ well             <chr> "B13", "B13", "B13", "B13", "B13", "B13", "B13", "B13…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 1…
    ## $ nuc_area         <dbl> 5144, 7474, 4040, 8707, 3760, 4622, 2924, 3316, 5427,…
    ## $ cyto_area        <dbl> 8053, 29211, 14156, 24134, 20069, 16717, 17771, 12435…
    ## $ cell_area        <dbl> 13197, 36685, 18196, 32841, 23829, 21339, 20695, 1575…
    ## $ nuc_marker_int   <dbl> 5032.96, 4549.93, 4691.14, 4136.36, 3911.96, 4743.30,…
    ## $ cyto_marker_int  <dbl> 9761.35, 5545.47, 5392.92, 5582.04, 3250.44, 5221.08,…
    ## $ ratio_marker_int <dbl> 0.515601, 0.820475, 0.869870, 0.741013, 1.203520, 0.9…

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

    ## Rows: 8,564
    ## Columns: 14
    ## $ screen           <chr> "191126-60x-TetON-siRNA-GRP78-PDI-plate1_20191126_122…
    ## $ plate            <chr> "Plate 1", "Plate 1", "Plate 1", "Plate 1", "Plate 1"…
    ## $ well             <chr> "B13", "B13", "B13", "B13", "B13", "B13", "B13", "B13…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 1…
    ## $ cell_line        <fct> Uninduced_siSUN2, Uninduced_siSUN2, Uninduced_siSUN2,…
    ## $ marker           <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78",…
    ## $ nuc_area         <dbl> 5144, 7474, 4040, 8707, 3760, 4622, 2924, 3316, 5427,…
    ## $ cyto_area        <dbl> 8053, 29211, 14156, 24134, 20069, 16717, 17771, 12435…
    ## $ cell_area        <dbl> 13197, 36685, 18196, 32841, 23829, 21339, 20695, 1575…
    ## $ nuc_marker_int   <dbl> 5032.96, 4549.93, 4691.14, 4136.36, 3911.96, 4743.30,…
    ## $ cyto_marker_int  <dbl> 9761.35, 5545.47, 5392.92, 5582.04, 3250.44, 5221.08,…
    ## $ ratio_marker_int <dbl> 0.515601, 0.820475, 0.869870, 0.741013, 1.203520, 0.9…
    ## $ sum_marker_int   <dbl> 14794.31, 10095.40, 10084.06, 9718.40, 7162.40, 9964.…

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
    ## $ screen                <chr> "191126-60x-TetON-siRNA-GRP78-PDI-plate1_2019112…
    ## $ well                  <chr> "B13", "B3", "B8", "C13", "C3", "C8", "D13", "D3…
    ## $ row                   <dbl> 2, 2, 2, 3, 3, 3, 4, 4, 4, 2, 2, 2, 3, 3, 3, 4, …
    ## $ column                <dbl> 13, 3, 8, 13, 3, 8, 13, 3, 8, 13, 3, 8, 13, 3, 8…
    ## $ cell_line             <fct> Uninduced_siSUN2, Uninduced_siCTRL, Uninduced_si…
    ## $ marker                <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GR…
    ## $ cell_n                <int> 193, 201, 165, 184, 182, 188, 166, 242, 209, 185…
    ## $ nuc_area_mean         <dbl> 5576.896, 5266.856, 5387.903, 5969.837, 5514.044…
    ## $ cyto_area_mean        <dbl> 17912.083, 18401.632, 20523.036, 19665.000, 1933…
    ## $ cell_area_mean        <dbl> 23488.98, 23668.49, 25910.94, 25634.84, 24852.73…
    ## $ nuc_marker_int_mean   <dbl> 3450.8893, 3400.9418, 3729.5036, 3507.5346, 3387…
    ## $ cyto_marker_int_mean  <dbl> 4780.1358, 4934.2266, 5117.5521, 4930.1832, 5185…
    ## $ ratio_marker_int_mean <dbl> 0.7591822, 0.7264506, 0.7521074, 0.7321775, 0.67…
    ## $ sum_marker_int_mean   <dbl> 8258.115, 8335.168, 8847.056, 8437.718, 8572.770…

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
    ## $ screen                <chr> "191126-60x-TetON-siRNA-GRP78-PDI-plate1_2019112…
    ## $ cell_line             <fct> Uninduced_siCTRL, Uninduced_siSUN1, Uninduced_si…
    ## $ marker                <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GR…
    ## $ cell_n                <dbl> 208.3333, 187.3333, 181.0000, 441.6667, 204.6667…
    ## $ nuc_area_mean         <dbl> 5389.389, 5410.146, 5844.479, 4667.544, 4948.468…
    ## $ cyto_area_mean        <dbl> 18175.365, 19535.747, 19315.809, 10414.747, 1650…
    ## $ cell_area_mean        <dbl> 23564.75, 24945.89, 25160.29, 15082.29, 21451.16…
    ## $ nuc_marker_int_mean   <dbl> 3464.5376, 3623.5486, 3447.3801, 2941.1068, 3592…
    ## $ cyto_marker_int_mean  <dbl> 5141.795, 5120.515, 4821.790, 3989.256, 4355.345…
    ## $ ratio_marker_int_mean <dbl> 0.7014815, 0.7299689, 0.7457416, 0.7595083, 0.86…
    ## $ sum_marker_int_mean   <dbl> 8606.332, 8747.527, 8283.456, 6929.136, 7948.312…

### Biological Replicates Level plot for the left plot in Fig.S4E

![](output/FigS4_E_left-1.png)<!-- -->

### Calculate Dunnett’s test for the continuous variables.

Define a custom function to run a Dunnett post-hoc test only on the Mean
marker intensity sum (Cyto + Nucleus), using the cell line as the
predictor variable, and fixing WT1 as the negative control. The output
of the Dunnett’s test is then rearranged to a tidy table to make it work
with `dplyr`.

``` r
calc_dunnett <- function(df){
  as.data.frame(as.table(DunnettTest(ratio_marker_int_mean ~ cell_line,
                          control = "Uninduced_siCTRL",
                          data = df)$Uninduced_siCTRL)) %>%
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

| marker | comparison                        |  diff | lwr.ci | upr.ci |  pval |
|:-------|:----------------------------------|------:|-------:|-------:|------:|
| GRP78  | Uninduced_siSUN1-Uninduced_siCTRL | 0.010 | -0.357 |  0.376 | 0.996 |
| GRP78  | Uninduced_siSUN2-Uninduced_siCTRL | 0.016 | -0.351 |  0.382 | 0.989 |

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
