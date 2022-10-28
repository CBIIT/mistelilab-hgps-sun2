Figure S1E and S1F: TetON-Lamin A / GRP78 and PDI
================
Sandra Vidak/Gianluca Pegoraro
October 28th 2022

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

    ## Rows: 24
    ## Columns: 4
    ## $ row       <dbl> 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, 3, 4, 8, 9, 10, 2, …
    ## $ column    <dbl> 3, 3, 3, 3, 3, 3, 5, 5, 5, 5, 5, 5, 14, 14, 14, 14, 14, 14, …
    ## $ marker    <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "PDI",…
    ## $ cell_line <fct> Uninduced, Uninduced, Uninduced, Induced, Induced, Induced, …

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

    ## Rows: 91,346
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

    ## Rows: 45,368
    ## Columns: 14
    ## $ screen           <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "171006-…
    ## $ plate            <chr> "Plate1", "Plate1", "Plate1", "Plate1", "Plate1", "Pl…
    ## $ well             <chr> "B14", "B14", "B14", "B14", "B14", "B14", "B14", "B14…
    ## $ row              <dbl> 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,…
    ## $ column           <dbl> 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 1…
    ## $ cell_line        <fct> Uninduced, Uninduced, Uninduced, Uninduced, Uninduced…
    ## $ marker           <chr> "GRP78", "GRP78", "GRP78", "GRP78", "GRP78", "GRP78",…
    ## $ nuc_area         <dbl> 2459, 1802, 1688, 1956, 1720, 1766, 1712, 2058, 2168,…
    ## $ cyto_area        <dbl> 5569, 1822, 1769, 13614, 9869, 11284, 4401, 6526, 559…
    ## $ cell_area        <dbl> 8028, 3624, 3457, 15570, 11589, 13050, 6113, 8584, 77…
    ## $ nuc_marker_int   <dbl> 2591.38, 1223.84, 1313.43, 1407.22, 1776.15, 2386.63,…
    ## $ cyto_marker_int  <dbl> 3412.72, 2498.72, 3226.09, 2910.42, 3375.52, 3737.21,…
    ## $ ratio_marker_int <dbl> 0.759330, 0.489789, 0.407128, 0.483511, 0.526184, 0.6…
    ## $ sum_marker_int   <dbl> 6004.10, 3722.56, 4539.52, 4317.64, 5151.67, 6123.84,…

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

    ## Rows: 48
    ## Columns: 14
    ## Groups: screen, well, row, column, cell_line [48]
    ## $ screen                <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "17…
    ## $ well                  <chr> "B14", "B16", "B3", "B5", "C14", "C16", "C3", "C…
    ## $ row                   <dbl> 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 8, 8, 8, 8, …
    ## $ column                <dbl> 14, 16, 3, 5, 14, 16, 3, 5, 14, 16, 3, 5, 14, 16…
    ## $ cell_line             <fct> Uninduced, Uninduced, Uninduced, Uninduced, Unin…
    ## $ marker                <chr> "GRP78", "PDI", "GRP78", "PDI", "GRP78", "PDI", …
    ## $ cell_n                <int> 1012, 918, 912, 937, 1039, 941, 966, 950, 1006, …
    ## $ nuc_area_mean         <dbl> 2072.559, 2023.108, 1961.951, 1966.489, 1951.420…
    ## $ cyto_area_mean        <dbl> 6014.120, 6766.755, 6709.349, 6424.189, 5889.753…
    ## $ cell_area_mean        <dbl> 8086.679, 8789.863, 8671.299, 8390.678, 7841.172…
    ## $ nuc_marker_int_mean   <dbl> 1832.2344, 1948.9771, 2070.3435, 1827.7416, 1877…
    ## $ cyto_marker_int_mean  <dbl> 2736.6367, 3746.7310, 3104.7562, 3521.1667, 2710…
    ## $ ratio_marker_int_mean <dbl> 0.6818954, 0.5281420, 0.6742553, 0.5260615, 0.70…
    ## $ sum_marker_int_mean   <dbl> 4568.8711, 5695.7081, 5175.0996, 5350.3126, 4589…

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

    ## Rows: 8
    ## Columns: 11
    ## Groups: screen, cell_line [4]
    ## $ screen                <chr> "171006-TetON-LA-GRP78-PDI_20171006_164850", "17…
    ## $ cell_line             <fct> Uninduced, Uninduced, Induced, Induced, Uninduce…
    ## $ marker                <chr> "GRP78", "PDI", "GRP78", "PDI", "GRP78", "PDI", …
    ## $ cell_n                <dbl> 985.5000, 916.1667, 1007.3333, 1010.5000, 956.16…
    ## $ nuc_area_mean         <dbl> 1995.510, 1988.292, 1698.177, 1666.531, 1874.642…
    ## $ cyto_area_mean        <dbl> 6285.398, 6695.843, 6371.813, 6416.049, 6108.565…
    ## $ cell_area_mean        <dbl> 8280.908, 8684.135, 8069.990, 8082.580, 7983.207…
    ## $ nuc_marker_int_mean   <dbl> 1911.2223, 1892.7044, 1595.8894, 1613.3233, 466.…
    ## $ cyto_marker_int_mean  <dbl> 2845.6897, 3667.0179, 2708.8865, 3944.9290, 612.…
    ## $ ratio_marker_int_mean <dbl> 0.6809448, 0.5238056, 0.6002293, 0.4186882, 0.78…
    ## $ sum_marker_int_mean   <dbl> 4757.1178, 5559.9563, 4304.7759, 5558.6187, 1078…

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
| GRP78  | Induced-Uninduced | -0.093 | -0.371 |  0.185 | 0.286 |
| PDI    | Induced-Uninduced | -0.182 | -0.820 |  0.456 | 0.345 |

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
