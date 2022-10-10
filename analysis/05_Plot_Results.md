Results Plots
================

- <a href="#cibersort" id="toc-cibersort">CIBERSORT</a>
  - <a href="#cibersort-results-online"
    id="toc-cibersort-results-online">CIBERSORT results online</a>
  - <a href="#cibersort-results-docker"
    id="toc-cibersort-results-docker">CIBERSORT results docker</a>
  - <a
    href="#separate-plot-per-subpopulation-using-online-version-of-cibersort"
    id="toc-separate-plot-per-subpopulation-using-online-version-of-cibersort">Separate
    plot per subpopulation, using online version of CIBERSORT</a>

# CIBERSORT

``` r
metadata <- read.csv(file = here("intermediate_files/metadata/reformat_samples_extra_info.csv"))
colnames(metadata)[1] <- "Mixture"
metadata$Mixture <- as.character(metadata$Mixture)
metadata$Condition <- ifelse(metadata$Condition == "Tumour", "CAF", "TAN")
```

## CIBERSORT results online

``` r
cibersort_plot_online
```

![](05_Plot_Results_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## Saving 7 x 5 in image

    ## [1] "saved to ... outfiles/cibersort_results_online_2022-10-10.png"

## CIBERSORT results docker

``` r
cibersort_plot_docker
```

![](05_Plot_Results_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## Saving 7 x 5 in image

    ## [1] "saved to ... outfiles/cibersort_results_docker_2022-10-10.png"

``` r
cibersort_grid
```

![](05_Plot_Results_files/figure-gfm/CIBERSORT%20plots-1.png)<!-- -->

    ## [1] "saved to ... outfiles/cibersort_plot_combined_batch_corrected_2022-10-10.png"

## Separate plot per subpopulation, using online version of CIBERSORT

``` r
cibersort_output_metadata <- full_join(cibersort_output, metadata, by = "Mixture")
ggplot(cibersort_output_metadata, 
                                aes(x = as.character(Mixture), y = S1, fill = Condition
)) +
  geom_col() + 
  ggtitle("CIBERSORT results S1") +   
  theme(plot.title = element_text(hjust = 0.5),  axis.text = element_text(size = 8, angle = 90)) +
  xlab("Mixture") + 
  ylab("Proportion") +
  theme(#panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0000001)) 
```

![](05_Plot_Results_files/figure-gfm/S1%20plot-1.png)<!-- -->

``` r
ggplot(cibersort_output_metadata, 
                                aes(x = as.character(Mixture), y = S3, fill = Condition
)) +
  geom_col() + 
  ggtitle("CIBERSORT results S3") +   
  theme(plot.title = element_text(hjust = 0.5),  axis.text = element_text(size = 8, angle = 90)) +
  xlab("Mixture") + 
  ylab("Proportion") +
  theme(#panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0000001)) 
```

![](05_Plot_Results_files/figure-gfm/S3%20plot-1.png)<!-- -->

``` r
ggplot(cibersort_output_metadata, 
                                aes(x = as.character(Mixture), y = S4, fill = Condition
)) +
  geom_col() + 
  ggtitle("CIBERSORT results S4") +   
  theme(plot.title = element_text(hjust = 0.5),  axis.text = element_text(size = 8, angle = 90)) +
  xlab("Mixture") + 
  ylab("Proportion") +
  theme(#panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0000001)) 
```

![](05_Plot_Results_files/figure-gfm/S4%20plot-1.png)<!-- -->
