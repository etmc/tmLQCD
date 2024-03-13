We use tmLQCD hierarchical timer and R to produce profiles of the HMC.

To collect time informations from the file `example_log.out` type

```
    Rscript make_profile.R example_log.out output_basename
```

The `make_profile.R` file will first collect the data and then render the file `profile.Rmd`, producing the pdf `output_basename.pdf`, where the name of the output file is set by passing something for `output_basename`.

Alternatively, the script can also be sourced from an interactive session:

```
> infile <- "example_log.out"
> outbase <- "output_basename"
> source("make_profile.R")
```

The R packages `dplyr, data.tree, stringr` and `microseq` are required for data extraction and the R packages `ggplot2, knitr, treemap, hadron` and `kableExtra` are required to render the Rmarkdown document.

The raw and summary data are stored in `profile.RData` which includes the `data.tree` structure `t_tree` which nicely shows the call stack.
