# gameofclones-data

Data on eco-evolutionary dynamics in pea aphids and parasitoids, 
plus scripts that analyze these data and use `gameofclones` package
to simulate the system.

Most of the code here was written by
[Lucas A. Nell](https://github.com/lucasnell) with some by 
[Anthony R. Ives](https://github.com/arives).
You can send questions to LAN unless a script explicitly says it was
written by ARI.



# Organization

Folder contents:

```
.
├── data-interm
├── data-raw
├── gameofclones-data.Rproj
├── LICENSE.md
├── plots
├── README.md
├── renv
├── renv.lock
└── scripts
```

The following folders should be present:

* `data-interm`: intermediate files with output from time-consuming simulations
* `data-raw`: raw data from the field, experiments, and preliminary assays
* `gameofclones-data.Rproj`: file saving this RStudio Project's preferences
* `LICENSE.md`: file containing the CC0 license for this repository
* `plots`: all the primary plots in the manuscript
* `README.md`: this file
* `renv`: files necessary to (optionally) use the R package `renv` to install
  the same R packages with the same versions we used.
* `renv.lock`: another file necessary for using `renv`
* `scripts`: the R scripts containing all the analyses and generating
  all the figures in the manuscript

The folders `data-interm`, `data-raw`, `plots`, and `scripts` have separate
`README.md` files that have more information on each.



# Replicating environment

I used R version 4.3.1 (platform: aarch64-apple-darwin20) for all my scripts.

This project uses the `renv` package, so if you want to use this, you must
first install it:

```r
install.packages("renv")
```

Then to install all the packages I used for these analyses, you can simply run
the following while having this project's main directory as your working
directory:

```r
renv::restore()
```


If you'd rather avoid `renv`, then you can install all the packages 
(in the versions I used) this way:

```r
pkgs <- c("tidyverse@2.0.0", "here@1.0.1", "viridisLite@0.4.1", "grid@4.3.0", 
          "scales@1.2.1", "patchwork@1.1.2", "ggtext@0.1.2", "lme4@1.1.33", 
          "numDeriv@2016.8.1.1", "parallel@4.3.0", "plot3D@1.4", 
          "sf@1.0.12", "ggmap@3.0.2", "rnaturalearth@0.3.2", 
          "rnaturalearthdata@0.1.0")
install.packages(pkgs)
install.packages("remotes")
remotes::install_github("lucasnell/gameofclones@v1.0.1")
```

Note that these only install the proper versions of the packages I manually 
installed, so dependencies might vary from what I used.


Similarly, but without version numbers at all:

```r
pkgs <- c("tidyverse", "here", "viridisLite", "grid", "scales", "patchwork", 
          "ggtext", "lme4", "numDeriv", "parallel", "plot3D", "sf", "ggmap", 
          "rnaturalearth", "rnaturalearthdata")
install.packages(pkgs)
install.packages("remotes")
remotes::install_github("lucasnell/gameofclones")
```



# Data from other sources


We used symbiont and parasitism data from
<https://doi.org/10.6084/m9.figshare.11828865.v1>.
We used the parasitism data unchanged.
For symbiont data, we ran the following to process the data:

```r
older_ham_df <- "<path>/Ives et al. 2020 Data Fig3A.csv" |>
    read_csv(col_types = cols()) |>
    mutate(season = case_when(is.na(date) ~ "fall",
                              late == 1 ~ "fall",
                              late == 0 ~ "spring")) |>
    select(year, season, date, field, clone, Hamiltonella) |>
    rename(ham = Hamiltonella) |>
    # Filling in this date manually based on info from Kerry Oliver:
    mutate(date = ifelse(year == 2012, "9/4/12", date),
           date = as.Date(date, format = "%m/%d/%y"))
```

We then merged it with our newer data and put this in the `data-raw`
folder as `hamiltonella-2012-2019.csv`.
