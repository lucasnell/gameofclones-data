# gameofclones-data

Data on eco-evolutionary dynamics in pea aphids and parasitoids, 
plus scripts that analyze these data and use `gameofclones` package
to simulate the system.


# Installing packages

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

Similarly, but without version numbers:

```r
pkgs <- c("tidyverse", "here", "viridisLite", "grid", "scales", "patchwork", 
          "ggtext", "lme4", "numDeriv", "parallel", "plot3D", "sf", "ggmap", 
          "rnaturalearth", "rnaturalearthdata")
install.packages(pkgs)
install.packages("remotes")
remotes::install_github("lucasnell/gameofclones")
```



