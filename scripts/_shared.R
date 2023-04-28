

#' This file is sourced by all R files inside `scripts` folder.
#' It loads all packages (except those used only once),
#' sets the default ggplot2 theme,
#' sets number of threads to use for multithreaded operations
#' defines function to only do multithreading on unix systems,
#' and defines the function used to save all plots.

suppressPackageStartupMessages({
    library(tidyverse)
    library(gameofclones)
    library(here)
    library(viridisLite)
    library(grid)
    library(scales)
    library(patchwork)
    library(ggtext)
    library(lme4)
    library(numDeriv)
    library(parallel)
})

#' The only other packages I used in any scripts are...
#'
#' `plot3D` (`figS8.R`)
#' `sf`, `ggmap`, `rnaturalearth`, and `rnaturalearthdata` (`field-data-maps.R`)
#'
#' See `README.md` in the parent directory for installing them all.
#'


#' Create new plot files?
#' `FALSE` allows you simply view them without writing new pdf files to disk.
#' `TRUE` will generate the plots in the `plots` directory.
write_plots <- FALSE


# Change to max number of threads you want to use:
.n_threads <- max(detectCores()-2L, 1L)

options(boot.ncpus = .n_threads)
options(mc.cores = .n_threads)
options(boot.parallel = ifelse(.Platform$OS.type != "windows", "multicore", "no"))


# colors for resistant, susceptible, and parasitoid wasps, respectively
col_pal <- list(r = viridis(100)[50],
                s = viridis(100)[95],
                w = viridis(100)[1])
# Just for two clones:
clone_pal <- c(col_pal$r, col_pal$s)
# Reduce opacity for parasitoid fill:
wasp_fill <- alpha(col_pal$w, 0.6)


#' Run mclapply on unix system, run lapply on anything else.
#'
#' @param X a vector (atomic or list) or an expressions vector.
#'     Other objects (including classed objects) will be coerced by `as.list`.
#' @param FUN the function to be applied to each element of `X`
#' @param ... Other arguments passed to mclapply (or lapply on Windows).
#'
safe_mclapply <- function(X, FUN, ...) {
    if (.Platform$OS.type == "unix") {
        out <- mclapply(X, FUN, ...)
    } else {
        out <- lapply(X, FUN, ...)
    }
    return(out)
}


# ggplot2 theme:
theme_set(theme_classic() %+replace%
              theme(strip.background = element_blank(),
                    strip.text = element_text(size = 11),
                    legend.background = element_blank(),
                    plot.title = element_text(size = 14, hjust = 0)))



#' Save a plot to a PDF file using `cairo_pdf`
#'
#' @param fn Filename of plot
#' @param p Plot object or function to create plot
#' @param w Width in inches
#' @param h Height in inches
#' @param seed Integer to seed RNG for consistent jittering.
#' @param ... Other arguments to `cairo_pdf`
#'
#' @return `NULL`
#' @export
#'
save_plot <- function(fn, p, w, h, seed = NULL, ...) {
    ext <- tail(strsplit(fn, "\\.")[[1]], 1)
    if (ext != "pdf") stop("ERROR: file name extension must be \".pdf\"")
    fn_dir <- dirname(fn)
    if (!dir.exists(fn_dir)) stop("ERROR: `", fn_dir, "` doesn't exist")
    if (!is.null(seed)) set.seed(seed)
    cairo_pdf(filename = fn, width = w, height = h, ...)
    if (is.function(p)) {
        p()
    } else {
        plot(p)
    }
    dev.off()
    invisible(NULL)
}





# Define clonal lines. Don't re-define these!
# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
# Resistant line: high resistance, low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      resistant = TRUE,
                      surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      surv_adult_apterous = "low",
                      repro_apterous = "low")
