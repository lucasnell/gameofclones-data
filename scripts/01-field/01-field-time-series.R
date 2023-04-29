
source("scripts/00-shared-all.R")

source("scripts/01-field/00-field-shared.R")



#'
#' Much of the data here are from https://doi.org/10.6084/m9.figshare.11828865.v1
#'
#' To run the scripts below, download this dataset, then rename
#' `Ives et al. 2020 Data Fig2_1.csv` to `parasitism-2001-2016.csv`
#' and
#' `Ives et al. 2020 Data Fig3A.csv` to `symbionts-2012-2017.csv`
#'
#' Then put both inside the `data-raw` folder.
#'





#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Shared objects ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------



#' Function to create time series plots and keep axes and font sizes consistent.
ts_p_maker <- function(.df, .y, .ylab, ...) {
    ggplot(.df, aes(plot_date, {{ .y }})) +
        geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
        list(...) +
        facet_wrap(~ year, ncol = 3, drop = FALSE) +
        scale_x_date(breaks = ymd(sprintf("2022-%02i-01", 5:11)),
                     labels = c("","Jun","","Aug", "","Oct",""),
                     limits = as.Date(c("2022-04-27", "2022-10-30"))) +
        scale_y_continuous(.ylab, limits = c(-0.2, 0.9), breaks = 0.4*0:2) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_text(color = "black", size = 8),
              axis.title.y = element_markdown(hjust = 0.5),
              legend.title = element_markdown(hjust = 0),
              legend.position = "none",
              strip.text = element_text(size = 9))
}






#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Hamiltonella dataset - read and organize ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------



older_ham_df <- here("data-raw/symbionts-2018-2019.csv") |>
    read_csv(col_types = cols()) |>
    base::`[`()



newer_ham_df <- here("data-raw/symbionts-2012-2017.csv") |>
    read_csv(col_types = cols()) |>
    mutate(season = case_when(is.na(date) ~ "fall",
                              late == 1 ~ "fall",
                              late == 0 ~ "spring")) |>
    select(year, season, date, field, clone, Hamiltonella) |>
    rename(ham = Hamiltonella) |>
    # Filling in this date manually based on info from Kerry Oliver:
    mutate(date = ifelse(year == 2012, "9/4/12", date),
           date = as.Date(date, format = "%m/%d/%y"))


ham_df <- bind_rows(newer_ham_df, older_ham_df) |>
    group_by(year, season, date, field) |>
    summarize(ham = mean(ham), n = n(), .groups = "drop") |>
    #' filter for fields where at least 10 aphids were surveyed:
    filter(n >= 10) |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer())) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(yday(date), origin = "2022-01-01"),
           field_id = interaction(field, year, drop = TRUE),
           year = factor(year, levels = 2011:2019))


#' If you want just some basic stats:
ham_df |>
    summarize(n = mean(n), mean = mean(ham), sd = sd(ham),
              min = min(ham), max = max(ham))




#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Make plots ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' The parasitism time series plot essentially replicates (and updates)
#' Ives et al. (2020, Nature E&E) paper.
#'
#' I'm not using lines here bc I think points illustrate my point better,
#' and because code I used to separate by cycle for years 2011--2016 doesn't
#' appear to work for 2017--2019.
#'


rr_rs_legend_title <- str_c("Relative fitness for<br>resistant aphids<br>",
                            "(*r<sub>r</sub>* / *r<sub>s</sub>*)")


par_ts_p <- par_df |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer())) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) |>
    add_rr_rs_fct() |>
    ts_p_maker(para, "Proportion parasitized",
               geom_segment(data = tibble(plot_date = yday(maps_dates) |>
                                              as.Date(origin = "2021-12-31"),
                                          year = year(maps_dates) |> factor(),
                                          para = -0.15),
                            aes(xend = plot_date, yend = -0.05),
                            linewidth = 0.5, linejoin = "mitre",
                            arrow = arrow(length = unit(0.1, "lines"), type = "closed")),
               geom_point(aes(color = rr_rs_fct, fill = rr_rs_fct),
                          size = 1, shape = 21, stroke = 0.5)) +
    scale_color_manual(rr_rs_legend_title,
                       values = rr_rs_pal$color) +
    scale_fill_manual(rr_rs_legend_title,
                      values = rr_rs_pal$fill2) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))


if (write_plots) {
    save_plot(here("plots/field-data/par-time.pdf"), par_ts_p,
              w = 4, h = 2.5)
} else {
    par_ts_p
}

par_ts_p_leg <- function() {
    legend <- (par_ts_p + theme(legend.position = "right")) |>
        (function(a.gplot){
            tmp <- ggplot_gtable(ggplot_build(a.gplot))
            leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
            legend <- tmp$grobs[[leg]]
            legend
        })()
    grid.newpage()
    grid.draw(legend)
}

if (write_plots) {
    save_plot(here("plots/field-data/par-time-legend.pdf"), par_ts_p_leg,
              w = 2, h = 2.5)
}


#' Same thing but without coloring by fitness and with lines connecting fields.
#' The stuff with `max_fc` is to keep the colors from all being similar when
#' not many fields were sampled in a given year
par_ts_wlines_p <- par_df |>
    split(~ year) |>
    map_dfr(~ mutate(.x, field_col = factor(field) |> as.integer(),
                     max_fc = max(field_col))) |>
    mutate(max_max_fc = max(max_fc)) |>
    split(~ year) |>
    map_dfr(\(x) {
        mfc <- x[["max_fc"]][[1]]
        mmfc <- x[["max_max_fc"]][[1]]
        x |>
            mutate(field_col = map_int(field_col,
                                        \(i) {
                                            ss <- seq(1, mmfc, length.out = mfc)
                                            as.integer(round(ss))[i]
                                        }))
    }) |>
    select(-ends_with("max_fc")) |>
    mutate(field_col = factor(field_col),
           # So they show as dates but can be plotted on same scale:
           plot_date = as.Date(day, origin = "2022-01-01")) |>
    add_rr_rs_fct() |>
    ts_p_maker(para, "Proportion parasitized",
               geom_line(aes(color = field_col), linewidth = 0.5),
               geom_point(aes(color = field_col, fill = field_col),
                          size = 0.5, shape = 21)) +
    scale_color_viridis_d(guide = "none") +
    scale_fill_viridis_d(guide = "none")

if (write_plots) {
    save_plot(here("plots/par-time-by-field.pdf"), par_ts_wlines_p,
              w = 6.5, h = 5)
} else {
    par_ts_wlines_p
}




ham_ts_p <- ham_df |>
    ts_p_maker(ham, "Proportion *H. defensa*",
               geom_line(aes(group = field_id), color = "gray60"),
               geom_jitter(size = 1, shape = 1, width = 3, height = 0))


if (write_plots) {
    save_plot(here("plots/field-data/ham-time.pdf"), ham_ts_p,
              w = 4, h = 2.5, seed = 380247925)
} else {
    ham_ts_p
}


