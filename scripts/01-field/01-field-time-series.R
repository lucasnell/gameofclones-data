
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
# Hamiltonella dataset - read and organize ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------



ham_df <- here("data-raw/hamiltonella-2012-2019.csv") |>
    read_csv(col_types = cols()) |>
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






#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Make main plot ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' The parasitism time series plot essentially replicates (and updates)
#' Ives et al. (2020, Nature E&E) paper.
#'
#'


ts_par_df <- obs_par_df |>
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
           plot_date = as.Date(day, origin = "2022-01-01"))

ts_par_mean_df <- ts_par_df |>
    group_by(year, obs) |>
    summarize(plot_date = as.Date(obs_day[1], origin = "2022-01-01"),
              para_sd = sd(para),
              para = mean(para),
              obs_n = obs_n[1],
              .groups = "drop") |>
    filter(obs_n >= 3)




#' The stuff with `max_fc` is to keep the colors from all being similar when
#' not many fields were sampled in a given year
par_ham_ts_p <- ts_par_df |>
    ggplot(aes(plot_date, para)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
    geom_hline(yintercept = par_rrrs0, color = "gray70",
               linewidth = 0.5, linetype = "22") +
    #' -------------
    #' Proportion resistant:
    geom_line(aes(group = field_col), linewidth = 0.5, color = "gray70") +
    geom_line(data = ts_par_mean_df,
              color = "black", linewidth = 0.75) +
    # stat_smooth(method = "gam", se = FALSE, color = "black", linewidth = 0.75) +
    #' -------------
    #' Proportion infected with Hamiltonella:
    geom_point(data = ham_df, aes(y = ham), size = 1, shape = 1,
               color = "#e69f00", stroke = 0.5) +
    geom_segment(data = ham_df |>
                     group_by(year, season) |>
                     summarize(x = mean(plot_date),
                               ham = mean(ham), .groups = "drop"),
                 aes(y = ham, yend = ham, x = x - 7, xend = x + 7),
                 linewidth = 1, color = "#e69f00") +
    #' -------------
    #' Arrows for dates on maps:
    geom_segment(data = tibble(plot_date = yday(maps_dates) |>
                                   as.Date(origin = "2021-12-31"),
                               year = year(maps_dates) |> factor(),
                               para = -0.2),
                 aes(xend = plot_date, yend = -0.05),
                 linewidth = 0.75, linejoin = "mitre", color = "black",
                 arrow = arrow(length = unit(0.15, "lines"), type = "closed")) +
    facet_wrap(~ year, ncol = 3, drop = FALSE) +
    scale_x_date(breaks = ymd(sprintf("2022-%02i-01", 5:11)),
                 labels = c("","Jun","","Aug", "","Oct",""),
                 limits = as.Date(c("2022-04-27", "2022-10-30"))) +
    scale_y_continuous(paste0("Proportion parasitized (lines)<br>",
                              "Proportion *H. defensa* (points)"),
                       breaks = 0.4*0:2) +
    coord_cartesian(ylim = c(-0.2, 0.9), expand = FALSE, clip = "off") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          axis.title.y = element_markdown(hjust = 0.5, lineheight = 1.1),
          legend.title = element_markdown(hjust = 0),
          legend.position = "none",
          strip.text = element_text(size = 9)) +
    scale_color_viridis_d(guide = "none") +
    scale_fill_viridis_d(guide = "none")
# par_ham_ts_p


if (write_plots) {
    save_plot(here("plots/01-field/field-mosaic/par-ham-time.pdf"),
              par_ham_ts_p, w = 5, h = 3)
} else {
    par_ts_wlines_p
}




#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# stdev plot ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

para_mod <- 0.5 #max(ts_par_mean_df$para_sd) / max(ts_par_mean_df$para)



par_sd_p <- ts_par_mean_df |>
    ggplot() +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
    geom_hline(data = ts_par_mean_df |>
                   group_by(year) |>
                   summarize(para_sd = sd(para)),
               aes(yintercept = para_sd), color = "red", linewidth = 1) +
    geom_line(aes(plot_date, para * para_mod), color = "gray70") +
    geom_point(aes(plot_date, para_sd)) +
    # geom_line(aes(plot_date, para_sd)) +
    scale_x_date(breaks = ymd(sprintf("2022-%02i-01", 5:11)),
                 labels = c("","Jun","","Aug", "","Oct",""),
                 limits = as.Date(c("2022-04-27", "2022-10-30"))) +
    scale_y_continuous("Standard deviation",
                       sec.axis = sec_axis(~ . / para_mod, "Parasitism")) +
    facet_wrap(~ year, nrow = 3) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = "black", size = 8),
          strip.text = element_text(size = 9))


if (write_plots) {
    save_plot(here("plots/01-field/field-stdevs.pdf"), par_sd_p, w = 5, h = 3)
} else {
    par_sd_p
}


