#'
#' Contains code that is shared among all `01-field` files,
#' including reading and cleaning up the data frame containing
#' parasitism data.
#'

#' Observation dates to use for maps.
#' We need to define these here to show them in the time series plots.
#' See "maps" section below for more.
maps_dates <- as.Date(c("2015-06-03", "2015-06-12", "2015-06-19",
                        "2013-08-01", "2013-08-09", "2013-08-19"))


# # Color palettes to use in plots (using binned rr_rs values):
# rr_rs_pal <- with(list(pal = inferno, inds = c(80, 60, 20)),
#                   list(color = c(pal(100)[inds], "gray40"),
#                        fill = c(pal(100)[inds], "white"),
#                        fill2 = c(pal(100, alpha = 0.5)[inds], "white")))



#' Calculating parasitism at which relative fitness for resistant aphids is 1.
#' This is based on internal data from `gameofclones` that were derived from
#' Ives et al. (2020).
par_rrrs0 <- with(list(
    # Index for which element in rr/rs lookup table is closest to one:
    idx = which(abs(gameofclones:::rr_rs_lookup$rr_rs - 1) ==
                    min(abs(gameofclones:::rr_rs_lookup$rr_rs - 1))),
    # Average proportion resistant from Ives et al. (2020):
    p_res = 0.48),
    # From this we can calculate what the parasitism is:
    ((1 - p_res) * gameofclones:::rr_rs_lookup$ps + p_res *
         gameofclones:::rr_rs_lookup$pr)[idx])


par_df <- list(
    here("data-raw/parasitism-2001-2016.csv") |>
        read_csv(col_types = cols()) |>
        select(field, Cycle, DateFormat, year, day, para, paraN) |>
        rename(cycle = Cycle, para_n = paraN, date = DateFormat) |>
        mutate(field = paste(field),
               # Field 18 is the same as field 110; the numbers weren't fixed
               # until 2015
               field = ifelse(field == "18", "110", field),
               # These were re-coded to stay as numbers, but I'm converting them
               # back to their original names:
               field = ifelse(field == "506.2", "506N", field),
               field = ifelse(field == "506.1", "506S", field),
               date = as.Date(date)) |>
        filter(!is.na(para), !is.na(para_n)),
    list.files(here("data-raw"), "parasitism-....\\.csv",
               full.names = TRUE) |>
        map_dfr(read_csv, show_col_types = FALSE) |>
        # I can't find this one in any of the Arlington maps...
        filter(Field != "N1902") |>
        # This one had clover for part of it, but where isn't clear
        filter(Field != "349 Clover") |>
        rename_with(tolower) |>
        select(field, cycle, date, starts_with("diss")) |>
        rename_with(function(.x) gsub("^diss_", "", gsub("\\+", "_", .x))) |>
        mutate(para = g_para + g_p_f + r_para + r_p_f,
               para_n = para + g_unpara + g_fungus + r_unpara + r_fungus,
               para = para / para_n) |>
        filter(!is.na(para)) |>
        select(field, cycle, date, para, para_n) |>
        # Because two date formats are used:
        mutate(date1 = as.Date(date, format = "%m/%d/%y"),
               date2 = as.Date(date, format = "%d-%b-%y"),
               # `structure(...` below is to keep `ifelse` from
               # changing to numeric
               date = structure(ifelse(is.na(date1), date2, date1),
                                class = class(date2))) |>
        select(-date1, -date2) |>
        mutate(field = paste(field),
               field = ifelse(field == "506SE", "506S", field),
               year = year(date),
               day = yday(date) - 1)  # <- previous years set Jan 1 as day 0
) |>
    do.call(what = bind_rows) |>
    # Not sure why, but there's a 2020 data point in the 2019 dataset:
    filter(year != 2020) |>
    # We need to have at least 10 aphids dissected:
    filter(para_n >= 10) |>
    mutate(cycle = floor(cycle),
           harvest = lead(cycle, default = tail(cycle, 1)) - cycle > 0 |
               # This is to force the first cell be `1`
               c(TRUE, rep(FALSE, n()-1)),
           year = factor(year, levels = sort(unique(year)))) |>
    #'
    #' I found that these dates have the same exact numbers for all fields
    #' on dates two days before.
    #' These must be duplicates, so I'm removing them.
    #'
    filter(!date %in% as.Date(c("2011-07-14", "2011-07-21", "2011-08-26",
                                "2012-08-17", "2012-07-11", "2012-06-14"))) |>
    # Add the relative fitness for resistance (r_r / r_s)
    mutate(rr_rs = rel_res_fitness(para))

#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Split parasitism into periods ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' To combine the field polygons with the parasitism data and to get average
#' parasitism by approximate date, I first need to prepare the parasitism
#' data for plotting by observation date, where the dates can be a
#' bit different:
#'
obs_par_df <- par_df |>
    select(-cycle, -harvest) |>
    arrange(date) |>
    mutate(obs = cut.Date(date, breaks = "3 days", labels = FALSE) |>
               as.numeric())

#'
#' These observation groups had multiple observations of the same fields,
#' so I split them up further.
#'
#' ```
#' obs_par_df |>
#'     group_by(obs) |>
#'     mutate(repeats = any(duplicated(field))) |>
#'     ungroup() |>
#'     filter(repeats) |>
#'     distinct(year, obs)
#'
#' # # A tibble: 5 × 2
#' #   year    obs
#' #   <fct> <dbl>
#' # 1 2014    371
#' # 2 2015    495
#' # 3 2016    609
#' # 4 2016    623
#' # 5 2016    637
#' ```
#'

obs_par_df <- obs_par_df |>
    group_by(obs) |>
    mutate(repeats = any(duplicated(field))) |>
    ungroup() |>
    mutate(tmpid = interaction(year, obs, drop = TRUE)) |>
    split(~ tmpid) |>
    map_dfr(function(.dd) {
        unq_dates <- length(unique(.dd$date))
        if (.dd$repeats[1]) {
            stopifnot(unq_dates > 1)
            for (i in 2:unq_dates) {
                frac <- (i - 1) / unq_dates
                .dd$obs[.dd$date == unique(.dd$date)[i]] <- .dd$obs[1] + frac
            }
        }
        return(.dd)
    }) |>
    select(-repeats, -tmpid) |>
    # Plus I'm going to manually combine these dates because it's clear that
    # half the fields were sampled one day, half the next.
    mutate(obs = ifelse(date %in% as.Date(c("2016-06-13", "2016-06-14")),
                        mean(obs[date %in% as.Date(c("2016-06-13", "2016-06-14"))]),
                        obs)) |>
    mutate(obs = factor(obs)) |>
    # I'm going to want to filter out groups that only sampled few fields:
    group_by(year) |>
    mutate(n_fields = length(unique(field))) |>
    group_by(year, obs) |>
    mutate(obs_n = n()) |>
    # To simplify plotting by having a single obs number by the obs
    mutate(obs_date = round(mean(date)),
           obs_day = mean(day)) |>
    ungroup()


#'
#' Plotting data through time where points are colored by obs.
#' Looks like the splitting worked well.
#'
# obs_par_df |>
#     mutate(plot_date = as.Date(day, origin = "2022-01-01")) |>
#     ggplot(aes(plot_date, para)) +
#     geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
#     geom_point(aes(color = obs), alpha = 0.5, size = 1) +
#     facet_wrap(~ year, ncol = 1) +
#     scale_color_manual(values = rep(c("#1b9e77", "#d95f02", "#7570b3"), 100),
#                        guide = "none") +
#     scale_x_date(date_breaks = "1 month", date_labels = "%b") +
#     scale_y_continuous("Proportion aphids parasitized", breaks = 0.4*0:2) +
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_text(color = "black"),
#           panel.grid.major.x = element_line(color = "gray80"),
#           strip.text = element_text(size = 9)) +
#     NULL
