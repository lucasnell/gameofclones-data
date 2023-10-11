#'
#' Contains code that is shared among all `01-field` files,
#' including reading and cleaning up the data frame containing
#' parasitism data.
#'


#' This should already be attached via `scripts/00-shared-all.R`,
#' but this is useful for testing:
if(! "tidyverse" %in% (.packages())) library(tidyverse)


#' Observation dates to use for maps.
#' We need to define these here to show them in the time series plots.
#' See "maps" section below for more.
maps_dates <- as.Date(c("2015-06-03", "2015-06-12", "2015-06-19",
                        "2013-08-01", "2013-08-09", "2013-08-19"))


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



#'
#' Read field monitoring data including parasitism and aphid abundances:
#'
field_df <- here("data-raw/parasitism-2011-2019.csv") |>
    read_csv(col_types = cols()) |>
    mutate(para = GrPara + RedPara,
           dissected = para + GrUnpara + RedUnpara) |>
    filter(!is.na(para)) |>
    rename_with(tolower) |>
    mutate(date = as.Date(paste(year, month, day, sep = "-")),
           day = yday(date) - 1L,
           log_aphids = log((sweepspeaaphids+1) / sweepcount),
           # cycle = as.integer(cycle),
           year = factor(year, levels = sort(unique(year)))) |>
    rename(alf_height = height_in) |>
    select(field, cycle, date, year, day, para, dissected,
           log_aphids, alf_height) |>
    arrange(year, field, day) |>
    #'
    #' I found that these dates have the same exact numbers for all fields
    #' on dates two days before.
    #' These must be duplicates, so I'm removing them.
    #'
    filter(!date %in% as.Date(c("2011-07-14", "2011-07-21", "2011-08-26"))) |>
    #' We need to have at least 10 aphids dissected:
    filter(dissected >= 10) |>
    #' Rename fields because some were updated:
    mutate(field = case_when(field == "4110" ~ "418",
                             field == "2110" ~ "218",
                             field == "110/110" ~ "110",
                             field == "635S" ~ "635",
                             TRUE ~ field))


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
obs_par_df <- field_df |>
    #' We don't have records of where these fields are, and this data frame
    #' will be used for maps:
    filter(field != "N1902", field != "349_Clover") |>
    select(-cycle, -log_aphids, -alf_height) |>
    arrange(date) |>
    mutate(obs = cut.Date(date, breaks = "3 days", labels = FALSE) |>
               as.numeric(),
           #' In these figures, parasitism is proportion:
           para = para / dissected)

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
