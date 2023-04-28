#'
#' Contains code that shared among all `scripts/field-*.R` files,
#' including reading and cleaning up the data frame containing
#' parasitism data.
#'

#' Observation dates to use for maps.
#' We need to define these here to show them in the time series plots.
#' See "maps" section below for more.
maps_dates <- as.Date(c("2015-06-03", "2015-06-12", "2015-06-19",
                        "2013-08-01", "2013-08-09", "2013-08-19"))


# Color palettes to use in plots (using binned rr_rs values):
rr_rs_pal <- with(list(pal = inferno, inds = c(80, 60, 20)),
                  list(color = c(pal(100)[inds], "gray40"),
                       fill = c(pal(100)[inds], "white"),
                       fill2 = c(pal(100, alpha = 0.5)[inds], "white")))

# Add factor that breaks rr_rs into bins for plotting:
add_rr_rs_fct <- function(.df) {
    .df |>
        mutate(rr_rs_fct = cut(rr_rs, c(0, 1, 1.05, 1.1, Inf),
                               labels = c("< 1", "1 – 1.05", "1.05 – 1.1", "> 1.1")),
               rr_rs_fct = factor(rr_rs_fct, levels = rev(levels(rr_rs_fct))))
}


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
