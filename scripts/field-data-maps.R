
source("scripts/_shared-all.R")

#' Contains code that shared among all `scripts/field-*.R` files,
#' including reading and cleaning up the data frame containing
#' parasitism data.
source("scripts/_shared-field.R")

library(sf)                 # st_* (e.g., st_read, st_transform, st_crs)
library(ggmap)              # get_googlemap
library(rnaturalearth)      # ne_countries
library(rnaturalearthdata)  # used for ne_countries



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
# Split parasitism into periods ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' To combine the field polygons with the parasitism data, I first need to
#' prepare the parasitism data for plotting by observation date, where
#' the dates can be a bit different:
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





#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Parasitism maps ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------




fields_sf <- st_read(here("data-raw/Arlington.geojson")) |>
    st_transform(st_crs(3857)) |>
    mutate(geometry = st_centroid(geometry)) |>
    rename(geom = geometry)


obs_fields_par <- obs_par_df |>
    filter(obs_date %in% maps_dates) |>
    # Make sure all fields are present in all dates within year:
    group_by(year, field) |>
    mutate(no = n()) |>
    ungroup() |>
    filter(no == 3) |>
    mutate(obs = obs |> fct_drop(),
           tmpid = 1:n()) |>
    split(~tmpid) |>
    map(function(.d) {
        stopifnot(nrow(.d) == 1)
        .f <- fields_sf |> filter(Name == .d$field)
        stopifnot(nrow(.f) == 1)
        .f$year <- .d$year
        .f$date <- .d$date
        .f$para <- .d$para
        .f$para_n <- .d$para_n
        .f$rr_rs <- .d$rr_rs
        .f$obs <- .d$obs
        .f$obs_date <- .d$obs_date
        return(.f)
    }) |>
    do.call(what = rbind) |>
    mutate(plot_date = factor(paste(obs_date),
                              levels = paste(sort(unique(obs_date))),
                              labels = format(sort(unique(obs_date)),
                                              "%e %b")),
           plot_date_by_yr = (as.integer(plot_date) - 1) %% 3 + 1,
           plot_date_by_yr = factor(plot_date_by_yr)) |>
    add_rr_rs_fct()


xy_lims <- st_bbox(obs_fields_par) |>
    as.list() |> as_tibble() |>
    mutate(xmin = xmin - 400,
           xmax = xmax + 400,
           ymin = ymin - 400,
           ymax = ymax + 400)


fields_par_p <- obs_fields_par |>
    ggplot() +
    geom_rect(xmin = xy_lims$xmin, xmax = xy_lims$xmax,
              ymin = xy_lims$ymin, ymax = xy_lims$ymax,
              fill = NA, color = "black", linewidth = 0.5) +
    geom_sf(aes(size = para, color = rr_rs_fct, fill = rr_rs_fct),
            shape = 21, stroke = 0.75) +
    scale_color_manual(NULL, guide = "none",
                       values = rr_rs_pal$color) +
    scale_fill_manual(NULL, guide = "none",
                      values = rr_rs_pal$fill) +
    scale_size("Proportion\nparasitized",
               limits = c(0, 0.85), range = c(0.5, 8),
               breaks = 0.2 * 0:4) +
    guides(size = guide_legend(override.aes = list(shape = 16))) +
    coord_sf(datum = st_crs(3857),
             xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
             ylim = as.numeric(xy_lims[c("ymin", "ymax")]),
             clip = "off") +
    facet_wrap(~ plot_date, nrow = 2) +
    theme_void() +
    theme(strip.text = element_text(size = 9, margin = margin(0,0,b=3,t=3)),
          legend.position = "none",
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent"))


if (write_plots) {
    save_plot(here("plots/field-data/par-map.pdf"), fields_par_p,
              w = 3.5, h = 2.5, bg = "transparent")
} else {
    fields_par_p
}





#' ----------------------
# Parasitism map separate legend
#' ----------------------


fields_par_p_leg <- function() {
    legend <- (fields_par_p + theme(legend.position = "right")) |>
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
    save_plot(here("plots/field-data/par-map-legend.pdf"), fields_par_p_leg,
              w = 2, h = 3.5)
}


#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# USA map inset ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------


usa_map_p <- ne_countries(country = "united states of america",
                          scale = "medium", returnclass = "sf") |>
    st_cast("POLYGON") |>
    (\(polies) {
        #' Find the polygon whose western edge is the furthest west without
        #' being in the eastern hemisphere. This should be the contiguous USA.
        bbx = lapply(1:nrow(polies),
                     \(i) st_bbox(polies[i,])[["xmax"]]) |>
            do.call(what = c)
        i <- which(bbx == max(bbx[bbx < 0]))
        return(polies[i,])
    })() |>
    st_transform(st_crs(3857)) |>
    ggplot() +
    geom_sf(size = 0.25) +
    geom_point(data = tibble(x = (xy_lims$xmin + xy_lims$xmax) / 2,
                             y = (xy_lims$ymin + xy_lims$ymax) / 2),
               aes(x, y), color = "black", shape = 20, size  = 4) +
    coord_sf(datum = st_crs(3857)) +
    theme_void()

if (write_plots) {
    save_plot(here("plots/field-data/par-map-usa-inset.pdf"), usa_map_p,
              w = 2.5, h = 1.5)
    # If you have pdfcrop installed:
    # system(paste("pdfcrop", here("plots/field-data/par-map-usa-inset.pdf")))
} else {
    usa_map_p
}



#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Satellite map ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------

#'
#' To use get_googlemap, you have to register a key like this:
#' `register_google(key = "mQkzTpiaLYjPqXQBotesgif3EfGL2dbrNVOrogg")`
#' (this is fake key provided in `register_google` docs)
#'

#' Accessed on 2023-04-26
arl_google <- get_googlemap(c(-89.35399, 43.30917), zoom = 13,
                            maptype = "satellite")

# from https://stackoverflow.com/a/50844502
# Define a function to fix the bbox to be in EPSG:3857
ggmap_bbox <- function(map) {
    if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
    # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector,
    # and set the names to what sf::st_bbox expects:
    map_bbox <- setNames(unlist(attr(map, "bb")),
                         c("ymin", "xmin", "ymax", "xmax"))

    # Coonvert the bbox to an sf polygon, transform it to 3857,
    # and convert back to a bbox (convoluted, but it works)
    bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))

    # Overwrite the bbox of the ggmap object with the transformed coordinates
    attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
    attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
    attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
    attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
    map
}


arl_goog_p <- arl_google |>
    ggmap_bbox() |>
    ggmap() +
    coord_sf(datum = st_crs(3857),
             xlim = as.numeric(xy_lims[c("xmin", "xmax")]),
             ylim = as.numeric(xy_lims[c("ymin", "ymax")]),
             clip = "off") +
    theme_void()


if (write_plots) {
     save_plot(here("plots/field-data/par-map-google-inset.pdf"), arl_goog_p,
               w = 1.5, h = 1.5)
} else {
    arl_goog_p
}




