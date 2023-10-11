
#'
#' Plots of parasitism at selected dates through space.
#'


source("scripts/00-shared-all.R")

source("scripts/01-field/00-field-shared.R")


library(sf)                 # st_* (e.g., st_read, st_transform, st_crs)
library(ggmap)              # get_googlemap
library(rnaturalearth)      # ne_countries
library(rnaturalearthdata)  # used for ne_countries




# Names of files produced here:
plots_out <- list(par_map = here("plots/01-field/field-mosaic/par-map.pdf"),
                  par_map_leg = here("plots/01-field/field-mosaic/par-map-legend.pdf"),
                  usa = here("plots/01-field/field-mosaic/par-map-usa-inset.pdf"),
                  terrain = here("plots/01-field/field-mosaic/par-map-google-inset.pdf"))



#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------
# Parasitism maps ----
#' ----------------------------------------------------------------------
#' ----------------------------------------------------------------------




fields_sf <- here("data-raw/Arlington.geojson") |>
    st_read() |>
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
        .f$dissected <- .d$dissected
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
    # Add factor that breaks rr_rs into bins for plotting:
    mutate(rr_rs_fct = factor(rr_rs >= 1, levels = c(FALSE, TRUE),
                              labels = c("<1", ">=1")))


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
    geom_sf(aes(size = para, fill = rr_rs_fct), color = "black",
            shape = 21, stroke = 0.75) +
    scale_fill_manual(NULL, guide = "none",
                      values = c("white", "#3DB7E9")) +
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
    save_plot(plots_out$par_map, fields_par_p,
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
    save_plot(plots_out$par_map_leg, fields_par_p_leg, w = 2, h = 3.5)
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
    save_plot(plot_out$usa, usa_map_p, w = 2.5, h = 1.5, bg = NA)
    # If you have pdfcrop installed:
    # system(paste("pdfcrop", here("plots/01-field/field-mosaic/par-map-usa-inset.pdf")))
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
     save_plot(plots_out$terrain, arl_goog_p, w = 1.5, h = 1.5)
} else {
    arl_goog_p
}




