


#' These simulations are to plan and create a priori hypotheses for the
#' eco-evo experiments.
#' I started with a model from a
#' [previous paper](http://doi.wiley.com/10.1890/13-1933.1)
#' from Tony and others, then I edited it to conform to our conditions
#' (most notably, it explicitly simulates alate (winged aphid) production
#' and plant death).
#' I won't be simulating plant death because simulations show that this only
#' adds complexity without changing outcomes.
#' I'll also have alate production be density-independent because previous
#' assays in our fields indicate density has little effect.
#'
#'
#' I'm simulating 1,000 days of interacting aphid and parasitoid wasp
#' populations.
#' There are two fields, one with wasps and another without.
#' Movement between fields must happen via alates, and
#' rates of alate production are density-independent.
#' Every day, 10% of all alates in each field disperse to the other field.
#' There are two aphid lines:
#' The susceptible line is susceptible to parasitism but has a higher
#' growth rate.
#' The resistant line is resistant to parasitism but has a lower growth rate.
#' Shifts in relative frequencies of susceptible and resistant lines
#' represents evolution in the aphid population.
#'


source("scripts/00-shared-all.R")


# Names of files produced here:
plots_out <- list(main = imap(c("no_dispersal", "dispersal"),
                              \(x, i) {
                                  paste0(here("plots/03-experiments/exper-and-sims/"),
                                         sprintf("sims-%i-%s.pdf", i, x))
                              }) |>
                      set_names(c("no_dispersal", "dispersal")),
                  wasps = here("plots/03-experiments/sims-exper-wasps-t0.pdf"),
                  p_res = here("plots/03-experiments/sims-exper-resist-t0.pdf"),
                  perturb = here("plots/03-experiments/sims-exper-perturbs.pdf"))




# To maintain the same y-axis limits:
wasp_mod <- 5.078307  # <-- should be max(disp_wasps$wasps) / max(log1p(disp_aphids$N))
max_N <- 7.785  # <-- should be ceiling(max(log(disp_aphids$N)) * 1e3) / 1e3
y_breaks <- log(10^(0:3))
y_labs <- 10^(0:3)




#'
#' Base function to do all the simulations.
#' Each section below builds on this one.
#'
do_base_sims <- function(...) {

    .args <- list(clonal_lines = c(line_s, line_r))
    .other_args <- list(...)
    if (length(.other_args) > 0) {
        stopifnot(!is.null(names(.other_args)))
        stopifnot(sum(duplicated(names(.other_args))) == 0)
        stopifnot(all(names(.other_args) %in% names(formals(sim_experiments))))
        for (n in names(.other_args)) .args[[n]] <- .other_args[[n]]
    }

    .sims <- do.call(sim_experiments, .args)

    .sims[["wasps"]] <- .sims[["wasps"]] |>
        select(-rep) |>
        mutate(field = factor(field, levels = 2:1,
                             labels = para_lvls))

    .sims[["aphids"]] <- .sims[["aphids"]] |>
        select(-rep, -plant) |>
        mutate(field = factor(field, levels = 2:1,
                             labels = para_lvls),
               line = factor(line, levels = c("resistant", "susceptible")))

    return(.sims)
}


#' Group aphids by type.
#' This never groups by line unless you're requesting mummies or all
#' parasitized aphids (which includes mummies).
#'
#' This is a separate function from `do_base_sims` so that output from
#' `do_base_sims` can be used to look at alates, etc.
#'
#' `.sim_aphids_df` should be the `"aphids"` field from `sim_gameofclones` output
#'
#' `.type` takes the following options:
#'     * `"living"`: all aphids including parasitized;
#'         this is better for plotting aphid abundances
#'     * `"unparasitized"`: unparasitized aphids;
#'         this is better for community matrices
#'     * `"mummies"`: just mummies
#'     * `"parasitized"`: mummies + parasitized aphids;
#'         this is better for community matrices
#'     * `"alates"`: alates
#'     * `"apterous"`: apterous
#'
group_aphids <- function(.sim_aphids_df, .type) {

    stopifnot(is.data.frame(.sim_aphids_df))

    needed_cols <- c("field", "time", "line", "type", "N")
    stopifnot(all(needed_cols %in% colnames(.sim_aphids_df)))

    # Other columns to potentially keep:
    other_cols <- colnames(.sim_aphids_df)
    other_cols <- other_cols[!other_cols %in% needed_cols]
    keep_other_cols <- c()
    if (length(other_cols) > 0) {
        for (.c in other_cols) {
            # Only keep them if they're identical for the entire data frame
            if (length(unique(.sim_aphids_df[[.c]])) == 1) {
                keep_other_cols <- c(keep_other_cols, .c)
            }
        }
    }

    if (.type == "living") {
        .out <- .sim_aphids_df |>
            filter(type != "mummy") |>
            group_by(field, time, line) |>
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "unparasitized") {
        .out <- .sim_aphids_df |>
            filter(type != "mummy", type != "parasitized") |>
            group_by(field, time, line) |>
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "mummies") {
        .out <- .sim_aphids_df |>
            filter(type == "mummy") |>
            group_by(field, time) |>
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "parasitized") {
        .out <- .sim_aphids_df |>
            filter(type == "mummy" | type == "parasitized") |>
            group_by(field, time) |>
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "alates") {
        .out <- .sim_aphids_df |>
            filter(type == "alate") |>
            group_by(field, time, line) |>
            summarize(N = sum(N), .groups = "drop")
    } else if (.type == "apterous") {
        .out <- .sim_aphids_df |>
            filter(type == "apterous") |>
            group_by(field, time, line) |>
            summarize(N = sum(N), .groups = "drop")
    } else {
        stop("ERROR: strange .type requested. Options are ",
             paste(c('"living"', '"mummies"', '"unparasitized"',
                     '"parasitized"', '"alates"', '"apterous"'),
                   collapse = ", "))
    }

    for (koc in keep_other_cols) .out[[koc]] <- .sim_aphids_df[[koc]][1]

    return(.out)

}

# ============================================================================*
# ============================================================================*

# Main simulations ----

# ============================================================================*
# ============================================================================*


main_sims <- rep(list(NA), 2)
names(main_sims) <- c("no aphid dispersal", "aphid dispersal")
main_sims[["no aphid dispersal"]] <- do_base_sims(alate_field_disp_p = 0)
main_sims[["aphid dispersal"]] <- do_base_sims()
for (d in names(main_sims)) {
    for (n in c("wasps", "aphids")) {
        main_sims[[d]][[n]] <- main_sims[[d]][[n]] |>
            mutate(disp = d)
    }
}; rm(d, n)


main_aphids <- map_dfr(main_sims, ~ group_aphids(.x[["aphids"]], "living")) |>
    mutate(disp = factor(disp, levels = names(main_sims)))
main_mummies <- map_dfr(main_sims, ~ group_aphids(.x[["aphids"]], "mummies")) |>
    mutate(disp = factor(disp, levels = names(main_sims)))
main_wasps <- map_dfr(main_sims, ~ .x[["wasps"]]) |>
    mutate(disp = factor(disp, levels = names(main_sims)))


main_p_list <- levels(main_aphids$disp) |>
    set_names(\(x) str_replace_all(x, " ", "_") |> str_remove("aphid_")) |>
    map(
    function(d) {
        # d = levels(main_aphids$disp)[1]
        # rm(mad, wad, d)
        mad <- main_aphids |>
            filter(disp == d)
        wad <- main_wasps |>
            filter(disp == d)
        mad |>
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N)) |>
            ggplot(aes(time, N)) +
            geom_line(data = wad |>
                          mutate(N = wasps / wasp_mod),
                      color = wasp_color) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(aes(color = line)) +
            scale_color_manual(values = clone_pal, guide = "none") +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Parasitoid abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days", limits = c(0, 250),
                               breaks = 0:5 * 50) +
            facet_grid( ~ field, scales = "fixed") +
            theme(strip.text = element_blank(),
                  axis.title = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size = 9),
                  panel.background = element_rect(fill = "transparent"),
                  plot.background = element_rect(fill = "transparent",
                                                 color = NA)) +
            coord_cartesian(clip = FALSE, ylim = c(0, max_N))
    })


stopifnot(all(names(main_p_list) %in% names(plots_out$main)))



if (write_plots) {
    for (n in names(main_p_list)) {
        save_plot(plots_out$main[[n]], main_p_list[[n]], 4, 1.25)
    }; rm(n)
} else wrap_plots(main_p_list, nrow = 1)







# ============================================================================*
# ============================================================================*

# Starting parasitoid abundance ----

# ============================================================================*
# ============================================================================*


#'
#' Do a call to `sim_gameofclones` that can vary by the proportion of alates
#' that move between fields every day (`alate_field_disp_p`).
#' How does dispersal affect maintenance of diversity in the resistance trait?
#'
do_wasp_sims <- function(.wasp_density_0, .max_t = 500) {

    stopifnot(is.numeric(.wasp_density_0) && length(.wasp_density_0) == 1 &&
                  .wasp_density_0 > 0)

    .sims <- do_base_sims(wasp_density_0 = c(.wasp_density_0, 0),
                          max_t = .max_t)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] |>
            mutate(wasps0 = .wasp_density_0)

    return(.sims)
}



wasp_sims <- safe_mclapply(c(1, 2, 3, 10), do_wasp_sims)


wasp_aphids <- map_dfr(wasp_sims, ~ group_aphids(.x[["aphids"]], "living")) |>
    mutate(wasps0 = factor(wasps0, labels = sprintf("starting parasitoid(s) = %i",
                                                      sort(unique(wasps0)))))
wasp_mummies <- map_dfr(wasp_sims, ~ group_aphids(.x[["aphids"]], "mummies")) |>
    mutate(wasps0 = factor(wasps0, labels = sprintf("starting parasitoid(s) = %i",
                                                      sort(unique(wasps0)))))
wasp_wasps <- map_dfr(wasp_sims, ~ .x[["wasps"]]) |>
    mutate(wasps0 = factor(wasps0, labels = sprintf("starting parasitoid(s) = %i",
                                                      sort(unique(wasps0)))))


wasp_p_list <- map(
    levels(wasp_aphids$wasps0),
    function(w0) {
        aw0 <- wasp_aphids |>
            filter(wasps0 == w0) |>
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N))
        ww0 <- wasp_wasps |>
            filter(wasps0 == w0)
        p <- aw0 |>
            filter(!is.na(N)) |>
            ggplot(aes(time, N)) +
            ggtitle(w0) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(data = ww0 |>
                          mutate(N = wasps / wasp_mod),
                      color = wasp_color) +
            geom_line(aes(color = line)) +
            scale_color_manual(NULL, values = clone_pal) +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Parasitoid abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap( ~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5)) +
            coord_cartesian(clip = FALSE, ylim = c(0, max_N))
        if (any(is.na(aw0$N))) {
            ext <- filter(aw0, is.na(N)) |>
                group_by(field) |>
                filter(time == min(time)) |>
                ungroup() |>
                mutate(N = 0)
            p <- p +
                geom_point(data = ext, aes(color = line), shape = 4, size = 2) +
                theme(legend.position = "none")
        }
        if (w0 == levels(wasp_aphids$wasps0)[1]) {
            p <- p +
                geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
                                        time = 500, N = log(50)),
                          aes(label = "parasitoids"), size = 9 / 2.8, hjust = 1, vjust = 0.5,
                          color = col_pal$w)
        }
        return(p)
    })



wasp_p <- wrap_plots(wasp_p_list, ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))


if (write_plots) {
    save_plot(plots_out$wasps, wasp_p, 10, 5)
} else wasp_p






# ============================================================================*
# ============================================================================*

# Stability - starts ----

# ============================================================================*
# ============================================================================*


#'
#' Vary the proportion of starting aphids that are resistant to parasitism.
#' How stable are results to starting conditions?
#'
do_stable_start_sims <- function(.p_res, .max_t = 5000) {

    if (round(sum(line_s$density_0[line_s$density_0 > 0]), 10) !=
        round(sum(line_r$density_0[line_r$density_0 > 0]), 10)) {
        stop("ERROR: starting densities in `line_s` and `line_r` should be ",
             "the same. (Let the later functions change this.)")
    }

    .n <- 2 * sum(line_s$density_0[line_s$density_0 > 0])

    .line_s <- line_s
    .line_r <- line_r

    # Make both densities sum to 1
    .line_s$density_0 <- .line_s$density_0 / sum(.line_s$density_0)
    .line_r$density_0 <- .line_r$density_0 /  sum(.line_r$density_0)

    # Now make them sum to correct abundances:
    .line_s$density_0 <- .line_s$density_0 * .n * (1 - .p_res)
    .line_r$density_0 <- .line_r$density_0 * .n * .p_res

    .sims <- do_base_sims(clonal_lines = c(.line_r, .line_s),
                          max_t = .max_t)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] |>
        mutate(p_res = .p_res)

    return(.sims)
}




p_res_lvls <- c(0.35, 0.4, 0.8, 0.85)

stable_start_sims <- safe_mclapply(p_res_lvls, do_stable_start_sims)


stable_start_aphids <- map_dfr(stable_start_sims,
                               ~ group_aphids(.x[["aphids"]], "living")) |>
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))
stable_start_mummies <- map_dfr(stable_start_sims,
                                ~ group_aphids(.x[["aphids"]], "mummies")) |>
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))
stable_start_wasps <- map_dfr(stable_start_sims, ~ .x[["wasps"]]) |>
    mutate(p_res = factor(p_res, labels = sprintf("%.0f%% starting resistance",
                                                  sort(unique(p_res)) * 100)))



stable_start_p_list <- map(
    levels(stable_start_aphids$p_res),
    function(pr) {
        ssa <- stable_start_aphids |>
            filter(p_res == pr, time <= 500) |>
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N))
        ssw <- stable_start_wasps |>
            filter(p_res == pr, time <= 500)
        p <- ssa |>
            filter(!is.na(N)) |>
            ggplot(aes(time, N)) +
            ggtitle(pr) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(data = ssw |>
                          mutate(N = wasps / wasp_mod),
                      color = wasp_color) +
            geom_line(aes(color = line)) +
            scale_color_manual(NULL, values = clone_pal) +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Parasitoid abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap(~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5)) +
            coord_cartesian(clip = FALSE, ylim = c(0, max_N))
        if (any(is.na(ssa$N))) {
            ext <- filter(ssa, is.na(N)) |>
                group_by(field) |>
                filter(time == min(time)) |>
                ungroup() |>
                mutate(N = 0)
            p <- p +
                geom_point(data = ext, aes(color = line), shape = 4, size = 2)
        } else {
            p <- p + theme(legend.position = "none")
        }
        if (pr == levels(stable_start_aphids$p_res)[1]) {
            lab_d <- tibble(field = factor(para_lvls[2], levels = para_lvls),
                         time = 500, N = log(70))
            p <- p +
                geom_text(data = lab_d, aes(label = "parasitoids"), size = 9 / 2.8,
                          hjust = 1, vjust = 1, color = col_pal$w)
        }
        return(p)
    })



stable_start_p <- wrap_plots(stable_start_p_list, ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

if (write_plots) {
    save_plot(plots_out$p_res, stable_start_p, 10, 5)
} else stable_start_p











# ============================================================================*
# ============================================================================*

# Stability - perturbs ----

# ============================================================================*
# ============================================================================*




#'
#' Vary the proportion of starting aphids that are resistant to parasitism.
#' How stable are results to starting conditions?
#'
do_stable_perturb_sims <- function(.perturb_df) {

    .max_t = max(.perturb_df$when) + 5e3

    .sims <- do_base_sims(max_t = .max_t,
                          perturb = .perturb_df)

    for (n in c("wasps", "aphids")) .sims[[n]] <- .sims[[n]] |>
        mutate(who = paste(unique(.perturb_df$who), collapse = " & "))

    return(.sims)
}

stable_perturb_sims <- list(tibble(when = 10e3, where = 1:2,
                                   who = "resistant", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "susceptible", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "wasps", how = 0.5),
                            tibble(when = 10e3, where = 1:2,
                                   who = "mummies", how = 0.5)) |>
    safe_mclapply(do_stable_perturb_sims)

pert_t_range <- c(9.8, 10.5) * 1e3

stable_perturb_aphids <- map_dfr(stable_perturb_sims,
                                 ~ group_aphids(.x[["aphids"]], "living")) |>
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps", "mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "adult parasitoids",
                                   "mummies & parasitized aphids"))) |>
    filter(time >= pert_t_range[1], time <= pert_t_range[2])
stable_perturb_mummies <- map_dfr(stable_perturb_sims,
                                  ~ group_aphids(.x[["aphids"]], "mummies")) |>
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps", "mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "adult parasitoids",
                                   "mummies & parasitized aphids"))) |>
    filter(time >= pert_t_range[1], time <= pert_t_range[2])
stable_perturb_wasps <- map_dfr(stable_perturb_sims, ~ .x[["wasps"]]) |>
    mutate(who = factor(who, levels = c("resistant", "susceptible",
                                        "wasps", "mummies"),
                        labels = c("resistant aphids", "susceptible aphids",
                                   "adult parasitoids",
                                   "mummies & parasitized aphids"))) |>
    filter(time >= pert_t_range[1], time <= pert_t_range[2])





stable_perturb_p_list <- map(
    levels(stable_perturb_aphids$who),
    function(w) {
        ssa <- stable_perturb_aphids |>
            filter(who == w)
        ssw <- stable_perturb_wasps |>
            filter(who == w)
        p <- ssa |>
            mutate(N = ifelse(N == 0, NA, N),
                   N = log(N)) |>
            ggplot(aes(time, N)) +
            ggtitle(w) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(data = ssw |>
                          mutate(N = wasps / wasp_mod),
                      color = wasp_color) +
            geom_line(aes(color = line)) +
            scale_color_manual(NULL, values = clone_pal) +
            scale_y_continuous("Aphid abundance",
                               breaks = y_breaks,
                               labels = y_labs,
                               sec.axis = sec_axis(~ . * wasp_mod,
                                                   "Parasitoid abundance",
                                                   breaks = 0:2 * 20)) +
            scale_x_continuous("Days") +
            facet_wrap(~ field, nrow = 1) +
            theme(strip.text = element_text(size = 10),
                  plot.title = element_text(size = 12, hjust = 0.5)) +
            coord_cartesian(clip = FALSE, ylim = c(0, max_N))
        if (w == levels(stable_perturb_aphids$who)[1]) {
            lab_d <- tibble(field = factor(para_lvls[2], levels = para_lvls),
                            time = 10.1e3, N = log(4))
            p <- p +
                geom_text(data = lab_d, aes(label = "parasitoids"), size = 9 / 2.8,
                          hjust = 0, vjust = 0, color = col_pal$w)
        }
        return(p)
    })




stable_perturb_p <- wrap_plots(stable_perturb_p_list, ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 14, face = "bold"))

if (write_plots) {
    save_plot(plots_out$perturb, stable_perturb_p, 10, 5)
} else stable_perturb_p



