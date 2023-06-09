
source("scripts/00-shared-all.R")



# What to name wasp and no wasp cages for plot:
cage_lvls <- c("parasitism cage" = "wasp", "no parasitism cage" = "no wasp")


exp_df <- here("data-raw/experiments-main.csv") |>
    read_csv(col_types = cols()) |>
    #'
    #' Rep 7 failed bc wasps got into the no-wasp cage, and we didn't
    #' adequately respond to this.
    #'
    filter(rep != 7) |>
    mutate(date = as.Date(date, format = "%m/%d/%Y"),
           start_date = as.Date(start_date, format = "%m/%d/%Y"),
           days = difftime(date, start_date, units = "days") |>
               as.integer(),
           # Adjusting for CES's different way of counting alate OUT
           # before "2021-08-05"
           alates_total_red = ifelse(date < as.Date("2021-08-05") &
                                       observer == "CES",
                                   alates_total_red * 2, alates_total_red),
           alates_total_green = ifelse(date < as.Date("2021-08-05") &
                                         observer == "CES",
                                   alates_total_green * 2,
                                   alates_total_green)) |>
    # Adjusting for these not being input for isolated treatment:
    mutate(across(starts_with("alates"), ~ ifelse(treatment == "ISOLATED" &
                                                      is.na(.), 0, .))) |>
    mutate(treatment = treatment |>
               tolower() |>
               factor(levels = c("dispersal", "isolated")),
           cage = cage |>
               tolower() |>
               factor(levels = c("no wasp", "wasp")),
           rep = factor(rep, levels = sort(unique(rep)))) |>
    select(everything(), -start_date, -red_line,
           -green_line) |>
    # To show when an experimental cage was terminated bc of an extinction:
    group_by(treatment, rep, cage) |>
    mutate(terminated = date == max(date[plant1_red >= 0])) |>
    ungroup() |>
    # Using `- 4` to prevent all cages that weren't sampled on the
    # latest date from returning TRUE here
    mutate(terminated = (terminated &
                             date < (max(date) - 4) &
                             days < (max(days) - 4)))



# "Pesky" wasps - wasps that made it into no-wasp cages:

pesky_df <- here("data-raw/experiments-pesky_wasps.csv") |>
    read_csv(col_types = cols()) |>
    select(rep, date, mummies, starts_with("adult")) |>
    #' On 3/16 and 3/22 I removed wasps/mummies but didn't record the number.
    filter(!is.na(mummies)) |>
    rename(females = `adult females`,
           males = `adult males`,
           unsexed = `adults unk.`) |>
    mutate(wasps = females + males + unsexed,
           date = as.Date(date, "%d-%b-%y")) |>
    select(rep, date, mummies, wasps, everything()) |>
    group_by(date, rep) |>
    summarize_all(sum) |>
    ungroup() |>
    # These are changed to make it join easier with the other wasp data:
    mutate(rep = factor(rep, levels = levels(exp_df$rep)),
           cage = factor("no wasp", levels = levels(exp_df$cage)),
           treatment = ifelse(rep %in% c(6,8,10,11),"dispersal","isolated") |>
               factor(levels = levels(exp_df$treatment)),
           start_date = map_chr(rep, ~ filter(exp_df, rep == .x)[["date"]] |>
                                    min() |> paste()) |> as.Date(),
           days = difftime(date, start_date, units = "days") |>
               as.integer()) |>
    select(-start_date) |>
    filter(date <= max(exp_df$date))




# ============================================================================*
# ============================================================================*

# Summarize by cage ----

# ============================================================================*
# ============================================================================*

# Cage-level wasp/mummy data:
wasp_cage_df <- exp_df |>
    # We don't want to include dates from pesky wasps bc it's redundant:
    filter(!interaction(rep, date, cage, drop = TRUE) %in%
               interaction(pesky_df$rep, pesky_df$date, pesky_df$cage,
                           drop = TRUE)) |>
    select(treatment, rep, cage, days, date, wasps, mummies, wasps_rm, mumm_rm) |>
    filter(wasps >= 0) |>
    mutate(pesky = 0) |>
    bind_rows(pesky_df |>
                  select(treatment, rep, cage, days, date, wasps, mummies) |>
                  mutate(wasps_rm = wasps, mumm_rm = mummies,
                         pesky = 1)) |>
    mutate(pesky = factor(pesky, levels = 0:1),
           cage = fct_recode(cage, !!!cage_lvls))



# Cage-level aphid abundance data:
aphid_cage_df <- exp_df |>
    mutate(red = rowSums(across(starts_with("plant") & ends_with("_red"))),
           green = rowSums(across(starts_with("plant") &
                                      ends_with("_green")))) |>
    select(-contains("plant"), -starts_with("wasps"), -mummies,
           -observer, -notes, -starts_with("alates")) |>
    pivot_longer(all_of(c("red", "green")), names_to = "line",
                 values_to = "aphids") |>
    mutate(line = case_when(line == "red" ~ "susceptible",
                            line == "green" ~ "resistant") |>
               factor(levels = c("resistant", "susceptible"))) |>
    filter(aphids >= 0) |>
    mutate(log_aphids = log1p(aphids),
           cage = fct_recode(cage, !!!cage_lvls)) |>
    select(treatment, rep, cage, days, line, aphids,
           log_aphids, terminated)



# Cage-level alate and replaced plant data:
alate_cage_df <- exp_df |>
    mutate(n_replaced = replaced_plants |>
               str_split(",") |>
               map_int(length)) |>
    select(treatment, rep, cage, days,
           starts_with("alates_"), n_replaced) |>
    pivot_longer(starts_with("alates_total_"), names_to = "line",
                 values_to = "alates_total") |>
    mutate(line = str_remove(line, "alates_total_"),
           alates_in = ifelse(line == "red", alates_in_red, alates_in_green),
           line = case_when(line == "red" ~ "susceptible",
                            line == "green" ~ "resistant") |>
               factor(levels = c("resistant", "susceptible")),
           cage = fct_recode(cage, !!!cage_lvls)) |>
    select(treatment, rep, cage, days, line,
           alates_total, alates_in, n_replaced)


# ============================================================================*
# ============================================================================*

# Main plot ----

# ============================================================================*
# ============================================================================*


# Adjust this if you'd prefer to use log-transformed or raw scale:
aphid_y <- c("aphids", "log_aphids")[2]

wasp_mod <- max(max(wasp_cage_df$wasps),
                max(alate_cage_df$alates_in, na.rm = TRUE)) /
    max(aphid_cage_df[[aphid_y]], na.rm = TRUE)

# maximum N used to define y axis limits and items that are near the top
# of plots:
max_N <- max(aphid_cage_df[[aphid_y]]) / 0.9


#' Function to plot experiments.
#' `r` is a string representing the rep (should be `%in% aphid_cage_df$rep`)
#' `ontop` is a logical for whether to have sub-plots (by cage) on top of each
#' other, rather than side by side.
#' Defaults to `FALSE` which is used in most cases, but I set it to `TRUE`
#' for the plot used to explain the contamination in rep 11.
#'
experiment_plotter <- function(r, ontop = FALSE, show_terminations = TRUE) {
    # r = "11"
    # rm(r, acd, lcd, wcd, p)
    acd <- aphid_cage_df |>
        filter(rep == r) |>
        rename(N = !!sym(aphid_y))
    lcd <- alate_cage_df |>
        filter(rep == r, alates_in > 0) |>
        mutate(N = alates_in / wasp_mod)
    wcd <- wasp_cage_df |>
        filter(rep == r) |>
        mutate(N = wasps / wasp_mod)
    #' Y breaks and labels will differ depending on whether we're
    #' using log1p(aphids) or just aphids.
    if (aphid_y == "aphids") {
        y_breaks <- 0:2 * 2e3
        y_labs <- paste0(0:2 * 2, "k")
    } else {
        y_breaks <- c(0, log1p(4 * 10^c(1,3)))
        y_labs <- c(0, 4 * 10^c(1,3))
    }
    p <- acd |>
        ggplot(aes(days, N, color = line))
    if (ontop) {
        p <- p +
            facet_grid(rows = vars(cage)) +
            theme(panel.spacing = unit(2, "lines"))
    } else {
        p <- p +
            facet_grid(cols = vars(cage))
    }
    if (show_terminations) {
        p <- p +
            # Vertical line(s) for early termination:
            geom_vline(data = acd |> filter(terminated),
                       aes(xintercept = days),
                       linewidth = 0.5, linetype = "22", color = "gray60")
    }
    p <- p +
        geom_area(data = wcd, fill = wasp_fill, color = NA) +
        geom_hline(yintercept = 0, color = "gray70") +
        # Main abundance lines:
        geom_line() +
        scale_color_manual(values = clone_pal, guide = "none") +
        scale_y_continuous(sec.axis = sec_axis(~ . * wasp_mod,
                                               breaks = 0:2 * 40),
                           breaks = y_breaks, labels = y_labs) +
        scale_x_continuous("Days", limits = c(0, 250),
                           breaks = 0:5 * 50) +
        theme(strip.text = element_blank(),
              axis.title = element_blank(),
              axis.text.x = element_blank(),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent",
                                             color = NA)) +
        coord_cartesian(clip = FALSE, ylim = c(0, max_N))
    return(p)
}


exp_p_list <- aphid_cage_df |>
    distinct(treatment, rep) |>
    arrange(desc(treatment), rep) |>
    getElement("rep") |>
    paste() |>
    set_names() |>
    map(experiment_plotter)


if (write_plots) {
    for (i in 1:length(exp_p_list)) {
        j <- as.integer(names(exp_p_list)[i])
        fn <- paste0(here("plots/03-experiments/exper-and-sims/"),
                     sprintf("%02i-rep%02i.pdf", i, j))
        save_plot(fn, exp_p_list[[i]], 4, 1.25)
    }; rm(i, j, fn)
} else {
    wrap_plots(exp_p_list, ncol = 1)
}




#' Plot used to explain the contamination in rep 11.
if (write_plots) {
    save_plot(here("plots/03-experiments/exper-rep11-contam-explain.pdf"),
          experiment_plotter("11", TRUE), 3, 3)
} else {
    experiment_plotter("11", TRUE)
}



# ============================================================================*
# ============================================================================*

# Wasp culling effects ----

# ============================================================================*
# ============================================================================*


#' We started culling adult wasps 2x / week on 26 July 2021.
#' We reduced this to 1x / week on 23 Aug 2021.
#' This only affected reps 6, 8, and 10 because the only other rep going
#' at that time (rep 9) had already lost all its aphids by the time we
#' started culling 2x / week.


cull_df <- exp_df |>
    filter(date <= as.Date("2021-08-23"), rep %in% c(6,8,10),
           cage == "wasp", days >= 7) |>
    select(rep, date, wasps, mummies)

#' Differences in peak wasp and mummy abundances before and after we
#' started culling:
cull_df |>
    mutate(cull = date >= as.Date("2021-07-26")) |>
    group_by(cull, rep) |>
    summarize(wasps = max(wasps),
              mummies = max(mummies),
              .groups = "drop") |>
    group_by(cull) |>
    summarize(wasps = mean(wasps),
              mummies = mean(mummies))



cull_df |>
    ggplot(aes(date, wasps)) +
    geom_line() +
    geom_point() +
    facet_wrap(~ rep, ncol = 1)






