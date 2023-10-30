
#'
#' Main experiment figures, plus supplemental figure for contamination in
#' rep 11 (trial "vii") and for dispersal pool sizes.
#'


source("scripts/00-shared-all.R")

# Names of files produced here:
plots_out <- list(main = c(9L, 12L, 13L, 6L, 8L, 10L, 11L) |>
                      imap_chr(\(j, i) {
                          paste0(here("plots/03-experiments/exper-and-sims/"),
                                 sprintf("%02i-rep%02i.pdf", i, j))
                      }) |>
                      set_names(c(9L, 12L, 13L, 6L, 8L, 10L, 11L)),
                  contam = here("plots/03-experiments/exper-rep11-contam-explain.pdf"),
                  dpool = here("plots/03-experiments/dispersal-pools.pdf"))


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
               map_int(~ length(.x[.x != "" & !is.na(.x)]))) |>
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
        # geom_area(data = wcd, fill = wasp_fill, color = NA) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(data = wcd, color = wasp_color) +
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

stopifnot(all(names(exp_p_list %in% names(plots_out$main))))

if (write_plots) {
    for (n in names(exp_p_list)) {
        save_plot(plots_out$main[[n]], exp_p_list[[n]], 4, 1.25)
    }; rm(i, j, fn)
} else {
    wrap_plots(exp_p_list, ncol = 1)
}




#' Plot used to explain the contamination in rep 11.
if (write_plots) {
    save_plot(plots_out$contam, experiment_plotter("11", TRUE), 3, 3)
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


here("data-raw/experiments-main.csv") |>
    read_csv(col_types = cols()) |>
    mutate(date = as.Date(date, format = "%m/%d/%Y"),
           start_date = as.Date(start_date, format = "%m/%d/%Y"),
           days = difftime(date, start_date, units = "days") |>
               as.integer(),
           cage = cage |>
               tolower() |>
               factor(levels = c("no wasp", "wasp")),
           cage = fct_recode(cage, !!!cage_lvls),
           rep = factor(rep, levels = sort(unique(rep)))) |>
    select(rep, cage, date, days, notes) |>
    filter(!is.na(notes)) |>
    filter(str_detect(notes, "wasp sex")) |>
    mutate(notes = notes |>
               str_remove(".*<wasp sex>") |>
               str_remove("</>.*") |>
               str_split(", ") |>
               map(\(x) {
                   males  <- x[str_detect(x, " male$")] |>
                       str_remove(" male$") |>
                       as.integer()
                   females  <- x[str_detect(x, " female$")] |>
                       str_remove(" female$") |>
                       as.integer()
                   fm <- c(females, males)
                   return(fm)
               }),
           female_ratio = notes |>
               map_dbl(\(x) {
                   x[1] / sum(x)
               }),
           n_sexed = notes |> map_int(sum)) |>
    select(-notes) |>
    # filter(n_sexed > 5) |>
    arrange(rep, days)


here("data-raw/experiments-pesky_wasps.csv") |>
    read_csv(col_types = cols()) |>
    select(rep, date, mummies, starts_with("adult")) |>
    #' On 3/16 and 3/22 I removed wasps/mummies but didn't record the number.
    filter(!is.na(mummies)) |>
    rename(females = `adult females`,
           males = `adult males`,
           unknown = `adults unk.`) |>
    mutate(date = as.Date(date, "%d-%b-%y"),
           rep = factor(rep, levels = levels(exp_df$rep))) |>
    group_by(rep, date) |>
    summarize_all(sum) |>
    ungroup() |>
    mutate(n_sexed = females + males,
           female_ratio = females / n_sexed) |>
    filter(n_sexed > 0) |>
    select(rep, date, female_ratio, n_sexed) |>
    # These are changed to make it join easier with the other wasp data:
    mutate(cage = factor("no parasitism cage", levels = names(cage_lvls)),
           start_date = map_chr(rep, ~ filter(exp_df, rep == .x)[["date"]] |>
                                    min() |> paste()) |> as.Date(),
           days = difftime(date, start_date, units = "days") |>
               as.integer()) |>
    select(rep, date, days, cage, female_ratio, n_sexed) |>
    filter(date <= max(exp_df$date))




# ============================================================================*
# ============================================================================*

# Dispersal pool sizes ----

# ============================================================================*
# ============================================================================*

dpool_exp_df <- left_join(alate_cage_df,
                       select(aphid_cage_df, -terminated, -log_aphids),
          by = c("treatment", "rep", "cage", "days", "line")) |>
    filter(treatment == "dispersal") |>
    mutate(rep_rmn = factor(rep, levels = c(6, 8, 10, 11),
                              labels = c("iv", "v", "vi", "vii"))) |>
    group_by(rep_rmn, days, line) |>
    summarize(dpool = sum(alates_total),
              aphids = sum(aphids), .groups = "drop")


#' Same as default simulations with dispersal, except only showing adult
#' aphid numbers:
dpool_sim_df <- sim_experiments(clonal_lines = c(line_s, line_r),
                              sep_adults = TRUE) |>
    getElement("aphids") |>
    select(-rep, -plant) |>
    mutate(field = factor(field, levels = 2:1,
                          labels = para_lvls),
           line = factor(line, levels = c("resistant", "susceptible"))) |>
    filter(!is.na(line)) |>
    group_by(time, line) |>
    summarize(aphids = sum(N),
              dpool = sum(N[type == "adult alate"]) * 0.1,
              .groups = "drop")




#' Include both experiments and simulation dispersal pools, where the latter
#' is corrected for being done daily by summing across all days between
#' consecutive samples for the associated experimental rep.
dpool_df <- dpool_exp_df |>
    arrange(rep_rmn, line, days) |>
    split(~ rep_rmn + line) |>
    map_dfr(\(dd) {
        get_sim_dp <- function(start_days, end_days, line) {
            lgl <- dpool_sim_df[["time"]] >= start_days &
                dpool_sim_df[["time"]] <= end_days &
                dpool_sim_df[["line"]] == line
            dp <- sum(dpool_sim_df[["dpool"]][lgl])
            return(dp)
        }
        end_days <- dd$days
        start_days <- c(0, head(end_days, -1) + 1)
        .line <- dd$line[[1]]
        dd[["dpool_sim"]] <- pmap_dbl(list(start_days, end_days),
                                      get_sim_dp, line = .line)
        dd[["aphids_sim"]] <- map_dbl(dd$days, \(d) {
            .lgl <- dpool_sim_df$time == d & dpool_sim_df$line == .line
            return(dpool_sim_df$aphids[.lgl])
        })
        return(dd)
    })




dpool_p <- dpool_df |>
    select(rep_rmn, days, line, starts_with("dpool")) |>
    pivot_longer(starts_with("dpool"), names_to = "source", values_to = "dpool") |>
    mutate(source = factor(source, levels = c("dpool", "dpool_sim"),
                           labels = c("Experiments", "Simulations"))) |>
    ggplot(aes(days, dpool)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(aes(color = line, linetype = source), linewidth = 0.75) +
    geom_text(data = dpool_df |>
                  distinct(rep_rmn) |>
                  mutate(days = 0, dpool = max(dpool_df$dpool, na.rm = TRUE)),
              aes(label = rep_rmn), hjust = 0, vjust = 1) +
    scale_color_manual(NULL, values = clone_pal) +
    scale_linetype_manual(NULL, values = c(1, 2)) +
    scale_y_continuous("Aphid dispersal pool", trans = "log1p",
                       breaks = c(0, 5^(1:3))) +
    scale_x_continuous("Days", limits = c(0, 250),
                       breaks = 0:5 * 50) +
    facet_grid(rep_rmn ~ .) +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          strip.text.y = element_blank(),
          legend.key.width = unit(2, "lines"))


if (write_plots) {
    save_plot(plots_out$dpool, dpool_p, 6.5, 4)
} else {
    dpool_p
}




dpool_df |>
    group_by(rep_rmn, days) |>
    summarize(aphids = sum(aphids),
              dpool = sum(dpool),
              dpool_sim = sum(dpool_sim),
              aphids_sim = sum(aphids_sim),
              .groups = "drop") |>
    filter(aphids > 0, aphids_sim > 0) |>
    mutate(p_disp = dpool / aphids,
           p_disp_sim = dpool_sim / aphids_sim) |>
    filter(days > 20) |>
    # group_by(line) |>
    summarize(p_disp = mean(p_disp),
              p_disp_sim = mean(p_disp_sim))

# # A tibble: 1 × 2
#   p_disp p_disp_sim
#    <dbl>      <dbl>
# 1 0.0118     0.00923


