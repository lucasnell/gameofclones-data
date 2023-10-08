
#'
#' These simulations are to simulate experiments with demographic stochasticity,
#' mortality of inoculations, and wasp harassment.
#'


source("scripts/00-shared-all.R")


# To maintain the same y-axis limits:
wasp_mod <- 5.078307
max_N <- 7.9
y_breaks <- log(10^(0:3))
y_labs <- 10^(0:3)



.initial_density <- 10
.aphid_demog_error <- TRUE
.wasp_demog_error <- TRUE
.wasp_badger_n <- 0.3
.n_reps <- 1000


# Susceptible line: no resistance, high population growth rate
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,.initial_density), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
# Resistant line: high resistance, low population growth rate
line_r <- clonal_line("resistant",
                      density_0 = cbind(c(0,0,0,0,.initial_density), rep(0, 5)),
                      resistant = TRUE, surv_paras = 0.57,
                      surv_juv_apterous = "low",
                      repro_apterous = "low",
                      surv_adult_apterous = "low")




set.seed(1740330095)
stoch_sims_nodisp <- sim_stochastic(c(line_s, line_r),
                                    alate_field_disp_p = 0,
                                    wasp_demog_error = .wasp_demog_error,
                                    aphid_demog_error = .aphid_demog_error,
                                    wasp_badger_n = .wasp_badger_n,
                                    n_reps = .n_reps,
                                    n_threads = .n_threads)

set.seed(159264885)
stoch_sims_disp <- sim_stochastic(c(line_s, line_r),
                                  wasp_demog_error = .wasp_demog_error,
                                  aphid_demog_error = .aphid_demog_error,
                                  wasp_badger_n = .wasp_badger_n,
                                  n_reps = .n_reps,
                                  n_threads = .n_threads)


stoch_time_series_plotter <- function(stoch_sims, .alpha = 0.1) {

    plot_reps <- sample.int(.n_reps, 100)

    stoch_sims$aphids <- stoch_sims$aphids |> filter(rep %in% plot_reps)
    stoch_sims$wasps <- stoch_sims$wasps |> filter(rep %in% plot_reps)

    mad <- stoch_sims$aphids |>
        filter(!is.na(line)) |>
        group_by(rep, field, time, line) |>
        summarize(N = sum(N), .groups = "drop") |>
        mutate(rep = factor(rep),
               field = factor(field, levels = 2:1,
                              labels = paste(c("no parasitism", "parasitism"),
                                             "patch")),
               logN = log(ifelse(N == 0, NA, N)),
               grp = interaction(rep, line))
    wad <- stoch_sims$wasps |>
        filter(field == 1) |>
        mutate(rep = factor(rep),
               field = factor(field, levels = 2:1,
                              labels = paste(c("no parasitism", "parasitism"),
                                             "patch")),
               logN = wasps / wasp_mod)

    mad |>
        ggplot(aes(time, logN)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(data = wad, aes(group = rep), color = wasp_color,
                  alpha = .alpha, linewidth = 0.5) +
        geom_line(aes(color = line, group = grp), na.rm = TRUE,
                  alpha = .alpha, linewidth = 0.5) +
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
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent", color = NA)) +
        coord_cartesian(clip = FALSE, ylim = c(0, max_N))
}


set.seed(2127406169)
stoch_time_series <- map(list(stoch_sims_nodisp, stoch_sims_disp),
                        stoch_time_series_plotter)
names(stoch_time_series) <- c("no_dispersal", "dispersal")


if (write_plots) {
    for (i in 1:length(stoch_time_series)) {
        fn <- paste0(here("plots/03-experiments/exper-and-sims/"),
                     sprintf("stoch-sims-%i-%s.pdf", i, names(stoch_time_series)[i]))
        save_plot(fn, stoch_time_series[[i]], 4, 1.25)
    }; rm(i, fn)
} else wrap_plots(stoch_time_series, nrow = 1)


stoch_sims_nodisp$wasps |>
    filter(time == max(time), field == 1) |>
    summarize(p_extinct = mean(wasps == 0))

stoch_sims_disp$wasps |>
    filter(time == max(time), field == 1) |>
    summarize(p_extinct = mean(wasps == 0))

stoch_sims_nodisp$aphids |>
    filter(type != "mummy") |>
    mutate(field = factor(field, levels = 2:1,
                              labels = paste(c("no parasitism", "parasitism"),
                                             "patch"))) |>
    group_by(field, rep, line) |>
    filter(time == max(time)) |>
    summarize(N = sum(N), .groups = "drop") |>
    pivot_wider(names_from = line, values_from = N) |>
    mutate(outcome = case_when(susceptible == 0 & resistant > 0 ~ "only_res",
                               susceptible > 0 & resistant == 0 ~ "only_sus",
                               susceptible > 0 & resistant > 0 ~ "both",
                               susceptible == 0 & resistant == 0 ~ "none",
                               TRUE ~ NA_character_)) |>
    group_by(field) |>
    summarize(only_res = mean(outcome == "only_res"),
              only_sus = mean(outcome == "only_sus"),
              both = mean(outcome == "both"),
              none = mean(outcome == "none"),
              ext_res = only_sus + none,
              ext_sus = only_res + none)

#'
#' Outcomes for no dispersal:
#'
safe_mclapply(1:.n_reps, \(i) {
    ad <- stoch_sims_nodisp$aphids[stoch_sims_nodisp$aphids$rep == i,] |>
        filter(time == max(time)) |>
        filter(type != "mummy") |>
        group_by(field, line) |>
        summarize(N = sum(N), .groups = "drop") |>
        pivot_wider(names_from = line, values_from = N) |>
        mutate(outcome = case_when(susceptible == 0 & resistant > 0 ~ "res",
                                   susceptible > 0 & resistant == 0 ~ "sus",
                                   susceptible > 0 & resistant > 0 ~ "both",
                                   susceptible == 0 & resistant == 0 ~ "none",
                                   TRUE ~ NA_character_)) |>
        select(field, outcome)
    wd <- stoch_sims_nodisp$wasps[stoch_sims_nodisp$wasps$rep == i,] |>
        filter(time == max(time), field == 1) |>
        getElement("wasps")
    if (wd > 0) {
        ad$outcome[ad$field == 1] <- paste0(ad$outcome[ad$field == 1], "_wasps")
    }
    return(ad)
}) |>
    do.call(what = bind_rows) |>
    mutate(field = factor(field, levels = 2:1,
                          labels = paste(c("no parasitism", "parasitism"),
                                         "patch"))) |>
    split(~ field) |>
    map(\(x) {
        table(x$outcome) / .n_reps
    })
# $`no parasitism patch`
#
# both  sus
# 0.76 0.24
#
# $`parasitism patch`
#
# both none  res
# 0.18 0.16 0.66



#'
#' Outcomes for dispersal:
#'
safe_mclapply(1:.n_reps, \(i) {
    ad <- stoch_sims_disp$aphids[stoch_sims_disp$aphids$rep == i,] |>
        filter(time == max(time)) |>
        filter(type != "mummy") |>
        group_by(field, line) |>
        summarize(N = sum(N), .groups = "drop") |>
        pivot_wider(names_from = line, values_from = N) |>
        mutate(outcome = case_when(susceptible == 0 & resistant > 0 ~ "res",
                                   susceptible > 0 & resistant == 0 ~ "sus",
                                   susceptible > 0 & resistant > 0 ~ "both",
                                   susceptible == 0 & resistant == 0 ~ "none",
                                   TRUE ~ NA_character_)) |>
        select(field, outcome)
    wd <- stoch_sims_disp$wasps[stoch_sims_disp$wasps$rep == i,] |>
        filter(time == max(time), field == 1) |>
        getElement("wasps")
    if (wd > 0) {
        ad$outcome[ad$field == 1] <- paste0(ad$outcome[ad$field == 1], "_wasps")
    }
    return(ad)
}) |>
    do.call(what = bind_rows) |>
    mutate(field = factor(field, levels = 2:1,
                          labels = paste(c("no parasitism", "parasitism"),
                                         "patch"))) |>
    split(~ field) |>
    map(\(x) {
        round(table(x$outcome) / .n_reps, 2)
    })

# $`no parasitism patch`
#
# both  sus
# 0.84 0.16
#
# $`parasitism patch`
#
# both both_wasps        sus  sus_wasps
# 0.11       0.51       0.00       0.38








#' In no-dispersal parasitoid cage, how often do resistant (and susceptible)
#' aphids go extinct?
nodisp_paras_outs <- safe_mclapply(1:.n_reps, \(i) {
    # This function returns abundance of resistant, then susceptible aphids
    # in the parasitism patch:
    dd <- stoch_sims_nodisp$aphids |>
        filter(time == max(time), rep == i, type != "mummy", field == 1)
    if (nrow(dd) == 0) return(c(0, 0))
    dd |>
        group_by(line) |>
        summarize(N = sum(N), .groups = "drop") |>
        getElement("N")
}) |>
    do.call(what = rbind)
mean(nodisp_paras_outs[,1] == 0 & nodisp_paras_outs[,2] == 0)
# [1] 0.16


stoch_hist_plotter <- function(stoch_sims) {
    # stoch_sims = stoch_sims_nodisp
    # rm(stoch_sims, w_df, p_df)
    w_df <- stoch_sims$wasps |>
        filter(time == max(time)) |>
        group_by(rep) |>
        summarize(wasps = sum(wasps), .groups = "drop")
    p_df <- stoch_sims$aphids |>
        filter(time == max(time)) |>
        filter(!is.na(line)) |>
        group_by(rep, line) |>
        summarize(N = sum(N), .groups = "drop") |>
        pivot_wider(names_from = line, values_from = N) |>
        mutate(rep = factor(rep),
               p = resistant / (resistant + susceptible),
               wasps = map_dbl(rep, \(r) w_df$wasps[w_df$rep == r]),
               has_wasps = factor(wasps > 0, levels = c(FALSE, TRUE),
                                  labels = c("wasps extinct", "wasps present")))
    p_df |>
        ggplot(aes(p, fill = has_wasps)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_histogram(binwidth = 0.025) +
        scale_fill_manual(values = c("#0077BB", "#EE7733"), guide = "none") +
        scale_y_continuous(breaks = c(0, 200, 400), labels = c(0, 0.2, 0.4)) +
        scale_x_continuous(breaks = 0:5 / 10)  +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent", color = NA)) +
        coord_cartesian(xlim = c(-0.02, 0.52),
                        ylim = c(-5, 540), clip = "off", expand = FALSE)
}

stoch_hists <- map(list(stoch_sims_nodisp, stoch_sims_disp), stoch_hist_plotter)
names(stoch_hists) <- c("no_dispersal", "dispersal")

if (write_plots) {
    for (i in 1:length(stoch_hists)) {
        fn <- paste0(here("plots/03-experiments/exper-and-sims/"),
                     sprintf("stoch-sims-hist-%i-%s.pdf", i, names(stoch_hists)[i]))
        save_plot(fn, stoch_hists[[i]], 3.68, 0.75)
    }; rm(i, fn)
} else wrap_plots(stoch_hists, nrow = 1)






# stoch_sims = stoch_sims_disp
# rm(stoch_sims, w_df, p_df)
w_df <- stoch_sims$wasps |>
    filter(time == max(time)) |>
    group_by(rep) |>
    summarize(wasps = sum(wasps), .groups = "drop")
p_df <- stoch_sims$aphids |>
    filter(time == max(time)) |>
    filter(!is.na(line)) |>
    group_by(rep, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    mutate(rep = factor(rep))
p_df |>
    ggplot(aes(log1p(N), color = line, fill = line)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_freqpoly(aes(color = line), binwidth = 0.5) +
    # geom_histogram(aes(fill = line), binwidth = 0.5,
    #                position = "dodge",
    #                color = NA) +
    scale_fill_manual(values = clone_pal, guide = "none") +
    scale_color_manual(values = clone_pal, guide = "none") +
    # scale_y_continuous(breaks = c(0, 200, 400), labels = c(0, 0.2, 0.4)) +
    # scale_x_continuous(breaks = 0:5 / 10)  +
    # theme(axis.title = element_blank(),
    #       axis.text.x = element_blank(),
    #       panel.background = element_rect(fill = "transparent"),
    #       plot.background = element_rect(fill = "transparent", color = NA)) +
    # coord_cartesian(xlim = c(-0.02, 0.52),
    #                 ylim = c(-5, 540), clip = "off", expand = FALSE) +
    NULL