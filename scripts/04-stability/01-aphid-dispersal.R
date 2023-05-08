
#'
#' This file makes figures that relate to how aphid dispersal
#' changes the stability and equilibrium states of the system using simulations
#' that match the experiments.
#' See `scripts/README.md` for the figure number(s).
#'


source("scripts/00-shared-all.R")

source("scripts/04-stability/00-stability-shared.R")



# Directory where plots produced here will be added:
plot_dir_out <- here("plots/04-stability/stable-sims")
if (!dir.exists(plot_dir_out) && write_plots) {
    dir.create(plot_dir_out, recursive = TRUE)
}
# Names of files produced here:
plots_out <- list(N = paste0(plot_dir_out, "/stable-sims-aphid_d-abundance.pdf"),
                  P = paste0(plot_dir_out, "/stable-sims-aphid_d-resistance.pdf"))
# Name of temporary results file produced here:
tmp_results <- here("data-interm/stable-sims-aphid_d.csv")




comm <- function(x, sim, new_starts = NULL, max_t = 1){
    if(is.null(new_starts)) new_starts <- sim$all_info[[1]]
    new_starts$N <- x
    restart.sim <- restart_experiment(sim, new_starts = new_starts,
                                      max_t = max_t)$all_info[[1]]
    return(restart.sim$N)
}
clone_eig <- function(sim){
    x <- sim$all_info[[1]]$N
    C <- jacobian(comm, x, sim=sim)
    eig <- eigen(C)
    value <- eig$value[abs(eig$values) == max(abs(eig$values))][1]
    magnitude <- Re(abs(eig$value[abs(eig$values) == max(abs(eig$values))][1]))
    return(list(value = value, magnitude = magnitude))
}


# shared end ----

max_t <- 1e4
tol <- 1e-3

# baseline
disp.list <- c(.01*(1:25))
disp.list <- c(.0001*(1:9), disp.list, .258 + .0001*(1:6), .26, .28)

pick.disp <- 15
delta.list <- exp(.25*(-20:20))
sim <- sim_experiments(clonal_lines = c(line_s, line_r),
                       alate_field_disp_p = disp.list[pick.disp],
                       extinct_N = 1e-10,
                       max_t = max_t, save_every = 1e0)
sim0 <- sim

# Re-run simulations?
rerun_sims <- !file.exists(tmp_results)



# NOTE: The starting conditions for each new disp value are based on the previous value of disp. This is necessary, because the stationary point is only locally stable. The most stable stationary point is near the middle of values of disp.list, so the 15th element was selected to start, and the range of disp.list is completed going in both directions from 15.

if (rerun_sims) {

    df <- data.frame(disp = rep(disp.list, each = length(delta.list)),
                     delta = rep(delta.list, times = length(disp.list)),
                     peak.resistant = NA, peak.susceptible = NA, peak.wasps = NA,
                     trough.resistant = NA, trough.susceptible = NA,
                     trough.wasps = NA, prop.delta = NA, converge = NA,
                     peak.prop = NA, trough.prop = NA)
    counter <- 0
    verbose <- FALSE

    for(case in 1:2){
        if(case == 1) {
            count.list <- pick.disp:length(disp.list)
        }else{
            sim <- sim0
            count.list <- (pick.disp-1):1
        }
        new.sim <- sim
        new.starts <- new.sim$all_info[[1]]
        for(count.disp in count.list){

            i.disp <- disp.list[count.disp]
            if (verbose) show(i.disp)

            old.sim <- new.sim
            new.sim <- restart_experiment(old.sim, new_starts = new.starts,
                                          alate_field_disp_p = i.disp, max_t = max_t)
            new.starts <- new.sim$all_info[[1]]

            S.resistant <- sum(new.starts$N[!is.na(new.starts$line) &
                                                new.starts$line == "resistant"])
            S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) &
                                                  new.starts$line == "susceptible"])
            test <- (S.resistant/(S.resistant + S.susceptible) > tol)
            if (!test) {
                for(i.delta in rev(delta.list)) {
                    test.starts <- old.sim$all_info[[1]]
                    pick <- (!is.na(test.starts$line)) &
                        (test.starts$line == "resistant")
                    test.starts[pick, "N"] <- i.delta * test.starts[pick, "N"]

                    test.sim <- restart_experiment(old.sim, new_starts = test.starts,
                                                   alate_field_disp_p = i.disp,
                                                   max_t = max_t)
                    test.starts <- test.sim$all_info[[1]]
                    S.resistant <- sum(test.starts$N[(!is.na(test.starts$line)) &
                                                         (test.starts$line == "resistant")])
                    S.susceptible <- sum(test.starts$N[(!is.na(test.starts$line)) &
                                                           (test.starts$line == "susceptible")])

                    if (verbose) show(c(i.disp, i.delta, S.resistant/
                                            (S.resistant + S.susceptible)))

                    new.test <- (S.resistant/(S.resistant + S.susceptible) > tol)
                    if(new.test) break
                }
                new.sim <- test.sim
                new.starts <- test.sim$all_info[[1]]
            }
            S.resistant <- sum(new.starts$N[!is.na(new.starts$line) &
                                                new.starts$line == "resistant"])
            S.susceptible <- sum(new.starts$N[!is.na(new.starts$line) &
                                                  new.starts$line == "susceptible"])

            # aphids
            d <- new.sim$aphids
            d <- d[d$type != "mummy",]
            d$line.type <- paste0(d$line,".", d$type)
            d <- d[,-c(1,5,6)]
            d <- spread(d, "line.type", "N")
            d <- d[d$time >= 9000,]

            d$resistant <- d$resistant.alate + d$resistant.apterous +
                d$resistant.parasitized
            d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
                d$susceptible.parasitized

            d <- aggregate(cbind(resistant, susceptible) ~ time, data=d, FUN = sum)
            d$prop <- d$resistant/(d$resistant + d$susceptible)

            df$peak.prop[df$disp == i.disp] <- max(d$prop)
            df$trough.prop[df$disp == i.disp] <- min(d$prop)
            df$peak.resistant[df$disp == i.disp] <- max(d$resistant)
            df$peak.susceptible[df$disp == i.disp] <- max(d$susceptible)
            df$trough.resistant[df$disp == i.disp] <- min(d$resistant)
            df$trough.susceptible[df$disp == i.disp] <- min(d$susceptible)

            # wasps
            d <- new.sim$wasps
            d <- d[d$time >= 9000,]

            d <- aggregate(wasps ~ time, data=d, FUN = sum)
            df$peak.wasps[df$disp == i.disp] <- max(d$wasps)
            df$trough.wasps[df$disp == i.disp] <- min(d$wasps)

            eig <- clone_eig(new.sim)$value
            df$eig[df$disp == i.disp] <- eig

            for(i.delta in delta.list) {
                df$converge[df$disp == i.disp & df$delta == i.delta] <-
                    clone_converge(new.sim, i.delta, max_t = max_t)
                df$prop.delta[df$disp == i.disp & df$delta == i.delta] <-
                    i.delta * sum(S.resistant) /
                    (i.delta * sum(S.resistant) + sum(S.susceptible))
            }
            if (verbose) show(df[df$disp == i.disp,])
        }
    }
    write.csv(df, tmp_results, row.names = FALSE)

} else {

    df <- read.csv(tmp_results)

}




# ============================================================================*
# pre-processing for figures ----
# ============================================================================*


w <- data.frame(disp = unique(df$disp), stable.equil1 = NA, stable.equil2 = NA, unstable.equil = NA, upper = NA)
for(i.disp in rev(w$disp)){
    ww <- df[df$disp == i.disp,]

    w$stable.equil1[w$disp == i.disp] <- ww$peak.prop[1]
    w$stable.equil2[w$disp == i.disp] <- ww$trough.prop[1]

    if(any(ww$converge[!is.na(ww$converge)])){
        www <- ww[ww$converge == TRUE,]
        w$unstable.equil[w$disp == i.disp] <- www$prop.delta[1]
        w$upper[w$disp == i.disp] <- www$prop.delta[nrow(www)]
    }
    w$peak.resistant[w$disp == i.disp] <- ww$peak.resistant[1]
    w$peak.susceptible[w$disp == i.disp] <- ww$peak.susceptible[1]
    w$peak.wasps[w$disp == i.disp] <- ww$peak.wasps[1]
    w$trough.resistant[w$disp == i.disp] <- ww$trough.resistant[1]
    w$trough.susceptible[w$disp == i.disp] <- ww$trough.susceptible[1]
    w$trough.wasps[w$disp == i.disp] <- ww$trough.wasps[1]

}
w <- w[w$disp >= 0.0005,]

w$trough.resistant <- log10(w$trough.resistant + .0001)
w$peak.resistant <- log10(w$peak.resistant + .0001)
w$trough.susceptible <- log10(w$trough.susceptible + .0001)
w$peak.susceptible <- log10(w$peak.susceptible + .0001)
w$trough.wasps <- log10(w$trough.wasps + .0001)
w$peak.wasps <- log10(w$peak.wasps + .0001)

ww <- w[!is.na(w$unstable.equil),]

ww.spline <- smooth.spline(x = ww$disp, y = ww$unstable.equil, df = 6)
ww$unstable.equil.spline <- predict(ww.spline, ww$disp)$y

www <- w[w$disp <= 0.2585,]

# This prevents line from continuing to the end for resistant aphids:
w$peak.resistant[tail(which(w$peak.resistant == log10(.0001)), -1)] <- NA
w$trough.resistant[tail(which(w$trough.resistant == log10(.0001)), -1)] <- NA




# ============================================================================*
# figures themselves ----
# ============================================================================*


# For equilibrium abundances ~ aphid dispersal

ss_aphid_d_abund <- function() {

    par(mai=c(0.9, 0.9, 0.1, 0.1))

    plot(peak.resistant ~ disp, data = w, typ="l", ylab = "Abundance",
         xlab = "Aphid dispersal", col=col_pal$r, lwd = 3, xlim = c(0,.28),
         ylim = c(-4,4), yaxt = "n")
    aty <- axTicks(2)
    y_labels <- sapply(aty,function(i) as.expression(bquote(10^ .(i))))
    axis(2,at=aty,labels=y_labels, las=1)
    lines(trough.resistant ~ disp, data = w, col=col_pal$r, lwd = 3)
    lines(peak.susceptible ~ disp, data = w, col=col_pal$s, lwd = 3)
    lines(trough.susceptible ~ disp, data = w, col=col_pal$s, lwd = 3)
    lines(peak.wasps ~ disp, data = w, col=col_pal$w, lwd = 3)
    lines(trough.wasps ~ disp, data = w, col=col_pal$w, lwd = 3)

    text(x = 0.05, y = 0, labels = "parasitoid", col=col_pal$w,
         adj = c(0, 1), font = 2)
    text(x = 0.25, y = -2, labels = "resistant", col=col_pal$r,
         adj = c(1, 0.5), font = 2)
    text(x = max(w$disp), y = 4, labels = "susceptible", col=col_pal$s,
         adj = c(1, 1), font = 2)

}

if (write_plots) {
    save_plot(plots_out$N, ss_aphid_d_abund, w = 6, h = 4)
} else {
    ss_aphid_d_abund()
}


# For equilibrium proportion resistance ~ aphid dispersal
ss_aphid_d_resist <- function() {
    # Threshold to manually prevent crossing over of lines at bifurcation
    cot <- 0.2582

    par(mai=c(0.5, 0.5, 0.1, 0.1))

    plot(stable.equil1 ~ disp, data = www[www$disp < cot,], typ="l",
         xlim = c(0,.28), ylim = c(0,1), # xaxt = "n",
         ylab = "", xlab = "")
         # ylab = "Proportion resistant", xlab = "Aphid dispersal")
    # axis(3)

    conf_bounds(x = ww$disp, y.lower = ww$upper, y.upper = rep(1,nrow(ww)),
                col="lightgray")
    conf_bounds(x = ww$disp, y.upper = ww$unstable.equil.spline,
                y.lower = rep(0,nrow(ww)), col="lightgray")
    conf_bounds(x = c(0.2584, .28), y.lower = c(0,0), y.upper = c(1,1),
                col="lightgray")

    lines(stable.equil1 ~ disp, data = www[www$disp < cot,],
          lwd = 3, col="dodgerblue")
    lines(stable.equil2 ~ disp, data = www[www$disp < cot,],
          lwd = 3, col="dodgerblue")
    points(stable.equil1 ~ disp, data = www[www$disp > .04 & www$disp <= .15,],
           col="dodgerblue")
    points(stable.equil2 ~ disp, data = www[www$disp > .04 & www$disp <= .15,],
           col="dodgerblue")
    lines(unstable.equil.spline ~ disp, data = ww, col = "black", lwd = 3, lty = "45")
}




if (write_plots) {
    save_plot(plots_out$P, ss_aphid_d_resist, w = 5.5, h = 3.6)
} else {
    ss_aphid_d_resist()
}






# ============================================================================*
# Figure showing what happens when we start with many resistant aphids ----
# This was written by LAN.
# ============================================================================*


# Make resistant aphids start at 99% of total population
line_s2 <- line_s
line_r2 <- line_r
# (all experimental simulations started with 64 total aphids)
line_s2$density_0[,1] <- line_s2$density_0[,1] |>
    (\(x){
        x / sum(x) * 64 * 0.01
    })()
line_r2$density_0[,1] <- line_r2$density_0[,1] |>
    (\(x){
        x / sum(x) * 64 * 0.99
    })()

ss_aphid_d_highr <- sim_experiments(clonal_lines = c(line_s2, line_r2),
                                   alate_field_disp_p = 0.1,
                                   extinct_N = 1e-5,
                                   max_t = max_t, save_every = 1)
para_lvls <- paste(c("no parasitism", "parasitism"), "patch")

ss_aphid_d_highr[["wasps"]] <- ss_aphid_d_highr[["wasps"]] |>
    select(-rep) |>
    mutate(field = factor(field, levels = 2:1,
                          labels = para_lvls))

ss_aphid_d_highr[["aphids"]] <- ss_aphid_d_highr[["aphids"]] |>
    mutate(field = factor(field, levels = 2:1,
                          labels = para_lvls),
           line = factor(line, levels = c("resistant", "susceptible"))) |>
    filter(type != "mummy") |>
    group_by(field, time, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    mutate(N = ifelse(N == 0, NA, N),
           N = log(N))


wasp_mod <- 5.078307  # <-- should be max(disp_wasps$wasps) / max(log1p(disp_aphids$N))
max_N <- 7.785  # <-- should be ceiling(max(log(disp_aphids$N)) * 1e3) / 1e3
y_breaks <- log(10^(0:3))
y_labs <- 10^(0:3)

ss_aphid_d_highr_p <- ss_aphid_d_highr[["aphids"]] |>
    filter(time <= 1000) |>
    ggplot(aes(time, N)) +
    geom_area(data = ss_aphid_d_highr[["wasps"]] |>
                  mutate(N = wasps / wasp_mod) |>
                  filter(time <= 1000),
              fill = wasp_fill, color = NA) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(aes(color = line))+
    geom_text(data = tibble(field = factor(para_lvls[2], levels = para_lvls),
                            time = 650, N = log(10)),
              aes(label = "wasps"), size = 9 / 2.8,
              hjust = 0.5, vjust = 0.5, color = col_pal$w) +
    geom_text(data = ss_aphid_d_highr[["aphids"]] |>
                  filter(time == 500, field == "no parasitism patch"),
              aes(label = line, color = line), size = 9 / 2.8,
              hjust = 0, vjust = 1, nudge_x = 50, nudge_y = -0.2) +
    scale_color_manual(values = clone_pal, guide = "none") +
    scale_y_continuous("Aphid abundance",
                       breaks = y_breaks,
                       labels = y_labs,
                       # limits = c(0)
                       sec.axis = sec_axis(~ . * wasp_mod,
                                           "Wasp abundance",
                                           breaks = 0:2 * 20)) +
    scale_x_continuous("Days") +
    facet_grid( ~ field, scales = "fixed") +
    theme(strip.text = element_text(size = 10)) +
    coord_cartesian(ylim = c(-0.1, NA))


if (write_plots) {
    save_plot("plots/04-stability/aphid-dispersal-high-resist.pdf", ss_aphid_d_highr_p, 6, 3)
} else ss_aphid_d_highr_p


