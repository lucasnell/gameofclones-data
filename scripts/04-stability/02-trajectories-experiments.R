
#' This file makes the following figure:
#' Example trajectories for the model parameterized for the experiment
#' in which parasitoids occur only in patch 1
#' See `scripts/README.md` for the figure number(s).


source("scripts/00-shared-all.R")
source("scripts/04-stability/00-stability-shared.R")


# File name for figure:
file_out <- here("plots/04-stability/trajectories-experiment.pdf")



clone_plot <- function(sim, delta.list, max_t = 500){

    max_t.converge = 1e4
    perturb = NULL

    d <- sim$aphids
    d$line.type <- paste0(d$line,".", d$type)
    d$line.type[d$line.type == "NA.mummy"] <- "mummy"
    d <- d[,-c(1,5,6)]
    d <- spread(d, "line.type", "N")
    d <- cbind(delta = 0, d, wasps = sim$wasps[,4], converge = TRUE)

    d <- d[d$time <= max_t,]

    for(i.delta in delta.list){
        dd <- cbind(clone_traj(sim = sim, delta = i.delta,
                               category = "resistant", max_t = max_t,
                               perturb = perturb),
                    converge = clone_converge(sim = sim, delta = i.delta,
                                              category = "resistant",
                                              max_t = max_t.converge,
                                              perturb = perturb))
        d <- rbind(d, dd)
    }

    d$resistant <- d$resistant.alate + d$resistant.apterous +
        d$resistant.parasitized
    d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
        d$susceptible.parasitized
    d$prop <- d$resistant/(d$resistant + d$susceptible)

    d$col <- "gray"
    d$col[d$converge] <- "dodgerblue"
    d$col[d$delta == 1] <- "black"

    d$lwd <- 1
    d$lwd[d$converge] <- 1.5
    d$lwd[d$delta == 1] <- 2

    par(mfrow=c(1, 3), mai = c(.6,.7,.4,.1))
    plot(prop ~ time, data=d[d$delta == 1 & d$field == 1,], typ="l", lwd = lwd,
         ylim=c(0,1), ylab = "Proportion resistant (para. patch)", xlab = "Time")
    mtext("A", side = 3, adj = 0, font = 2, cex = 1.5, line = 1)
    for (i.delta in delta.list) {
        lines(prop ~ time, data=d[d$delta == i.delta & d$field == 1,],
              col=col, lwd = lwd)
    }
    lines(prop ~ time, data=d[d$delta == 1 & d$field == 1,], lwd = lwd)

    plot(prop ~ time, data=d[d$delta == 1 & d$field == 2,], typ="l", lwd = lwd,
         ylim=c(0,1), ylab = "Proportion resistant (no-para. patch)", xlab = "Time")
    mtext("B", side = 3, adj = 0, font = 2, cex = 1.5, line = 1)
    for (i.delta in delta.list) {
        lines(prop ~ time, data=d[d$delta == i.delta & d$field == 2,],
              col=col, lwd = lwd)
    }
    lines(prop ~ time, data=d[d$delta == 1 & d$field == 2,], lwd = lwd)

    plot(wasps ~ time, data=d[d$delta == 1 & d$field == 1,], typ="l", lwd = lwd,
         ylim = c(0, max(d$wasps)), ylab = "Parasitoid abundance (para. patch)", xlab = "Time")
    mtext("C", side = 3, adj = 0, font = 2, cex = 1.5, line = 1)
    for(i.delta in delta.list) {
        lines(wasps ~ time, data=d[d$delta == i.delta & d$field == 1,],
              col=col, lwd = lwd)
    }
    lines(wasps ~ time, data=d[d$delta == 1 & d$field == 1,], lwd = lwd)
}


# shared end ----


# domain of attraction
max_t <- 1e4
tol <- 1e-3

# baseline
disp.list <- c(.01*(1:25))
disp.list <- c(.0001*(1:9), disp.list, .258 + .0001*(1:6), .26, .28)

pick.disp <- 25
delta.list <- exp(.25*(-20:20))
sim <- sim_experiments(clonal_lines = c(line_s, line_r),
                       alate_field_disp_p = disp.list[pick.disp],
                       extinct_N = 1e-10,
                       max_t = max_t, save_every = 1e0)

traj_exp_p <- function() clone_plot(sim, delta.list)

if (write_plots) {
    save_plot(file_out, traj_exp_p, w = 8, h = 3)
} else {
    traj_exp_p()
}

