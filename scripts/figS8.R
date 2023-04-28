
source("scripts/_shared.R")
source("scripts/_shared-stable.R")


library(plot3D)


# File name for figure:
file_out <- here("plots/stationary-points.pdf")




# shared end ----


max_t <- 1e3
disp_p <- .258
delta.list <- c(50, 40,30, 2, .85, .8)
sim <- sim_experiments(clonal_lines = c(line_s, line_r),
                       alate_field_disp_p = .15,
                       extinct_N = 1e-10,
                       max_t = max_t, save_every = 1e2)
sim <- restart_experiment(sim,
                          alate_field_disp_p = disp_p,
                          new_starts = sim$all_info[[1]]
)


d <- clone_traj(sim, delta = delta.list[1], max_t = max_t)
for(i.delta in delta.list[2:length(delta.list)]) {
    d <- rbind(d, clone_traj(sim, delta = i.delta, max_t = max_t))
}
d$resistant <- d$resistant.alate + d$resistant.apterous +
    d$resistant.parasitized
d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
    d$susceptible.parasitized

d <- aggregate(cbind(resistant, susceptible, wasps) ~ time + delta, data=d,
               FUN = sum)
d <- d[d$time > 100,]
d$prop <- d$resistant/(d$susceptible + d$resistant)



stat_points_p <- function() {

    par(mfrow=c(1,3))

    dd <- d[d$delta == delta.list[1],]
    lines3D(dd$resistant, dd$susceptible, dd$wasps, colvar = dd$wasps,
            colkey = FALSE, theta = 0, phi = 0, xlab = "Resistant",
            ylab = "Susceptible", zlab = "Wasps")
    mtext("A", side = 3, adj = 0, font = 2, cex = 1.5, line = 2)
    for(i.delta in delta.list[2:length(delta.list)]) {
        dd <- d[d$delta == i.delta,]
        lines3D(dd$resistant, dd$susceptible, dd$wasps, colvar = dd$wasps,
                colkey = FALSE, add = TRUE)
    }

    dd <- d[d$delta == delta.list[1],]
    lines2D(dd$resistant, dd$susceptible, colvar = dd$wasps, colkey = FALSE,
            xlab = "Resistant", ylab = "Susceptible")
    mtext("B", side = 3, adj = 0, font = 2, cex = 1.5, line = 2)
    points2D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], add=TRUE)
    for(i.delta in delta.list[2:length(delta.list)]) {
        dd <- d[d$delta == i.delta,]
        lines2D(dd$resistant, dd$susceptible, add = TRUE, colvar = dd$wasps,
                colkey = FALSE)
        points2D(dd$resistant[nrow(dd)], dd$susceptible[nrow(dd)], add=TRUE)
    }

    dd <- d[d$delta == delta.list[1],]
    lines2D(dd$time, dd$prop, colvar = dd$wasps, colkey = FALSE, xlab = "Time",
            ylab = "Proportion resistant")
    mtext("C", side = 3, adj = 0, font = 2, cex = 1.5, line = 2)
    for(i.delta in delta.list[2:length(delta.list)]) {
        dd <- d[d$delta == i.delta,]
        lines2D(dd$time, dd$prop, add = TRUE, colvar = dd$wasps, colkey = FALSE)
    }

    invisible(NULL)
}


if (write_plots) {
    save_plot(file_out, stat_points_p, w = 10, h = 4)
} else {
    stat_points_p()
}


