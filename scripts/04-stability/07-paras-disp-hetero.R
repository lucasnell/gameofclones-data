
#'
#' Stability in relation to parasitoid dispersal heterogeneity.
#'

source("scripts/00-shared-all.R")
source("scripts/04-stability/00-stability-shared.R")

# Directory where plots produced here will be added:
plot_dir_out <- here("plots/04-stability")
if (!dir.exists(plot_dir_out) && write_plots) {
    dir.create(plot_dir_out, recursive = TRUE)
}

# Names of files produced here:
plots_out <- list(dynamics = paste0(plot_dir_out, "/example-field-dynamics.pdf"),
                  resist = paste0(plot_dir_out, "/stable-sims/field-hetero-resist.pdf"),
                  abunds = paste0(plot_dir_out, "/field-hetero-abunds.pdf"),
                  peak_trough = paste0(plot_dir_out, "/field-peaks-troughs.pdf"))
# Name of temporary results file produced here:
tmp_results <- list(df = here("data-interm/paras-disp-hetero.csv"),
                    traj = here("data-interm/paras-disp-hetero-traj.csv"))



glmer_cont <- glmerControl(check.conv.singular =
                               .makeCC(action = "ignore",  tol = 1e-4))


#' ###################################################################
# varying sd and wasp_disp_m1 ----
#' ###################################################################

n_fields <- 28L
wasp_density_0 <- rep(0.2, n_fields)

delta.list <- exp(0.5*(-10:10))

# wasp dispersal
wasp.var.list <- 0.1*(0:20)
pick.wasp.var <- 11L

set.seed(0)
#sd.wasp.field.attract <- 0.5844
sd.wasp.field.attract <- 0.5844 * wasp.var.list[pick.wasp.var]
wasp_field_attract <- sort(exp(-rnorm(0, sd = sd.wasp.field.attract, n = n_fields)))
wasp_field_attract <- wasp_field_attract/sum(wasp_field_attract)
wasp_field_attract

wasp_disp_m0 <- 0.3
#wasp_disp_m1 <- 0.34906
wasp_disp_m1 <- 0.34906 * wasp.var.list[pick.wasp.var]

harvesting.length <- 28
day.interval <- harvesting.length/n_fields
# max_t needs to be a multiple of the harvesting length
max_t <- harvesting.length*500

extinct_N <- 1e-10

min.surv <- .01
max.surv <- .04
n.events <- round(max_t/day.interval)

perturb <- data.frame(when=rep(day.interval*(1:n.events), each=3),
                      where=rep(1:n_fields, each=3),
                      who=c("resistant","susceptible","mummies"),
                      how=runif(3*n.events*n_fields, min=min.surv, max=max.surv))
# kill all mummies
perturb$how[perturb$who == "mummies"] <- 0

perturb <- perturb[1:(3*n.events),]

sim <- sim_experiments(clonal_lines = c(line_s, line_r),
						n_fields = n_fields,
						wasp_density_0 = wasp_density_0,
						wasp_disp_m0 = wasp_disp_m0,
						wasp_disp_m1 = wasp_disp_m1,
						perturb = perturb,
						extinct_N = extinct_N,
                        wasp_field_attract = wasp_field_attract,
						max_t = max_t)
sim <- rm_tibs(sim)

example_field_dynamics <- function() {

    grid.newpage()

    par(mfrow=c(2,2), mai=c(0.7, 0.8, 0.4, 0.1))

    w <- sim$aphids[sim$aphids$time > harvesting.length*400,]
    ww <- sim$wasps[sim$wasps$time > harvesting.length*400,]

    w$time0 <- w$time - min(w$time)
    ww$time0 <- ww$time - min(ww$time)

    lab_xy <- list(t = 0.95, b = 0.4, l = 0.01, r = 0.55)
    y_limits <- c(0.001, 10000)

    #' ===============
    #' Part A: field 1
    #'
    i.field <- 1
    lgl <- w$field == i.field & w$type == "apterous"
    plot(N ~ time0, data=w[lgl & w$line == "resistant",], typ="l",
         col=col_pal$r, log="y", xlab="", ylim = y_limits)
    title(sprintf("field %i", i.field), cex.main = 1, font.main = 1)
    lines(N ~ time0, data=w[lgl & w$line == "susceptible",],
          col=col_pal$s, lwd = 2)
    lines(wasps ~ time0, data=ww[ww$field == i.field,], col=col_pal$w, lwd = 2)
    text(x = min(w$time0), y = y_limits[2] / 2, labels = "susceptible", col=col_pal$s,
         adj = c(0, 1), font = 2)
    mtext("A", line = 1, adj = 0, font = 2, cex = 1.5)

    #' ===============
    #' Part B: field 14
    #'
    i.field <- n_fields %/% 2L
    lgl <- w$field == i.field & w$type == "apterous"
    # par(mai=c(.4,.4,.6,.1))
    plot(N ~ time0, data=w[lgl & w$line == "resistant",], typ="l",
         col=col_pal$r, log="y", xlab="", ylim = y_limits, lwd = 2)
    title(sprintf("field %i", i.field), cex.main = 1, font.main = 1)
    lines(N ~ time0, data=w[lgl & w$line == "susceptible",],
          col=col_pal$s, lwd = 2)
    lines(wasps ~ time0, data=ww[ww$field == i.field,], col=col_pal$w, lwd = 2)
    text(x = min(w$time0), y = y_limits[2] / 2, labels = "resistant", col=col_pal$r,
         adj = c(0, 1), font = 2)
    text(x = min(w$time0), y = y_limits[1] * 5, labels = "parasitoid", col=col_pal$w,
         adj = c(0, 0), font = 2)
    mtext("B", line = 1, adj = 0, font = 2, cex = 1.5)

    #' ===============
    #' Part C: field 28
    #'
    i.field <- n_fields
    lgl <- w$field == i.field & w$type == "apterous"
    # par(mai=c(.9,.4,.1,.1))
    plot(N ~ time0, data=w[lgl & w$line == "resistant",], typ="l", col=col_pal$r,
         log="y", xlab="Days", ylim = y_limits, lwd = 2)
    title(sprintf("field %i", i.field), cex.main = 1, font.main = 1)
    lines(N ~ time0, data=w[lgl & w$line == "susceptible",],
          col=col_pal$s, lwd = 2)
    lines(wasps ~ time0, data=ww[ww$field == i.field,], col=col_pal$w, lwd = 2)
    mtext("C", line = 1, adj = 0, font = 2, cex = 1.5)

    #' ===============
    #' Part D: simulated data histogram
    #'
    set.seed(0)

    w <- sim$aphids[sim$aphids$time > harvesting.length*400,]
    w$time0 <- w$time - min(w$time)

    h <- aggregate(N ~ time0 + field, data = w, FUN = sum)
    h$parasitized <- aggregate(N ~ time0 + field,
                               data = w[w$type == "parasitized",], FUN = sum)$N
    # head(h, 40)

    sample.times <- sort(sample(unique(h$time0),size = 100))
    h <- h[is.element(h$time0, sample.times),]
    h <- h[order(h$time0),]

    sigma <- data.frame(time = sample.times)
    counter <- 0
    for(time in sample.times){
        hh <- h[h$time0 == time,]
        hh <- hh[is.element(hh$field, sample(hh$field, 8)),]
        hh$prob <- hh$parasitized/hh$N
        hh$para <- rbinom(8, size = 100, prob = hh$prob)
        hh$obs <- 1:8
        if(all(hh$para == 0)) next
        counter <- counter + 1
        z <- glmer(cbind(para, 100 - para) ~ 1 + (1|obs), data=hh,
                   family = binomial, control = glmer_cont)
        z0 <- glm(cbind(para, 100 - para) ~ 1, data=hh, family = binomial)
        sigma$sd[counter] <- as.numeric(VarCorr(z))^.5
        sigma$p.value[counter] <- pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F)
        sigma$mean.prob[counter] <- mean(hh$prob)
        sigma$mean.N[counter] <- mean(hh$N)
    }
    sigma <- sigma[1:counter,]
    # show(c(mean(sigma$mean.prob),mean(sigma$mean.N)))
    sigma <- sigma[sigma$mean.prob > 0.009,]

    # par(mai=c(.9,.4,.1,.1))
    f <- ceiling(5*max(sigma$sd))
    hist(sigma$sd, col = "blue", breaks = .2*(0:f), main = "",
         xlab = expression(sigma))
    hist(sigma$sd[sigma$p.value > 0.05], add = T,
         col=adjustcolor("white", alpha.f = 0.5), breaks = .2*(0:f))
    mtext("D", line = 1, adj = 0, font = 2, cex = 1.5)

    invisible(NULL)
}


if (write_plots) {
    save_plot(plots_out$dynamics, example_field_dynamics, w = 8, h = 8)
} else {
    example_field_dynamics()
}



#' #####################
# supp. figure to describe peaks and troughs ----
#' #####################


# Dynamics to use for peaks and troughs:
pt_dyn_df <- sim$aphids |>
    as_tibble() |>
    filter(!is.na(line),
           # field %in% round(seq(1, 28, length.out = 3)),
           time > harvesting.length*400,
           time < harvesting.length*500) |>
    group_by(time, field, line) |>
    summarize(N = sum(N), .groups = "drop") |>
    pivot_wider(names_from = "line", values_from = "N") |>
    mutate(time = time - min(time),
           prop = resistant / (resistant + susceptible),
           field = factor(field))

pt_pts_df <- pt_dyn_df |>
    group_by(field) |>
    filter(prop %in% range(prop)) |>
    mutate(type = ifelse(prop == min(prop), "trough", "peak")) |>
    group_by(field, prop) |>
    # Make sure there's only one per field and proportion:
    filter(time == min(time)) |>
    ungroup() |>
    select(time, field, type, prop)

pt_summ_df <- pt_pts_df |>
    group_by(type) |>
    summarize(hi = max(prop), prop = min(prop), .groups = "drop") |>
    mutate(time = max(pt_dyn_df$time) * 1.05)



peak_trough_p <- pt_dyn_df |>
    ggplot(aes(time, prop)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_hline(yintercept = 1, color = "gray70") +
    geom_line(aes(color = field), linewidth = 0.75, alpha = 0.25) +
    geom_point(data = pt_pts_df, aes(color = field, shape = type),
               size = 3) +
    geom_errorbar(data = pt_summ_df, aes(ymin = prop, ymax = hi),
                  width = 40, color = safe_pals$main[2], linewidth = 1) +
    geom_text(data = pt_summ_df, aes(y = (prop + hi) / 2, label = type),
              color = safe_pals$main[2], hjust = 0, nudge_x = 50,
              fontface = "bold", size = 12 / 2.8) +
    scale_color_viridis_d(guide = "none") +
    scale_shape_manual(values = c(17, 15), guide = "none") +
    scale_y_continuous("Proportion resistant", limits = c(0, 1)) +
    scale_x_continuous("Days", limits = c(0, max(pt_dyn_df$time) * 1.125),
                       breaks = seq(0, max(pt_dyn_df$time) %/% 500 * 500, 500))


if (write_plots) {
    save_plot(plots_out$peak_trough, peak_trough_p, w = 6.5, h = 5)
} else {
    peak_trough_p
}





#' #####################
# phase portrait varying sd and wasp_disp_m1 ----
#' #####################

sim0 <- sim

if (! file.exists(tmp_results$df) | ! file.exists(tmp_results$traj)) {
    # Takes ~1.5 min
    df <- data.frame(wasp.var = rep(wasp.var.list, each = length(delta.list)),
                     delta = rep(delta.list, times = length(wasp.var.list)),
                     prop.delta = NA, converge = NA,
                     peak.resistant1 = NA, peak.resistant2 = NA,
                     peak.resistant3 = NA,
                     trough.resistant1 = NA, trough.resistant2 = NA,
                     trough.resistant3 = NA,
                     peak.susceptible1 = NA, peak.susceptible2 = NA,
                     peak.susceptible3 = NA,
                     trough.susceptible1 = NA, trough.susceptible2 = NA,
                     trough.susceptible3 = NA,
                     peak.wasps1 = NA, peak.wasps2 = NA, peak.wasps3 = NA,
                     trough.wasps1 = NA, trough.wasps2 = NA, trough.wasps3 = NA,
                     peak.prop1 = NA, peak.prop2 = NA,
                     trough.prop1 = NA, trough.prop2 = NA)

    pick.field.list <- c(1,14,28)
    traj.df <- data.frame(wasp.var = wasp.var.list)
    traj.list <- list(traj.df,traj.df,traj.df)
    for(case in 1:2){
        if(case == 1) {
            count.list <- pick.wasp.var:length(wasp.var.list)
        }else{
            sim <- sim0
            count.list <- (pick.wasp.var-1):1
        }
        #' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #' >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        #' Below is to create the version of the plot from Tony's supplemental
        #' info he emailed. The original version had the second line commented
        #' instead of the first. Waiting on Tony's reply to confirm what he
        #' recommends here.
        #'
        # new.sim <- sim # use previous iteration
        new.sim <- sim0 # use same initial starting conditions
        #'
        #' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        #' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        new.starts <- new.sim$all_info[[1]]
        for(count.wasp.var in count.list){
            i.wasp.var <- wasp.var.list[count.wasp.var]
            # show(i.wasp.var)

            #~~~~~~~~~~~~~~~~~
            set.seed(0)
            sd.wasp.field.attract <- 0.5844 * i.wasp.var
            wasp_field_attract <- sort(exp(-rnorm(0, sd = sd.wasp.field.attract,
                                                  n = n_fields)))
            wasp_field_attract <- wasp_field_attract/sum(wasp_field_attract)

            wasp_disp_m1 <- 0.34906 * i.wasp.var

            #~~~~~~~~~~~~~~~~~

            old.sim <- new.sim
            new.sim <- restart_experiment(sims_obj = old.sim,
                                          new_starts = new.starts,
                                          max_t = max_t,
                                          perturb = perturb,
                                          wasp_field_attract = wasp_field_attract,
                                          wasp_disp_m1 = wasp_disp_m1)
            new.starts <- new.sim$all_info[[1]]

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
            d <- d[d$time >= harvesting.length*400,]
            d <- d[d$time %% day.interval == 0,]

            d$resistant <- d$resistant.alate + d$resistant.apterous +
                d$resistant.parasitized
            d$susceptible <- d$susceptible.alate + d$susceptible.apterous +
                d$susceptible.parasitized

            d <- aggregate(cbind(resistant, susceptible) ~ time + field,
                           data=d, FUN = sum)
            d$prop <- d$resistant/(d$resistant + d$susceptible)

            # wasps
            dw <- new.sim$wasps
            dw <- dw[dw$time >= harvesting.length*400,]
            dw <- dw[dw$time %% day.interval == 0,]

            dw <- aggregate(wasps ~ time + field, data=dw, FUN = sum)

            d <- merge(d, dw)

            max.matrix <- matrix(NA,n_fields,4)
            min.matrix <- matrix(NA,n_fields,4)
            for(i.field in 1:n_fields){
                max.matrix[i.field,1] <- max(d$resistant[d$field == i.field])
                min.matrix[i.field,1] <- min(d$resistant[d$field == i.field])

                max.matrix[i.field,2] <- max(d$susceptible[d$field == i.field])
                min.matrix[i.field,2] <- min(d$susceptible[d$field == i.field])

                max.matrix[i.field,3] <- max(d$wasps[d$field == i.field])
                min.matrix[i.field,3] <- min(d$wasps[d$field == i.field])

                max.matrix[i.field,4] <- max(d$prop[d$field == i.field])
                min.matrix[i.field,4] <- min(d$prop[d$field == i.field])
            }

            for(i in 1:3){
                lgl <- traj.list[[i]]$wasp.var == i.wasp.var
                traj.list[[i]]$peak.resistant[lgl] = max.matrix[pick.field.list[i],1]
                traj.list[[i]]$peak.susceptible[lgl] = max.matrix[pick.field.list[i],2]
                traj.list[[i]]$peak.wasps[lgl] = max.matrix[pick.field.list[i],3]
                traj.list[[i]]$trough.resistant[lgl] = min.matrix[pick.field.list[i],1]
                traj.list[[i]]$trough.susceptible[lgl] = min.matrix[pick.field.list[i],2]
                traj.list[[i]]$trough.wasps[lgl] = min.matrix[pick.field.list[i],3]
            }
            df$peak.prop1[df$wasp.var == i.wasp.var] <- max(max.matrix[,4])
            df$peak.prop2[df$wasp.var == i.wasp.var] <- min(max.matrix[,4])
            df$trough.prop1[df$wasp.var == i.wasp.var] <- max(min.matrix[,4])
            df$trough.prop2[df$wasp.var == i.wasp.var] <- min(min.matrix[,4])

            # traj abundances for time %% harvesting.length == 0
            max.matrix <- matrix(NA,n_fields,4)
            min.matrix <- matrix(NA,n_fields,4)
            for(i.field in 1:n_fields){
                lgl <- d$field == i.field & d$time %% harvesting.length == 0
                max.matrix[i.field,1] <- max(d$resistant[lgl])
                min.matrix[i.field,1] <- min(d$resistant[lgl])

                max.matrix[i.field,2] <- max(d$susceptible[lgl])
                min.matrix[i.field,2] <- min(d$susceptible[lgl])

                max.matrix[i.field,3] <- max(d$wasps[lgl])
                min.matrix[i.field,3] <- min(d$wasps[lgl])

                max.matrix[i.field,4] <- max(d$prop[lgl])
                min.matrix[i.field,4] <- min(d$prop[lgl])
            }

            # show(cbind(min.matrix, max.matrix))

            for(i in 1:3){
                lgl <- traj.list[[i]]$wasp.var == i.wasp.var
                traj.list[[i]]$peak.resistant.int[lgl] <-
                    max.matrix[pick.field.list[i],1]
                traj.list[[i]]$peak.susceptible.int[lgl] <-
                    max.matrix[pick.field.list[i],2]
                traj.list[[i]]$peak.wasps.int[lgl] <-
                    max.matrix[pick.field.list[i],3]
                traj.list[[i]]$trough.resistant.int[lgl] <-
                    min.matrix[pick.field.list[i],1]
                traj.list[[i]]$trough.susceptible.int[lgl] <-
                    min.matrix[pick.field.list[i],2]
                traj.list[[i]]$trough.wasps.int[lgl] <-
                    min.matrix[pick.field.list[i],3]
            }

            df$peak.prop1.int[df$wasp.var == i.wasp.var] <- max(max.matrix[,4])
            df$peak.prop2.int[df$wasp.var == i.wasp.var] <- min(max.matrix[,4])
            df$trough.prop1.int[df$wasp.var == i.wasp.var] <- max(min.matrix[,4])
            df$trough.prop2.int[df$wasp.var == i.wasp.var] <- min(min.matrix[,4])

            w <- new.sim$aphids[new.sim$aphids$time %% harvesting.length == 0 &
                                    new.sim$aphids$time > harvesting.length*400,]
            ww <- new.sim$wasps[new.sim$wasps$time %% harvesting.length == 0 &
                                    new.sim$wasps$time > harvesting.length*400,]
            # par(mfrow=c(8,1), mai=c(.3,.5,.1,.1))
            # for(i.field in 1:7){
            #     lgl <- w$field == i.field & w$type == "apterous"
            # 	plot(N ~ time, data=w[lgl & w$line == "resistant",],
            # 	     typ="l", col=col_pal$r, log="y", xlab="", ylim = c(.01,10000))
            # 	lines(N ~ time, data=w[lgl & w$line == "susceptible",])
            # 	lines(wasps ~ time, data=ww[ww$field == i.field,], col=col_pal$w)
            # }

            www <- w[w$line == "susceptible" & w$type == "apterous", 2:3]
            www$susceptible <- w$N[w$line == "susceptible" & w$type == "apterous"] +
                w$N[w$line == "susceptible" & w$type == "alate"]
            www$resistant <- w$N[w$line == "resistant" & w$type == "apterous"] +
                w$N[w$line == "resistant" & w$type == "alate"]
            www$parasitized <- w$N[w$line=="susceptible" & w$type=="parasitized"] +
                w$N[w$line == "resistant" & w$type == "parasitized"]
            www$aphids <- www$susceptible + www$resistant
            www$para <- www$parasitized/(www$parasitized + www$aphids)

            sigma <- aggregate(para ~ time, data = www,
                               FUN = function(x) sd(log(x/(1-x))))
            names(sigma)[2] <- "sd"
            # if(any(!is.nan(sigma$sd))) hist(sigma$sd)
            df$mean.sd[df$wasp.var == i.wasp.var] <- mean(sigma$sd)
            df$sd.sd[df$wasp.var == i.wasp.var] <- sd(sigma$sd)
        }
    }

    write_csv(df, tmp_results$df)

    # Simplify `traj.list` to save as csv:
    traj.df <- lapply(1:length(traj.list), \(i) {
        z <- traj.list[[i]]
        z$index <- i
        return(z)
    }) |>
        do.call(what = rbind)
    write_csv(traj.df, tmp_results$traj)

} else {
    #' #############################
    df <- read.csv(tmp_results$df)
    traj.list <- read.csv(tmp_results$traj) |>
        split(~ index) |>
        lapply(\(x) {
            x$index <- NULL
            rownames(x) <- 1:nrow(x)
            return(x)
        }) |>
        unname()
}





#' ======================================================================
#' ======================================================================

# For equilibrium proportion resistance ~ parasitoid dispersal heterogeneity

field_hetero_resist <- function() {

    w <- data.frame(wasp.var = unique(df$wasp.var))
    for (i.wasp.var in w$wasp.var){

        ww <- df[df$wasp.var == i.wasp.var,]

        w$peak.prop1[w$wasp.var == i.wasp.var] <- ww$peak.prop1[1]
        w$peak.prop2[w$wasp.var == i.wasp.var] <- ww$peak.prop2[1]
        w$trough.prop1[w$wasp.var == i.wasp.var] <- ww$trough.prop1[1]
        w$trough.prop2[w$wasp.var == i.wasp.var] <- ww$trough.prop2[1]

        w$peak.prop1.int[w$wasp.var == i.wasp.var] <- ww$peak.prop1.int[1]
        w$peak.prop2.int[w$wasp.var == i.wasp.var] <- ww$peak.prop2.int[1]
        w$trough.prop1.int[w$wasp.var == i.wasp.var] <- ww$trough.prop1.int[1]
        w$trough.prop2.int[w$wasp.var == i.wasp.var] <- ww$trough.prop2.int[1]

        w$mean.sd[w$wasp.var == i.wasp.var] <- ww$mean.sd[1]

        www <- ww[ww$converge == T,]
        if(nrow(www) > 0 & nrow(www) < length(delta.list)){
            # show(www)
            if (www$prop.delta[1] < delta.list[round(length(delta.list)/2)]) {
                w$lower[w$wasp.var == i.wasp.var] <- www$prop.delta[1]
            }
            if(www$prop.delta[nrow(www)] > delta.list[round(length(delta.list)/2)]){
                w$upper[w$waspvar == i.wasp.var] <- www$prop.delta[nrow(www)]
            }
        }
    }

	par(mfrow = c(1,1), mai=c(0.5, 0.5, 0.1, 0.1))

	plot(peak.prop1 ~ wasp.var, data = w, typ="l", ylim = c(0,1),
	     # ylab = "Proportion resistant",
	     # xlab = expression(Parasitoid~dispersal~heterogeneity~gamma))
	     xlab = "", ylab = "")

	env_col <- safe_pals$main[2]

	conf_bounds(x = w$wasp.var, y.lower = (w$peak.prop2 > 0),
	            y.upper = rep(1,length(w$peak.prop2)), col="lightgray")
	conf_bounds(x = w$wasp.var, y.lower = w$peak.prop2,
	            y.upper = w$peak.prop1, col = alpha(env_col, 0.3))
	conf_bounds(x = w$wasp.var, y.lower = w$trough.prop2,
	            y.upper = w$trough.prop1, col = alpha(env_col, 0.3))

	lines(peak.prop1 ~ wasp.var, data = w, col=env_col, lwd = 2)
	lines(peak.prop2 ~ wasp.var, data = w, col=env_col, lwd = 2)
	lines(trough.prop1 ~ wasp.var, data = w, col=env_col, lwd = 2)
	lines(trough.prop2 ~ wasp.var, data = w, col=env_col, lwd = 2)

	# variation in parasitism:
	# lines(mean.sd/2 ~ wasp.var, data = w, lty = 1, col="red", lwd = 2)
}


if (write_plots) {
    save_plot(plots_out$resist, field_hetero_resist, w = 5.5, h = 3.6)
} else {
    field_hetero_resist()
}





#' ======================================================================
#' ======================================================================


# For equilibrium abundance ~ parasitoid dispersal heterogeneity


field_hetero_abunds <- function() {

	par(mfrow = c(2,3), mai=c(0.6, 0.6, 0.3, 0.05))

	for (i in 1:3) {
		ww <- traj.list[[i]]
		ww[,-1] <- log10(ww[,-1]+.0001)

		ylim = c(min(c(ww$peak.resistant, ww$peak.susceptible, ww$peak.wasps)),
		         max(c(ww$peak.resistant, ww$peak.susceptible, ww$peak.wasps)))
		plot(peak.resistant ~ wasp.var, data = ww, typ="l",
		     ylab = ifelse(i == 1, "Abundance", ""),
		     xlab = "", col=col_pal$r, yaxt = "n", ylim = ylim)
		aty <- axTicks(2)
		y_labels <- sapply(aty, function(i) as.expression(bquote(10^ .(i))))
		axis(2,at=aty,labels=y_labels, las=1)

		conf_bounds(x = ww$wasp.var, y.upper = ww$peak.resistant,
		            y.lower = ww$trough.resistant, col=alpha(col_pal$r, 0.3))
		# lines(peak.resistant ~ wasp.var, data = ww, col=col_pal$r)
		lines(trough.resistant ~ wasp.var, data = ww, col=col_pal$r)

		conf_bounds(x = ww$wasp.var, y.upper = ww$peak.susceptible,
		            y.lower = ww$trough.susceptible, col=alpha(col_pal$s, 0.3))
		lines(peak.susceptible ~ wasp.var, data = ww, col=col_pal$s)
		lines(trough.susceptible ~ wasp.var, data = ww, col=col_pal$s)

		conf_bounds(x = ww$wasp.var, y.upper = ww$peak.wasps,
		            y.lower = ww$trough.wasps, col=alpha(col_pal$w, 0.3))
		lines(peak.wasps ~ wasp.var, data = ww, col=col_pal$w)
		lines(trough.wasps ~ wasp.var, data = ww, col=col_pal$w)
		mtext(LETTERS[i], line = 0.5, adj = 0, font = 2, cex = 1)

	}
	for (i in 1:3) {
		ww <- traj.list[[i]]
		ww[,-1] <- log10(ww[,-1]+.0001)

		ylim = c(min(c(ww$peak.resistant, ww$peak.susceptible, ww$peak.wasps)),
		         max(c(ww$peak.resistant, ww$peak.susceptible, ww$peak.wasps)))
		plot(peak.resistant.int ~ wasp.var, data = ww, typ="l",
		     ylab = ifelse(i == 1, "Abundance", ""),
		     xlab = expression(Parasitoid~dispersal~heterogeneity~gamma),
		     col=col_pal$r, yaxt = "n", ylim = ylim)
		aty <- axTicks(2)
		y_labels <- sapply(aty, function(i) as.expression(bquote(10^ .(i))))
		axis(2,at=aty,labels=y_labels, las=1)

		conf_bounds(x = ww$wasp.var, y.upper = ww$peak.resistant.int,
		            y.lower = ww$trough.resistant.int, col=alpha(col_pal$r, 0.3))
		# lines(peak.resistant.int ~ wasp.var, data = ww, col=col_pal$r)
		lines(trough.resistant.int ~ wasp.var, data = ww, col=col_pal$r)

		conf_bounds(x = ww$wasp.var, y.upper = ww$peak.susceptible.int,
		            y.lower = ww$trough.susceptible.int, col=alpha(col_pal$s, 0.3))
		lines(peak.susceptible.int ~ wasp.var, data = ww, col=col_pal$s)
		lines(trough.susceptible.int ~ wasp.var, data = ww, col=col_pal$s)

		conf_bounds(x = ww$wasp.var, y.upper = ww$peak.wasps.int,
		            y.lower = ww$trough.wasps.int, col=alpha(col_pal$w, 0.3))
		lines(peak.wasps.int ~ wasp.var, data = ww, col=col_pal$w)
		lines(trough.wasps.int ~ wasp.var, data = ww, col=col_pal$w)
		mtext(LETTERS[i+3], line = 0.5, adj = 0, font = 2, cex = 1)

		if (i == 3) {
		    text(x = 0, y = 2.5, labels = "susceptible", col=col_pal$s,
		         adj = c(0, 0), font = 2)
		    text(x = 1, y = 2, labels = "resistant", col=col_pal$r,
		         adj = c(0, 0), font = 2)
		    text(x = 1, y = -2, labels = "parasitoid", col=col_pal$w,
		         adj = c(0, 1), font = 2)
		}
	}

	invisible(NULL)
}



if (write_plots) {
    save_plot(plots_out$abunds, field_hetero_abunds, w = 8, h = 6)
} else {
    field_hetero_abunds()
}




