
#'
#' Stability in relation to parasitoid dispersal heterogeneity, but
#' in this file, aphid-density-dependent parasitoid dispersal is set to zero.
#'

source("scripts/00-shared-all.R")
source("scripts/04-stability/00-stability-shared.R")

# Directory where plots produced here will be added:
plot_dir_out <- here("plots/04-stability")
if (!dir.exists(plot_dir_out) && write_plots) {
    dir.create(plot_dir_out, recursive = TRUE)
}

# Name of file produced here:
plot_out <- here("plots/04-stability/field-hetero-abunds-no-aphid.pdf")

# Name of temporary results file produced here:
tmp_results <- list(df = here("data-interm/paras-disp-hetero-no-aphid.csv"),
                    traj = here("data-interm/paras-disp-hetero-no-aphid-traj.csv"))



glmer_cont <- glmerControl(check.conv.singular =
                               .makeCC(action = "ignore",  tol = 1e-4))


#' ###################################################################
# varying sd only ----
#' ###################################################################

n_fields <- 28
wasp_density_0 <- rep(0.2, n_fields)

delta.list <- exp(0.5*(-10:10))

# wasp dispersal
wasp.var.list <- 0.1*(0:20)
# <<<< Below was changed by LAN bc using `11L` causes >>>>
# <<<< resistant aphids and parasitoids to go extinct >>>>
pick.wasp.var <- 12L

set.seed(0)
#sd.wasp.field.attract <- 0.5844
sd.wasp.field.attract <- 0.5844 * wasp.var.list[pick.wasp.var]
wasp_field_attract <- sort(exp(-rnorm(0, sd = sd.wasp.field.attract, n = n_fields)))
wasp_field_attract <- wasp_field_attract/sum(wasp_field_attract)
# wasp_field_attract

wasp_disp_m0 <- 0.3
#wasp_disp_m1 <- 0.34906
wasp_disp_m1 <- 0


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

sim0 <- sim



#' #####################
# phase portrait varying sd only
#' #####################


if (! file.exists(tmp_results$df) | ! file.exists(tmp_results$traj)) {

    # Takes ~1.25 min

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

            wasp_disp_m1 <- 0

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

            d$resistant <- d$resistant.alate + d$resistant.apterous + d$resistant.parasitized
            d$susceptible <- d$susceptible.alate + d$susceptible.apterous + d$susceptible.parasitized

            d <- aggregate(cbind(resistant, susceptible) ~ time + field, data=d, FUN = sum)
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
                traj.list[[i]]$peak.resistant[lgl] <- max.matrix[pick.field.list[i],1]
                traj.list[[i]]$peak.susceptible[lgl] <- max.matrix[pick.field.list[i],2]
                traj.list[[i]]$peak.wasps[lgl] <- max.matrix[pick.field.list[i],3]
                traj.list[[i]]$trough.resistant[lgl] <- min.matrix[pick.field.list[i],1]
                traj.list[[i]]$trough.susceptible[lgl] <- min.matrix[pick.field.list[i],2]
                traj.list[[i]]$trough.wasps[lgl] <- min.matrix[pick.field.list[i],3]
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
                traj.list[[i]]$peak.resistant.int[lgl] <- max.matrix[pick.field.list[i],1]
                traj.list[[i]]$peak.susceptible.int[lgl] <- max.matrix[pick.field.list[i],2]
                traj.list[[i]]$peak.wasps.int[lgl] <- max.matrix[pick.field.list[i],3]
                traj.list[[i]]$trough.resistant.int[lgl] <- min.matrix[pick.field.list[i],1]
                traj.list[[i]]$trough.susceptible.int[lgl] <- min.matrix[pick.field.list[i],2]
                traj.list[[i]]$trough.wasps.int[lgl] <- min.matrix[pick.field.list[i],3]
            }

            df$peak.prop1.int[df$wasp.var == i.wasp.var] <- max(max.matrix[,4])
            df$peak.prop2.int[df$wasp.var == i.wasp.var] <- min(max.matrix[,4])
            df$trough.prop1.int[df$wasp.var == i.wasp.var] <- max(min.matrix[,4])
            df$trough.prop2.int[df$wasp.var == i.wasp.var] <- min(min.matrix[,4])

            # for(i.delta in delta.list) {
            #     lgl <- df$wasp.field.attract == i.wasp.var & df$delta == i.delta
            #     df$converge[lgl] <- clone_wasp_converge(new.sim, i.delta,
            #                                             max_t = max_t,
            #                                             perturb = perturb,
            #                                             tol = extinct_N,
            #                                             day.interval = day.interval)
            #     df$prop.delta[lgl] <- i.delta*sum(S.resistant)/
            #         (i.delta*sum(S.resistant) + sum(S.susceptible))
            # }
            # show(df[df$wasp.field.attract == i.wasp.var,])

            # par(mfrow=c(8,1), mai=c(.3,.5,.1,.1))
            lgl <- new.sim$aphids$time %% harvesting.length == 0 &
                new.sim$aphids$time > harvesting.length*400
            w <- new.sim$aphids[lgl,]
            lgl <- new.sim$wasps$time %% harvesting.length == 0 &
                new.sim$wasps$time > harvesting.length*400
            ww <- new.sim$wasps[lgl,]

            www <- w[w$line == "susceptible" & w$type == "apterous", 2:3]
            www$susceptible <- w$N[w$line == "susceptible" & w$type == "apterous"] +
                w$N[w$line == "susceptible" & w$type == "alate"]
            www$resistant <- w$N[w$line == "resistant" & w$type == "apterous"] +
                w$N[w$line == "resistant" & w$type == "alate"]
            www$parasitized <- w$N[w$line == "susceptible" & w$type == "parasitized"] +
                w$N[w$line == "resistant" & w$type == "parasitized"]
            www$aphids <- www$susceptible + www$resistant
            www$para <- www$parasitized/(www$parasitized + www$aphids)

            sigma <- aggregate(para ~ time, data = www,
                               FUN = function(x) sd(log(x/(1-x))))
            names(sigma)[2] <- "sd"
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

#' ======================================================================
#' ======================================================================


# For equilibrium abundance ~ parasitoid dispersal heterogeneity (no aphid)


field_hetero_abunds_aphids0 <- function() {

    par(mfrow = c(2,3), mai=c(0.6, 0.6, 0.3, 0.05))

	ylim <- c(-4, 4)

    for(i in 1:3){
		ww <- traj.list[[i]]
		ww[,-1] <- log10(ww[,-1]+.0001)

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
	for(i in 1:3){
		ww <- traj.list[[i]]
		ww[,-1] <- log10(ww[,-1]+.0001)

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
		    text(x = 2, y = 2, labels = "resistant", col=col_pal$r,
		         adj = c(1, 0), font = 2)
		    text(x = 1.35, y = -3.9, labels = "parasitoid", col=col_pal$w,
		         adj = c(0, 0), font = 2)
		}
	}

	invisible(NULL)
}



if (write_plots) {
    save_plot(plot_out, field_hetero_abunds_aphids0, w = 8, h = 6)
} else {
    field_hetero_abunds_aphids0()
}


