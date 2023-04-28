
clone_converge <- function(sim, delta, category = "resistant", max_t = 1e4,
                           tol = 1e-10, perturb = NULL){

    new.starts <- sim$all_info[[1]]
    pick <- !is.na(new.starts$line) & (new.starts$line == category)
    new.starts[pick, "N"] <- delta * new.starts[pick, "N"]

    new.sim <- restart_experiment(sim, new_starts = new.starts, max_t = max_t,
                                  perturb = perturb)
    ns.df <- new.sim$all_info[[1]]

    S.resistant <- sum(ns.df$N[!is.na(ns.df$line) & ns.df$line == "resistant"])
    S.susceptible <- sum(ns.df$N[!is.na(ns.df$line) &
                                     ns.df$line == "susceptible"])
    test <- (S.resistant/(S.resistant + S.susceptible) > tol)
    return(test)
}


conf_bounds <- function(x, y.lower, y.upper, col="lightgray"){
    polygon(c(x, rev(x)), c(y.upper, rev(y.lower)), col=col, border=NA)
}

clone_traj <- function(sim, delta, category = "resistant", max_t = 400,
                       perturb = NULL){

    new.starts <- sim$all_info[[1]]
    pick <- !is.na(new.starts$line) & (new.starts$line == category)
    new.starts[pick, "N"] <- delta * new.starts[pick, "N"]

    new.sim <- restart_experiment(sim, new_starts = new.starts, max_t = max_t,
                                  perturb = perturb)
    new.sim <- rm_tibs(new.sim)
    d <- new.sim$aphids
    d$line.type <- paste0(d$line,".", d$type)
    d$line.type[d$line.type == "NA.mummy"] <- "mummy"
    d <- d[,-c(1,5,6)]
    d <- spread(d, "line.type", "N")
    d <- cbind(delta = delta, d, wasps = new.sim$wasps[,4])
    return(d)
}


