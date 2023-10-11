#'
#' This file analyses observational field data to describe heterogeneity in
#' parasitism.
#' See `scripts/README.md` for the table and figure number(s).
#'



source("scripts/00-shared-all.R")
source("scripts/01-field/00-field-shared.R")

# Names of files produced here:
plots_out <- list(sigma_hist = "plots/01-field/sigma-histograms.pdf",
                  obs_sd     = "plots/01-field/obs-sigma-scatter.pdf",
                  examples   = "plots/01-field/example-parasitism.pdf")


#' ############################################################################


analysis_df <- field_df |>
    filter(dissected >= 50) |>
    arrange(date) |>
    group_by(date) |>
    mutate(n_fields = n()) |>
    ungroup()

# number of dates by min. fields:
n_dates_mf <- data.frame(min_fields = c(5L, 8L))
n_dates_mf$n_dates <- sapply(n_dates_mf$min_fields, \(mf) {
    analysis_df |>
        filter(n_fields >= mf) |>
        distinct(date) |>
        nrow()
})
n_dates_mf
#   min_fields n_dates
# 1          5     100
# 2          8      29


#' ############################################################################

d <- as.data.frame(analysis_df)

field.table <- table(d$date)

glmer_cont <- glmerControl(check.conv.singular =
                               .makeCC(action = "ignore",  tol = 1e-4))


X <- data.frame(min_fields = rep(n_dates_mf$min_fields, n_dates_mf$n_dates * 3L),
                date = lapply(n_dates_mf$min_fields, \(mf) {
                    x <- analysis_df$date[analysis_df$n_fields >= mf]
                    return(rep(sort(unique(x)), 3))
                }) |>
                    do.call(what = c),
                model = c(rep(1:3, each = n_dates_mf$n_dates[1]),
                          rep(1:3, each = n_dates_mf$n_dates[2])),
                dev = NA, p.value = NA,
                obs.sd = NA, mean.dissected = NA, mean.para = NA,
                mean.para_prop = NA,
                log_aphids = NA, log_aphids.P = NA,
                height = NA, height.P = NA)

# Takes ~10 sec
set.seed(1795238152)
for (i.min_fields in n_dates_mf$min_fields) {

    dates <- as.Date(names(field.table)[field.table>=i.min_fields])

    for (i.model in 1:3) {

        for(i.date in dates){

            i.inds <- X$min_fields == i.min_fields & X$date == i.date &
                X$model == i.model

            dd <- d[d$date == i.date,]
            dd$obs <- as.factor(1:nrow(dd))

            if (i.model == 1) condition <- sum(dd$para>0) >= 2
            if (i.model == 2) {
                condition <- sum(!is.na(dd$log_aphids)) >= i.min_fields &
                    sum(dd$para>0) >= 2
            }
            if (i.model == 3) {
                condition <- sum(!is.na(dd$log_aphids)) >= i.min_fields &
                    sum(!is.na(dd$alf_height)) >= i.min_fields &
                    sum(dd$para>0) >= 2
            }

            if(condition){
                if(i.model==1){
                    z <- glmer(cbind(para, dissected - para) ~ 1 + (1|obs),
                               data=dd, family = binomial, nAGQ = 2,
                               control = glmer_cont)
                    z0 <- glm(cbind(para, dissected - para) ~ 1, data=dd,
                              family = binomial)
                }
                if(i.model==2){
                    z <- glmer(cbind(para, dissected - para) ~ log_aphids +
                                   (1|obs),
                               data=dd, family = binomial, nAGQ = 2,
                               control = glmer_cont)
                    z0 <- glm(cbind(para, dissected - para) ~ log_aphids,
                              data=dd, family = binomial)

                    X$log_aphids[i.inds] <- summary(z)[
                        "coefficients"]$coefficients[2,1]
                    X$log_aphids.P[i.inds] <- summary(z)[
                        "coefficients"]$coefficients[2,4]
                }
                if(i.model==3){
                    z <- glmer(cbind(para, dissected - para) ~ log_aphids +
                                   alf_height + (1|obs),
                               data=dd, family = binomial, nAGQ = 2,
                               control = glmer_cont)
                    z0 <- glm(cbind(para, dissected - para) ~ log_aphids +
                                  alf_height,
                              data=dd, family = binomial)

                    X$log_aphids[i.inds] <- summary(z)[
                        "coefficients"]$coefficients[2,1]
                    X$log_aphids.P[i.inds] <- summary(z)[
                        "coefficients"]$coefficients[2,4]

                    X$height[i.inds] <- summary(z)[
                        "coefficients"]$coefficients[3,1]
                    X$height.P[i.inds] <- summary(z)[
                        "coefficients"]$coefficients[3,4]
                }

                X$dev[i.inds] <- -deviance(z) + deviance(z0)
                X$p.value[i.inds] <- pchisq(-deviance(z) + deviance(z0),
                                                      df = 1, lower.tail = F)
                X$obs.sd[i.inds] <- as.numeric(VarCorr(z))^.5
                X$mean.dissected[i.inds] <- mean(dd$dissected)
                X$mean.para[i.inds] <- mean(dd$para)
                X$mean.para_prop[i.inds] <- mean(dd$para/dd$dissected)

                # show(as.numeric(VarCorr(z))^.5)
                # show(dd[,c("para","dissected")])
            }
        }
    }
}

# de-clutter environment:
rm(dd, z, z0, condition, dates, i.date, i.inds, i.min_fields, i.model)

XX <- X[!is.na(X$dev),]


summ_df <- XX |>
    as_tibble() |>
    # So that taking means below gets proportion that are significant:
    mutate(p.value = p.value < 0.05,
           log_aphids.P = log_aphids.P < 0.05,
           height.P = height.P < 0.05,
           # so it shows up the same as in table:
           min_fields = factor(min_fields, levels = c(8, 5))) |>
    group_by(min_fields, model) |>
    summarize(pos_aphids = mean(log_aphids > 0),
              pos_aphids_sig = mean(log_aphids[log_aphids.P] > 0),
        across(c(p.value, mean.dissected:height.P), mean),
              .groups = "drop") |>
    select(min_fields, model, p.value, log_aphids.P, pos_aphids, pos_aphids_sig,
           height.P, everything())

#' This is info for supplemental table:
summ_df[,1:7]
# # A tibble: 6 × 7
#   min_fields model p.value log_aphids.P pos_aphids pos_aphids_sig height.P
#   <fct>      <int>   <dbl>        <dbl>      <dbl>          <dbl>    <dbl>
# 1 8              1  0.8621      NA         NA             NA      NA
# 2 8              2  0.7857       0.1786     0.7143         1      NA
# 3 8              3  0.75         0.1786     0.7143         0.8     0.07143
# 4 5              1  0.7835      NA         NA             NA      NA
# 5 5              2  0.6495       0.1856     0.6186         0.5556 NA
# 6 5              3  0.5          0.2128     0.5638         0.55    0.1915


#' Other potentially useful info:
summ_df[,c(1:2,8:ncol(summ_df))]

# # A tibble: 6 × 7
#   min_fields model mean.dissected mean.para mean.para_prop log_aphids    height
#   <fct>      <int>          <dbl>     <dbl>          <dbl>      <dbl>     <dbl>
# 1 8              1          88.91     12.41         0.1414   NA       NA
# 2 8              2          88.77     12.75         0.1453    0.1642  NA
# 3 8              3          88.77     12.75         0.1453    0.1551  -0.02428
# 4 5              1          92.60     17.02         0.1806   NA       NA
# 5 5              2          92.60     17.02         0.1806    0.07674 NA
# 6 5              3          92.81     16.63         0.1754    0.03916 -0.006202











#' ############################################################################



sigma_histogram <- function() {
    expression.list <- c(expression(sigma),
                         expression(sigma[resid]~(aphids)),
                         expression(sigma[resid]~(aphids~and~plants)))
    par(mfrow=c(3,2), mai=c(0.6, 0.6, 0.4, 0.1))
    for(i.type in 1:3){
        for (.min_fields in c(8L, 5L)) {
            xx <- XX[XX$min_fields == .min_fields & XX$model == i.type,]
            hist(xx$obs.sd, col = "blue", breaks = .2*(0:7),
                 main = "", xlab = expression.list[i.type])
            if (i.type == 1) {
                mtext(sprintf("min. fields = %i", .min_fields),
                      side = 3, line = 1)
            }
            hist(xx$obs.sd[xx$p.value > 0.05], add = TRUE,
                 col=adjustcolor("white", alpha.f = 0.5), breaks = .2*(0:7))
        }
    }
}


if (write_plots) {
    save_plot(plots_out$sigma_hist, sigma_histogram, w = 8, h = 6)
} else {
    sigma_histogram()
}





obs_sd_scatter <- function() {
    expression.list <- c(expression(sigma),
                         expression(sigma[resid]~(aphids)),
                         expression(sigma[resid]~(aphids~and~plants)))
    # Only doing this for min. 8 fields:
    .min_fields <- 8
    xx <- XX[XX$min_fields == .min_fields,]
    # common dates for all model types:
    comm_dates <- Reduce(intersect, lapply(1:3, \(i) xx$date[xx$model == i])) |>
        as.Date()
    xx <- xx[xx$date %in% comm_dates,]
    par(mfrow=c(2,1), mai=c(1, 1, .1, .1))
    plot(xx$obs.sd[xx$model == 1], xx$obs.sd[xx$model == 2],
         xlab = expression.list[1], ylab = expression.list[2])
    lines(c(0,10),c(0,10))
    plot(xx$obs.sd[xx$model == 1], xx$obs.sd[xx$model == 3],
         xlab = expression.list[1], ylab = expression.list[3])
    lines(c(0,10),c(0,10))
}



if (write_plots) {
    save_plot(plots_out$obs_sd, obs_sd_scatter, w = 4, h = 7.5)
} else {
    obs_sd_scatter()
}




field_examples <- function() {
    # Only doing this for min. 8 fields:
    .min_fields <- 8
    xx <- XX[XX$min_fields == .min_fields & XX$model == 1,]
    xx <- xx[order(xx$obs.sd),]
    date.list <- xx$date[c(1, round(nrow(xx)/2), nrow(xx))]
    # print(date.list)
    # [1] "2014-06-02" "2015-06-12" "2017-08-15"
    par(mfrow=c(3,1), mai = c(0.4, 0.7, 0.2, 0.2))
    nrep <- 10^5
    ylim.list <- c(.15, .5, .18)
    for(i in 1:3){
        dd <- d[d$date == date.list[i],]
        dd$freq <- dd$para/dd$dissected
        dd <- dd[order(dd$freq, decreasing = F),]
        barplot(dd$freq, main="", ylim=c(0,ylim.list[i]), width = .9)
        mtext(bquote(sigma ~ "=" ~ .(round(
            xx$obs.sd[xx$date == date.list[i]], digits = 2))))
        z0 <- glm(cbind(para, dissected - para) ~ 1, data=dd, family = binomial)
        prob <- 1/(1+exp(-coef(z0)))
        w <- matrix(NA,nrep,nrow(dd))
        for(i.rep in 1:nrep) {
            w[i.rep,] <- sort(rbinom(n = nrow(dd), size=dd$dissected,
                                     prob = prob), decreasing = FALSE)
        }
        for(ii in 1:ncol(w)) w[,ii] <- sort(w[,ii])
        y <- colMeans(w)/dd$dissected
        x <- -.5 + 1.09*(1:ncol(w))
        y.upper <- w[.975*nrep,]/dd$dissected
        y.lower <- w[.025*nrep,]/dd$dissected
        points(x, y)
        # error bars:
        arrows(x0=x, x1=x, y0=y.lower, y1=y.upper, code=3, angle=90,
               0.05, col="black")
    }
    mtext(side = 1, "Rank by parasitism", padj = 2)
    mtext(side = 2, "Parasitism", padj = -3.5, adj = 3.7)
}



if (write_plots) {
    save_plot(plots_out$examples, field_examples, w = 6, h = 6)
} else {
    field_examples()
}


