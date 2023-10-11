
#'
#' This file analyzes experimental field data to describe heterogeneity in
#' parasitoid abundance.
#' See `scripts/README.md` for the table and figure number(s).
#'



source("scripts/00-shared-all.R")
source("scripts/01-field/00-field-shared.R")


# Names of file produced here:
plot_out <- here("plots/01-field/sigma-histograms-disturbance.pdf")



d <- read.csv(here("data-raw/disturbance.csv"))
# colnames(d)

# For consistency with monitoring analysis:
d$log_aphids <- d$logAphids
d$logAphids <- NULL
d$para <- d$total_para
d$total_para <- NULL


# select days with at least 9 fields sampled (i.e., complete experiment)
field.table <- table(d$Day_of_sampling)

min.fields <- 9
days.9 <- names(field.table)[field.table>=min.fields]
# days.9

d <- d[d$Day_of_sampling %in% days.9,]
d$obs <- 1:nrow(d)

glmer_cont <- glmerControl(optimizer = "bobyqa",
                           check.conv.singular = .makeCC(action = "ignore",
                                                         tol = 1e-4))


#' ========================================================================
#' ========================================================================


# analysis of wasp dispersion using log_aphids
# wasps excluding trt
# you can ignore warning; doesn't affect outcome
z <- glmer(wasp ~ 1 + as.factor(Day_of_sampling) + log_aphids + (1|obs) +
               offset(log(sweeps.wasps)), data=d, family = poisson)
z0 <- glm(wasp ~ 1 + as.factor(Day_of_sampling) + log_aphids +
              offset(log(sweeps.wasps)), data=d, family = poisson)
# summary(z)
# c(deviance(z), deviance(z0), pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F))
# Anova(z)

#' -----------
#' Information for table:
summary(z)[["coefficients"]][c("(Intercept)", "log_aphids"), 1:2]
#               Estimate Std. Error
# (Intercept) -2.2560406 0.22178224
# log_aphids   0.3490601 0.09238651
as.data.frame(Anova(z))
#                               Chisq Df   Pr(>Chisq)
# as.factor(Day_of_sampling) 98.26491 12 1.216901e-15
# log_aphids                 14.27524  1 1.579287e-04
summary(z)[["varcor"]]
# Groups Name        Std.Dev.
# obs    (Intercept) 0.58438
# LRT on glmm vs glm:
pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F)
# [1] 5.124567e-74
#' -----------



#' ========================================================================
#' ========================================================================

# wasps including trt
z <- glmer(wasp ~ 1 + as.factor(Day_of_sampling) + log_aphids + trt + (1|obs) +
               offset(log(sweeps.wasps)), data=d, family = poisson,
           control = glmer_cont)
z0 <- glm(wasp ~ 1 + as.factor(Day_of_sampling) + log_aphids + trt +
              offset(log(sweeps.wasps)), data=d, family = poisson)
# summary(z)
# summary(z0)
# c(deviance(z), deviance(z0), pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F))
# Anova(z)

#' ^^^^^^^^^^^^
#' This model was not included bc treatment was not significant.
#' ^^^^^^^^^^^^





#' ========================================================================
#' ========================================================================

# analysis of parasitism from dissections

# similar to above, you can ignore warning; addressing this actually
# reduces the log likelihood
z <- glmer(cbind(para, dissected - para) ~ 1 +
               as.factor(Day_of_sampling) + trt + log_aphids + (1|obs),
           data=d, family = binomial)
z0 <- glm(cbind(para, dissected - para) ~ 1 +
              as.factor(Day_of_sampling) + trt + log_aphids,
          data=d, family = binomial)

# summary(z)
# c(deviance(z), deviance(z0), pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F))
# Anova(z)

#' -----------
#' Information for table:
summary(z)[["coefficients"]][c("(Intercept)", "log_aphids", "trtInsecticide", "trtStrips"), 1:2]
#                   Estimate  Std. Error
# (Intercept)     1.14899853 0.003273049
# log_aphids     -0.17590798 0.003466949
# trtInsecticide  0.12624248 0.003271015
# trtStrips      -0.01225341 0.003271555
as.data.frame(Anova(z))
#                                   Chisq Df Pr(>Chisq)
# as.factor(Day_of_sampling) 15589548.103 12          0
# trt                            1503.736  2          0
# log_aphids                     2574.401  1          0
summary(z)[["varcor"]]
# Groups Name        Std.Dev.
# obs    (Intercept) 0.36841
# LRT on glmm vs glm:
pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F)
# [1] 1.288354e-23

#' -----------





#' ========================================================================
#' ========================================================================

# analyses of parasitism on separate dates

X <- data.frame(day = rep(unique(d$Day_of_sampling), 2),
                model = rep(1:2, each = length(unique(d$Day_of_sampling))),
                dev = NA, p.value = NA,
                obs.sd = NA, mean.dissected = NA, mean.para = NA,
                mean.para_prop = NA, log_aphids = NA, log_aphids.P = NA)
# Takes ~0.5 sec
set.seed(1695781797)
for(i.day in unique(d$Day_of_sampling)){

    dd <- d[d$Day_of_sampling == i.day,]
    dd$obs <- as.factor(1:nrow(dd))

    for (i.model in 1:2){

        i.inds <- X$model == i.model & X$day == i.day

        if(i.model==1){
            z <- glmer(cbind(para, dissected - para) ~ 1 + (1|obs), data=dd,
                       family = binomial, control = glmer_cont)
            z0 <- glm(cbind(para, dissected - para) ~ 1, data=dd,
                      family = binomial)
        }
        if(i.model==2){
            z <- glmer(cbind(para, dissected - para) ~ log_aphids + (1|obs),
                       data=dd, family = binomial, control = glmer_cont)
            z0 <- glm(cbind(para, dissected - para) ~ log_aphids, data=dd,
                      family = binomial)

            X$log_aphids[i.inds] <- summary(z)["coefficients"]$coefficients[2,1]
            X$log_aphids.P[i.inds] <- summary(z)["coefficients"]$coefficients[2,4]
        }

        X$dev[i.inds] <- -deviance(z) + deviance(z0)
        X$p.value[i.inds] <- pchisq(-deviance(z) + deviance(z0), df = 1, lower.tail = F)
        X$obs.sd[i.inds] <- as.numeric(VarCorr(z))^.5
        X$mean.dissected[i.inds] <- mean(dd$dissected)
        X$mean.para[i.inds] <- mean(dd$para)
        X$mean.para_prop[i.inds] <- mean(dd$para/dd$dissected)

        # show(as.numeric(VarCorr(z))^.5)
        # show(dd[,c("para","dissected")])
    }
}

XX <- X[!is.na(X$dev),]



#' ========================================================================
#' ========================================================================


summ_df <- XX |>
    as_tibble() |>
    # So that taking means below gets proportion that are significant:
    mutate(p.value = p.value < 0.05,
           log_aphids.P = log_aphids.P < 0.05) |>
    group_by(model) |>
    summarize(pos_aphids = mean(log_aphids > 0),
              pos_aphids_sig = mean(log_aphids[log_aphids.P] > 0),
              across(c(p.value, mean.dissected:log_aphids.P), mean),
              .groups = "drop") |>
    select(model, p.value, log_aphids.P, pos_aphids, pos_aphids_sig,
           everything())

# For results table:
summ_df[,1:5]
# # A tibble: 2 × 5
#   model p.value log_aphids.P pos_aphids pos_aphids_sig
#   <int>   <dbl>        <dbl>      <dbl>          <dbl>
# 1     1  0.4615      NA         NA                  NA
# 2     2  0.3846       0.2308     0.4615              0


# Other info
summ_df[,c(1, 6:ncol(summ_df))]
# # A tibble: 2 × 5
#   model mean.dissected mean.para mean.para_prop log_aphids
#   <int>          <dbl>     <dbl>          <dbl>      <dbl>
# 1     1          45.43     7.068         0.1493    NA
# 2     2          45.43     7.068         0.1493    -0.1317


disturb_hist <- function() {
    expression.list <- c(expression(sigma), expression(sigma[resid]~(aphids)),
                         expression(sigma[resid]~(aphids~and~wasps)),
                         expression(sigma[resid]~(aphids~and~ladybeetles)))
    par(mfrow=c(2,1), mai=c(.9, .9, .1, .1))
    for(i.model in 1:2){
        hist(XX$obs.sd[XX$model == i.model], col = "blue", breaks = .2*(0:7),
             main = "", xlab = expression.list[i.model])
        hist(XX$obs.sd[XX$model == i.model & XX$p.value > 0.05], add = T,
             col=adjustcolor("white", alpha.f = 0.5), breaks = .2*(0:7))
    }
}


if (write_plots) {
    save_plot(plot_out, disturb_hist, w = 6, h = 6)
} else {
    disturb_hist()
}
