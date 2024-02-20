## Fit models

dir.create("./figures/fit/", showWarnings = FALSE)


## Base model ----------------------------------------------

## Define priors
mean(redd_yr$abund_ln)
pr = c(prior(student_t(3, 3, 2.5), class = "Intercept"),
       prior(student_t(3, 0, 2.5), class = "ar", lb = -1, ub = 1),
       prior(student_t(3, 0, 2.5), class = "b"),
       prior(student_t(3, 0, 2.5), class = "sd", lb = 0),
       prior(student_t(3, 0, 2.5), class = "sigma", lb = 0))


## Fit model
fit_brm = brm(abund_ln ~ treatment + stage + treatment:stage + (1 | stream_fac) + (1 | year_fac) +
                    ar(year, gr = stream_fac, p = 1),
                  data = redd_yr,
                  prior = pr,
                  iter = 2000,
                  warmup = 500,
                  chains = 4, cores = 4,
                  control = list(adapt_delta = 0.90),
                  seed = 4242)
fit_brm = add_criterion(fit_brm, c("bayes_R2"))
save(fit_brm, file = "./outputs/fit_brm.RData")


## Model checks
prior_summary(fit_brm)
check_hmc_diagnostics(fit_brm$fit)
rhat_highest(fit_brm$fit)
neff_lowest(fit_brm$fit)
hist(fit_brm$criteria$bayes_R2)
pp_check(fit_brm, type = "dens_overlay", ndraws = 50)
pp_check(fit_brm, type = "scatter_avg", ndraws = 10)


## Model summaries
summary(fit_brm)
ranef(fit_brm)
fixef(fit_brm)

ce = conditional_effects(fit_brm)
plot(ce, ask = FALSE)

em = emmeans(fit_brm, pairwise ~ stage | treatment)
plot(pairs(em))



## Fitted values
fit_dat = copy(redd_yr)
fit_brm_fitted = fitted(fit_brm)
fit_dat$fitted = fit_brm_fitted[ , "Estimate"]
fit_dat$upper95 =  fit_brm_fitted[ , "Q2.5"]
fit_dat$lower95 =  fit_brm_fitted[ , "Q97.5"]
g = ggplot(fit_dat) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund_ln, color = treatment)) +
    geom_line(aes(x = year, y = fitted)) +
    geom_point(aes(x = year, y = fitted), size = 1) +
    geom_ribbon(aes(x = year, ymin = lower95, ymax = upper95), alpha = 0.20) +
    labs(x = "Year", y = "log Abundance", color = "Treatment", title = "Fitted values") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit/fit_brm_fitted.jpg", width = 8, height = 4)


## Treatment means by stream
nd_stream = redd_yr_stream[ , .(stream, treatment, stage, year_min, year_max)]
setnames(nd_stream, "stream", "stream_fac")
pred_mat = posterior_epred(fit_brm, newdata = nd_stream,
                           re_formula = ~ (1 | stream_fac),
                           incl_autocor = FALSE)
pred_df = epred_df(pred_mat, nd_stream)
pred_s = pred_df[ , .(year_min = unique(year_min),
                      year_max = unique(year_max),
                      median = median(sample),
                      lower95 = quantile(sample, probs = 0.025),
                      upper95 = quantile(sample, probs = 0.975)),
                 by = .(stream_fac, treatment, stage)]
setnames(pred_s, "stream_fac", "stream")

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund_ln, color = treatment), alpha = 0.8) +
    geom_segment(data = pred_s, linewidth = 1,
                 aes(x = year_min, xend = year_max, y = median, yend = median,
                     color = treatment)) +
    geom_rect(data = pred_s,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = lower95,
                  ymax = upper95, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "log Abundance", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit/fit_brm_stream_mean.jpg", width = 8, height = 4)



## Pairwise comparisons
nd = expand.grid(stage = c("before", "during", "after"),
                 treatment = c("control", "supplemented"))
pred_mat = posterior_epred(fit_brm, newdata = nd, re_formula = NA, incl_autocor = FALSE)

x = data.table(c_before_during = pred_mat[ , 1] - pred_mat[ , 2],
               c_during_after  = pred_mat[ , 2] - pred_mat[ , 3],
               c_before_after  = pred_mat[ , 1] - pred_mat[ , 3],
               s_before_during = pred_mat[ , 4] - pred_mat[ , 5],
               s_during_after  = pred_mat[ , 5] - pred_mat[ , 6],
               s_before_after  = pred_mat[ , 4] - pred_mat[ , 6])
m = melt(x, measure.vars = names(x))
m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
m[ , comparison := gsub("^[cs]_", "", variable)]
m[ , comparison := gsub("_", " - ", comparison)]

s = m[ , .(median = median(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975)), by = .(comparison, treatment)]

g = ggplot(m) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
    geom_segment(data = s, linewidth = 1,
                 aes(x = comparison, xend = comparison, y = lower95, yend = upper95, color = comparison)) +
    geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
    labs(x = "Contrast", y = "Estimate") +
    scale_color_manual(values = M1[3:5]) +
    scale_fill_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit/fit_brm_pairs.jpg", width = 6, height = 4)
