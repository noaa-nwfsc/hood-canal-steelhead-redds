## Fit spawn timing models

dir.create("./figures/fit-spawn/", showWarnings = FALSE)


## Spawn timing: base model --------------------------------

## Define priors
mean(redd_yr$spawn_doy, na.rm = TRUE)
pr_spt = c(prior(student_t(3, 95, 2.5), class = "Intercept"),
           prior(student_t(3, 0, 2.5), class = "ar", lb = -1, ub = 1),
           prior(student_t(3, 0, 2.5), class = "b"),
           prior(student_t(3, 0, 2.5), class = "sd", lb = 0),
           prior(student_t(3, 0, 2.5), class = "sigma", lb = 0))


## Fit model
fit_spt = brm(spawn_doy ~ treatment + stage + treatment:stage + (1 | stream_fac) + (1 | year_fac) +
                    ar(year, gr = stream_fac, p = 1),
                  data = redd_yr[!is.na(spawn_doy), ],
                  prior = pr_spt,
                  iter = 2000,
                  warmup = 500,
                  chains = 4, cores = 4,
                  control = list(adapt_delta = 0.99),
                  seed = 4242)
fit_spt = add_criterion(fit_spt, c("bayes_R2", "loo"))
save(fit_spt, file = "./outputs/fit_spt.RData")


## Model checks
prior_summary(fit_spt)
check_hmc_diagnostics(fit_spt$fit)
rhat_highest(fit_spt$fit)
neff_lowest(fit_spt$fit)
hist(fit_spt$criteria$bayes_R2)
pp_check(fit_spt, type = "dens_overlay", ndraws = 50)
pp_check(fit_spt, type = "scatter_avg", ndraws = 10)
mcmc_trace(fit_spt)


## Model summaries
summary(fit_spt)
ranef(fit_spt)
fixef(fit_spt)

ce_spt = conditional_effects(fit_spt)
plot(ce_spt, ask = FALSE)

em_spt = emmeans(fit_spt, pairwise ~ stage | treatment)
plot(pairs(em_spt))



## Fitted values
fit_dat_spt = copy(redd_yr)
fit_dat_spt = fit_dat_spt[!is.na(spawn_doy), ]
fit_spt_fitted = fitted(fit_spt)
fit_dat_spt$fitted = fit_spt_fitted[ , "Estimate"]
fit_dat_spt$upper95 =  fit_spt_fitted[ , "Q2.5"]
fit_dat_spt$lower95 =  fit_spt_fitted[ , "Q97.5"]
g = ggplot(fit_dat_spt) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = spawn_doy, color = treatment)) +
    geom_line(aes(x = year, y = fitted)) +
    geom_point(aes(x = year, y = fitted), size = 1) +
    geom_ribbon(aes(x = year, ymin = lower95, ymax = upper95), alpha = 0.20) +
    labs(x = "Year", y = "Median Spawn Day of Year", color = "Treatment") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn/fit_spt_fitted.jpg", width = 8, height = 4)


## Treatment means by stream
nd_stream = redd_yr_stream[ , .(stream, treatment, stage, year_min, year_max)]
setnames(nd_stream, "stream", "stream_fac")
pred_mat = posterior_epred(fit_spt, newdata = nd_stream,
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

g = ggplot(redd_yr[!is.na(spawn_doy), ]) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = spawn_doy, color = treatment), alpha = 0.8) +
    geom_segment(data = pred_s, linewidth = 1,
                 aes(x = year_min, xend = year_max, y = median, yend = median,
                     color = treatment)) +
    geom_rect(data = pred_s,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = lower95,
                  ymax = upper95, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "Median Spawn Day of Year", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn/fit_spt_stream_mean.jpg", width = 8, height = 4)


## Stage means
nd = expand.grid(stage = c("before", "during", "after"),
                 treatment = c("control", "supplemented"))
pred_mat = posterior_epred(fit_spt, newdata = nd, re_formula = NA, incl_autocor = FALSE)
x = data.table(c_before = pred_mat[ , 1],
               c_during = pred_mat[ , 2],
               c_after  = pred_mat[ , 3],
               s_before = pred_mat[ , 4],
               s_during = pred_mat[ , 5],
               s_after  = pred_mat[ , 6])
m = melt(x, measure.vars = names(x))
m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
m[ , comparison := gsub("^[cs]_", "", variable)]
m[ , comparison := gsub("_", " - ", comparison)]
m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
m[ , comparison := factor(comparison, levels = c("before", "during", "after"))]

s = m[ , .(median = median(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975),
           perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]

g = ggplot(m) +
    geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
    geom_segment(data = s, linewidth = 1,
                 aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
                     color = comparison)) +
    geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
    labs(x = "Stage", y = "Median Spawn Day of Year") +
    scale_color_manual(values = M1[3:5]) +
    scale_fill_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn/fit_spt_stage_mean.jpg", width = 6, height = 4)


## Pairwise comparisons
nd = expand.grid(stage = c("before", "during", "after"),
                 treatment = c("control", "supplemented"))
pred_mat = posterior_epred(fit_spt, newdata = nd, re_formula = NA, incl_autocor = FALSE)
sig = as.data.frame(fit_spt)$sigma

x = data.table(c_during_before = pred_mat[ , 2] - pred_mat[ , 1],
               c_after_during  = pred_mat[ , 3] - pred_mat[ , 2],
               c_after_before  = pred_mat[ , 3] - pred_mat[ , 1],
               s_during_before = pred_mat[ , 5] - pred_mat[ , 4],
               s_after_during  = pred_mat[ , 6] - pred_mat[ , 5],
               s_after_before  = pred_mat[ , 6] - pred_mat[ , 4])
m = melt(x, measure.vars = names(x))
m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
m[ , comparison := gsub("^[cs]_", "", variable)]
m[ , comparison := gsub("_", " - ", comparison)]
m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
m[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]

s = m[ , .(median = median(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975),
           perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
s[ , y := c(rep(0.01, 3), rep(0, 3))]

g = ggplot(m) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
    geom_segment(data = s, linewidth = 1,
                 aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
                     color = comparison)) +
    geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
    geom_text(data = s, aes(x = comparison, y = 14, label = perc_pos_lab), size = 3, color = "grey25") +
    labs(x = "Comparison", y = "Difference in Median Spawn Day") +
    scale_color_manual(values = M1[3:5]) +
    scale_fill_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn/fit_spt_pairs.jpg", width = 6, height = 4)


ove = m[ , .(overlap = overlapping::overlap(x = list(value[treatment == "Control"],
                                                     value[treatment == "Supplemented"]))$OV),
        by = .(comparison)]
ove[ , overlap_lab := paste0("Overlap = ", formatC(overlap * 100, digits = 1, format = "f"), "%")]
leg = data.table(comparison = rep("during - before", 2), x = c(5, -8), y = c(0.1, 0.1),
                 treatment = c("Control", "Supplemented"))
leg[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
g = ggplot(m) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_density(aes(x = value, fill = treatment), color = NA, alpha = 0.2) +
    geom_linerange(data = s, linewidth = 0.5,
                 aes(y = y, xmin = lower95, xmax = upper95, color = treatment)) +
    geom_point(data = s, aes(x = median, y = y, color = treatment), size = 1) +
    labs(x = "Difference in Median Spawn Day", y = "Posterior density") +
    geom_text(data = ove, aes(x = -15, y = 0.2, label = overlap_lab), size = 3, color = "grey25", hjust = 0) +
    geom_text(data = leg, aes(x = x, y = y, label = treatment, color = treatment), size = 3) +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ comparison, ncol = 1) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn/fit_spt_pairs_overlap.jpg", width = 4, height = 5)


## Year effect
ref_mat = ranef(fit_spt)[[2]][,,1]
ref_spt = data.table(year = row.names(ref_mat),
                     estimate = ref_mat[ , "Estimate"],
                     lower95 = ref_mat[ , "Q2.5"],
                     upper95 = ref_mat[ , "Q97.5"])
g = ggplot(ref_spt) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_linerange(linewidth = 0.5, aes(x = year, ymin = lower95, ymax = upper95)) +
    geom_point(aes(x = year, y = estimate), size = 1) +
    labs(x = "Year", y = "Effect") +
    theme_simple()
print(g)
ggsave("./figures/fit-spawn/fit_spt_yr_effect.jpg", width = 6, height = 4)


## R-squared 
r2 = as.data.frame(bayes_R2(fit_spt, summary = FALSE))
g = ggplot(r2) +
    geom_histogram(aes(x = R2), bins = 20, color = "white",
                   fill = "grey30", alpha = 0.5, na.rm = TRUE) +
    labs(x = bquote(R^2), y = "Posterior density") +
    theme_simple()
print(g)
ggsave("./figures/fit-spawn/fit_spt_r2.jpg", width = 5, height = 4)



## PP checks
pp = posterior_predict(fit_spt, draw_ids = 1:30)
lst = vector("list", nrow(pp))
for(i in 1:nrow(pp)) {
    dt = data.table(yrep = pp[i, ], n = i)
    lst[[i]] = dt
}
ppc = rbindlist(lst)
#
g = ggplot(ppc) +
    geom_density(aes(x = yrep, group = n), color = "grey75", linewidth = 0.2) +
    geom_density(data = redd_yr[!is.na(spawn_doy)], aes(x = spawn_doy), color = "black", linewidth = 1) +
    labs(x = "Median Spawn Day of Year", y = "Density") +
    theme_simple()
print(g)
ggsave("./figures/fit-spawn/fit_spt_ppc.jpg", width = 5, height = 4)
