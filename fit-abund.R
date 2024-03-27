## Fit abundance models

dir.create("./figures/fit-abund/", showWarnings = FALSE)


## Abundance: base model -----------------------------------

## Define priors
mean(redd_yr$abund_ln)
pr_abn = c(prior(student_t(3, 3, 2.5), class = "Intercept"),
           prior(student_t(3, 0, 2.5), class = "ar", lb = -1, ub = 1),
           prior(student_t(3, 0, 2.5), class = "b"),
           prior(student_t(3, 0, 2.5), class = "sd", lb = 0),
           prior(student_t(3, 0, 2.5), class = "sigma", lb = 0))


## Fit model
fit_abn = brm(abund_ln ~ treatment + stage + treatment:stage + (1 | stream_fac) + (1 | year_fac) +
                    ar(year, gr = stream_fac, p = 1),
                  data = redd_yr,
                  prior = pr_abn,
                  iter = 2000,
                  warmup = 500,
                  chains = 4, cores = 4,
                  control = list(adapt_delta = 0.99),
                  seed = 4242)
fit_abn = add_criterion(fit_abn, c("bayes_R2"))
save(fit_abn, file = "./outputs/fit_abn.RData")


## Model checks
prior_summary(fit_abn)
check_hmc_diagnostics(fit_abn$fit)
rhat_highest(fit_abn$fit)
neff_lowest(fit_abn$fit)
hist(fit_abn$criteria$bayes_R2)
pp_check(fit_abn, type = "dens_overlay", ndraws = 50)
pp_check(fit_abn, type = "scatter_avg", ndraws = 10)


## Model summaries
summary(fit_abn)
ranef(fit_abn)
fixef(fit_abn)

ce_abn = conditional_effects(fit_abn)
plot(ce_abn, ask = FALSE)

em_abn = emmeans(fit_abn, pairwise ~ stage | treatment)
plot(pairs(em_abn))



## Fitted values
fit_dat_abn = copy(redd_yr)
fit_abn_fitted = fitted(fit_abn)
fit_dat_abn$fitted = fit_abn_fitted[ , "Estimate"]
fit_dat_abn$upper95 =  fit_abn_fitted[ , "Q2.5"]
fit_dat_abn$lower95 =  fit_abn_fitted[ , "Q97.5"]
g = ggplot(fit_dat_abn) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund_ln, color = treatment)) +
    geom_line(aes(x = year, y = fitted)) +
    geom_point(aes(x = year, y = fitted), size = 1) +
    geom_ribbon(aes(x = year, ymin = lower95, ymax = upper95), alpha = 0.20) +
    labs(x = "Year", y = "log Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-abund/fit_abn_fitted.jpg", width = 8, height = 4)


## Treatment means by stream
nd_stream = redd_yr_stream[ , .(stream, treatment, stage, year_min, year_max)]
setnames(nd_stream, "stream", "stream_fac")
pred_mat = posterior_epred(fit_abn, newdata = nd_stream,
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
ggsave("./figures/fit-abund/fit_abn_stream_mean.jpg", width = 8, height = 4)


## Stage means
nd = expand.grid(stage = c("before", "during", "after"),
                 treatment = c("control", "supplemented"))
pred_mat = posterior_epred(fit_abn, newdata = nd, re_formula = NA, incl_autocor = FALSE)
sig = as.data.frame(fit_abn)$sigma
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
m[ , value_exp := exp(value) + exp(..sig^2 / 2)]
m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
m[ , comparison := factor(comparison, levels = c("before", "during", "after"))]

## log scale
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
    labs(x = "Stage", y = "log Abundance") +
    scale_color_manual(values = M1[3:5]) +
    scale_fill_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-abund/fit_abn_stage_mean.jpg", width = 6, height = 4)


## Pairwise comparisons
nd = expand.grid(stage = c("before", "during", "after"),
                 treatment = c("control", "supplemented"))
pred_mat = posterior_epred(fit_abn, newdata = nd, re_formula = NA, incl_autocor = FALSE)
sig = as.data.frame(fit_abn)$sigma

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
m[ , value_exp := exp(value) + exp(..sig^2 / 2)]
m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
m[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]

## log scale
s = m[ , .(median = median(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975),
           perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
s[ , y := c(rep(0.06, 3), rep(0, 3))]

g = ggplot(m) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
    geom_segment(data = s, linewidth = 1,
                 aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
                     color = comparison)) +
    geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
    geom_text(data = s, aes(x = comparison, y = 2.3, label = perc_pos_lab), size = 3, color = "grey25") +
    labs(x = "Comparison", y = "Difference in log abundance") +
    scale_color_manual(values = M1[3:5]) +
    scale_fill_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-abund/fit_abn_pairs.jpg", width = 6, height = 4)


ove = m[ , .(overlap = overlapping::overlap(x = list(value[treatment == "Control"],
                                                     value[treatment == "Supplemented"]))$OV),
        by = .(comparison)]
ove[ , overlap_lab := paste0("Overlap = ", formatC(overlap * 100, digits = 1, format = "f"), "%")]
leg = data.table(comparison = rep("during - before", 2), x = c(-0.95, 2.1), y = c(0.9, 0.9),
                 treatment = c("Control", "Supplemented"))
leg[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
g = ggplot(m) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_density(aes(x = value, fill = treatment), color = NA, alpha = 0.2) +
    geom_linerange(data = s, linewidth = 0.5,
                 aes(y = y, xmin = lower95, xmax = upper95, color = treatment)) +
    geom_point(data = s, aes(x = median, y = y, color = treatment), size = 1) +
    labs(x = "Difference in log abundance", y = "Posterior density") +
    geom_text(data = ove, aes(x = -2.75, y = 1.3, label = overlap_lab), size = 3, color = "grey25", hjust = 0) +
    geom_text(data = leg, aes(x = x, y = y, label = treatment, color = treatment), size = 3) +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ comparison, ncol = 1) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-abund/fit_abn_pairs_overlap.jpg", width = 4, height = 5)


## original scale
x = data.table(c_during_before = (exp(pred_mat[ , 2]) + exp((sig*sig) / 2)) - (exp(pred_mat[ , 1]) + exp((sig*sig) / 2)),
               c_after_during  = (exp(pred_mat[ , 3]) + exp((sig*sig) / 2)) - (exp(pred_mat[ , 2]) + exp((sig*sig) / 2)),
               c_after_before  = (exp(pred_mat[ , 3]) + exp((sig*sig) / 2)) - (exp(pred_mat[ , 1]) + exp((sig*sig) / 2)),
               s_during_before = (exp(pred_mat[ , 5]) + exp((sig*sig) / 2)) - (exp(pred_mat[ , 4]) + exp((sig*sig) / 2)),
               s_after_during  = (exp(pred_mat[ , 6]) + exp((sig*sig) / 2)) - (exp(pred_mat[ , 5]) + exp((sig*sig) / 2)),
               s_after_before  = (exp(pred_mat[ , 6]) + exp((sig*sig) / 2)) - (exp(pred_mat[ , 4]) + exp((sig*sig) / 2)))
m = melt(x, measure.vars = names(x))
m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
m[ , comparison := gsub("^[cs]_", "", variable)]
m[ , comparison := gsub("_", " - ", comparison)]
m[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
# m[ , value_exp := exp(value) + exp(..sig^2 / 2)]

s = m[ , .(median = median(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975)), by = .(comparison, treatment)]

g = ggplot(m) +
    geom_hline(yintercept = 1, color = "grey50", linetype = 2) +
    # geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
    geom_segment(data = s, linewidth = 1,
                 aes(x = comparison, xend = comparison, y = lower95, yend = upper95, color = comparison)) +
    geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
    labs(x = "Contrast", y = "Difference in log abundance") +
    # ylim(-300, 300) +
    scale_color_manual(values = M1[3:5]) +
    scale_fill_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-abund/fit_abn_pairs_exp.jpg", width = 6, height = 4)


## Year effect
ref_mat = ranef(fit_abn)[[2]][,,1]
ref_abn = data.table(year = row.names(ref_mat),
                     estimate = ref_mat[ , "Estimate"],
                     lower95 = ref_mat[ , "Q2.5"],
                     upper95 = ref_mat[ , "Q97.5"])
g = ggplot(ref_abn) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_linerange(linewidth = 0.5, aes(x = year, ymin = lower95, ymax = upper95)) +
    geom_point(aes(x = year, y = estimate), size = 1) +
    labs(x = "Year", y = "Effect") +
    theme_simple()
print(g)
ggsave("./figures/fit-abund/fit_abn_yr_effect.jpg", width = 6, height = 4)


## R-squared 
r2 <- as.data.frame(bayes_R2(fit_abn, summary = FALSE))
g <- ggplot(r2) +
    geom_histogram(aes(x = R2), bins = 20, color = "white",
                   fill = "grey30", alpha = 0.5, na.rm = TRUE) +
    labs(x = bquote(R^2), y = "Posterior density") +
    theme_simple()
print(g)
ggsave("./figures/fit-abund/fit_abn_r2.jpg", width = 5, height = 4)



## PP checks
pp = posterior_predict(fit_abn, draw_ids = 1:30)
lst = vector("list", nrow(pp))
for(i in 1:nrow(pp)) {
    dt = data.table(yrep = pp[i, ], n = i)
    lst[[i]] = dt
}
ppc = rbindlist(lst)
#
g = ggplot(ppc) +
    geom_density(aes(x = yrep, group = n), color = "grey75", linewidth = 0.2) +
    geom_density(data = redd_yr, aes(x = abund_ln), color = "black", linewidth = 1) +
    labs(x = "log Abundance", y = "Density") +
    theme_simple()
print(g)
ggsave("./figures/fit-abund/fit_abn_ppc.jpg", width = 5, height = 4)
