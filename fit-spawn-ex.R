## Fit spawn timing models -- expanded data


dir.create("./figures/fit-spawn-ex/", showWarnings = FALSE)


## Expand redd data to for spawn timing --------------------
## One row for each redd observed
redd_sub = redd[redds > 0, ]
redd_ex = vector("list", nrow(redd_sub))
for(i in 1:nrow(redd_sub)) {
    redd_i = redd_sub[i, ]
    n = redd_i$redds
    if(n == 1) {
        redd_ex[[i]] = redd_i
    } else {
        lst = vector("list", n)
        for(j in 1:n) lst[[j]] = redd_i
        redd_ex[[i]] = rbindlist(lst)
    }
}
redd_ex = rbindlist(redd_ex)
redd_ex = redd_ex[order(stream, year, date), ]
redd_ex[ , redds := 1]
redd_ex[ , redds_sum := cumsum(redds), by = .(stream, year)]
redd_ex[ , year_fac := as.factor(year)]
redd_ex[ , stream_fac := as.factor(stream)]
redd_ex[ , doy_anom := doy - mean(doy), .(stream)]
redd_ex[ , stage := factor(stage, levels = c("before", "during", "after"))]
save(redd_ex, file = "./outputs/redd_ex.RData")


## Summarize by year
redd_ex_yr = redd_ex[ , .(spawn_date = median(date)), by = .(stream, year, treatment)]
redd_ex_yr[ , spawn_doy := as.numeric(format(spawn_date, "%j"))]
redd_ex_yr = redd_ex_yr[order(stream, year), ]
save(redd_ex_yr, file = "./outputs/redd_ex_yr.RData")


## verify calc matches previous method
x = redd_yr[!is.na(spawn_doy), ]
all.equal(x$spawn_doy, redd_ex_yr$spawn_doy)  ## should equal TRUE


g = ggplot(redd_ex) +
    geom_histogram(aes(doy), bins = 10) +
    labs(x = "Day of year", y = "Count") +
    facet_wrap( ~ stream, scales = "free") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn-ex/spawn_doy_hist.jpg", width = 6, height = 6)


g = ggplot(redd_ex) +
    geom_point(aes(x = year, y = doy, color = treatment), alpha = 0.25) +
    geom_point(data = redd_ex_yr, aes(x = year, y = spawn_doy), color = "grey20", size = 2) +
    geom_line(data = redd_ex_yr, aes(x = year, y = spawn_doy), color = "grey20") +
    labs(x = "Year", y = "Day of year", color = "") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn-ex/spawn_doy_dot.jpg", width = 7, height = 6)


g = ggplot(redd_ex) +
    geom_jitter(aes(x = stage, y = doy, color = stage, fill = stage), width = 0.2, alpha = 0.1) +
    geom_violin(aes(x = stage, y = doy, color = stage, fill = stage), alpha = 0.25, adjust = 1.1) +
    labs(x = "Year", y = "Day of year", color = "", fill = "") +
    scale_fill_manual(values = M1[3:5]) +
    scale_color_manual(values = M1[3:5]) +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn-ex/spawn_doy_violin.jpg", width = 7, height = 6)


dat = redd_ex[ , .(doy_var = var(doy), doy_anom_var = var(doy_anom)), by = .(stage, treatment)]
g = ggplot(dat) +
    geom_col(aes(x = stage, y = doy_anom_var, fill = stage, color = stage), alpha = 0.25) +
    labs(x = "Year", y = "Day of year", color = "") +
    facet_wrap( ~ treatment) +
    labs(x = "", y = "Spawn day variance") +
    scale_fill_manual(values = M1[3:5]) +
    scale_color_manual(values = M1[3:5]) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn-ex/spawn_doy_var.jpg", width = 7, height = 6)


dat = redd_ex[ , .(doy_var = var(doy), doy_anom_var = var(doy_anom)), by = .(stage, treatment, stream)]
g = ggplot(dat) +
    geom_col(aes(x = stage, y = doy_anom_var, fill = stage, color = stage), alpha = 0.25) +
    labs(x = "Year", y = "Day of year", color = "") +
    facet_wrap( ~ stream) +
    labs(x = "", y = "Spawn day variance") +
    scale_fill_manual(values = M1[3:5]) +
    scale_color_manual(values = M1[3:5]) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn-ex/spawn_doy_var_stream.jpg", width = 7, height = 6)



## F-test -- changing variance -----------------------------
cb = redd_ex[treatment == "control" & stage == "before", ]
ca = redd_ex[treatment == "control" & stage == "after", ]
var.test(cb$doy, ca$doy)             ## sig
var.test(cb$doy_anom, ca$doy_anom)   ## no sig

sb = redd_ex[treatment == "supplemented" & stage == "before", ]
sa = redd_ex[treatment == "supplemented" & stage == "after", ]
var.test(sb$doy, sa$doy)             ## sig
var.test(sb$doy_anom, sa$doy_anom)   ## sig



## Bootstrap test -- changing variance ---------------------
## Test statistic: var before / var after

set.seed(4242)
boot_dat_lst = vector("list", 1000)
s = split(redd_ex, by = c("stream", "year"))
for(i in 1:1000) {
    s_lst = vector("list", length(s))
    for(j in seq_along(s)) {
        s_j = s[[j]]
        ind = sample(1:nrow(s_j), nrow(s_j), replace = TRUE)
        dat_j = s_j[ind]
        dat_j[ , iter := i]
        s_lst[[j]] = dat_j
    }
    boot_dat_lst[[i]] = rbindlist(s_lst)
}
spt_boot = rbindlist(boot_dat_lst)
spt_boot[ , doy_anom := doy - mean(doy), .(iter, stream)]
save(spt_boot, file = "./outputs/spt_boot.RData")


## Summarize by stream
spt_boot_stream = spt_boot[ , .(var = var(doy), var_anom = var(doy_anom)),
                           by = .(iter, treatment, stage, stream)]
save(spt_boot_stream, file = "./outputs/spt_boot_stream.RData")

spt_boot_stream_s = spt_boot_stream[ , .(mean_var = mean(var),
                                         lower95_var = quantile(var, probs = 0.025),
                                         upper95_var = quantile(var, probs = 0.975),
                                         mean_var_anom = mean(var_anom),
                                         lower95_var_anom = quantile(var_anom, probs = 0.025),
                                         upper95_var_anom = quantile(var_anom, probs = 0.975)),
                                    by = .(treatment, stage, stream)]
save(spt_boot_stream_s, file = "./outputs/spt_boot_stream_s.RData")


spt_boot_stream_ratio = dcast(spt_boot_stream, iter + treatment + stream ~ stage, value.var = c("var", "var_anom"))
spt_boot_stream_ratio[ , var_ratio := var_before / var_after]
spt_boot_stream_ratio[ , var_anom_ratio := var_anom_before / var_anom_after]
save(spt_boot_stream_ratio, file = "./outputs/spt_boot_stream_ratio.RData")

spt_boot_stream_ratio_s = spt_boot_stream_ratio[ , .(var_ratio = median(var_ratio),
                                                     lower95 = quantile(var_ratio, probs = 0.025),
                                                     upper95 = quantile(var_ratio, probs = 0.975)),
                                                by = .(treatment, stream)]
save(spt_boot_stream_ratio_s, file = "./outputs/spt_boot_stream_ratio_s.RData")

g = ggplot(spt_boot_stream) +
    geom_violin(aes(x = stage, y = var_anom, fill = treatment, color = treatment), alpha = 0.25) +
    geom_segment(data = spt_boot_stream_s, linewidth = 1,
                 aes(x = stage, xend = stage, y = lower95_var_anom, yend = upper95_var_anom,
                     color = treatment)) +
    geom_point(data = spt_boot_var_s,
               aes(x = stage, y = mean_var_anom, color = treatment), size = 2) +
    scale_fill_manual(values = M1) +
    scale_color_manual(values = M1) +
    labs(x = "", y = "Spawn day variance") +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn-ex/boot_var_stream.jpg", width = 7, height = 6)

g = ggplot(spt_boot_stream_ratio) +
    geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = var_anom_ratio, color = treatment, fill = treatment),
                   bins = 30, position = "identity", alpha = 0.25) +
    labs(x = "Variance ratio (var before / var after)", y = "Count") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn-ex/boot_var_stream_ratio_hist.jpg", width = 7, height = 6)

g = ggplot(spt_boot_stream_ratio_s) +
    geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
    geom_point(aes(x = var_ratio, y = stream, color = treatment)) +
    geom_linerange(linewidth = 0.5,
                 aes(y = stream, xmin = lower95, xmax = upper95, color = treatment)) +
    labs(x = "Variance ratio (var before / var after)", y = "", color = "") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn-ex/boot_var_stream_ratio_dot.jpg", width = 7, height = 6)




## Summarize by stage
spt_boot_stage = spt_boot[ , .(var = var(doy), var_anom = var(doy_anom)),
                          by = .(iter, treatment, stage)]
save(spt_boot_stream, file = "./outputs/spt_boot_stage.RData")

spt_boot_stage_s = spt_boot_stage[ , .(mean_var = mean(var),
                                       lower95_var = quantile(var, probs = 0.025),
                                       upper95_var = quantile(var, probs = 0.975),
                                       mean_var_anom = mean(var_anom),
                                       lower95_var_anom = quantile(var_anom, probs = 0.025),
                                       upper95_var_anom = quantile(var_anom, probs = 0.975)),
                                  by = .(treatment, stage)]
save(spt_boot_stage_s, file = "./outputs/spt_boot_stage_s.RData")


spt_boot_stage_ratio = dcast(spt_boot_stage, iter + treatment ~ stage, value.var = c("var", "var_anom"))
spt_boot_stage_ratio[ , var_ratio := var_before / var_after]
spt_boot_stage_ratio[ , var_anom_ratio := var_anom_before / var_anom_after]
save(spt_boot_stage_ratio, file = "./outputs/spt_boot_stage_ratio.RData")

spt_boot_stage_ratio_s = spt_boot_stage_ratio[ , .(var_ratio = median(var_ratio),
                                                   lower95 = quantile(var_ratio, probs = 0.025),
                                                   upper95 = quantile(var_ratio, probs = 0.975)),
                                              by = .(treatment)]
save(spt_boot_stage_ratio_s, file = "./outputs/spt_boot_stage_ratio_s.RData")

g = ggplot(spt_boot_stage) +
    geom_violin(aes(x = stage, y = var_anom, fill = stage, color = stage), alpha = 0.25) +
    geom_segment(data = spt_boot_stage_s, linewidth = 1,
                 aes(x = stage, xend = stage, y = lower95_var_anom, yend = upper95_var_anom,
                     color = stage)) +
    geom_point(data = spt_boot_stage_s,
               aes(x = stage, y = mean_var_anom, color = stage), size = 2) +
    scale_fill_manual(values = M1[3:5]) +
    scale_color_manual(values = M1[3:5]) +
    labs(x = "", y = "Spawn day variance") +
    facet_wrap( ~ treatment) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn-ex/boot_var_stage.jpg", width = 7, height = 6)

g = ggplot(spt_boot_stage_ratio) +
    geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
    geom_histogram(aes(x = var_anom_ratio, color = treatment, fill = treatment),
                   bins = 30, position = "identity", alpha = 0.25) +
    labs(x = "Variance ratio (var before / var after)", y = "Count") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/fit-spawn-ex/boot_var_stage_ratio_hist.jpg", width = 7, height = 6)

g = ggplot(spt_boot_stage_ratio_s) +
    geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
    geom_point(aes(x = var_ratio, y = treatment, color = treatment)) +
    geom_linerange(linewidth = 0.5,
                 aes(y = treatment, xmin = lower95, xmax = upper95, color = treatment)) +
    labs(x = "Variance ratio (var before / var after)", y = "", color = "") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/fit-spawn-ex/boot_var_stage_ratio_dot.jpg", width = 4, height = 6)



## Fit model w/ constant residual variance -----------------
## mean(redd_ex$doy, na.rm = TRUE)
## pr_stx1 = c(prior(student_t(3, 122, 25), class = "Intercept"),
##             prior(student_t(3, 0, 25), class = "b"),
##             prior(student_t(3, 0, 25), class = "sd", lb = 0),
##             prior(student_t(3, 0, 25), class = "sigma", lb = 0))
##
## fm_stx1 = bf(doy ~ treatment + stage + treatment:stage + (1 | stream_fac) + (1 | year_fac))
## get_prior(fm_stx1, data = redd_ex)
## fit_stx1 = brm(fm_stx1,
##                data = redd_ex,
##                prior = pr_stx1,
##                iter = 2000,
##                warmup = 500,
##                chains = 4, cores = 4,
##                control = list(adapt_delta = 0.999, max_treedepth = 13),
##                seed = 4242)
## fit_stx1 = add_criterion(fit_stx1, c("bayes_R2", "loo"))
## save(fit_stx1, file = "./outputs/fit_stx1.RData")
##
##
## ## Model checks
## prior_summary(fit_stx1)
## check_hmc_diagnostics(fit_stx1$fit)
## rhat_highest(fit_stx1$fit)
## neff_lowest(fit_stx1$fit)
## hist(fit_stx1$criteria$bayes_R2)
## pp_check(fit_stx1, type = "dens_overlay", ndraws = 50)
## pp_check(fit_stx1, type = "scatter_avg", ndraws = 10)
## mcmc_trace(fit_stx1)
##
##
## ## Conditional effects
## ce_stx1 = conditional_effects(fit_stx1)
## plot(ce_stx1, ask = FALSE)
##
##
## ## Fitted values
## fit_dat_stx1 = copy(redd_ex)
## fit_stx1_fitted = fitted(fit_stx1)
## fit_dat_stx1$fitted = fit_stx1_fitted[ , "Estimate"]
## fit_dat_stx1$upper95 =  fit_stx1_fitted[ , "Q2.5"]
## fit_dat_stx1$lower95 =  fit_stx1_fitted[ , "Q97.5"]
##
## g = ggplot(fit_dat_stx1) +
##     geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
##     geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
##     geom_point(aes(x = year, y = doy, color = treatment), alpha = 0.25) +
##     geom_line(aes(x = year, y = fitted)) +
##     geom_point(aes(x = year, y = fitted), size = 1) +
##     geom_ribbon(aes(x = year, ymin = lower95, ymax = upper95), alpha = 0.40) +
##     labs(x = "Year", y = "Spawn day", color = "Treatment") +
##     scale_color_manual(values = M1) +
##     facet_wrap( ~ stream, scale = "free_y") +
##     theme_simple(grid = TRUE)
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_fitted.jpg", width = 8, height = 4)
##
##
## ## Treatment means by stream
## nd_stream = redd_yr_stream[ , .(stream, treatment, stage, year_min, year_max)]
## setnames(nd_stream, "stream", "stream_fac")
## pred_mat = posterior_epred(fit_stx1, newdata = nd_stream,
##                            re_formula = ~ (1 | stream_fac),
##                            incl_autocor = FALSE)
## pred_df = epred_df(pred_mat, nd_stream)
## pred_s = pred_df[ , .(year_min = unique(year_min),
##                       year_max = unique(year_max),
##                       median = median(sample),
##                       lower95 = quantile(sample, probs = 0.025),
##                       upper95 = quantile(sample, probs = 0.975)),
##                  by = .(stream_fac, treatment, stage)]
## setnames(pred_s, "stream_fac", "stream")
##
## g = ggplot(redd_ex) +
##     geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
##     geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
##     geom_point(aes(x = year, y = doy, color = treatment), alpha = 0.1) +
##     geom_segment(data = pred_s, linewidth = 1,
##                  aes(x = year_min, xend = year_max, y = median, yend = median,
##                      color = treatment)) +
##     geom_rect(data = pred_s,
##               aes(xmin = year_min,
##                   xmax = year_max,
##                   ymin = lower95,
##                   ymax = upper95, fill = treatment), alpha = 0.25) +
##     labs(x = "Year", y = "Spawn day", color = "Treatment", fill = "Treatment") +
##     scale_color_manual(values = M1) +
##     scale_fill_manual(values = M1) +
##     facet_wrap( ~ stream, scale = "free_y") +
##     theme_simple(grid = TRUE)
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_stream_mean.jpg", width = 8, height = 4)
##
##
## ## Stage means
## nd = expand.grid(stage = c("before", "during", "after"),
##                  treatment = c("control", "supplemented"))
## pred_mat = posterior_epred(fit_stx1, newdata = nd, re_formula = NA, incl_autocor = FALSE)
## x = data.table(c_before = pred_mat[ , 1],
##                c_during = pred_mat[ , 2],
##                c_after  = pred_mat[ , 3],
##                s_before = pred_mat[ , 4],
##                s_during = pred_mat[ , 5],
##                s_after  = pred_mat[ , 6])
## m = melt(x, measure.vars = names(x))
## m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
## m[ , comparison := gsub("^[cs]_", "", variable)]
## m[ , comparison := gsub("_", " - ", comparison)]
## m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
## m[ , comparison := factor(comparison, levels = c("before", "during", "after"))]
##
## s = m[ , .(median = median(value),
##            lower95 = quantile(value, probs = 0.025),
##            upper95 = quantile(value, probs = 0.975),
##            perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
## s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
##
## g = ggplot(m) +
##     geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
##     geom_segment(data = s, linewidth = 1,
##                  aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
##                      color = comparison)) +
##     geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
##     labs(x = "Stage", y = "Spawn day") +
##     scale_color_manual(values = M1[3:5]) +
##     scale_fill_manual(values = M1[3:5]) +
##     facet_wrap( ~ treatment) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_stage_mean.jpg", width = 6, height = 4)
##
##
## ## Pairwise comparisons
## nd = expand.grid(stage = c("before", "during", "after"),
##                  treatment = c("control", "supplemented"))
## pred_mat = posterior_epred(fit_stx1, newdata = nd, re_formula = NA, incl_autocor = FALSE)
## sig = as.data.frame(fit_stx1)$sigma
##
## x = data.table(c_during_before = pred_mat[ , 2] - pred_mat[ , 1],
##                c_after_during  = pred_mat[ , 3] - pred_mat[ , 2],
##                c_after_before  = pred_mat[ , 3] - pred_mat[ , 1],
##                s_during_before = pred_mat[ , 5] - pred_mat[ , 4],
##                s_after_during  = pred_mat[ , 6] - pred_mat[ , 5],
##                s_after_before  = pred_mat[ , 6] - pred_mat[ , 4])
## m = melt(x, measure.vars = names(x))
## m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
## m[ , comparison := gsub("^[cs]_", "", variable)]
## m[ , comparison := gsub("_", " - ", comparison)]
## m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
## m[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
##
## s = m[ , .(median = median(value),
##            lower95 = quantile(value, probs = 0.025),
##            upper95 = quantile(value, probs = 0.975),
##            perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
## s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
## s[ , y := c(rep(0.01, 3), rep(0, 3))]
##
## g = ggplot(m) +
##     geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
##     geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
##     geom_segment(data = s, linewidth = 1,
##                  aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
##                      color = comparison)) +
##     geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
##     geom_text(data = s, aes(x = comparison, y = 14, label = perc_pos_lab), size = 3, color = "grey25") +
##     labs(x = "Comparison", y = "Difference in mean spawn day") +
##     scale_color_manual(values = M1[3:5]) +
##     scale_fill_manual(values = M1[3:5]) +
##     facet_wrap( ~ treatment) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_pairs.jpg", width = 6, height = 4)
##
##
## ove = m[ , .(overlap = overlapping::overlap(x = list(value[treatment == "Control"],
##                                                      value[treatment == "Supplemented"]))$OV),
##         by = .(comparison)]
## ove[ , overlap_lab := paste0("Overlap = ", formatC(overlap * 100, digits = 1, format = "f"), "%")]
## leg = data.table(comparison = rep("during - before", 2), x = c(5, -13), y = c(0.1, 0.1),
##                  treatment = c("Control", "Supplemented"))
## leg[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
## g = ggplot(m) +
##     geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
##     geom_density(aes(x = value, fill = treatment), color = NA, alpha = 0.2) +
##     geom_linerange(data = s, linewidth = 0.5,
##                  aes(y = y, xmin = lower95, xmax = upper95, color = treatment)) +
##     geom_point(data = s, aes(x = median, y = y, color = treatment), size = 1) +
##     labs(x = "Difference in mean spawn day", y = "Posterior density") +
##     geom_text(data = ove, aes(x = -20, y = 0.13, label = overlap_lab), size = 3, color = "grey25", hjust = 0) +
##     geom_text(data = leg, aes(x = x, y = y, label = treatment, color = treatment), size = 3) +
##     scale_color_manual(values = M1) +
##     scale_fill_manual(values = M1) +
##     facet_wrap( ~ comparison, ncol = 1) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_pairs_overlap.jpg", width = 4, height = 5)
##
##
## ## Year effect
## ref_mat = ranef(fit_stx1)[[2]][,,1]
## ref_spt = data.table(year = row.names(ref_mat),
##                      estimate = ref_mat[ , "Estimate"],
##                      lower95 = ref_mat[ , "Q2.5"],
##                      upper95 = ref_mat[ , "Q97.5"])
## g = ggplot(ref_spt) +
##     geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
##     geom_linerange(linewidth = 0.5, aes(x = year, ymin = lower95, ymax = upper95)) +
##     geom_point(aes(x = year, y = estimate), size = 1) +
##     labs(x = "Year", y = "Effect") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_yr_effect.jpg", width = 6, height = 4)
##
##
## ## R-squared
## r2 = as.data.frame(fit_stx1$criteria$bayes_R2)
## g = ggplot(r2) +
##     geom_histogram(aes(x = R2), bins = 20, color = "white",
##                    fill = "grey30", alpha = 0.5, na.rm = TRUE) +
##     labs(x = bquote(R^2), y = "Posterior density") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_r2.jpg", width = 5, height = 4)
##
##
## ## PP checks
## pp = posterior_predict(fit_stx1, draw_ids = 1:30)
## lst = vector("list", nrow(pp))
## for(i in 1:nrow(pp)) {
##     dt = data.table(yrep = pp[i, ], n = i)
##     lst[[i]] = dt
## }
## ppc = rbindlist(lst)
## #
## g = ggplot(ppc) +
##     geom_density(aes(x = yrep, group = n), color = "grey75", linewidth = 0.2) +
##     geom_density(data = redd_ex, aes(x = doy), color = "black", linewidth = 1, adjust = 1.1) +
##     labs(x = "Spawn day", y = "Density") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx1_ppc.jpg", width = 5, height = 4)



## Fit model w/ varying residual variances -----------------
## mean(redd_ex$doy, na.rm = TRUE)
## pr_stx2 = c(prior(student_t(3, 122, 25), class = "Intercept"),
##             prior(student_t(3, 0, 5), class = "Intercept", dpar = "sigma"),
##             prior(student_t(3, 0, 25), class = "b"),
##             prior(student_t(3, 0, 5), class = "b", dpar = "sigma"),
##             prior(student_t(3, 0, 25), class = "sd", lb = 0))
##
## fm_stx2 = bf(doy ~ treatment + stage + treatment:stage + (1 | stream_fac) + (1 | year_fac),
##              sigma ~ treatment + stage + treatment:stage)
## get_prior(fm_stx2, data = redd_ex)
## fit_stx2 = brm(fm_stx2,
##                data = redd_ex,
##                prior = pr_stx2,
##                iter = 2000,
##                warmup = 500,
##                chains = 4, cores = 4,
##                control = list(adapt_delta = 0.99, max_treedepth = 13),
##                seed = 4242)
## fit_stx2 = add_criterion(fit_stx2, c("bayes_R2", "loo"))
## save(fit_stx2, file = "./outputs/fit_stx2.RData")
##
##
## ## Model checks
## prior_summary(fit_stx2)
## check_hmc_diagnostics(fit_stx2$fit)
## rhat_highest(fit_stx2$fit)
## neff_lowest(fit_stx2$fit)
## hist(fit_stx2$criteria$bayes_R2)
## pp_check(fit_stx2, type = "dens_overlay", ndraws = 50)
## pp_check(fit_stx2, type = "scatter_avg", ndraws = 10)
## mcmc_trace(fit_stx2)
##
##
## ## Conditional effects
## ce_stx2 = conditional_effects(fit_stx2)
## plot(ce_stx2, ask = FALSE)
##
##
## ## Fitted values
## fit_dat_stx2 = copy(redd_ex)
## fit_stx2_fitted = fitted(fit_stx2)
## fit_dat_stx2$fitted = fit_stx2_fitted[ , "Estimate"]
## fit_dat_stx2$upper95 =  fit_stx2_fitted[ , "Q2.5"]
## fit_dat_stx2$lower95 =  fit_stx2_fitted[ , "Q97.5"]
##
## g = ggplot(fit_dat_stx2) +
##     geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
##     geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
##     geom_point(aes(x = year, y = doy, color = treatment), alpha = 0.25) +
##     geom_line(aes(x = year, y = fitted)) +
##     geom_point(aes(x = year, y = fitted), size = 1) +
##     geom_ribbon(aes(x = year, ymin = lower95, ymax = upper95), alpha = 0.40) +
##     labs(x = "Year", y = "Spawn day", color = "Treatment") +
##     scale_color_manual(values = M1) +
##     facet_wrap( ~ stream, scale = "free_y") +
##     theme_simple(grid = TRUE)
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_fitted.jpg", width = 8, height = 4)
##
##
## ## Treatment means by stream
## nd_stream = redd_yr_stream[ , .(stream, treatment, stage, year_min, year_max)]
## setnames(nd_stream, "stream", "stream_fac")
## pred_mat = posterior_epred(fit_stx2, newdata = nd_stream,
##                            re_formula = ~ (1 | stream_fac),
##                            incl_autocor = FALSE)
## pred_df = epred_df(pred_mat, nd_stream)
## pred_s = pred_df[ , .(year_min = unique(year_min),
##                       year_max = unique(year_max),
##                       median = median(sample),
##                       lower95 = quantile(sample, probs = 0.025),
##                       upper95 = quantile(sample, probs = 0.975)),
##                  by = .(stream_fac, treatment, stage)]
## setnames(pred_s, "stream_fac", "stream")
##
## g = ggplot(redd_ex) +
##     geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
##     geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
##     geom_point(aes(x = year, y = doy, color = treatment), alpha = 0.1) +
##     geom_segment(data = pred_s, linewidth = 1,
##                  aes(x = year_min, xend = year_max, y = median, yend = median,
##                      color = treatment)) +
##     geom_rect(data = pred_s,
##               aes(xmin = year_min,
##                   xmax = year_max,
##                   ymin = lower95,
##                   ymax = upper95, fill = treatment), alpha = 0.25) +
##     labs(x = "Year", y = "Spawn day", color = "Treatment", fill = "Treatment") +
##     scale_color_manual(values = M1) +
##     scale_fill_manual(values = M1) +
##     facet_wrap( ~ stream, scale = "free_y") +
##     theme_simple(grid = TRUE)
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_stream_mean.jpg", width = 8, height = 4)
##
##
## ## Stage means
## nd = expand.grid(stage = c("before", "during", "after"),
##                  treatment = c("control", "supplemented"))
## pred_mat = posterior_epred(fit_stx2, newdata = nd, re_formula = NA, incl_autocor = FALSE)
## x = data.table(c_before = pred_mat[ , 1],
##                c_during = pred_mat[ , 2],
##                c_after  = pred_mat[ , 3],
##                s_before = pred_mat[ , 4],
##                s_during = pred_mat[ , 5],
##                s_after  = pred_mat[ , 6])
## m = melt(x, measure.vars = names(x))
## m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
## m[ , comparison := gsub("^[cs]_", "", variable)]
## m[ , comparison := gsub("_", " - ", comparison)]
## m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
## m[ , comparison := factor(comparison, levels = c("before", "during", "after"))]
##
## s = m[ , .(median = median(value),
##            lower95 = quantile(value, probs = 0.025),
##            upper95 = quantile(value, probs = 0.975),
##            perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
## s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
##
## g = ggplot(m) +
##     geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
##     geom_segment(data = s, linewidth = 1,
##                  aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
##                      color = comparison)) +
##     geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
##     labs(x = "Stage", y = "Spawn day") +
##     scale_color_manual(values = M1[3:5]) +
##     scale_fill_manual(values = M1[3:5]) +
##     facet_wrap( ~ treatment) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_stage_mean.jpg", width = 6, height = 4)
##
##
## ## Pairwise comparisons
## nd = expand.grid(stage = c("before", "during", "after"),
##                  treatment = c("control", "supplemented"))
## pred_mat = posterior_epred(fit_stx2, newdata = nd, re_formula = NA, incl_autocor = FALSE)
## sig = as.data.frame(fit_stx2)$sigma
##
## x = data.table(c_during_before = pred_mat[ , 2] - pred_mat[ , 1],
##                c_after_during  = pred_mat[ , 3] - pred_mat[ , 2],
##                c_after_before  = pred_mat[ , 3] - pred_mat[ , 1],
##                s_during_before = pred_mat[ , 5] - pred_mat[ , 4],
##                s_after_during  = pred_mat[ , 6] - pred_mat[ , 5],
##                s_after_before  = pred_mat[ , 6] - pred_mat[ , 4])
## m = melt(x, measure.vars = names(x))
## m[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
## m[ , comparison := gsub("^[cs]_", "", variable)]
## m[ , comparison := gsub("_", " - ", comparison)]
## m[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
## m[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
##
## s = m[ , .(median = median(value),
##            lower95 = quantile(value, probs = 0.025),
##            upper95 = quantile(value, probs = 0.975),
##            perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
## s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
## s[ , y := c(rep(0.01, 3), rep(0, 3))]
##
## g = ggplot(m) +
##     geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
##     geom_violin(aes(x = comparison, y = value, fill = comparison), color = NA, alpha = 0.2) +
##     geom_segment(data = s, linewidth = 1,
##                  aes(x = comparison, xend = comparison, y = lower95, yend = upper95,
##                      color = comparison)) +
##     geom_point(data = s, aes(x = comparison, y = median, color = comparison), size = 2) +
##     geom_text(data = s, aes(x = comparison, y = 14, label = perc_pos_lab), size = 3, color = "grey25") +
##     labs(x = "Comparison", y = "Difference in mean spawn day") +
##     scale_color_manual(values = M1[3:5]) +
##     scale_fill_manual(values = M1[3:5]) +
##     facet_wrap( ~ treatment) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_pairs.jpg", width = 6, height = 4)
##
##
## ove = m[ , .(overlap = overlapping::overlap(x = list(value[treatment == "Control"],
##                                                      value[treatment == "Supplemented"]))$OV),
##         by = .(comparison)]
## ove[ , overlap_lab := paste0("Overlap = ", formatC(overlap * 100, digits = 1, format = "f"), "%")]
## leg = data.table(comparison = rep("during - before", 2), x = c(8, -18), y = c(0.07, 0.07),
##                  treatment = c("Control", "Supplemented"))
## leg[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
## g = ggplot(m) +
##     geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
##     geom_density(aes(x = value, fill = treatment), color = NA, alpha = 0.2) +
##     geom_linerange(data = s, linewidth = 0.5,
##                  aes(y = y, xmin = lower95, xmax = upper95, color = treatment)) +
##     geom_point(data = s, aes(x = median, y = y, color = treatment), size = 1) +
##     labs(x = "Difference in mean spawn day", y = "Posterior density") +
##     geom_text(data = ove, aes(x = -30, y = 0.095, label = overlap_lab), size = 3, color = "grey25", hjust = 0) +
##     geom_text(data = leg, aes(x = x, y = y, label = treatment, color = treatment), size = 3) +
##     scale_color_manual(values = M1) +
##     scale_fill_manual(values = M1) +
##     facet_wrap( ~ comparison, ncol = 1) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_pairs_overlap.jpg", width = 4, height = 5)
##
##
## ## Year effect
## ref_mat = ranef(fit_stx2)[[2]][,,1]
## ref_spt = data.table(year = row.names(ref_mat),
##                      estimate = ref_mat[ , "Estimate"],
##                      lower95 = ref_mat[ , "Q2.5"],
##                      upper95 = ref_mat[ , "Q97.5"])
## g = ggplot(ref_spt) +
##     geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
##     geom_linerange(linewidth = 0.5, aes(x = year, ymin = lower95, ymax = upper95)) +
##     geom_point(aes(x = year, y = estimate), size = 1) +
##     labs(x = "Year", y = "Effect") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_yr_effect.jpg", width = 6, height = 4)
##
##
## ## R-squared
## r2 = as.data.frame(fit_stx2$criteria$bayes_R2)
## g = ggplot(r2) +
##     geom_histogram(aes(x = R2), bins = 20, color = "white",
##                    fill = "grey30", alpha = 0.5, na.rm = TRUE) +
##     labs(x = bquote(R^2), y = "Posterior density") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_r2.jpg", width = 5, height = 4)
##
##
## ## PP checks
## pp = posterior_predict(fit_stx2, draw_ids = 1:30)
## lst = vector("list", nrow(pp))
## for(i in 1:nrow(pp)) {
##     dt = data.table(yrep = pp[i, ], n = i)
##     lst[[i]] = dt
## }
## ppc = rbindlist(lst)
## #
## g = ggplot(ppc) +
##     geom_density(aes(x = yrep, group = n), color = "grey75", linewidth = 0.2) +
##     geom_density(data = redd_ex, aes(x = doy), color = "black", linewidth = 1, adjust = 1.1) +
##     labs(x = "Spawn day", y = "Density") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_ppc.jpg", width = 5, height = 4)
##
##
## ## Sigma changes
## post_stx2 = as.matrix(fit_stx2)
## head(post_stx2[ , 1:15])
## sig_c_b = post_stx2[ , "b_sigma_Intercept"]
## sig_c_d = sig_c_b + post_stx2[ , "b_sigma_stageduring"]
## sig_c_a = sig_c_b + post_stx2[ , "b_sigma_stageafter"]
## #
## sig_t_b = sig_c_b + post_stx2[ , "b_sigma_treatmentsupplemented"]
## sig_t_d = sig_c_d + post_stx2[ , "b_sigma_treatmentsupplemented"] +
##     post_stx2[ , "b_sigma_treatmentsupplemented:stageduring"]
## sig_t_a = sig_c_a + post_stx2[ , "b_sigma_treatmentsupplemented"] +
##     post_stx2[ , "b_sigma_treatmentsupplemented:stageafter"]
##
## sig = c(sig_c_b,
##         sig_c_d,
##         sig_c_a,
##         sig_t_b,
##         sig_t_d,
##         sig_t_a)
## treat = c(rep("control", 6000 * 3), rep("supplemented", 6000 * 3))
## stage = c(rep("before", 6000),
##           rep("during", 6000),
##           rep("after", 6000),
##           rep("before", 6000),
##           rep("during", 6000),
##           rep("after", 6000))
## post_stx2_sig = data.table(treatment = treat, stage = stage, value = sig)
## post_stx2_sig[ , stage := factor(stage, levels = c("before", "during", "after"))]
##
## g = ggplot(post_stx2_sig) +
##     geom_density(aes(x = exp(value), color = stage)) +
##     scale_color_manual(values = M1) +
##     labs(y = "Posterior density", x = "Residual SD", color = "") +
##     facet_wrap( ~ treatment) +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_sigma_dens.jpg", width = 6, height = 4)
##
##
## post_stx2_sig_s = post_stx2_sig[ , .(median = median(value),
##                                      lower95 = quantile(value, probs = 0.025),
##                                      upper95 = quantile(value, probs = 0.975)), by = .(treatment, stage)]
## g = ggplot(post_stx2_sig) +
##     geom_violin(aes(x = stage, y = exp(value), fill = stage), color = NA, alpha = 0.2) +
##     geom_segment(data = post_stx2_sig_s, linewidth = 1,
##                  aes(x = stage, xend = stage, y = exp(lower95), yend = exp(upper95),
##                      color = stage)) +
##     geom_point(data = post_stx2_sig_s, aes(x = stage, y = exp(median), color = stage), size = 2) +
##     labs(x = "", y = "Residual SD") +
##     scale_color_manual(values = M1[3:5]) +
##     scale_fill_manual(values = M1[3:5]) +
##     facet_wrap( ~ treatment) +
##     theme_simple(grid = TRUE) +
##     theme(legend.position = "none")
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_sigma_violin.jpg", width = 6, height = 4)
##
##
## ## Sigma ratio (before / after)
## sig_c = data.table(treatment = "control", sd_ratio = exp(sig_c_b) / exp(sig_c_a))
## sig_t = data.table(treatment = "supplemented", sd_ratio = exp(sig_t_b) / exp(sig_t_a))
## post_stx2_sig_ratio = rbind(sig_c, sig_t)
## g = ggplot(post_stx2_sig_ratio) +
##     geom_vline(xintercept = 1, color = "grey50", linetype = 2) +
##     geom_histogram(aes(x = sd_ratio, color = treatment, fill = treatment),
##                    bins = 30, position = "identity", alpha = 0.25) +
##     scale_color_manual(values = M1) +
##     scale_fill_manual(values = M1) +
##     labs(y = "Posterior samples", x = "Residual SD ratio (before / after)", color = "", fill = "") +
##     theme_simple()
## print(g)
## ggsave("./figures/fit-spawn-ex/fit_stx2_sigma_ratio_hist.jpg", width = 6, height = 4)



## Model comparison ----------------------------------------
## loo_compare(fit_stx1, fit_stx2)
## loo(fit_stx1)
## loo(fit_stx2)
