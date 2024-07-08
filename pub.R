## Publication figures


dir.create("./figures/pub/", showWarnings = FALSE)


## Numbers -------------------------------------------------


## Abundance means
tab_a = redd_yr[, .(treatment = unique(treatment),
                    mean_abund = round(mean(abund), 0),
                    min_abund = min(abund),
                    max_abund = max(abund),
                    year_min = year[which.min(abund)],
                    year_max = year[which.max(abund)]), by = .(stream)]

## Abundance means
tab_s = redd_yr[, .(treatment = unique(treatment),
                    mean_spawn = round(mean(spawn_doy, na.rm = TRUE), 0),
                    min_spawn = min(spawn_doy, na.rm = TRUE),
                    max_spawn = max(spawn_doy, na.rm = TRUE),
                    year_min = year[which.min(spawn_doy)],
                    year_max = year[which.max(spawn_doy)]), by = .(stream)]
max(tab_s$mean_spawn) - min(tab_s$mean_spawn)


## Abundance raw -------------------------------------------
dat = copy(redd_yr)
dat[ , stream_lab := as.character(stream)]
dat[ , stream_lab := gsub(" + Hatchery Cr", "", stream_lab, fixed = TRUE)]
dat[ , stream_lab := gsub(", SF + Vance Cr", "", stream_lab, fixed = TRUE)]
unique(dat$stream_lab)

or = dat[ , .(mean = mean(abund)), .(treatment, stream_lab)]
or = or[order(treatment, mean), ]
levs = or$stream_lab[c(1, 5, 2, 6, 3, 7, 4)]
dat[ , stream_lab := factor(stream_lab, levels = c(levs))]

g = ggplot(dat) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2, linewidth = 0.25) +
    aes(x = year, y = abund_ln, color = treatment) +
    geom_point(size = 1) +
    geom_line(linewidth = 0.5) +
    ylim(0, 7) +
    labs(x = "Year", y = "log Abundance", color = "") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream_lab, ncol = 2) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(legend.position = c(0.75, 0.12))
print(g)
ggsave("./figures/pub/abundance_ln_ts.jpg", width = 4, height = 4)



## Pairwise: abundance -------------------------------------
nd = expand.grid(stage = c("before", "during", "after"),
                 treatment = c("control", "supplemented"))
pred_mat = posterior_epred(fit_abn, newdata = nd, re_formula = NA, incl_autocor = FALSE)
sig = as.data.frame(fit_abn)$sigma

## log scale
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
#
s = m[ , .(mean = mean(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975),
           perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
s[ , y := c(rep(0.06, 3), rep(0, 3))]

s_exp = m[ , .(mean = mean(value_exp),
               lower95 = quantile(value_exp, probs = 0.025),
               upper95 = quantile(value_exp, probs = 0.975),
               perc_pos = sum(value_exp > 0) / .N), by = .(comparison, treatment)]

## original scale
x_exp = data.table(c_during_before = (exp(pred_mat[ , 2]) + exp(sig*sig / 2)) -
                       (exp(pred_mat[ , 1]) + exp(sig*sig / 2)),
                   c_after_during  = (exp(pred_mat[ , 3]) + exp(sig*sig / 2)) -
                       (exp(pred_mat[ , 2]) + exp(sig*sig / 2)),
                   c_after_before  = (exp(pred_mat[ , 3]) + exp(sig*sig / 2)) -
                       (exp(pred_mat[ , 1]) + exp(sig*sig / 2)),
                   s_during_before = (exp(pred_mat[ , 5]) + exp(sig*sig / 2)) -
                       (exp(pred_mat[ , 4]) + exp(sig*sig / 2)),
                   s_after_during  = (exp(pred_mat[ , 6]) + exp(sig*sig / 2)) -
                       (exp(pred_mat[ , 5]) + exp(sig*sig / 2)),
                   s_after_before  = (exp(pred_mat[ , 6]) + exp(sig*sig / 2)) -
                       (exp(pred_mat[ , 4]) + exp(sig*sig / 2)))
m_exp = melt(x_exp, measure.vars = names(x_exp))
m_exp[ , treatment := ifelse(grepl("^c_", variable), "control", "supplemented")]
m_exp[ , comparison := gsub("^[cs]_", "", variable)]
m_exp[ , comparison := gsub("_", " - ", comparison)]
m_exp[ , treatment := ifelse(treatment == "control", "Control", "Supplemented")]
m_exp[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
#
s_exp = m_exp[ , .(mean = mean(value),
                   lower95 = quantile(value, probs = 0.025),
                   upper95 = quantile(value, probs = 0.975)), by = .(comparison, treatment)]
s_exp_tab = copy(s_exp)
s_exp_tab[ , mean := formatC(mean, digits = 1, format = "f")]
s_exp_tab[ , lower95 := formatC(lower95, digits = 1, format = "f")]
s_exp_tab[ , upper95 := formatC(upper95, digits = 1, format = "f")]
fwrite(s_exp_tab, "./figures/pub/model_abund.csv")


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
    geom_point(data = s, aes(x = mean, y = y, color = treatment), size = 1) +
    labs(x = "Difference in log abundance", y = "Posterior density") +
    geom_text(data = ove, aes(x = -2.75, y = 1.3, label = overlap_lab), size = 3, color = "grey25", hjust = 0) +
    geom_text(data = leg, aes(x = x, y = y, label = treatment, color = treatment), size = 3) +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ comparison, ncol = 1) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/pub/fit_abn_pairs_overlap.jpg", width = 4, height = 5)



## Spawn timing raw ----------------------------------------
dat = copy(redd_yr)
dat[ , stream_lab := as.character(stream)]
dat[ , stream_lab := gsub(" + Hatchery Cr", "", stream_lab, fixed = TRUE)]
dat[ , stream_lab := gsub(", SF + Vance Cr", "", stream_lab, fixed = TRUE)]
unique(dat$stream_lab)

or = dat[ , .(mean = mean(abund)), .(treatment, stream_lab)]
or = or[order(treatment, mean), ]
levs = or$stream_lab[c(1, 5, 2, 6, 3, 7, 4)]
dat[ , stream_lab := factor(stream_lab, levels = c(levs))]

g = ggplot(dat) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2, linewidth = 0.25) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2, linewidth = 0.25) +
    aes(x = year, y = spawn_doy, color = treatment) +
    geom_point(size = 1) +
    geom_line(linewidth = 0.5) +
    labs(x = "Year", y = "Median spawn day of year", color = "") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream_lab, ncol = 2) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(legend.position = c(0.75, 0.12))
print(g)
ggsave("./figures/pub/spawn_ts.jpg", width = 4, height = 4)



## Spawn timing bar ----------------------------------------
dat = copy(redd_yr)
dat[ , stream_lab := as.character(stream)]
dat[ , stream_lab := gsub(" + Hatchery Cr", "", stream_lab, fixed = TRUE)]
dat[ , stream_lab := gsub(", SF + Vance Cr", "", stream_lab, fixed = TRUE)]
unique(dat$stream_lab)
#
or = dat[ , .(mean = mean(spawn_doy, na.rm = TRUE)), .(stream_lab)]
or = or[order(mean), ]
levs = or$stream_lab[c(1, 3, 4, 5, 2, 6, 7)]
dat[ , stream_lab := factor(stream_lab, levels = c(levs))]
dat[ , year_fac := factor(year, levels = 2023:2007)]
#
xint = seq(25, 150, 25)
M1_col = alpha(M1[3:5], 0.6)
M1_fill = alpha(M1[3:5], 0.2)
M1_solid = M1[3:5]
xlim = c(0, 155)
brks = seq(0, 150, 25)
expd = c(0, 0)
#
leg = data.table(stream_lab = rep("Big Beef Creek", 3),
                 stage = c("Before", "During", "After"),
                 year_fac = c("2008", "2015", "2022"),
                 x = 83)
leg[ , stage := factor(stage, levels = unique(stage))]
#
ylab = c("2007" =  "2007",
         "2008" =  "",
         "2009" =  "",
         "2010" =  "",
         "2011" =  "2011",
         "2012" =  "",
         "2013" =  "",
         "2014" =  "",
         "2015" =  "2015",
         "2016" =  "",
         "2017" =  "",
         "2018" =  "",
         "2019" =  "2019",
         "2020" =  "",
         "2021" =  "",
         "2022" =  "",
         "2023" =  "2023")
#
lst = vector("list", length(levs))
for(i in seq_along(lst)) {
    if(i == 7) xlab = "Median spawn day of year" else xlab = NULL
    g = ggplot(dat[stream_lab == levs[i], ]) +
        geom_vline(xintercept = xint, color = "grey85", linetype = 1) +
        geom_col(aes(y = year_fac, x = spawn_doy, fill = stage, color = stage), na.rm = TRUE) +
        scale_color_manual(values = M1_col) +
        scale_fill_manual(values = M1_fill) +
        scale_x_continuous(breaks = brks, limits = xlim, expand = expd) +
        scale_y_discrete(labels = ylab) +
        labs(x = xlab, y = NULL, subtitle = levs[i]) +
        theme_simple(base_size = 9) +
        theme(legend.position = "none") +
        theme(plot.title = element_text(hjust = 1))
    # if(i %in% 1:4) g = g + ylab("Control")
    # if(i %in% 5:7) g = g + ylab("Supplemented")
    if(i == 1) g = g + ggtitle("Control")
    if(i == 5) g = g + ggtitle("Supplemented")
    if(i == 1) g = g + geom_text(data = leg, aes(x = x, y = year_fac, label = stage),
                                 color = M1_solid, hjust = 0, size = 3)
    if(i == 4) bot = 10 else bot = 5
    g = g + theme(plot.margin = unit(c(0, 0, bot, 0), "pt"))
    if(!i %in% c(4, 7)) g = g + theme(axis.text.x = element_blank())
    print(g)
    lst[[i]] = g
}
#
g = lst[[1]] + lst[[2]] + lst[[3]] + lst[[4]] + lst[[5]] + lst[[6]] + lst[[7]] + plot_layout(ncol = 1,
 axis_titles = "collect")
print(g)
ggsave("./figures/pub/spawn_bar.jpg", width = 4, height = 8)



## Spawn tile ----------------------------------------------
dat = copy(redd)
dat[ , stream_lab := as.character(stream)]
dat[ , stream_lab := gsub(" + Hatchery Cr", "", stream_lab, fixed = TRUE)]
dat[ , stream_lab := gsub(", SF + Vance Cr", "", stream_lab, fixed = TRUE)]

or = redd_yr[ , .(mean = mean(spawn_doy, na.rm = TRUE)), .(treatment, stream)]
or = or[order(treatment, mean), ]
levs = as.character(or$stream)
levs = gsub(" + Hatchery Cr", "", levs, fixed = TRUE)
levs = gsub(", SF + Vance Cr", "", levs, fixed = TRUE)
dat[ , stream_lab := factor(stream_lab, levels = levs)]

g = ggplot(dat) +
    geom_tile(aes(x = doy, y = year, fill = log(redds))) +
    geom_hline(yintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_hline(yintercept = yrs_break[2], color = "grey50", linetype = 2) +
    labs(x = "Day of year", y = "Year") +
    scale_fill_viridis_c(na.value = "grey85") +
    xlim(10, 200) +
    facet_wrap( ~ stream_lab, ncol = 1) +
    theme_simple()
print(g)
ggsave("./figures/pub/spawn_doy_tile.jpg", width = 10, height = 10)


## Pairwise: spawn timing ----------------------------------

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

s = m[ , .(mean = mean(value),
           lower95 = quantile(value, probs = 0.025),
           upper95 = quantile(value, probs = 0.975),
           perc_pos = sum(value > 0) / .N), by = .(comparison, treatment)]
s[ , perc_pos_lab := paste0(formatC(perc_pos * 100, digits = 2, format = "f"), "%")]
s[ , y := c(rep(0.01, 3), rep(0, 3))]

s_tab = copy(s)
s_tab[ , perc_pos := NULL]
s_tab[ , perc_pos_lab := NULL]
s_tab[ , y := NULL]
s_tab[ , mean := formatC(mean, digits = 1, format = "f")]
s_tab[ , lower95 := formatC(lower95, digits = 1, format = "f")]
s_tab[ , upper95 := formatC(upper95, digits = 1, format = "f")]
fwrite(s_tab, "./figures/pub/model_spawn.csv")

ove = m[ , .(overlap = overlapping::overlap(x = list(value[treatment == "Control"],
                                                     value[treatment == "Supplemented"]))$OV),
        by = .(comparison)]
ove[ , overlap_lab := paste0("Overlap = ", formatC(overlap * 100, digits = 1, format = "f"), "%")]
leg = data.table(comparison = rep("during - before", 2), x = c(5.5, -8), y = c(0.1, 0.1),
                 treatment = c("Control", "Supplemented"))
leg[ , comparison := factor(comparison, levels = c("during - before", "after - during", "after - before"))]
g = ggplot(m) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_density(aes(x = value, fill = treatment), color = NA, alpha = 0.2) +
    geom_linerange(data = s, linewidth = 0.5,
                 aes(y = y, xmin = lower95, xmax = upper95, color = treatment)) +
    geom_point(data = s, aes(x = mean, y = y, color = treatment), size = 1) +
    labs(x = "Difference in Median Spawn Day", y = "Posterior density") +
    geom_text(data = ove, aes(x = -15, y = 0.2, label = overlap_lab), size = 3, color = "grey25", hjust = 0) +
    geom_text(data = leg, aes(x = x, y = y, label = treatment, color = treatment), size = 3) +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ comparison, ncol = 1) +
    theme_simple(grid = TRUE) +
    theme(legend.position = "none")
print(g)
ggsave("./figures/pub/fit_spt_pairs_overlap.jpg", width = 4, height = 5)



## Year effects --------------------------------------------

## Abundance
ref_mat = ranef(fit_abn)[[2]][,,1]
ref_abn = data.table(year = row.names(ref_mat),
                     estimate = ref_mat[ , "Estimate"],
                     lower95 = ref_mat[ , "Q2.5"],
                     upper95 = ref_mat[ , "Q97.5"])
ref_abn[ , year := as.numeric(year)]
g = ggplot(ref_abn) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_linerange(linewidth = 0.5, aes(x = year, ymin = lower95, ymax = upper95)) +
    geom_point(aes(x = year, y = estimate), size = 1) +
    labs(x = "Year", y = "Effect") +
    theme_simple()
print(g)
ggsave("./figures/pub/fit_abn_yr_effect.jpg", width = 6, height = 4)


## Spawn timing
ref_mat = ranef(fit_spt)[[2]][,,1]
ref_spt = data.table(year = row.names(ref_mat),
                     estimate = ref_mat[ , "Estimate"],
                     lower95 = ref_mat[ , "Q2.5"],
                     upper95 = ref_mat[ , "Q97.5"])
ref_spt[ , year := as.numeric(year)]
g = ggplot(ref_spt) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_linerange(linewidth = 0.5, aes(x = year, ymin = lower95, ymax = upper95)) +
    geom_point(aes(x = year, y = estimate), size = 1) +
    labs(x = "Year", y = "Effect") +
    theme_simple()
print(g)
ggsave("./figures/pub/fit_spt_yr_effect.jpg", width = 6, height = 4)



## Stream effects ------------------------------------------

levs = c("Big Beef Creek",
         "Union River",
         "Little Quilcene River",
         "Tahuya River",
         "Dewatto River",
         "Duckabush River",
         "Skokomish River")

## Abundance
ref_mat = ranef(fit_abn)[[1]][,,1]
ref_abn = data.table(stream = row.names(ref_mat),
                     estimate = ref_mat[ , "Estimate"],
                     lower95 = ref_mat[ , "Q2.5"],
                     upper95 = ref_mat[ , "Q97.5"])
ref_abn[ , stream := gsub(" + Hatchery Cr", "", stream, fixed = TRUE)]
ref_abn[ , stream := gsub(", SF + Vance Cr", "", stream, fixed = TRUE)]
ref_abn[ , stream := factor(stream, levels = levs)]
g = ggplot(ref_abn) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_linerange(linewidth = 0.5, aes(y = stream, xmin = lower95, xmax = upper95)) +
    geom_point(aes(x = estimate, y = stream), size = 1) +
    labs(x = "Effect", y = "") +
    theme_simple()
print(g)
ggsave("./figures/pub/fit_abn_stream_effect.jpg", width = 4, height = 4)


## Spawn timing
ref_mat = ranef(fit_spt)[[1]][,,1]
ref_spt = data.table(stream = row.names(ref_mat),
                     estimate = ref_mat[ , "Estimate"],
                     lower95 = ref_mat[ , "Q2.5"],
                     upper95 = ref_mat[ , "Q97.5"])
ref_spt[ , stream := gsub(" + Hatchery Cr", "", stream, fixed = TRUE)]
ref_spt[ , stream := gsub(", SF + Vance Cr", "", stream, fixed = TRUE)]
ref_spt[ , stream := factor(stream, levels = levs)]
g = ggplot(ref_spt) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_linerange(linewidth = 0.5, aes(y = stream, xmin = lower95, xmax = upper95)) +
    geom_point(aes(x = estimate, y = stream), size = 1) +
    labs(x = "Effect", y = "") +
    theme_simple()
print(g)
ggsave("./figures/pub/fit_spt_stream_effect.jpg", width = 4, height = 4)

