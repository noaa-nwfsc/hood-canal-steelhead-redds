## Data

dir.create("./figures/data/", showWarnings = FALSE)


## Read raw data -------------------------------------------
redd_raw = fread("./data/hood_canal_redd_data.csv")

## Save outupts
save(redd_raw, file = "./outputs/redd_raw.RData")



## Clean raw data ------------------------------------------
redd_full = copy(redd_raw)
redd_full[ , date := as.Date(date)]
redd_full[ , doy := as.numeric(format(date, "%j"))]
class(redd_full$date)
sum(is.na(redd_full$date))

## Barry indicated this survey should be removed
redd_full = redd_full[!(stream == "Big Beef Creek" & year == 2011 & is.na(rm_length))]

## Fix missing survey lengths based on Barry's corrections
ind = which(redd_full$stream == "Tahuya River" & is.na(redd_full$rm_length))
redd_full$rm_upper[ind] = 2.6
redd_full$rm_lower[ind] = 1.0
redd_full$rm_length[ind] = 1.6

ind = which(redd_full$stream == "Skokomish River, SF" & is.na(redd_full$rm_length))
redd_full$rm_upper[ind] = 7.8
redd_full$rm_lower[ind] = 6.8
redd_full$rm_length[ind] = 1.0

## Upper and lower rm were switched
ind = which(redd_full$rm_length < 0)
rmlo = redd_full$rm_lower[ind]
rmup = redd_full$rm_upper[ind]
redd_full$rm_lower[ind] = rmup
redd_full$rm_upper[ind] = rmlo
redd_full$rm_length[ind] = rmlo - rmup

## Handle survey length for "spot check survey": set to 0.01
redd_full[ , rm_length := ifelse(rm_length == 0 & survey_type == "SPOT", 0.01, rm_length)]
redd_full[ , rkm_length := rm_length * 1.60934]

## Add a combined Skokomish SF + Vance Cr. stream
sk_va = redd_full[stream %in% c("Skokomish River, SF", "Vance Creek (SF Skok trib)")]
sk_va[ , stream := "Skokomish River, SF + Vance Cr"]
redd_full = rbind(redd_full, sk_va)

## Add a combined Duckabush + Hatchery Cr
unique(redd_full$stream)
du_hc = redd_full[stream %in% c("Duckabush River", "Hatchery Creek (Duck trib)")]
du_hc[ , stream := "Duckabush River + Hatchery Cr"]
redd_full = rbind(redd_full, du_hc)


## Calc cumulative sum of redds -- order by date first
redd_sp = split(redd_full, by = c("stream", "year"))
redd_lst = lapply(redd_sp, function(m) {
    ord = m[order(date), ]
    ord[ , redds_sum := cumsum(redds)]
    return(ord)
})
redd_full = rbindlist(redd_lst)


## Add redds / rkm
# redd_full[ , redds_rkm := redds / rkm_length]
# redd_full[ , redds_rkm_sum := cumsum(redds_rkm), by = .(stream, year)]


## Add a treatment indicator
## Supplemented:
##    Dewatto River
##    South Fork Skokomish River + Vance
##    Duckabush River
## Controls:
##    Tahuya River
##    Big Beef Creek
##    Little Quilcene River
unique(redd_full$stream)
redd_full[ , treatment := NA]
contr = c("Tahuya River", "Big Beef Creek", "Little Quilcene River", "Union River")
redd_full[ , treatment := ifelse(stream %in% contr, "control", treatment)]
suppl = c("Duckabush River", "Duckabush River + Hatchery Cr", "Dewatto River", "Skokomish River, SF",
          "Vance Creek (SF Skok trib)", "Skokomish River, SF + Vance Cr")
redd_full[ , treatment := ifelse(stream %in% suppl, "supplemented", treatment)]


## Add an experimental stage indicator
##   before = pre-supplementation       2007-2010
##   during = during supplementation    2011-2019
##   after  = post-supplementation      2020-2023
redd_full[ , stage := "before"]
redd_full[ , stage := ifelse(year > 2010, "during", stage)]
redd_full[ , stage := ifelse(year > 2019, "after", stage)]


## Save outupts
save(redd_full, file = "./outputs/redd_full.RData")



## Filter data ---------------------------------------------
##  Limit 'stream' to only include those listed in MS
##
inc = c("Dewatto River",
        "Duckabush River + Hatchery Cr",
        "Skokomish River, SF + Vance Cr",
        "Big Beef Creek",
        "Little Quilcene River",
        "Tahuya River",
        "Union River")
redd = redd_full[stream %in% inc, ]
redd[ , stream := factor(stream, levels = inc)]


g = ggplot(redd) +
    geom_point(aes(x = date, y = redds)) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)

g = ggplot(redd) +
    geom_point(aes(x = date, y = redds_sum)) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)

g = ggplot(redd) +
    geom_histogram(aes(x = redds), bins = 10) +
    facet_wrap( ~ stream, scale = "free_x") +
    theme_simple(grid = TRUE)
print(g)

g = ggplot(redd) +
    geom_boxplot(aes(x = as.factor(year), y = redds)) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)


## Save outputs
save(redd, file = "./outputs/redd.RData")



## Summarize by year ---------------------------------------
## Summarize quantities:
##   - Annual abundance
##   - Mean spawn date
redd_sp = split(redd, by = c("stream", "year"), drop = TRUE)
redd_lst = vector("list", length(redd_sp))
for(i in seq_along(redd_lst)) {
    m = redd_sp[[i]]
    if(sum(m$redds) > 0) {
        ## Calc median spawn date: rep the date for each redd and take median
        sub = m[redds > 0, ]
        lst = vector("list", nrow(sub))
        lst_sd = vector("list", nrow(sub))
        for(j in 1:nrow(sub)) lst[[j]] = rep(sub$date[j], sub$redds[j])
        for(j in 1:nrow(sub)) lst_sd[[j]] = rep(sub$doy[j], sub$redds[j])
        msd = median(do.call("c", lst))
        # browser()
        sdsd = sd(do.call("c", lst_sd))
    } else {
        ## systems that had no redds for a year are set to NA
        msd = as.Date(NA)
        sdsd = as.numeric(NA)
    }
    dt = data.table(stream = unique(m$stream),
                    year = unique(m$year),
                    treatment = unique(m$treatment),
                    stage = unique(m$stage),
                    n_surveys = nrow(m),
                    abund = max(m$redds_sum),
                    # abund_rkm = max(m$redds_rkm_sum),
                    spawn_date = msd,
                    spawn_doy_sd = sdsd)
    redd_lst[[i]] = dt
}
redd_yr = rbindlist(redd_lst)
redd_yr[ , spawn_doy := as.numeric(format(spawn_date, "%j"))]


## Use weir count of females to get Big Beef Creek redd abundance in 2007
## 16 females passed the weir and assume 1.23 redds / female
bbc_2007 = data.table(stream = "Big Beef Creek",
                      year = 2007,
                      treatment = "control",
                      stage = "before",
                      n_surveys = 0,
                      abund = round(16 * 1.23, digits = 0),
                      # abund_rkm = round(16 * 1.23, digits = 0) / (5.7 * 1.60934),
                      spawn_date = as.Date(NA),
                      spawn_doy_sd = NA,
                      spawn_doy = NA)
redd_yr = rbind(redd_yr, bbc_2007)
redd_yr = redd_yr[order(stream, year), ]


## Check oddities
# redd_raw[Stream == "Big Beef Creek" & RunYear == 2007, ]  ## no data
# redd[stream == "Big Beef Creek" & year == 2022, ] ## zero redds
unique(redd_yr$stream)


redd_yr[ , rkm_surveyed := NA]
redd_yr[ , rkm_surveyed := ifelse(stream == "Dewatto River", 4.8, rkm_surveyed)]
redd_yr[ , rkm_surveyed := ifelse(stream == "Duckabush River + Hatchery Cr", 6.1, rkm_surveyed)]
redd_yr[ , rkm_surveyed := ifelse(stream == "Skokomish River, SF + Vance Cr", 39.8, rkm_surveyed)]
redd_yr[ , rkm_surveyed := ifelse(stream == "Big Beef Creek", 11.4, rkm_surveyed)]
redd_yr[ , rkm_surveyed := ifelse(stream == "Little Quilcene River", 8.5, rkm_surveyed)]
redd_yr[ , rkm_surveyed := ifelse(stream == "Tahuya River", 17.7, rkm_surveyed)]
redd_yr[ , rkm_surveyed := ifelse(stream == "Union River", 6.4, rkm_surveyed)]
redd_yr[ , abund_rkm := abund / rkm_surveyed]


redd_yr[ , year_fac := as.factor(year)]
redd_yr[ , stream_fac := as.factor(stream)]
redd_yr[ , stage := factor(stage, levels = c("before", "during", "after"))]
redd_yr[ , abund_stnd := scale(abund), by = .(stream)]
redd_yr[ , treatment := as.factor(treatment)]
redd_yr[ , abund_ln := log(abund + 1)]
redd_yr[ , abund_ln_stnd := scale(abund_ln), by = .(stream)]
redd_yr[ , abund_rkm_ln := log(abund_rkm + 1)]


redd_yr_stream = redd_yr[ , .(abund_mean = mean(abund),
                              abund_sd = sd(abund),
                              abund_rkm_mean = mean(abund_rkm),
                              abund_rkm_sd = sd(abund_rkm),
                              abund_ln_mean = mean(abund_ln),
                              abund_ln_sd = sd(abund_ln),
                              abund_rkm_ln_mean = mean(abund_rkm_ln),
                              abund_rkm_ln_sd = sd(abund_rkm_ln),
                              abund_stnd_mean = mean(abund_stnd),
                              abund_stnd_sd = sd(abund_stnd),
                              spawn_doy_mean = mean(spawn_doy, na.rm = TRUE),
                              spawn_doy_sd = sd(spawn_doy, na.rm = TRUE),
                              year_min = min(year),
                              year_max = max(year)), by = .(stream, treatment, stage)]
redd_yr_stream[ , stage := factor(stage, levels = c("before", "during", "after"))]

redd_yr_inter = redd_yr[ , .(abund_mean = mean(abund),
                             abund_sd = sd(abund),
                             abund_rkm_mean = mean(abund_rkm),
                             abund_rkm_sd = sd(abund_rkm),
                             abund_ln_mean = mean(abund_ln),
                             abund_ln_sd = sd(abund_ln),
                             abund_rkm_ln_mean = mean(abund_rkm_ln),
                             abund_rkm_ln_sd = sd(abund_rkm_ln),
                             abund_stnd_mean = mean(abund_stnd),
                             abund_stnd_sd = sd(abund_stnd)), by = .(treatment, stage)]
redd_yr_inter[ , stage := factor(stage, levels = c("before", "during", "after"))]

b1 = min(redd_yr$year[redd_yr$stage == "during"]) - 0.5
b2 = max(redd_yr$year[redd_yr$stage == "during"]) + 0.5
yrs_break = c(b1, b2)

## Save outputs
save(redd_yr, file = "./outputs/redd_yr.RData")
save(redd_yr_stream, file = "./outputs/redd_yr_stream.RData")
save(redd_yr_inter, file = "./outputs/redd_yr_inter.RData")
save(yrs_break, file = "./outputs/yrs_break.RData")



## Plots: abundance ts -------------------------------------

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    aes(x = year, y = abund, color = treatment) +
    geom_point() +
    geom_line() +
    labs(x = "Year", y = "Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_ts.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    aes(x = year, y = abund_rkm, color = treatment) +
    geom_point() +
    geom_line() +
    labs(x = "Year", y = "Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_rkm_ts.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    aes(x = year, y = abund_ln, color = treatment) +
    geom_point() +
    geom_line() +
    labs(x = "Year", y = "log Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_ln_ts.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    aes(x = year, y = abund_ln, color = treatment, groups = stream) +
    geom_point() +
    geom_line() +
    labs(x = "Year", y = "log Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_ln_ts_single.jpg", width = 6, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    aes(x = year, y = abund_stnd, color = treatment, groups = stream) +
    geom_point() +
    geom_line() +
    labs(x = "Year", y = "Abundance") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_stnd_ts.jpg", width = 8, height = 4)



## Plots: abundance ts + group means -----------------------
g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund, color = treatment), alpha = 0.8) +
    geom_segment(data = redd_yr_stream, linewidth = 1,
                 aes(x = year_min, xend = year_max, y = abund_mean, yend = abund_mean,
                     color = treatment)) +
    geom_rect(data = redd_yr_stream,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = abund_mean - abund_sd,
                  ymax = abund_mean + abund_sd, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "Abundance", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_avg_stage.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund_stnd, color = treatment), alpha = 0.8) +
    geom_segment(data = redd_yr_stream, linewidth = 1,
                 aes(x = year_min, xend = year_max, y = abund_stnd_mean, yend = abund_stnd_mean,
                     color = treatment)) +
    geom_rect(data = redd_yr_stream,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = abund_stnd_mean - abund_stnd_sd,
                  ymax = abund_stnd_mean + abund_stnd_sd, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "Abundance standardized", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_stnd_avg_stage.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund_ln, color = treatment), alpha = 0.8) +
    geom_segment(data = redd_yr_stream, linewidth = 1,
                 aes(x = year_min, xend = year_max, y = abund_ln_mean, yend = abund_ln_mean,
                     color = treatment)) +
    geom_rect(data = redd_yr_stream,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = abund_ln_mean - abund_ln_sd,
                  ymax = abund_ln_mean + abund_ln_sd, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "log Abundance", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream, scale = "free_y") +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_ln_avg_stage.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = abund_rkm_ln, color = treatment), alpha = 0.8) +
    geom_segment(data = redd_yr_stream, linewidth = 1,
                 aes(x = year_min, xend = year_max, y = abund_rkm_ln_mean, yend = abund_rkm_ln_mean,
                     color = treatment)) +
    geom_rect(data = redd_yr_stream,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = abund_rkm_ln_mean - abund_rkm_ln_sd,
                  ymax = abund_rkm_ln_mean + abund_rkm_ln_sd, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "log Abundance / rkm", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_rkm_ln_avg_stage.jpg", width = 8, height = 4)



## Plots: interaction --------------------------------------
g = ggplot(redd_yr_inter) +
    geom_point(aes(x = stage, y = abund_mean, color = treatment)) +
    geom_line(aes(x = stage, y = abund_mean, color = treatment, group = treatment)) +
    labs(x = "", y = "Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_inter.jpg", width = 6, height = 4)

g = ggplot(redd_yr_inter) +
    geom_point(aes(x = stage, y = abund_ln_mean, color = treatment)) +
    geom_line(aes(x = stage, y = abund_ln_mean, color = treatment, group = treatment)) +
    labs(x = "", y = "log Abundance", color = "Treatment") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/abundance_ln_inter.jpg", width = 6, height = 4)



## Plots: spawn timing -------------------------------------
g = ggplot(redd_yr) +
    aes(x = year, y = spawn_doy, color = treatment) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    labs(x = "Year", y = "Median spawn day", color = "Treatment") +
    scale_color_manual(values = M1) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/spawn_doy_ts.jpg", width = 8, height = 4)

g = ggplot(redd_yr) +
    geom_vline(xintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_vline(xintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(aes(x = year, y = spawn_doy, color = treatment), na.rm = TRUE) +
    geom_segment(data = redd_yr_stream, linewidth = 1, na.rm = TRUE,
                 aes(x = year_min, xend = year_max, y = spawn_doy_mean, yend = spawn_doy_mean,
                     color = treatment)) +
    geom_rect(data = redd_yr_stream,
              aes(xmin = year_min,
                  xmax = year_max,
                  ymin = spawn_doy_mean - spawn_doy_sd,
                  ymax = spawn_doy_mean + spawn_doy_sd, fill = treatment), alpha = 0.25) +
    labs(x = "Year", y = "Median spawn day", color = "Treatment", fill = "Treatment") +
    scale_color_manual(values = M1) +
    scale_fill_manual(values = M1) +
    facet_wrap( ~ stream) +
    theme_simple(grid = TRUE)
print(g)
ggsave("./figures/data/spawn_doy_avg_stage.jpg", width = 8, height = 4)

g = ggplot(redd[redds > 0, ]) +
    geom_hline(yintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_hline(yintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_point(data = redd[redds == 0, ], aes(x = doy, y = year), shape = 4, size = 0.2,
               color = "grey60") +
    geom_point(aes(x = doy, y = year, size = log(redds), color = log(redds)), alpha = 0.75, shape = 16) +
    scale_size_area(max_size = 3) +
    labs(x = "Day of year", y = "Year") +
    scale_color_viridis_c() +
    facet_wrap( ~ stream, ncol = 2) +
    theme_simple()
print(g)
ggsave("./figures/data/spawn_doy_bubble.jpg", width = 10, height = 10)

g = ggplot(redd) +
    geom_hline(yintercept = yrs_break[1], color = "grey50", linetype = 2) +
    geom_hline(yintercept = yrs_break[2], color = "grey50", linetype = 2) +
    geom_tile(aes(x = doy, y = year, fill = log(redds))) +
    labs(x = "Day of year", y = "Year") +
    scale_fill_viridis_c() +
    xlim(10, 200) +
    facet_wrap( ~ stream, ncol = 1) +
    theme_simple()
print(g)
ggsave("./figures/data/spawn_doy_tile.jpg", width = 10, height = 10)


## AR1 estimates -------------------------------------------
ar1 = redd_yr[ , .(phi = acf(abund, plot = FALSE)$acf[2]), by = .(stream)]
mean(ar1$phi)



## Response distribution -----------------------------------
g1 = ggplot(redd_yr) +
    geom_density(aes(x = abund), adjust = 1.2) +
    labs(x = "Abundance", y = "Density") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g1)
#
g2 = ggplot(redd_yr) +
    geom_density(aes(x = abund_ln), adjust = 1.2) +
    labs(x = "log Abundance", y = "Density") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g2)
#
g = g1 + g2
print(g)
ggsave("./figures/data/abundance_dist.jpg", width = 8, height = 4)


g1 = ggplot(redd_yr) +
    geom_density(aes(x = abund_rkm), adjust = 1.2) +
    labs(x = "Abundance / rkm", y = "Density") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g1)
#
g2 = ggplot(redd_yr) +
    geom_density(aes(x = abund_rkm_ln), adjust = 1.2) +
    labs(x = "log Abundance / rkm", y = "Density") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g2)
#
g = g1 + g2
print(g)
ggsave("./figures/data/abundance_rkm_dist.jpg", width = 8, height = 4)


g1 = ggplot(redd_yr[!is.na(spawn_doy), ]) +
    geom_density(aes(x = spawn_doy), adjust = 1.2) +
    labs(x = "Median Spawn Day of Year", y = "Density") +
    scale_color_manual(values = M1) +
    theme_simple(grid = TRUE)
print(g1)
ggsave("./figures/data/spawn_doy_density.jpg", width = 8, height = 4)
