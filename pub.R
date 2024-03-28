## Publication figures


dir.create("./figures/pub/", showWarnings = FALSE)


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
