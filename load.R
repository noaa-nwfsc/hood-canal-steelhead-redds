## Load packages


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(emmeans))
suppressPackageStartupMessages(library(bayesplot))
library(overlapping)
library(patchwork)
library(ggsimple) ## see: https://github.com/michaelmalick/ggsimple

dir.create("./figures/", showWarnings = FALSE)
dir.create("./outputs/", showWarnings = FALSE)

source("./functions.R")

M1 = c("#239DD7", "#CD2343", "#C28500", "#8E41BB","#000000", "#949494")

## Load output data
for(i in list.files(path = "./outputs/", pattern = "*.RData$")) {
    load(paste("./outputs/", i, sep = ""))
}
