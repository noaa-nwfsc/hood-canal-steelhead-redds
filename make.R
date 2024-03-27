## Reproduce project
##
## ** make.R last run on 2024-03-27 **

unlink("./figures", recursive = TRUE)
unlink("./outputs", recursive = TRUE)

rm(list = ls())

cat("Sourcing load.R ...", "\n");      source("./load.R")
cat("Sourcing data.R ...", "\n");      source("./data.R")
cat("Sourcing fit-abund.R ...", "\n"); source("./fit-abund.R")
cat("Sourcing fit-spawn.R ...", "\n"); source("./fit-spawn.R")
cat("Done!", "\n")

writeLines(capture.output(sessionInfo()), "session.txt")
