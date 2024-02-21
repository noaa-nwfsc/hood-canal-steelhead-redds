## Reproduce project
##
## ** make.R last run on 2024-02-21 **

unlink("./figures", recursive = TRUE)
unlink("./outputs", recursive = TRUE)

rm(list = ls())

cat("Sourcing load.R ...", "\n");     source("./load.R")
cat("Sourcing data.R ...", "\n");     source("./data.R")
cat("Sourcing fit.R ...", "\n");      source("./fit.R")
cat("Done!", "\n")

writeLines(capture.output(sessionInfo()), "session.txt")
