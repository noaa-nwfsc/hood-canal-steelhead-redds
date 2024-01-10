## Functions for analysis


## epred_df ------------------------------------------------
epred_df = function(epred, newdata) {
    ## Convert epred output to a dataframe
    ##
    ## epred = output from brms::posterior_epred
    ## newdata = data used as input to brms::posterior_epred

    epred = as.data.table(epred)
    newdata = as.data.table(newdata)

    lst = vector("list", nrow(newdata))
    for(i in seq_along(lst)) {
        dt = as.data.table(newdata[i, ])
        pr = epred[ , ..i]
        names(pr) = "sample"
        lst[[i]] = cbind(dt, pr)
    }
    out = rbindlist(lst)

    return(out)
}



## diagnostic checks ---------------------------------------
rhat_highest_dfa <- function(dfa, k = 4, pars = c("sigma", "x\\[", "Z\\[")) {

    rhat <- dfa$monitor[which(grepl(paste(pars,
        collapse = "|"), rownames(dfa$monitor)) == TRUE),
        "Rhat"]
    rhat.max <- rev(sort(rhat)[(length(rhat) - k):length(rhat)])
    return(rhat.max)
}

neff_lowest_dfa <- function(dfa, k = 4, pars = c("sigma", "x\\[", "Z\\[")) {
    neff <- dfa$monitor[which(grepl(paste(pars,
        collapse = "|"), rownames(dfa$monitor)) == TRUE),
        "n_eff"]
    neff.min <- sort(neff)[1:k]
    return(neff.min)
}

rhat_highest <- function(stanfit, k = 4, pars) {
    rhat <- get_rhat(stanfit, pars = pars)
    rhat.max <- rev(sort(rhat)[(length(rhat) - k):length(rhat)])
    return(rhat.max)
}

neff_lowest <- function(stanfit, k = 4, pars) {
    neff <- get_neff(stanfit, pars = pars)
    neff.min <- sort(neff)[1:k]
    return(neff.min)
}

pairs_lowest <- function(stanfit, k = 4, pars) {
    n <- get_neff(stanfit, pars = pars)
    n.min <- names(sort(n))[1:k]
    pairs(stanfit, pars = n.min)
}

get_rhat <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    na.omit(summary(stanfit, pars = pars)$summary[ , "Rhat"])
}

get_neff <- function(stanfit, pars) {
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    summary(stanfit, pars = pars)$summary[ , "n_eff"]
}

total_draws <- function(stanfit) {
    ## N chains * N draws  -- post warmup
    if(!is(stanfit, "stanfit")) {
        stop("Input not of class stanfit")
    }
    dim(stanfit)[1] * dim(stanfit)[2]
}


