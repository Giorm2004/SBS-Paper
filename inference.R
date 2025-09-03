source('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/ptd_sfs.R')

comp_kl <- function(obs, ex){
    vec <- sapply(1:length(obs), function(i) ifelse(obs[i]==0, 0,obs[i]*log(obs[i]/ex[i])))
    sum(vec)}
comp_like <- function(obs, ex){
    vec <- sapply(1:length(obs), function(i) -1 * obs[i] * log(ex[i]))
    sum(vec)
}

comp_2kl <- function(obs, ex){
abs(obs[1] * log(obs[1]/ex[1])) + abs(obs[2] * log(obs[2]/ex[2]))
}

time_series_sfs_sbs <- function(k, p_vec, ph, rho){
as.vector(sapply(p_vec, function(i) norm_sfs_sbs(k, i, ph, rho)))
}
time_series_sfs_over <- function(k, p_vec, eq, rho){
as.vector(sapply(p_vec, function(i) norm_sfs_over(k, i, eq, rho)))
}

get_complike_sbs <- function(ph, obs, p_vec, k, rho){
comp_like(obs, time_series_sfs_sbs(k, p_vec, ph, rho))
}
get_complike_over <- function(eq, obs, p_vec, k, rho){
    comp_like(obs, time_series_sfs_over(k, p_vec, eq, rho))
}
get_kl_sbs <- function(ph, obs, p_vec, k, rho){
comp_kl(obs, time_series_sfs_sbs(k, p_vec, ph, rho))
}
get_kl_over <- function(eq, obs, p_vec, k, rho){
    comp_kl(obs, time_series_sfs_over(k, p_vec, eq, rho))
}

get_2kl_sbs <- function(ph, obs, p_vec, k, rho){
comp_2kl(obs, time_series_sfs_sbs(k, p_vec, ph, rho))
}
get_2kl_over <- function(eq, obs, p_vec, k, rho){
    comp_2kl(obs, time_series_sfs_over(k, p_vec, eq, rho))
}
get_mpl_sbs <- function(obs, p_vec, k, rho){
    optim(0.4, get_complike_sbs, obs = obs, p_vec = p_vec, k = k, rho = rho, method="Brent", lower=0, upper=0.45, control = list(reltol = 1e-3))
}

get_mpl_over <- function(obs, p_vec, k, rho){
optim(0.5, get_complike_over, obs = obs, p_vec = p_vec, k = k, rho = rho, method="Brent", lower=0.1, upper=0.9, control = list(reltol = 1e-3))
}

best_kl_sbs <- function(obs, p_vec, k, rho){
    optim(0.4, get_kl_sbs, obs = obs, p_vec = p_vec, k = k, rho = rho, method="Brent", lower=0, upper=0.45, control = list(reltol = 1e-3))
}

best_kl_over <- function(obs, p_vec, k, rho){
optim(0.5, get_kl_over, obs = obs, p_vec = p_vec, k = k, rho = rho, method="Brent", lower=0.1, upper=0.9, control = list(reltol = 1e-3))
}

best_2kl_sbs <- function(obs, p, k, rho){
    optim(0.4, get_2kl_sbs, obs = obs, p_vec = p, k = k, rho = rho, method="Brent", lower=0, upper=0.45, control = list(reltol = 1e-3))
}

best_2kl_over <- function(obs, p, k, rho){
    optim(0.4, get_2kl_over, obs = obs, p_vec = p, k = k, rho = rho, method="Brent", lower=0.1, upper=0.9, control = list(reltol = 1e-3))
}