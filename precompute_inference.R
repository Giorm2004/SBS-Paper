s_20 <- read.table('/home/nzx3cc/sbs_precompute/sbs_sfss_20.out', header=TRUE)
s_25 <- read.table('/home/nzx3cc/sbs_precompute/sbs_sfss_25.out', header=TRUE)
o_20 <- read.table('/home/nzx3cc/over_precompute/over_sfss_20.out', header=TRUE)
o_25 <- read.table('/home/nzx3cc/over_precompute/over_sfss_25.out', header=TRUE)
source('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/ptd_sfs.R')

window_mpl <- function(obs, ex){
   ex = unlist(ex)
   obs = unlist(obs)
   vec = sapply(c(1:length(obs)), function(i) -1 * obs[i] * log(ex[i]))
   sum(unlist(vec))
}

lookup_mpl <- function(ph, tab, obs_vec){
    #assumes windows are always sorted close to far
   ph1 <- round(ph, 2)
  # print(ph1)
   rows = tab[tab$ph == ph1, ]
   #print("rows:")
  # print(rows)
   ex = rows[,-c(1:4)]
   #print("ex:")
   #print(ex)
   vec = sapply(c(1:(dim(obs_vec)[1])), function(i) window_mpl(obs_vec[i, ], unlist(ex[i, ])))
   sum(vec)
}

best_lookup_mpl <- function(obs, tab, p, rho_vec){
    #assumes windows are always sorted close to far
    p1 = round(p, 2)
    subtab = tab[tab$p == p1 & tab$rho %in% rho_vec,]
    #print(subtab)
    max = max(subtab$ph)
    min = min(subtab$ph)
    optim(0.4, lookup_mpl, obs = obs, tab = subtab, method="Brent", lower=min, upper=max, control = list(reltol = 1e-3))
}

neut_mpl <- function(obs_vec, k){
    neut = nnsfs(k)
    vec = sapply(c(1:(dim(obs_vec)[1])), function(i) window_mpl(obs_vec[i, ], neut))
    sum(vec)
}