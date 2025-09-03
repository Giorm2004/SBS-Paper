source('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/precompute_inference.R')
source('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/inference.R')
args <- commandArgs(trailingOnly = TRUE)
filename = toString(args[1])
k = as.numeric(args[2])
p = as.numeric(args[3])
rho = as.numeric(args[4:length(args)])

# following unlist should be removed, as will cause problems for multiple windows
obs <- read.table(filename)
if (k == 20){
    s = s_20
    o = o_20
}
if (k == 25){
    s = s_25
    o = o_25
}
sbs = best_lookup_mpl(obs, s, p, rho)
over = best_lookup_mpl(obs, o, p, rho)
neut = neut_mpl(obs, k)
#print(sbs$value)
#print(over$value)
#print(neut)
obs <- unlist(obs)
sbs2 = optim(round(sbs$par,2), get_complike_sbs, obs = obs, p_vec = p, k = k, rho = rho, method="Brent", lower=(sbs$par-0.0075), upper=(sbs$par + 0.0075), control = list(reltol = 1e-2))
over2 = optim(round(over$par,2), get_complike_over, obs = obs, p_vec = p, k = k, rho = rho, method="Brent", lower=(over$par-0.0075), upper=(over$par + 0.0075), control = list(reltol = 1e-2))
print(sbs2$value)
print(over2$value)
print(neut)