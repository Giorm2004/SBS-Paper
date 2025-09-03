source('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/ptd_sfs.R')

args <- commandArgs(trailingOnly = TRUE)
k = as.numeric(args[1])
p = as.numeric(args[2])
ph = as.numeric(args[3])
rho = as.numeric(args[4])

tab = t(norm_sfs_sbs(k, p, ph, rho))
write.table(tab, file = paste0("sbs_sfs_", k, "_",p ,"_", ph, "_", rho, ".out" ), row.names=FALSE, col.names=FALSE)

