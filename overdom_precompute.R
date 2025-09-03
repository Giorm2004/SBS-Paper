source('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/ptd_sfs.R')

args <- commandArgs(trailingOnly = TRUE)
k = as.numeric(args[1])
p = as.numeric(args[2])
peq = as.numeric(args[3])
rho = as.numeric(args[4])

tab = t(norm_sfs_over(k, p, peq, rho))
write.table(tab, file = paste0("over_sfs_", k, "_",p ,"_", peq, "_", rho, ".out" ), row.names=FALSE, col.names=FALSE)

