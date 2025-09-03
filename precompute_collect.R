sbstab <- read.table('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/sbs_sfs_control.txt', col.names = c("k", "p", "ph", "rho"), row.names=NULL)
overtab <- read.table('/scratch/nzx3cc/nzx3cc/sbs_paper/scripts/over_sfs_control.txt', col.names = c("k", "p", "ph", "rho"), row.names=NULL)
get_sbs_tab <- function(){
    setwd('/home/nzx3cc/sbs_precompute')
    t_20= data.frame(t(sapply(seq(1,dim(sbstab)[1], 2), function(i) as.numeric(read.table(paste0("sbs_sfs_", sbstab[i,1], "_",sbstab[i,2] ,"_", sbstab[i,3], "_", sbstab[i,4], ".out" ), header=FALSE)))))
    params_20 = data.frame(sbstab[seq(1,dim(sbstab)[1], 2),])
    tab_20 = do.call(cbind, list(params_20, t_20))
    t_25= data.frame(t(sapply(seq(2,dim(sbstab)[1], 2), function(i) as.numeric(read.table(paste0("sbs_sfs_", sbstab[i,1], "_",sbstab[i,2] ,"_", sbstab[i,3], "_", sbstab[i,4], ".out" ), header=FALSE)))))
    params_25 = data.frame(sbstab[seq(2,dim(sbstab)[1], 2),])
    tab_25 = do.call(cbind, list(params_25, t_25))
    write.table(tab_20, file="sbs_sfss_20.out", row.names=FALSE, col.names=TRUE)
    write.table(tab_25, file="sbs_sfss_25.out", row.names=FALSE, col.names=TRUE)
}

get_over_tab <- function(){
    setwd('/home/nzx3cc/over_precompute')
    t_20= data.frame(t(sapply(seq(1,dim(overtab)[1], 2), function(i) as.numeric(read.table(paste0("over_sfs_", overtab[i,1], "_",overtab[i,2] ,"_", overtab[i,3], "_", overtab[i,4], ".out" ), header=FALSE)))))
    params_20 = data.frame(overtab[seq(1,dim(overtab)[1], 2),])
    tab_20 = do.call(cbind, list(params_20, t_20))
    t_25= data.frame(t(sapply(seq(2,dim(overtab)[1], 2), function(i) as.numeric(read.table(paste0("over_sfs_", overtab[i,1], "_",overtab[i,2] ,"_", overtab[i,3], "_", overtab[i,4], ".out" ), header=FALSE)))))
    params_25 = data.frame(overtab[seq(2,dim(overtab)[1], 2),])
    tab_25 = do.call(cbind, list(params_25, t_25))
    write.table(tab_20, file="over_sfss_20.out", row.names=FALSE, col.names=TRUE)
    write.table(tab_25, file="over_sfss_25.out", row.names=FALSE, col.names=TRUE)
}

