library(stringr)
get_res <- function(i){
    line <- readLines(paste0('mpl_', i, '.out'))
    foundset<-line[which(str_sub(line, 1, 1) == '[' )]
    char_mat <- str_split(foundset, ' ', simplify=TRUE)
    sapply(c(1:3), function(x) as.numeric(char_mat[x, 2]))
}

full_res <- function(vec){
    data.frame(t(sapply(vec, get_res)))
}

get_res_kl <- function(i){
    line <- readLines(paste0('kl_', i, '.out'))
    foundset<-line[which(str_sub(line, 1, 1) == '[' )]
    char_mat <- str_split(foundset, ' ', simplify=TRUE)
    sapply(c(1:3), function(x) as.numeric(char_mat[x, 2]))
}

full_res_kl <- function(vec){
    data.frame(t(sapply(vec, get_res_kl)))
}


classes <- function(tab){
    tot = dim(tab)[1]
    sbs = dim(tab[tab$X1 < tab$X2 & tab$X1 < tab$X3, ])[1]
    over = dim(tab[tab$X2 < tab$X1 & tab$X2 < tab$X3, ])[1]
    neut = tot - sbs - over
    c(sbs, over, neut)
}

classes_bal_only <- function(tab){
    tot = dim(tab)[1]
    sbs = dim(tab[tab$X1 < tab$X2, ])[1]
    over = tot - sbs
    c(sbs, over)
}