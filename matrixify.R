#####build matrixify function####
matrixify=function(stat_list,ivs){
stat_matrix <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix[k] <- as.double()
for (i in 1:length(stat_list)){stat_matrix[i,] <- stat_list[[i]]}
for (n in 1:length(ivs)) row.names(stat_matrix)[n] <- ivs[n]
return(stat_matrix)
}