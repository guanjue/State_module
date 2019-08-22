### get parameters
args = commandArgs(trailingOnly=TRUE)

gene_state_seqmat1 = args[1]#'G1Ehigh.G1E.gene.withbinID.statemat.500.txt'
gene_state_seqmat2 = args[2]#'G1Ehigh.ER4.gene.withbinID.statemat.500.txt'
gene_state_tranmat1 = args[3]#'G1Ehigh.G1E.gene.withbinID.tranmat.500.txt'
gene_state_tranmat2 = args[4]#'G1Ehigh.ER4.gene.withbinID.tranmat.500.txt'
bin_length = as.numeric(args[5])#50
step = as.numeric(args[6])#5


get_transition_mat_seq = function(x, start, end, bin_length, step){
transition_vec_all = c()
s1_0 = x
seq_len = length(s1_0)
for (h in seq(1, seq_len, step)){
if ((h+bin_length-1)<=seq_len){
### get tmp seq
s1 = s1_0[h:(h+bin_length-1)]
transition_vec = rep(0, end+1)
### get transit mat
for (i in start:end){
count = sum(s1==i)
transition_vec[i+1] = count
}
###
transition_vec_all = rbind(transition_vec_all, c(transition_vec))
}
}
return(transition_vec_all)
}


get_transition_mat_seq_pip = function(state_mat1, bin_length, step){
### get state number range
start = min(state_mat1)
end = max(state_mat1)
print(start)
print(end)
### initialize matrix
transit_mat_pool = matrix(0, ncol=(end+1), nrow=10000000)
### get seq
nrow_place1 = 1
for (i in 1:dim(state_mat1)[1]){
#for (i in 1:1000){
if (i%%10==0){
print(i)
#print(dim(transit_mat_pool))
}
transit_mat_tmp = get_transition_mat_seq(state_mat1[i,], start, end, bin_length, step)
nrow_len = dim(transit_mat_tmp)[1]
nrow_place2 = nrow_place1 + nrow_len -1
transit_mat_pool[nrow_place1:nrow_place2,] = transit_mat_tmp
nrow_place1 = nrow_place1 + nrow_len
}
### save used matrix
transit_mat_pool_filled = transit_mat_pool[1:nrow_place2,]
transit_mat_pool0 = transit_mat_pool_filled
return(transit_mat_pool0)
}


### get state seq
state_mat1 = as.matrix(read.table(gene_state_seqmat1, sep='\t', header=F))
state_mat2 = as.matrix(read.table(gene_state_seqmat2, sep='\t', header=F))
print(dim(state_mat1))
print(dim(state_mat2))
### get transit_mat seq
transit_mat_pool1 = get_transition_mat_seq_pip(state_mat1, bin_length, step)
transit_mat_pool2 = get_transition_mat_seq_pip(state_mat2, bin_length, step)

write.table(transit_mat_pool1, gene_state_tranmat1, sep='\t', quote=F, col.names=F, row.names=F)
write.table(transit_mat_pool2, gene_state_tranmat2, sep='\t', quote=F, col.names=F, row.names=F)

