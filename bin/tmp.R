### get parameters
args = commandArgs(trailingOnly=TRUE)

state_file = args[1]
gene_withid = args[2]
outputfile_mat = args[3]
outputfile_transit_mat = args[4]
outputfile_wgtransit_mat = args[5]
win_num_start = as.numeric(args[6])
win_num_end = as.numeric(args[7])
type = args[8]
bin_length = as.numeric(args[9])

state_file = 'G1E.state.txt'
gene_withid = 'G1E_exp.gene.withbinID.txt'
outputfile_mat='G1E_exp.gene.withbinID.statemat.txt'
outputfile_transit_mat='G1E_exp.gene.withbinID.transmat.txt'
win_num_start = 200
win_num_end = 0
type = 'u'
bin_length = 50
step = 5
#time Rscript get_transition_matrix.R G1E.state.txt G1E_exp.gene.withbinID.txt G1E_exp.gene.withbinID.statemat.txt G1E_exp.gene.withbinID.transmat.txt G1E_wg.transmat.txt



get_transition_mat_seq = function(x, start, end, bin_length, step){
transition_vec_all = c()
#transition_mat = matrix(0, ncol=end-start+1, nrow=end-start+1)
s1_0 = x[-length(x)]
s2_0 = x[-1]
seq_len = length(s1_0)
for (h in seq(1, seq_len, step)){
if ((h+bin_length-1)<=seq_len){
### get tmp seq
s1 = s1_0[h:(h+bin_length-1)]
s2 = s2_0[h:(h+bin_length-1)]
transition_mat = matrix(0, ncol=end+1, nrow=end+1)
### get transit mat
for (i in start:end){
for (j in start:end){
count = sum((s1==i) & (s2==j))
transition_mat[i+1,j+1] = count
}
}
###
transition_vec_all = rbind(transition_vec_all, c(transition_mat))
}
}
return(transition_vec_all)
}

s = read.table(state_file, header=F)
ss = s[,2]
d1 = read.table(gene_withid, header=F)
d1 = d1[d1[,6]!='.',]
ids = as.numeric(as.character(d1[,6]))
d1[,6] = ids

state_mat = c()
bins_num = dim(s)[1]
#win_num = 50

for (i in 1:dim(d1)[1]){
info = d1[i,]
if (i%%1000==0){
print(i)
}
id_tmp = as.numeric(info[6])
if (((id_tmp-win_num_start)>0) & ((id_tmp-win_num_start)<bins_num)){
if (type=='ud'){
start = id_tmp-win_num_start
end = id_tmp+win_num_start
} else if (type=='u'){
start = id_tmp-win_num_start
end = id_tmp-win_num_end
} else if (type=='d'){
start = id_tmp+win_num_end
end = id_tmp+win_num_start
}
used_id = c(start:end)
if (info[5]=='+'){
state_mat = rbind(state_mat, ss[used_id])
} else{
state_mat = rbind(state_mat, ss[rev(used_id)])
}
}
}

#write.table(state_mat, outputfile_mat, sep='\t', quote=F, col.names=F, row.names=F)
dim(state_mat)

### get transition matrix
all_state = as.numeric(rownames(table(state_mat)))
start = min(all_state)
end = max(all_state)

#wg_transit_mat = get_transition_mat(ss, start, end)
#write.table(wg_transit_mat, outputfile_transit_mat, sep='\t', quote=F, col.names=T, row.names=T)

transit_mat_pool = c()
for (i in 1:dim(state_mat)[1]){
#for (i in 1:1000){
if (i%%100==0){
print(i)
print(dim(transit_mat_pool))
}
transit_mat_tmp = get_transition_mat_seq(state_mat[1,], start, end, bin_length, step)
transit_mat_pool = rbind(transit_mat_pool, transit_mat_tmp)
}

### clustering positions
all_counts_clust = kmeans(t(transit_mat_pool), centers=10)
all_counts_clust_clusters_id = all_counts_clust$cluster
all_state_kmeans = as.numeric(rownames(table(all_counts_clust_clusters_id)))
start_kmeans = min(all_state_kmeans)
end_kmeans = max(all_state_kmeans)

### write transition matrices
for (i in start_kmeans:end_kmeans){
transit_mat_pool_tmp = transit_mat_pool[all_counts_clust_clusters_id==i,]
transit_mat_pool_tmp = colSums(transit_mat_pool_tmp)
transit_mat_tmp = matrix(transit_mat_pool_tmp, ncol=sqrt(length(transit_mat_pool_tmp)), nrow=sqrt(length(transit_mat_pool_tmp)))
outputfile_transit_mat_tmp = paste('d12/', outputfile_transit_mat, '.kmeans.K', i, '.txt', sep='')
write.table(transit_mat_tmp, outputfile_transit_mat_tmp, sep='\t', quote=F, col.names=T, row.names=T)
}






