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
#state_file = 'G1E.state.txt'
#gene_withid = 'G1E_exp.gene.withbinID.txt'
#outputfile_mat='G1E_exp.gene.withbinID.statemat.txt'
#outputfile_transit_mat='G1E_exp.gene.withbinID.transmat.txt'
#time Rscript get_transition_matrix.R G1E.state.txt G1E_exp.gene.withbinID.txt G1E_exp.gene.withbinID.statemat.txt G1E_exp.gene.withbinID.transmat.txt G1E_wg.transmat.txt

get_transition_mat = function(x, start, end){
	transition_mat = matrix(0, ncol=end-start+1, nrow=end-start+1)
	s1 = x[-length(x)]
	s2 = x[-1]
	for (i in start:end){
	for (j in start:end){
	count = sum((s1==i) & (s2==j))
	transition_mat[i+1,j+1] = count
	}
	}
	return(transition_mat)
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

write.table(state_mat, outputfile_mat, sep='\t', quote=F, col.names=F, row.names=F)

### get transition matrix
all_state = as.numeric(rownames(table(state_mat)))
start = min(all_state)
end = max(all_state)

#wg_transit_mat = get_transition_mat(ss, start, end)
#write.table(wg_transit_mat, outputfile_transit_mat, sep='\t', quote=F, col.names=T, row.names=T)

transit_mat = get_transition_mat(state_mat[1,], start, end)

for (i in 2:dim(state_mat)[1]){
if (i%%1000==0){
print(i)
}
transit_mat_tmp = get_transition_mat(state_mat[i,], start, end)
transit_mat = transit_mat+transit_mat_tmp
}

colnames(transit_mat) = all_state
rownames(transit_mat) = all_state
write.table(transit_mat, outputfile_transit_mat, sep='\t', quote=F, col.names=T, row.names=T)




