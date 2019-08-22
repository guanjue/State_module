### get parameters
args = commandArgs(trailingOnly=TRUE)

state_file = args[1]
gene_withid = args[2]
win_num = as.numeric(args[3])
outputfile_mat = args[4]

#state_file = 'G1E.state.txt'
#gene_withid = 'G1E_exp.gene.withbinID.txt'
#outputfile_mat='G1E_exp.gene.withbinID.statemat.txt'
#outputfile_transit_mat='G1E_exp.gene.withbinID.transmat.txt'
#time Rscript get_transition_matrix.R G1E.state.txt G1E_exp.gene.withbinID.txt G1E_exp.gene.withbinID.statemat.txt G1E_exp.gene.withbinID.transmat.txt G1E_wg.transmat.txt

s = read.table(state_file, header=F)
ss = s[,2]
d1 = read.table(gene_withid, header=F)
d1 = d1[d1[,6]!='.',]
ids = as.numeric(as.character(d1[,6]))
d1[,6] = ids

state_mat = c()
bins_num = dim(s)[1]

for (i in 1:dim(d1)[1]){
info = d1[i,]
if (i%%1000==0){
print(i)
}
id_tmp = as.numeric(info[6])
if (((id_tmp-win_num)>0) & ((id_tmp-win_num)<bins_num)){
start = id_tmp-win_num
end = id_tmp+win_num
used_id = c(start:end)
ss_seq = ss[used_id]
if (sum(is.na(ss_seq))==0){
if (info[5]=='+'){
state_mat = rbind(state_mat, ss_seq)
} else{
state_mat = rbind(state_mat, rev(ss_seq))
}
}
}
}

write.table(state_mat, outputfile_mat, sep='\t', quote=F, col.names=F, row.names=F)

