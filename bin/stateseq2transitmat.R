### get parameters
args = commandArgs(trailingOnly=TRUE)

state_seq_file = args[1]
outputfile_transit_mat = args[2]
startseq = as.numeric(args[3])
endseq = as.numeric(args[4])

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
	#print(i)
	for (j in start:end){
	count = sum((s1==i) & (s2==j))
	transition_mat[i+1,j+1] = count
	}
	}
	return(transition_mat)
}

state_mat = read.table(state_seq_file, header=F)[,startseq:endseq]

### get transition matrix
all_state = as.numeric(rownames(table(state_mat[,1:5])))
start = min(all_state)
end = max(all_state)

transit_mat = get_transition_mat(state_mat[1,], start, end)

for (i in 2:dim(state_mat)[1]){
if (i%%100==0){
print(i)
}
transit_mat_tmp = get_transition_mat(state_mat[i,], start, end)
transit_mat = transit_mat+transit_mat_tmp
}

colnames(transit_mat) = all_state
rownames(transit_mat) = all_state
write.table(transit_mat, outputfile_transit_mat, sep='\t', quote=F, col.names=T, row.names=T)




