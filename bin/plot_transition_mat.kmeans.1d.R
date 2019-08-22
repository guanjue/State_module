library(pheatmap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

gene_state_tranmat1 = args[1]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat2 = args[2]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
kmeans_decide_k = args[3]#'kmeans_ss.pairwise.1d.pdf'
K_try = as.numeric(args[4])
add_smallnum = as.numeric(args[5])


order_used = c(16,11,18,15,2,0,3,4,6,5,20,7,13,9,19,10,24,17,25,8,1,21,23,14,22,12)+1

k=20
d12_dif_all = c()
pdf('G1E_vs_ER4.motifs.1d.K.pdf', width=7.15, height=7)
for (pos_i in 1:k){
print(pos_i)
d1 = read.table(paste('kmeans_pairwise_motif_1d/kmeans_motif_1d.kmeans.K', pos_i, '.1.txt', sep=''), header=F)
d2 = read.table(paste('kmeans_pairwise_motif_1d/kmeans_motif_1d.kmeans.K', pos_i, '.2.txt', sep=''), header=F)

small_num = 100
d12_dif = as.matrix(log2((d1+small_num)/(d2/sum(d2)*sum(d1)+small_num)))
d12_dif_all = cbind(d12_dif_all, d12_dif)
d12_dif = cbind(d12_dif, d12_dif)
print(max(d12_dif))
print(min(d12_dif))

#diag(d12_dif) = 0
colnames(d12_dif) = 1:2
rownames(d12_dif) = 0:(dim(d12_dif)[1]-1)

corlo_lim_scale = 3.5#max(max(abs(d12_dif)))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))
d12_dif[d12_dif>corlo_lim_scale] = corlo_lim_scale

pheatmap(d12_dif[order_used,],cluster_rows=F,cluster_cols=F, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)

}

dev.off()

corlo_lim_scale = 3.5#max(max(abs((d12_dif_all[order_used,]))))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.2)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

colnames(d12_dif_all) = 1:k
rownames(d12_dif_all) = (1:dim(d12_dif_all)[1])-1
pdf('G1E_vs_ER4.motifs.1d.all.K.pdf', width=7.15, height=7)
pheatmap((d12_dif_all[order_used,]),cluster_rows=F,cluster_cols=T, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()



