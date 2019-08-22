### get parameters


win_num_start = 200
win_num_end = 0
type = 'u'
bin_length = 50
step = 5
add_smallnum = 1


gene_state_tranmat1 = 'G1Ehigh.G1E.gene.withbinID.tranmat.500.txt'
gene_state_tranmat2 = 'G1Ehigh.ER4.gene.withbinID.tranmat.500.txt'

pca_mat = 'G1Ehigh.G1E_vs_ER4.pca.mat.txt'
kmeans_decide_k = 'kmeans_ss.pairwise.pdf'
output_folder = 'kmeans_pairwise_motif'
outputfile_transit_mat = 'kmeans_motif'

### get transit_mat seq

transit_mat_pool1 = read.table(gene_state_tranmat1, sep='\t', header=F)
transit_mat_pool2 = read.table(gene_state_tranmat2, sep='\t', header=F)

transit_mat_pool_log2fc = log2((transit_mat_pool1+add_smallnum) / (transit_mat_pool2+add_smallnum))

transit_mat_pool_log2fc_pca0 = read.table(pca_mat, sep='\t', header=F)

set.seed(2019)
k=15
used_dim = 20
all_counts_clust = kmeans((transit_mat_pool_log2fc_pca0[,1:used_dim]), centers=k)
all_counts_clust_clusters_id = all_counts_clust$cluster
all_state_kmeans = as.numeric(rownames(table(all_counts_clust_clusters_id)))
start_kmeans = min(all_state_kmeans)
end_kmeans = max(all_state_kmeans)

library(pheatmap)
png('test.heatmap.png')
corlo_lim_scale = max(max(abs(transit_mat_pool_log2fc_pca0)))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue','white', 'red'))(n = length(breaksList_scale))
#my_colorbar_scale=colorRampPalette(c('blue','white', 'red'))(n = 100)

transit_mat_pool_log2fc_plot = transit_mat_pool_log2fc_pca0[order(all_counts_clust_clusters_id),1:used_dim][seq(1,dim(transit_mat_pool_log2fc)[1], 10),]
#d12_dif[d12_dif>corlo_lim_scale] = corlo_lim_scale
#pheatmap(transit_mat_pool_log2fc_pca$x[order(all_counts_clust_clusters_id),1:10],cluster_rows=F,cluster_cols=F,color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
#pheatmap(transit_mat_pool_log2fc_pca$x[order(all_counts_clust_clusters_id),1:10],cluster_rows=F,cluster_cols=F,color=my_colorbar_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
pheatmap(transit_mat_pool_log2fc_plot,color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows=F,cluster_cols=F,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()

### write transition matrices
dir.create(output_folder)
for (i in 1:k){
transit_mat_pool_tmp = transit_mat_pool_log2fc[all_counts_clust_clusters_id==i,]
transit_mat_pool_tmp = colSums(transit_mat_pool_tmp)
transit_mat_tmp = matrix(transit_mat_pool_tmp, ncol=sqrt(length(transit_mat_pool_tmp)), nrow=sqrt(length(transit_mat_pool_tmp)))
outputfile_transit_mat_tmp = paste(output_folder, '/', outputfile_transit_mat, '.kmeans.K', i, '.txt', sep='')
write.table(transit_mat_tmp, outputfile_transit_mat_tmp, sep='\t', quote=F, col.names=T, row.names=T)
}





library(pheatmap)

order_used = c(16,11,18,15,2,0,3,4,6,5,20,7,13,9,19,10,24,17,25,8,1,21,23,14,22,12)+1

k=15
pdf('G1E_vs_ER4.motifs.K.pdf', width=7.15, height=7)

for (pos_i in 1:k){
print(pos_i)
d1 = read.table(paste('kmeans_pairwise_motif/kmeans_motif.kmeans.K', pos_i, '.txt', sep=''))
d2 = read.table('G1E_noexp.gene.withbinID.transmat.50u.txt')

small_num = 1000
d12_dif = d1
#log2((d1+small_num)/(d2/sum(d2)*sum(d1)+small_num))

#diag(d12_dif) = 0
colnames(d12_dif) = 0:(dim(d12_dif)[1]-1)
rownames(d12_dif) = 0:(dim(d12_dif)[1]-1)
#order_used = c(16,11,18,15,2,0,3,4,6,5,20,7,13,9,19,10,24,17,25,8,1,21,23,14,22,12)+1

corlo_lim_scale = max(max(abs(d12_dif)))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))
d12_dif[d12_dif>corlo_lim_scale] = corlo_lim_scale

pheatmap(d12_dif[order_used,order_used],cluster_rows=T,cluster_cols=T, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)

}

dev.off()



