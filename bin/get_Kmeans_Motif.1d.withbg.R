### get parameters
args = commandArgs(trailingOnly=TRUE)

gene_state_tranmat1 = args[1]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat2 = args[2]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat_bg1 = args[3]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat_bg2 = args[4]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'

output_folder = args[5]#'kmeans_pairwise_motif_1d'
outputfile_transit_mat = args[6]#'kmeans_motif_1d'
heatmap_name = args[7]#'test.heatmap.1d.png'
add_smallnum = as.numeric(args[8])
k = as.numeric(args[9])

### get transit_mat seq
transit_mat_pool1 = read.table(gene_state_tranmat1, sep='\t', header=F)
transit_mat_pool2 = read.table(gene_state_tranmat2, sep='\t', header=F)
transit_mat_bg1 = read.table(gene_state_tranmat_bg1, sep='\t', header=F)
transit_mat_bg2 = read.table(gene_state_tranmat_bg2, sep='\t', header=F)

transit_mat_pool_bg1 = rbind(transit_mat_pool1, transit_mat_pool2, transit_mat_bg1, transit_mat_bg2)
transit_mat_pool_bg2 = rbind(transit_mat_pool2, transit_mat_pool1, transit_mat_bg2, transit_mat_bg1)


### get transit_mat seq
transit_mat_pool_high_log2fc = log2((transit_mat_pool1+add_smallnum) / (transit_mat_pool2+add_smallnum))
transit_mat_pool_low_log2fc = log2((transit_mat_pool2+add_smallnum) / (transit_mat_pool1+add_smallnum))
transit_mat_bg1_log2fc = log2((transit_mat_bg1+add_smallnum) / (transit_mat_bg2+add_smallnum))
transit_mat_bg2_log2fc = log2((transit_mat_bg2+add_smallnum) / (transit_mat_bg1+add_smallnum))

transit_mat_pool_bg_log2fc = rbind(transit_mat_pool_high_log2fc, transit_mat_pool_low_log2fc, transit_mat_bg1_log2fc, transit_mat_bg2_log2fc)

set.seed(2019)

all_counts_clust = kmeans((transit_mat_pool_bg_log2fc), centers=k)
all_counts_clust_clusters_id = all_counts_clust$cluster
all_state_kmeans = as.numeric(rownames(table(all_counts_clust_clusters_id)))
start_kmeans = min(all_state_kmeans)
end_kmeans = max(all_state_kmeans)

library(pheatmap)
png(heatmap_name)
corlo_lim_scale = max(max(abs(transit_mat_pool_bg_log2fc)))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue','white', 'red'))(n = length(breaksList_scale))
transit_mat_pool_log2fc_plot = transit_mat_pool_bg_log2fc[order(all_counts_clust_clusters_id),][seq(1,dim(transit_mat_pool_bg_log2fc)[1], 3),]
pheatmap(transit_mat_pool_log2fc_plot,color=my_colorbar_scale, breaks = breaksList_scale, cluster_rows=F,cluster_cols=F,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()

### write transition matrices
dir.create(output_folder)
for (i in 1:k){
print(i)
### 1st
transit_mat_pool_tmp = (transit_mat_pool_bg1)[all_counts_clust_clusters_id==i,]
print(dim(transit_mat_pool_tmp))
transit_mat_pool_tmp = colSums(transit_mat_pool_tmp)
outputfile_transit_mat_tmp = paste(output_folder, '/', outputfile_transit_mat, '.kmeans.K', i, '.1.txt', sep='')
write.table(transit_mat_pool_tmp, outputfile_transit_mat_tmp, sep='\t', quote=F, col.names=F, row.names=F)
### 2nd
transit_mat_pool_tmp = (transit_mat_pool_bg2)[all_counts_clust_clusters_id==i,]
print(dim(transit_mat_pool_tmp))
transit_mat_pool_tmp = colSums(transit_mat_pool_tmp)
outputfile_transit_mat_tmp = paste(output_folder, '/', outputfile_transit_mat, '.kmeans.K', i, '.2.txt', sep='')
write.table(transit_mat_pool_tmp, outputfile_transit_mat_tmp, sep='\t', quote=F, col.names=F, row.names=F)
}


counts_mat_high = table(all_counts_clust_clusters_id[1:dim(transit_mat_pool_high_log2fc)[1]])
counts_mat_low = table(all_counts_clust_clusters_id[(1+dim(transit_mat_pool_high_log2fc)[1]) : (dim(transit_mat_pool_high_log2fc)[1]+dim(transit_mat_pool_low_log2fc)[1])])
counts_mat_bg1 = table(all_counts_clust_clusters_id[(1+dim(transit_mat_pool_high_log2fc)[1]+dim(transit_mat_pool_low_log2fc)[1]) : (dim(transit_mat_pool_high_log2fc)[1]+dim(transit_mat_pool_low_log2fc)[1]+dim(transit_mat_bg1_log2fc)[1])])
counts_mat_bg2 = table(all_counts_clust_clusters_id[(1+dim(transit_mat_pool_high_log2fc)[1]+dim(transit_mat_pool_low_log2fc)[1]+dim(transit_mat_bg1_log2fc)[1]) : (dim(transit_mat_pool_high_log2fc)[1]+dim(transit_mat_pool_low_log2fc)[1]+dim(transit_mat_bg1_log2fc)[1]+dim(transit_mat_bg2_log2fc)[1])])

counts_mat = cbind(counts_mat_high, counts_mat_low, counts_mat_bg1, counts_mat_bg2)

write.table(counts_mat, paste(outputfile_transit_mat, '.counts.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')




