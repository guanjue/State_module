library(pheatmap)

d1 = read.table('G1E_exp.gene.withbinID.transmat.txt')
d2 = read.table('G1E_noexp.gene.withbinID.transmat.txt')
d3 = read.table('G1E_wg.transmat.txt')

small_num = 500
d12 = log2((d1+d2+small_num)/(d3/sum(d3)*sum(d1+d2)+small_num))
d12_dif = log2((d1+small_num)/(d2/sum(d2)*sum(d1)+small_num))

#diag(d12) = 1
#diag(d12_dif) = 0
colnames(d12) = rownames(d12)
colnames(d12_dif) = rownames(d12_dif)


corlo_lim_scale = max(max(abs(d12_dif)))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

pdf('d12_dif.pdf', width=7.15, height=7)
pheatmap(d12_dif,cluster_rows=F,cluster_cols=F, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()

corlo_lim_scale = max(max(abs(d12)))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

pdf('d12.pdf', width=7.15, height=7)
pheatmap(d12,cluster_rows=F,cluster_cols=F,color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()

