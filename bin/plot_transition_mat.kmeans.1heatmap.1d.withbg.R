library(pheatmap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

input_folder = args[1]
input_file_name = args[2]
count_mat = args[3]
ouput_name = args[4]#'G1E_vs_ER4.motifs.1d.all.K.pdf'
small_num = as.numeric(args[5])#100
k = as.numeric(args[6])#100
kcount = as.numeric(args[7])
order_used = c(16,11,18,15,2,0,3,4,6,5,20,7,13,9,19,10,24,17,25,8,1,21,23,14,22,12)+1


d12_dif_all = c()

for (pos_i in 1:k){
print(pos_i)
d1 = read.table(paste(input_folder, '/', input_file_name, '.kmeans.K', pos_i, '.1.txt', sep=''), header=F)
d2 = read.table(paste(input_folder, '/', input_file_name, '.kmeans.K', pos_i, '.2.txt', sep=''), header=F)
d12_dif = as.matrix(log2((d1+small_num)/(d2/sum(d2)*sum(d1)+small_num)))
d12_dif_all = cbind(d12_dif_all, d12_dif)
}

bg = read.table(count_mat, header=F)
bg_colmeans = colMeans(bg)
bg = t(apply(bg, 1, function(x) x/bg_colmeans*bg_colmeans[1]))
bg[,3] = bg[,3]/2+bg[,4]/2
bg = bg[,1:3]
print(bg)
small_num_bg = 0
bglog2fc = t(apply(bg, 1, function(x) log2((x[1]+small_num_bg)/(x[2:length(x)]+small_num_bg))))

print(dim(bglog2fc))

rownames(bglog2fc) = paste(c(1:dim(bglog2fc)[1]), 'k', sep='_')
set.seed(2019)
bg_kmeans = kmeans(bglog2fc, centers=kcount)
bg_kmeans_clusters_id = bg_kmeans$cluster


corlo_lim_scale = 3.5#max(max(abs((d12_dif_all[order_used,]))))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.01)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

colnames(d12_dif_all) = 1:k
rownames(d12_dif_all) = (1:dim(d12_dif_all)[1])-1
pdf(ouput_name, width=7, height=7)
pheatmap(t(d12_dif_all[order_used,order(rowMeans(bglog2fc), decreasing=T)]),cluster_rows=F,cluster_cols=T, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=F,show_colnames=TRUE,annotation_names_row=F,annotation_names_col=TRUE)
dev.off()


print(summary(bglog2fc))
#print(bglog2fc)
corlo_lim_scale = max(abs(bglog2fc))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.01)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

pdf(paste(ouput_name, '.count.pdf', sep=''), width=2, height=7)
pheatmap((bglog2fc[order(rowMeans(bglog2fc), decreasing=T),]),cluster_rows=F,cluster_cols=F, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()


