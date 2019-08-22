library(pheatmap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

input_folder = args[1]
input_file_name = args[2]
ouput_name = args[3]#'G1E_vs_ER4.motifs.1d.all.K.pdf'
small_num = as.numeric(args[4])#100
k = as.numeric(args[5])#100
#order_used = c(16,11,18,15,2,0,3,4,6,5,20,7,13,9,19,10,24,17,25,8,1,21,23,14,22,12)+1


d12_dif_all = c()

for (pos_i in 1:k){
print(pos_i)
d1 = read.table(paste(input_folder, '/', input_file_name, '.kmeans.K', pos_i, '.1.txt', sep=''), header=F)
d2 = read.table(paste(input_folder, '/', input_file_name, '.kmeans.K', pos_i, '.2.txt', sep=''), header=F)
d12_dif = as.matrix(log2((d1+small_num)/(d2/sum(d2)*sum(d1)+small_num)))
d12_dif_all = cbind(d12_dif_all, d12_dif)
}


corlo_lim_scale = 3.5#max(max(abs((d12_dif_all[order_used,]))))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.1)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

colnames(d12_dif_all) = 1:k
rownames(d12_dif_all) = (1:dim(d12_dif_all)[1])-1
pdf(ouput_name, width=7.15, height=7)
pheatmap((d12_dif_all[,]),cluster_rows=T,cluster_cols=T, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()



