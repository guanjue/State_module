### get parameters
args = commandArgs(trailingOnly=TRUE)

gene_state_tranmat1 = args[1]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat2 = args[2]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat_bg1 = args[3]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat_bg2 = args[4]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
para_file = args[5]#../ideasM2/ideas_result2/ideas_m2.para.modified.para

kmeans_decide_k = args[6]#'kmeans_ss.pairwise.1d.pdf'
K_try = as.numeric(args[7])
add_smallnum = as.numeric(args[8])

### get transit_mat seq
transit_mat_pool1 = read.table(gene_state_tranmat1, sep='\t', header=F)
transit_mat_pool2 = read.table(gene_state_tranmat2, sep='\t', header=F)
transit_mat_bg1 = read.table(gene_state_tranmat_bg1, sep='\t', header=F)
transit_mat_bg2 = read.table(gene_state_tranmat_bg2, sep='\t', header=F)
para_signal = read.table(para_file, header=T)

### get para mean signal
k = dim(para_signal)[2]
p = (sqrt(9+8*(k-1))-3)/2
para_all_signal = para_signal[,2:p]
para_mean_signal = para_all_signal / para_signal[,1]
para_mean_sq_signal = rowSums(para_mean_signal^2)

### transit_mat_signal
transit_mat_pool1_signal = t(t(transit_mat_pool1) * para_mean_sq_signal)
transit_mat_pool2_signal = t(t(transit_mat_pool2) * para_mean_sq_signal)
transit_mat_bg1_signal = t(t(transit_mat_bg1) * para_mean_sq_signal)
transit_mat_bg2_signal = t(t(transit_mat_bg2) * para_mean_sq_signal)

### log2fc
transit_mat_pool_high_log2fc = log2((transit_mat_pool1_signal+add_smallnum) / (transit_mat_pool2_signal+add_smallnum))
transit_mat_pool_low_log2fc = log2((transit_mat_pool2_signal+add_smallnum) / (transit_mat_pool1_signal+add_smallnum))
transit_mat_bg1_log2fc = log2((transit_mat_bg1_signal+add_smallnum) / (transit_mat_bg2_signal+add_smallnum))
transit_mat_bg2_log2fc = log2((transit_mat_bg2_signal+add_smallnum) / (transit_mat_bg1_signal+add_smallnum))

transit_mat_pool_bg_log2fc = rbind(transit_mat_pool_high_log2fc, transit_mat_pool_low_log2fc, transit_mat_bg1_log2fc, transit_mat_bg2_log2fc)

### decide K
set.seed(2019)
colmeans = colMeans(transit_mat_pool_bg_log2fc)
ss0 = sum(apply(transit_mat_pool_bg_log2fc, 1, function(x) sum((x-colmeans)^2)))
withinss_all = c(ss0)
betweenss_all = c(ss0)
for (k in 2:K_try){
	print(k)
	all_counts_clust = kmeans(transit_mat_pool_bg_log2fc, centers=k)
	withinss_all = c(withinss_all, all_counts_clust$tot.withinss)
	betweenss_all = c(betweenss_all, all_counts_clust$betweenss)
}
betweenss_all[1] = 0
pdf(kmeans_decide_k)
plot(1:K_try, withinss_all, type='p', pch=16, col='red', ylim = c(0, max(withinss_all)))
lines(1:K_try, withinss_all, lwd=1.5, col='red')
points(1:K_try, betweenss_all, pch=16, col='blue')
lines(1:K_try, betweenss_all, lwd=1.5, col='blue')
log2_ratio = log2((withinss_all+1)/(betweenss_all+1))
plot(2:K_try, log2_ratio[-1], type='p', pch=16, col='blue')
lines(2:K_try, log2_ratio[-1], lwd=1.5, col='blue')
plot(3:K_try, log2_ratio[-c(1,length(log2_ratio))]-log2_ratio[-c(1,2)], type='p', pch=16, col='blue')
lines(3:K_try, log2_ratio[-c(1,length(log2_ratio))]-log2_ratio[-c(1,2)], lwd=1.5, col='blue')
dev.off()



