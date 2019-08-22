### get parameters
args = commandArgs(trailingOnly=TRUE)

gene_state_tranmat1 = args[1]#'G1Ehigh.G1E.gene.withbinID.tranmat1d.500.txt'
gene_state_tranmat2 = args[2]#'G1Ehigh.ER4.gene.withbinID.tranmat1d.500.txt'
kmeans_decide_k = args[3]#'kmeans_ss.pairwise.1d.pdf'
K_try = as.numeric(args[4])
add_smallnum = as.numeric(args[5])

### get transit_mat seq
transit_mat_pool1 = read.table(gene_state_tranmat1, sep='\t', header=F)
transit_mat_pool2 = read.table(gene_state_tranmat2, sep='\t', header=F)

### log2fc
transit_mat_pool_log2fc = log2((transit_mat_pool1+add_smallnum) / (transit_mat_pool2+add_smallnum))

### decide K
set.seed(2019)
colmeans_pca = colMeans(transit_mat_pool_log2fc)
ss0 = sum(apply(transit_mat_pool_log2fc, 1, function(x) sum((x-colmeans_pca)^2)))
withinss_all = c(ss0)
betweenss_all = c(ss0)
for (k in 2:K_try){
	print(k)
	all_counts_clust = kmeans(transit_mat_pool_log2fc, centers=k)
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


