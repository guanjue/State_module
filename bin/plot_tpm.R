library(LSD)

d = read.table('rnaTPMall_withcoordinates.0.txt', header=T)

small_num = 1e-1
ER4 = log2(d[,26]/2+d[,27]/2+small_num)
G1E = log2(d[,28]/2+d[,29]/2+small_num)

lim = max(cbind(G1E, ER4))

pdf('G1E_vs_ER4.TPM.pdf')
heatscatter(ER4, G1E, xlim=c(log2(small_num), lim), ylim=c(log2(small_num), lim))
abline(0, 1, lwd=1.5)
dev.off()

pdf('G1E.TPMhist.pdf')
hist(G1E, breaks=50)
abline(v=7.5, lwd=1.5)
box()
dev.off()

G1E_exp = d[G1E>=7.5,1:6]
G1E_noexp = d[G1E<7.5,1:6]

write.table(G1E_exp, 'G1E_exp.gene.txt', quote=F,sep='\t', col.names=F, row.names=F)
write.table(G1E_noexp, 'G1E_noexp.gene.txt', quote=F,sep='\t', col.names=F, row.names=F)

