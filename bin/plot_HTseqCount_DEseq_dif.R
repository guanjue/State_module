library(LSD)
#install.packages('ggplot2', dependencies=TRUE, INSTALL_opts = c('--no-lock'))
library(DESeq2)


### get parameters
args = commandArgs(trailingOnly=TRUE)

count_table = args[1]#'rnaHtseqCountsall_withcoordinates.0.txt'
ct1 = args[2]#'G1E'
ct2 = args[3]#'ER4'
log2fc_thresh = as.numeric(args[4])#1
p_thresh = as.numeric(args[5])#0.05
ct1_rep1 = as.numeric(args[6])#28
ct1_rep2 = as.numeric(args[7])#29
ct2_rep1 = as.numeric(args[8])#26
ct2_rep2 = as.numeric(args[9])#27

d = read.table(count_table, header=T, sep='\t')

d_G1E_ER4 = d[,c(ct2_rep1, ct2_rep2, ct1_rep1, ct1_rep2)]

coldata = cbind(c(ct2, ct2, ct1, ct1))
colnames(coldata) = c('condition')

dds <- DESeqDataSetFromMatrix(countData = d_G1E_ER4,
                              colData = coldata,
                              design= ~ condition)



dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_G1E_vs_ER4")
# or to shrink log fold changes association with condition:
#
#resLFC <- lfcShrink(dds, coef="condition_G1E_vs_ER4", type="apeglm")
#resLFC

res <- results(dds, contrast=c("condition",ct1,ct2))

png(paste(ct1 ,'_vs_', ct2, '.MAplot.png', sep=''))
plotMA(res, ylim=c(-2,2))
dev.off()


G1E_high_exp = d[(!is.na(res$padj)) & (res$log2FoldChange>=log2fc_thresh) & (res$padj<p_thresh),1:6]
ER4_high_exp = d[(!is.na(res$padj)) & (res$log2FoldChange<=(-log2fc_thresh)) & (res$padj<p_thresh),1:6]
nodif = d[(!is.na(res$padj)) & (abs(res$log2FoldChange)<=(log2fc_thresh)) & (res$padj>=p_thresh),1:6]

write.table(G1E_high_exp, paste(ct1, '_vs_', ct2, '.', ct1, 'high.gene.txt', sep=''), quote=F,sep='\t', col.names=F, row.names=F)
write.table(ER4_high_exp, paste(ct1, '_vs_', ct2, '.', ct2, 'high.gene.txt', sep=''), quote=F,sep='\t', col.names=F, row.names=F)
write.table(nodif, paste(ct1, '_vs_', ct2, '.nodif.gene.txt', sep=''), quote=F,sep='\t', col.names=F, row.names=F)


