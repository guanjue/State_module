library(pheatmap)

### get parameters
args = commandArgs(trailingOnly=TRUE)

input_folder = args[1]
input_file_name = args[2]
count_mat = args[3]
para_file = args[4]#../ideasM2/ideas_result2/ideas_m2.para.modified.para

ouput_name = args[5]#'G1E_vs_ER4.motifs.1d.all.K.pdf'
small_num = as.numeric(args[6])#100
k = as.numeric(args[7])#100
kcount = as.numeric(args[8])
#order_used = c(16,11,18,15,2,0,3,4,6,5,20,7,13,9,19,10,24,17,25,8,1,21,23,14,22,12)+1

### read para file & para mean signal
para_signal = read.table(para_file, header=T)
k0 = dim(para_signal)[2]
p = (sqrt(9+8*(k0-1))-3)/2
para_all_signal = para_signal[,2:p]
para_mean_signal = para_all_signal / para_signal[,1]
para_mean_sq_signal = rowSums(para_mean_signal^2)

### get number of each state in each K
d12_dif_all = c()
for (pos_i in 1:k){
print(pos_i)
d1 = read.table(paste(input_folder, '/', input_file_name, '.kmeans.K', pos_i, '.1.txt', sep=''), header=F)
d2 = read.table(paste(input_folder, '/', input_file_name, '.kmeans.K', pos_i, '.2.txt', sep=''), header=F)
d12_dif = as.matrix(log2((d1+small_num)/(d2/sum(d2)*sum(d1)+small_num)))
d12_dif_all = cbind(d12_dif_all, d12_dif)
}

### get state new orders
dist <- dist(d12_dif_all, method = "euclidean") # distance matrix
fit <- hclust(dist, method="ward.D2") 
hclust_order = fit$order

### get count mat
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
pheatmap(t(d12_dif_all[hclust_order,order(rowMeans(bglog2fc), decreasing=T)]),cluster_rows=F,cluster_cols=F, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=F,show_colnames=TRUE,annotation_names_row=F,annotation_names_col=TRUE)
dev.off()

state_order_for_tomtom = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)+1
write.table(t(d12_dif_all[state_order_for_tomtom,order(rowMeans(bglog2fc), decreasing=T)]), paste(ouput_name, '.heatmap.signal.log2fc.txt', sep=''), quote=F, sep='\t', col.names=F, row.names=F)
write.table(t(2^d12_dif_all[state_order_for_tomtom,order(rowMeans(bglog2fc), decreasing=T)]), paste(ouput_name, '.heatmap.signal.txt', sep=''), quote=F, sep='\t', col.names=F, row.names=F)
d12_dif_all_write = t(d12_dif_all[state_order_for_tomtom,order(rowMeans(bglog2fc), decreasing=T)])
d12_dif_all_write = t(apply(d12_dif_all_write, 1, function(x) x-min(x)))
write.table(d12_dif_all_write, paste(ouput_name, '.heatmap.signal.log2fc.pos.txt', sep=''), quote=F, sep='\t', col.names=F, row.names=F)


print(summary(bglog2fc))
#print(bglog2fc)
corlo_lim_scale = max(abs(bglog2fc))
breaksList_scale = seq(-corlo_lim_scale, corlo_lim_scale, by = 0.01)
my_colorbar_scale=colorRampPalette(c('blue', 'white', 'red'))(n = length(breaksList_scale))

pdf(paste(ouput_name, '.count.pdf', sep=''), width=2, height=7)
pheatmap((bglog2fc[order(rowMeans(bglog2fc), decreasing=T),]),cluster_rows=F,cluster_cols=F, color=my_colorbar_scale, breaks = breaksList_scale, show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=TRUE,annotation_names_col=TRUE)
dev.off()



### get state reordered heatmap
createHeatmap_sort_eRP = function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), show=TRUE,fout=NULL, sortstate=TRUE,scale=FALSE)
{	x=read.table(parafile,comment="!",header=T);
	k=dim(x)[2];
	l=dim(x)[1];
	p=(sqrt(9+8*(k-1))-3)/2;
	m=as.matrix(x[,1+1:p]/x[,1]);
	colnames(m) = colnames(x)[1+1:p];
	marks=colnames(m);
	rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");

if(sortstate)
{
#o=hclust(dist(m),method="ward.D2")$order;
o=hclust_order
m=m[hclust_order,];
if(length(statecolor) != 0)
{	statecolor=statecolor[o,];
}
}
om=m;
if(scale)
{	m = t((t(m) - apply(m,2,min))/(apply(m,2,max)-apply(m,2,min)+1e-10));
}

	if(length(fout)!=0)
	{	pdf(fout);	}	
	par(mar=c(6,1,1,6));
	rg=range(m);
	colors=0:100/100*(rg[2]-rg[1])+rg[1];
        my_palette=colorRampPalette(cols)(n=100);
	defpalette=palette(my_palette);

if(show)
{
	plot(NA,NA,xlim=c(0,p+0.7),ylim=c(0,l),xaxt="n",yaxt="n",xlab=NA,ylab=NA,frame.plot=F);
	axis(1,at=1:p-0.5,labels=colnames(m),las=2);
	axis(4,at=1:l-0.5,labels=rownames(m),las=2);
	rect(rep(1:p-1,l),rep(1:l-1,each=p),rep(1:p,l),rep(1:l,each=p),col=round((t(m)-rg[1])/(rg[2]-rg[1])*100));#,border=NA);
}
if(scale)
{	m = om;
}

	if(length(statecolor)==0)
	{	if(length(markcolor)==0)
		{	markcolor=t(col2rgb(terrain.colors(ceiling(p))[1:p]));	
			for(i in 1:length(marks))
			{	if(regexpr("h3k4me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(255,0,0);	}
				if(regexpr("h3k4me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,100,0);	}
				if(regexpr("h3k4me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,250,0);	}
				if(regexpr("h3k36me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,0);	}
				if(regexpr("h2a",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,150);	}
				if(regexpr("dnase",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,200,200);	}
				if(regexpr("atac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,50,150);	}
				if(regexpr("h3k9ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,0,200);	}
				if(regexpr("h3k9me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(100,100,100);	}
				if(regexpr("h3k27ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,150,0);	}
				if(regexpr("h3k27me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,0,225);	}
				if(regexpr("h3k79me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,200);	}
				if(regexpr("h4k20me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(50,200,50);	}
				if(regexpr("ctcf",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,250);	}
			}
		}
		statecolor=array(stateColor(m,markcolor),dim=c(dim(m)[1],2));
	}
	if(show)
	{	rect(rep(p+0.2,l),1:l-0.8,rep(p+0.8,l),1:l-0.2,col=statecolor[,2]);
	}
	if(sortstate)	statecolor[o,]=statecolor;

	palette(defpalette);
	if(length(fout)!=0)
	{	dev.off();	}
	return(statecolor);
}

stateColor<-function(statemean, markcolor=NULL)
{	
	if(length(markcolor)==0)
	{	markcolor=rep("",dim(statemean)[2]);
		markcolor[order(apply(statemean,2,sd),decreasing=T)]=hsv((1:dim(statemean)[2]-1)/dim(statemean)[2],1,1)
		markcolor=t(col2rgb(markcolor));
	}

	rg=apply(statemean,1,range);
	mm=NULL;
	for(i in 1:dim(statemean)[1])
	{	mm=rbind(mm,(statemean[i,]-rg[1,i]+1e-10)/(rg[2,i]-rg[1,i]+1e-10));
	}
	mm = mm^5; 
	if(dim(mm)[2]>1) mm = mm / (apply(mm, 1, sum)+1e-10);
	mycol=mm%*%markcolor;
	s=apply(statemean,1,max);
	s=(s-min(s))/(max(s)-min(s)+1e-10);
#s=s^0.5;
	
mycol=round(255-(255-mycol)*s/0.5);
mycol[mycol<0]=0;
rt=paste(mycol[,1],mycol[,2],mycol[,3],sep=",");
h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
rt=cbind(rt,h);
return(rt);

	h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
	h[,2]=h[,2]*s;
	#h[,3]=1-(1-h[,3])*s^0.5;
	h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
	rt=cbind(apply(t(col2rgb(h)),1,function(x){paste(x,collapse=",")}),h);
	
	return(rt);
}
pdf(paste(ouput_name, '.state_heatmap.pdf', sep=''), width=2.5, height=5.5)
createHeatmap_sort_eRP(para_file, sortstate=TRUE)
dev.off()



