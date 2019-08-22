
cd /storage/home/g/gzx103/group/projects/vision_human/state_module

cp /storage/home/g/gzx103/group/projects/vision_human/ideasM2/windowsNoBlackForIdeas.bed ./

cp /storage/home/g/gzx103/group/projects/vision/rna/rnaTPMall_withcoordinates.0.txt ./


### get bed files
time cat windowsNoBlackForIdeas.bed | awk -F ' ' -v OFS='\t' '{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n > windowsNoBlack.sort.bed
time tail -n+2 rnaTPMall_withcoordinates.0.txt | awk -F '\t' -v OFS='\t' '{if ($5=="protein_coding" && $6=="+") print $1,$2,$2+1,$4"_"$6; else if ($5=="protein_coding" && $6=="-") print $1,$3-1,$3,$4"_"$6}' | sort -k1,1 -k2,2n > gene_id.pc.sort.bed

### get TSS bin ID
time bedtools map -a gene_id.pc.sort.bed -b windowsNoBlack.sort.bed -c 4 -o concat > gene_id_binID.pc.sort.bed



ct1=G1E
ct2=ER4
### get G1E vs ER4 state
tail -n+2 ../ideasM2/ideas_result2/ideas_m2.state | awk -F ' ' -v OFS='\t' '{print $1,$5}' > $ct1'.state.txt'
tail -n+2 ../ideasM2/ideas_result2/ideas_m2.state | awk -F ' ' -v OFS='\t' '{print $1,$22}' > $ct2'.state.txt'


ct1=G1E
ct2=ER4
ct1rep1=28
ct1rep2=29
ct2rep1=26
ct2rep2=27
k=30


### get G1E exp & noexp genes
time Rscript plot_HTseqCount_DEseq_dif.R rnaHtseqCountsall_withcoordinates.0.txt $ct1 $ct2 1 0.05 $ct1rep1 $ct1rep2 $ct2rep1 $ct2rep2
###
### get tss with binID
cat $ct1'_vs_'$ct2'.'$ct1'high.gene.txt' | awk -F ' ' -v OFS='\t' '{if ($5=="protein_coding") print $1,$2,$3,$4"_"$6}' > $ct1'_vs_'$ct2'.'$ct1'high.gene.named.txt'
cat $ct1'_vs_'$ct2'.'$ct2'high.gene.txt' | awk -F ' ' -v OFS='\t' '{if ($5=="protein_coding") print $1,$2,$3,$4"_"$6}' > $ct1'_vs_'$ct2'.'$ct2'high.gene.named.txt'
cat $ct1'_vs_'$ct2'.nodif.gene.txt' | awk -F ' ' -v OFS='\t' '{if ($5=="protein_coding") print $1,$2,$3,$4"_"$6}' > $ct1'_vs_'$ct2'.nodif.gene.named.txt'
### get bin ID
python vlookup.py -t gene_id_binID.pc.sort.bed -m 4 -s $ct1'_vs_'$ct2'.'$ct1'high.gene.named.txt' -n 4 -o $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt' 
cat $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt' | awk -F '_' -v OFS='\t' '{print $1,$2}' > $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt.tmp' && mv $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt.tmp' $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt' 
python vlookup.py -t gene_id_binID.pc.sort.bed -m 4 -s $ct1'_vs_'$ct2'.'$ct2'high.gene.named.txt' -n 4 -o $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt' 
cat $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt' | awk -F '_' -v OFS='\t' '{print $1,$2}' > $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt.tmp' && mv $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt.tmp' $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt' 
python vlookup.py -t gene_id_binID.pc.sort.bed -m 4 -s $ct1'_vs_'$ct2'.nodif.gene.named.txt' -n 4 -o $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt' 
cat $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt' | awk -F '_' -v OFS='\t' '{print $1,$2}' > $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt.tmp' && mv $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt.tmp' $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt' 
###
### get state seq
time Rscript get_state_seq.R $ct1'.state.txt' $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt' 500 $ct1'_vs_'$ct2'.'$ct1'high.'$ct1'.gene.withbinID.statemat.500.txt'
time Rscript get_state_seq.R $ct2'.state.txt' $ct1'_vs_'$ct2'.'$ct1'high.gene.withbinID.txt' 500 $ct1'_vs_'$ct2'.'$ct1'high.'$ct2'.gene.withbinID.statemat.500.txt'
time Rscript get_state_seq.R $ct1'.state.txt' $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt' 500 $ct1'_vs_'$ct2'.'$ct2'high.'$ct1'.gene.withbinID.statemat.500.txt'
time Rscript get_state_seq.R $ct2'.state.txt' $ct1'_vs_'$ct2'.'$ct2'high.gene.withbinID.txt' 500 $ct1'_vs_'$ct2'.'$ct2'high.'$ct2'.gene.withbinID.statemat.500.txt'
time Rscript get_state_seq.R $ct1'.state.txt' $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt' 500 $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.statemat.500.txt'
time Rscript get_state_seq.R $ct2'.state.txt' $ct1'_vs_'$ct2'.nodif.gene.withbinID.txt' 500 $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.statemat.500.txt'
###
###
###
### get G1E & ER4 transit matrix
cat $ct1'_vs_'$ct2'.'$ct1'high.'$ct1'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct2'.gene.withbinID.statemat.500.txt' > 'high.'$ct1'_'$ct2'.gene.withbinID.statemat.500.txt'
cat $ct1'_vs_'$ct2'.'$ct1'high.'$ct2'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct1'.gene.withbinID.statemat.500.txt' > 'low.'$ct1'_'$ct2'.gene.withbinID.statemat.500.txt'
###
time Rscript get_transition_matrix_pairwise.1d.R 'high.'$ct1'_'$ct2'.gene.withbinID.statemat.500.txt' 'low.'$ct1'_'$ct2'.gene.withbinID.statemat.500.txt' 'high.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' 'low.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' 50 5
time Rscript get_transition_matrix_pairwise.1d.R $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' 50 5
### get Kmeans state motifs
time Rscript get_Kmeans_decide_K.1d.withbg.R 'high.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' 'low.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' 'high.'$ct2'_vs_'$ct1'.kmeans_ss.pairwise.1d.pdf' 50 1
time Rscript get_Kmeans_Motif.1d.withbg.R 'high.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' 'low.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' 'high_'$ct2'_vs_'$ct1'_kmeans_pairwise_motif_1d' 'high_'$ct2'_vs_'$ct1'_kmeans_motif_1d' 'high_'$ct2'_vs_'$ct1'.heatmap.1d.png' 1 $k
time Rscript plot_transition_mat.kmeans.1heatmap.1d.withbg.signal.R 'high_'$ct2'_vs_'$ct1'_kmeans_pairwise_motif_1d' 'high_'$ct2'_vs_'$ct1'_kmeans_motif_1d' 'high_'$ct2'_vs_'$ct1'_kmeans_motif_1d.counts.txt' '../ideasM2/ideas_result2/ideas_m2.para.modified.para' 'high_'$ct2'_vs_'$ct1'.motifs.1d.all.K.pdf' 100 $k 3











k=30
time Rscript get_Kmeans_Motif.1d.withbg.signal.R 'high.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' 'low.'$ct1'_'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' '../ideasM2/ideas_result2/ideas_m2.para.modified.para' 'high_'$ct2'_vs_'$ct1'_kmeans_pairwise_motif_1d_signal' 'high_'$ct2'_vs_'$ct1'_kmeans_motif_1d_signal' 'high_'$ct2'_vs_'$ct1'.heatmap.1d.signal.png' 1 $k
time Rscript plot_transition_mat.kmeans.1heatmap.1d.withbg.signal.R 'high_'$ct2'_vs_'$ct1'_kmeans_pairwise_motif_1d_signal' 'high_'$ct2'_vs_'$ct1'_kmeans_motif_1d_signal' 'high_'$ct2'_vs_'$ct1'_kmeans_motif_1d_signal.counts.txt' '../ideasM2/ideas_result2/ideas_m2.para.modified.para' 'high_'$ct2'_vs_'$ct1'.motifs.1d.all.K.signal.pdf' 100 $k 3











### get G1E transit matrix
time Rscript get_transition_matrix_pairwise.1d.R $ct1'_vs_'$ct2'.'$ct1'high.'$ct1'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.'$ct1'high.'$ct2'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.'$ct1'high.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct1'high.'$ct2'.gene.withbinID.tranmat1d.500.txt' 50 5
time Rscript get_Kmeans_decide_K.1d.R $ct1'_vs_'$ct2'.'$ct1'high.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct1'high.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.kmeans_ss.pairwise.1d.pdf' 50 1
time Rscript get_Kmeans_Motif.1d.R $ct1'high.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct1'high.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct1'_vs_'$ct2'_kmeans_pairwise_motif_1d' $ct1'_vs_'$ct2'_kmeans_motif_1d' $ct1'_vs_'$ct2'.heatmap.1d.png' 1 20
time Rscript plot_transition_mat.kmeans.1heatmap.1d.R $ct1'_vs_'$ct2'_kmeans_pairwise_motif_1d' $ct1'_vs_'$ct2'_kmeans_motif_1d' $ct1'_vs_'$ct2'.motifs.1d.all.K.pdf' 100 20

### get ER4 transit matrix
time Rscript get_transition_matrix_pairwise.1d.R $ct1'_vs_'$ct2'.'$ct2'high.'$ct2'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct1'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct1'.gene.withbinID.tranmat1d.500.txt' 50 5
time Rscript get_Kmeans_decide_K.1d.R $ct1'_vs_'$ct2'.'$ct2'high.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct2'_vs_'$ct1'.kmeans_ss.pairwise.1d.pdf' 50 1
time Rscript get_Kmeans_Motif.1d.R $ct1'_vs_'$ct2'.'$ct2'high.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.'$ct2'high.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct2'_vs_'$ct1'_kmeans_pairwise_motif_1d' $ct2'_vs_'$ct1'_kmeans_motif_1d' $ct2'_vs_'$ct1'.heatmap.1d.png' 1 20
time Rscript plot_transition_mat.kmeans.1heatmap.1d.R $ct2'_vs_'$ct1'_kmeans_pairwise_motif_1d' $ct2'_vs_'$ct1'_kmeans_motif_1d' $ct2'_vs_'$ct1'.motifs.1d.all.K.pdf' 100 20







### get fasta
time python vec2fasta.py -i G1Ehigh.G1E.gene.withbinID.statemat.500.txt -o G1Ehigh.G1E.gene.withbinID.statemat.500.fa
time python vec2fasta.py -i G1Ehigh.ER4.gene.withbinID.statemat.500.txt -o G1Ehigh.ER4.gene.withbinID.statemat.500.fa
time python vec2fasta.py -i ER4high.G1E.gene.withbinID.statemat.500.txt -o ER4high.G1E.gene.withbinID.statemat.500.fa
time python vec2fasta.py -i ER4high.ER4.gene.withbinID.statemat.500.txt -o ER4high.ER4.gene.withbinID.statemat.500.fa







### get nodif transit matrix
time Rscript get_transition_matrix_pairwise.1d.R $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.statemat.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' 50 5
time Rscript get_Kmeans_decide_K.1d.R $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.kmeans_ss.pairwise.1d.pdf' 50 1
time Rscript get_Kmeans_Motif.1d.R $ct1'_vs_'$ct2'.nodif.'$ct2'.gene.withbinID.tranmat1d.500.txt' $ct1'_vs_'$ct2'.nodif.'$ct1'.gene.withbinID.tranmat1d.500.txt' $ct2'_vs_'$ct1'_nodif_kmeans_pairwise_motif_1d' $ct2'_vs_'$ct1'_nodif_kmeans_motif_1d' $ct2'_vs_'$ct1'.nodif.heatmap.1d.png' 1 20
time Rscript plot_transition_mat.kmeans.1heatmap.1d.R $ct2'_vs_'$ct1'_nodif_kmeans_pairwise_motif_1d' $ct2'_vs_'$ct1'_nodif_kmeans_motif_1d' $ct2'_vs_'$ct1'.nodif.motifs.1d.all.K.pdf' 100 20



### get mat
time Rscript get_transition_matrix.R G1E.state.txt G1E_exp.gene.withbinID.txt G1E_exp.gene.withbinID.statemat.txt G1E_exp.gene.withbinID.transmat.txt G1E_wg.transmat.txt 50 u
time Rscript get_transition_matrix.R G1E.state.txt G1E_noexp.gene.withbinID.txt G1E_noexp.gene.withbinID.statemat.txt G1E_noexp.gene.withbinID.transmat.txt G1E_wg.transmat.txt 50 u

time Rscript plot_transition_mat.R

time Rscript get_transition_matrix.R G1E.state.txt G1E_exp.gene.withbinID.txt G1E_exp.gene.withbinID.statemat.50u.txt G1E_exp.gene.withbinID.transmat.50u.txt G1E_wg.transmat.txt 50 u
time Rscript get_transition_matrix.R G1E.state.txt G1E_noexp.gene.withbinID.txt G1E_noexp.gene.withbinID.statemat.50u.txt G1E_noexp.gene.withbinID.transmat.50u.txt G1E_wg.transmat.txt 50 u

for i in {1..10}
do
echo $i
pos_start=$(($i * 50))
pos_end=$(($i - 1 * 50))
time Rscript get_transition_matrix.R G1E.state.txt G1E_exp.gene.withbinID.txt b12/G1E_exp.gene.withbinID.statemat.50u.txt b12/G1E_exp.gene.withbinID.transmat.50u.txt b12/G1E_wg.transmat.txt $pos_start $pos_end u
done


### get fasta
python vec2fasta.py -i G1E_exp.gene.withbinID.statemat.txt -o G1E_exp.gene.withbinID.statemat.fa
python vec2fasta.py -i G1E_noexp.gene.withbinID.statemat.txt -o G1E_noexp.gene.withbinID.statemat.fa
python vec2fasta.py -i G1E_exp.gene.withbinID.statemat.500.txt -o G1E_exp.gene.withbinID.statemat.500.fa
python vec2fasta.py -i G1E_noexp.gene.withbinID.statemat.500.txt -o G1E_noexp.gene.withbinID.statemat.500.fa



