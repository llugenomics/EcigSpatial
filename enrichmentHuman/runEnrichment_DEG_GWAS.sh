#! /bin/bash


   genesDEG=/genomics/chiragProjects/bgiStereoseq/analysisGwas/3months_snRNAseq/genesDEG.bed
genesNonDEG=/genomics/chiragProjects/bgiStereoseq/analysisGwas/3months_snRNAseq/genesNonDEG.bed


##################################################################################################################################################################################################
##################################################################################################################################################################################################
output1=outputOverlapDEG_GWAS
rm -rf $output1

# Uplifted RAT GWAS
gwasDIR=/genomics/chiragProjects/bgiStereoseq/analysisGwas/dataset/pshycEncode/gwas

for snv in $(ls $gwasDIR/upLifted_rn7AgeFirstBirth_Mills_meta $gwasDIR/upLifted_rn7CigarettesPerDay_meta_Koskeridis $gwasDIR/upLifted_rn7MDD_PGC_meta $gwasDIR/upLifted_rn7Bipolar_PGC_meta $gwasDIR/upLifted_rn7EverSmoked_meta_Karlsson $gwasDIR/upLifted_rn7SCZ_meta_PGC)
do
	out=$(echo $snv | sed 's/\//\t/g' | awk '{ print $NF }' | sed 's/upLifted_rn7//g' | sed 's/_meta//g')

	# Focus analyses only SNVs that overlap genes of interest

	count1=$(cat $snv |  intersectBed -wa -wb -b stdin -a $genesDEG    | cut -f5 |sort | uniq |wc -l)
	count2=$(cat $snv |  intersectBed      -v -b stdin -a $genesDEG    | cut -f5 |sort | uniq |wc -l)
	count3=$(cat $snv |  intersectBed -wa -wb -b stdin -a $genesNonDEG | cut -f5 |sort | uniq |wc -l)
	count4=$(cat $snv |  intersectBed      -v -b stdin -a $genesNonDEG | cut -f5 |sort | uniq |wc -l)
	echo $out $count1 $count2 $count3 $count4 | awk 'BEGIN {OFS="\t"} {out=$1; for ( i=2; i<= NF; i++) out=out"\t"$i; print out }' 
        echo $out $count1 $count2 $count3 $count4 | awk 'BEGIN {OFS="\t"} {out=$1; for ( i=2; i<= NF; i++) out=out"\t"$i; print out }' | sed -e 's/_filtered//g' >> $output1
	
done



