
DEG=/genomics/chiragProjects/bgiStereoseq/analysisGwas/3months_snRNAseq/orthologousDEG
NoDEG=/genomics/chiragProjects/bgiStereoseq/analysisGwas/3months_snRNAseq/orthologousNonDEG

count=count
rm -rf $count

for name in $(ls deg*csv); 
do 
	out=$(echo $name | sed 's/.csv//g'); 
	# Filter DEGs with p-adjusted 0.05  
	cat $name | sed 's/,/\t/g' | awk '{ print $0 }' | awk '{ if ($(NF-1) < 0.05) { print $0 }}' > out$out;

	join  -t $'\t' <(cat $DEG | cut -f2 |sort -k1,1) <(cat out$out |sort -k1,1) > overalp$out
	
	count1=$(join      <(cat $DEG | cut -f2 |sort -k1,1) <(cat out$out |sort -k1,1) |wc -l)
	count2=$(join -v 1 <(cat $DEG | cut -f2 |sort -k1,1) <(cat out$out |sort -k1,1) |wc -l)
	count3=$(join      <(cat $NoDEG | cut -f2 |sort -k1,1) <(cat out$out |sort -k1,1) |wc -l)
	count4=$(join -v 1 <(cat $NoDEG | cut -f2 |sort -k1,1) <(cat out$out |sort -k1,1) |wc -l)
	rm -rf out$out;

	echo $name $count1 $count2 $count3 $count4 | awk '{ out=$1; for (i=2; i<=NF; i++) out=out"\t"$i; print out }' >> $count

done




