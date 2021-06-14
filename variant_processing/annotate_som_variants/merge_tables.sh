#!/bin/bash


dir=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/scripts/annotate_som_variants/output/
out_txt=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/scripts/annotate_som_variants/mut_profile_som.txt.gz

counter=0
for in_txt in $(find $dir -type f -name '*txt.gz'); do
	counter=$((counter+1))

	sample_name=$(basename $(dirname $in_txt))
	echo -ne "Processing [$counter]: $sample_name\r"

	## Write header using first file
	if [[ $counter -eq 1 ]]; then
		zcat $in_txt | head -n 1 | 
		awk '{print "sample""\t"$0}' | 
		gzip -c > $out_txt
	fi

	zcat $in_txt | tail -n +2 | 
	#grep BRCA | 
	awk -v sample_name="$sample_name" '{print sample_name"\t"$0}' | 
	gzip -c >> $out_txt

done

echo -e '\n'


