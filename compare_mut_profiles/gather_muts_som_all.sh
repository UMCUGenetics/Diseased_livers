#!/bin/bash

base_dir=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/

manifest_path=$base_dir/metadata/manifest.txt.gz
out_file=$base_dir/analysis3/compare_mut_profiles/muts_som_all.txt.gz

echo -e 'sample\tchrom\tpos\tref\talt' | gzip -c > $out_file

header=$(zcat $manifest_path | head -n1)
counter=0
zcat $manifest_path | tail -n +2 | while read $header; do
	counter=$((counter+1))

	echo -e "[$counter] $sample_name"

	zcat $som_vcf | grep -v '^#' | grep 'PASS' |
	awk -v sample="$sample_name" '{print sample"\t"$1"\t"$2"\t"$4"\t"$5}' |
	gzip -c >> $out_file

	#if [[ $counter -eq 1 ]]; then break; fi
done
