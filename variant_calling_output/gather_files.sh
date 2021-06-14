#!/bin/bash

## Paths ================================
base_dir=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/
wd=$base_dir/analysis3/variant_calling_output/

manifest=$wd/manifest.ss.txt.gz
out_dir=$wd/data/; mkdir -p $out_dir

## Main ================================
header=$(zcat $manifest | head -n1)
counter=0
zcat $manifest | tail -n+2 | while read $header; do

	counter=$((counter+1))

	echo "[$counter]: $sample_name"
	sample_dir=$out_dir/$sample_name/; mkdir -p $sample_dir
	done_file=$sample_dir/cp.done

	if [[ ! -f $done_file ]]; then
		cp $som_vcf $sv_vcf $cnv_tsv $sample_dir
	fi

	#if [[ $counter -eq 1 ]]; then break; fi
done

