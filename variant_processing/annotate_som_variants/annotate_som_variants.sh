#!/bin/bash

base_dir=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/

manifest_path=$base_dir/metadata/manifest.txt.gz
parent_out_dir=$base_dir/scripts/annotate_som_variants/output/; mkdir -p $parent_out_dir

RSCRIPT=$base_dir/scripts/annotate_som_variants/annotate_som_variants.R

header=$(zcat $manifest_path | head -n1)
counter=0
zcat $manifest_path | tail -n +2 | while read $header; do
	counter=$((counter+1))

	echo -e "\n######### [$counter] $sample_name #########"

	out_dir=$parent_out_dir/$sample_name/; mkdir -p $out_dir

	job_dir=$out_dir/jobs/; mkdir -p $job_dir
	job_file=$job_dir/${sample_name}_annSomVar.job
	done_file=${job_file}.done

echo "#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --job-name=$(basename $job_file)
#SBATCH --output=${job_file}.o
guixr load-profile ~/.guix-profile --<<EOF
Rscript $RSCRIPT $som_vcf $out_dir $sample_name && touch $done_file
EOF
" > $job_file
	
	if [[ ! -f $done_file ]]; then
		sbatch $job_file
	else
		echo "Skipping: $(basename $job_file)"
	fi

	#if [[ $counter -eq 1 ]]; then break; fi
done
