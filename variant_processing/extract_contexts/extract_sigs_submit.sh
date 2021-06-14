#!/bin/bash

manifest=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/metadata/manifest.txt.gz

wd=/hpc/cuppen/projects/P0003_FOOTPRINTS/WGS_clones/processed/Liver_footprint/scripts/extract_sigs/
rscript=$wd/extract_sigs.R
job_dir=$wd/jobs/; mkdir -p $job_dir

header=$(zcat $manifest | head -n1)
counter=0
zcat $manifest | tail -n+2 | while read $header; do
	counter=$((counter+1))
	job_file=$job_dir/${sample_name}.job
	
	echo "#!/bin/bash
#SBATCH --job-name=extract_sigs_${sample_name}
#SBATCH --output=${job_file}.o
#SBATCH --time=01:00:00
#SBATCH --mem=10G

guixr load-profile ~/.guix-profile/ --<<EOF
Rscript $rscript $som_vcf $sv_vcf $sample_name && touch ${job_file}.done
EOF
" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch $job_file
	fi

	#if [[ $counter -eq 1 ]]; then break; fi
done
