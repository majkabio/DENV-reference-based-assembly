#!/bin/bash

##################################################################
#
#	v2.2 update:	added a column in alignment stats that contains average depth
#					reference file is required
#
#	Modifications from v2:
#	-- removed prinseq step; to use prinseq fastq from canu_metagenomics pipeline
#
# 	Name of program: ref-based-assmebly_v2.sh
# 	Created by: Majo Galarion
# 	Date created: 21 February 2024
# 	Revisions:
#	do medaka polishing of consensus sequence
#	updated way to create a neat stat summary table
#	
#	
# 	This script is a pipeline to do a reference-based assembly of
#	datasets from a non-tiling PCR pipeline
#
##################################################################

# Always run this script inside directory of your input files!
# Always use a custom reference in same directory as input fastq

### Important dependencies:

# Define tools and how you call them in your current system
porechop="porechop-runner.py"
prinseq="prinseq-lite.pl"
minimap2="minimap2"
samtools="samtools"
bcftools="bcftools"
medaka=". ~/medaka/venv/bin/activate"

# Takes a non-compressed, concatenated FASTQ file


usage() {
	echo "";
	echo "Reference-based assembly for metagenomics dataset"
	echo "";
	echo "Usage $0 -i reads.fastq -r reference.fasta" 1>&2;
	echo "";
	echo "Optional parameters:";
	echo "-t: number of threads (default: 8)";
	#echo "-r: Reference genome (default: DENV-1 RefSeq)";
	echo "-M: Medaka polishing model (default: r941_min_hac_g507m)";
	echo "-m: minimum read length (default: 200)";
	echo "-x: maximum read length (default: none)";
	echo "-d: minimum read depth (default: 20)";
	echo "-q: minimum read Q-score (default: 9)";
	echo ""
	exit 1;
}

if [ $# -eq 0 ]; then usage; exit 1; fi


### These are the default parameters
threads="8";
#ref="/home3/2509094g/DENV-ref/DENV-1_refseq.fasta";
medakaModel="r941_min_hac_g507";
minLength=200;
maxLength="none";
minDepth=20;
minQual=9;

### Read these arguments when specified in command
while getopts ":i:t:r:b:M:m:x:d:q:" opt
do
	case "${opt}" in
		i) fq=${OPTARG} ;;
		t) threads=${OPTARG} ;;
		r) ref=${OPTARG} ;;
		M) medakaModel=${OPTARG} ;;
		m) minLength=${OPTARG} ;;
		x) maxLength=${OPTARG} ;;
		d) minDepth=${OPTARG} ;;
		q) minQual=${OPTARG} ;;
		?) echo "Option -${OPTARG} requires an argument" >&2
		exit 1
		;;
	esac
done
shift $((OPTIND-1))

### Define name of sample to be used as prefix
name=${fq%.fastq.gz}
name=${name%.fq.gz}
name=${name%.fq}
name=${name%.fastq}

### Create output directory
mkdir ref_based_assembly_${name}

### Output arguments specified to terminal stdout
echo "";
echo -e "\033[31m Reference-based assembly for metagenomics dataset\033[0m";
echo -e "\033[31m by Majo Galarion\033[0m";
echo "";
echo -e "Input file: \033[33m${fq}\033[0m";
echo -e "Number of threads: \033[33m${threads}\033[0m";
echo -e "Reference file: \033[33m${ref}\033[0m";
echo -e "Medaka model: \033[33m${medakaModel}\033[0m";
echo -e "Min read length: \033[33m${minLength}\033[0m";
echo -e "Max read length: \033[33m${maxLength}\033[0m";
echo -e "Minimum depth: \033[33m${minDepth}\033[0m";
echo -e "Minimum read Q-score: \033[33m${minQual}\033[0m"


### Size selection with prinseq-lite
#echo "";
#echo -e "\033[32m Step 1/7: Size selection with prinseq-lite...\033[0m"
#echo "";
	#${prinseq} -fastq ${fq} -min_len ${minLength} -min_qual_mean ${minQual} -out_format "3" -out_good ref_based_assembly_${name}/${name}.prinseq
	#rm -f *_prinseq_bad_*

### Get read length distribution on raw reads before and after porechop and prinseq
	#cat ${fq} | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ref_based_assembly_${name}/${fq}_read-len-dist.csv
	#cat ref_based_assembly_${name}/${name}.prinseq.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > ref_based_assembly_${name}/${name}.prinseq_read-len-dist.csv

### Reference-based assembly
echo "";
echo -e "\033[32m Step 1/6: Reference-based assembly with Minimap2...\033[0m"
echo "";
	${minimap2} -t $threads -ax map-ont ${ref} ${fq} | ${samtools} view -bS - | ${samtools} sort -o ref_based_assembly_${name}/${name}.srt.bam
	${samtools} index ref_based_assembly_${name}/${name}.srt.bam
	${samtools} idxstats ref_based_assembly_${name}/${name}.srt.bam | sort -n -r -k3 > ref_based_assembly_${name}/${name}_alignment_sumstats.txt

### Generate VCF file via bcftools mpileup and filter positions with specified minimum depth
echo "";
echo -e "\033[32m Step 2/6: Generate VCF file...\033[0m"
echo "";
	${bcftools} mpileup --threads ${threads} -Ou -d 1000000 -f ${ref} ref_based_assembly_${name}/${name}.srt.bam | ${bcftools} call -c -Oz -o ref_based_assembly_${name}/${name}.vcf.gz --ploidy 1
	${bcftools} index ref_based_assembly_${name}/${name}.vcf.gz


### Generate consensus using the VCF file and mask positions with N if depth is below minimum depth
echo "";
echo -e "\033[32m Step 3/6: Generate draft consensus...\033[0m"
echo "";
	${bcftools} consensus -f ${ref} ref_based_assembly_${name}/${name}.vcf.gz > ref_based_assembly_${name}/${name}.draft.consensus1.fasta

### Extract mapped-only, clipped reads from original BAM file and convert to FASTQ file
echo "";
echo -e "\033[32m Step 4/7: Extract mapped-only reads...\033[0m"
echo "";
	${samtools} view -@ 8 -b -F 256 ref_based_assembly_${name}/${name}.srt.bam | ${samtools} fastq > ref_based_assembly_${name}/${name}.mapped.fastq

### Run medaka to generate a polished consensus
echo "";
echo -e "\033[32m Step 5/6: Polishing with medaka with model ${medakaModel}...\033[0m"
echo "";
	${medaka}
	medaka_consensus -i ref_based_assembly_${name}/${name}.mapped.fastq -d ref_based_assembly_${name}/${name}.draft.consensus1.fasta -o ref_based_assembly_${name}/medaka-consensus_log -t ${threads} -m ${medakaModel}
	deactivate
	
	echo ">${name}_medaka" > ref_based_assembly_${name}/${name}.medaka.consensus.fasta
	grep -v ">" ref_based_assembly_${name}/medaka-consensus_log/consensus.fasta >> ref_based_assembly_${name}/${name}.medaka.consensus.fasta
	
	echo ">${name}_draft" > ref_based_assembly_${name}/${name}.draft.consensus.fasta
	grep -v ">" ref_based_assembly_${name}/${name}.draft.consensus1.fasta >> ref_based_assembly_${name}/${name}.draft.consensus.fasta
	rm ref_based_assembly_${name}/${name}.draft.consensus1.fasta

### Mask low depth positions with N
echo ">${name}_medaka_masked" > ref_based_assembly_${name}/${name}.medaka.consensus.fasta
samtools depth -d 0 -a ref_based_assembly_${name}/medaka-consensus_log/calls_to_draft.bam > ref_based_assembly_${name}/depth
paste ref_based_assembly_${name}/depth <(grep -v ">" ref_based_assembly_${name}/medaka-consensus_log/consensus.fasta | tr -d '\n' | grep -o . | \
	cat -n -) | \
	awk -v minD=${minDepth} '{if($2==$4){if ($3>=minD) printf "%c", $5; else printf "%c","N";} else print "N"}'| \
	fold -w 70 >> ref_based_assembly_${name}/${name}.medaka.consensus.fasta 
rm ref_based_assembly_${name}/depth



### Get stats for plotting
echo "";
echo -e "\033[32m Step 6/6: Getting alignment stats for plotting...\033[0m"
echo "";
	mkdir ref_based_assembly_${name}/alignment_stats
	mv ref_based_assembly_${name}/${name}_alignment_sumstats.txt ref_based_assembly_${name}/alignment_stats

	sed -i '5d' ref_based_assembly_${name}/alignment_stats/${name}_alignment_sumstats.txt
	cut -f1,2,3 ref_based_assembly_${name}/alignment_stats/${name}_alignment_sumstats.txt > ref_based_assembly_${name}/alignment_stats/${name}_alignment_summary_stats.txt
	rm ref_based_assembly_${name}/alignment_stats/${name}_alignment_sumstats.txt


grep ">" ${ref} | cut -c2- | while read -r line;
do
	# get summary stats of fastq files
	ref_length=$(grep ${line} ref_based_assembly_${name}/alignment_stats/${name}_alignment_summary_stats.txt | cut -f2)
	reads=$(cat ${fq} | wc -l | awk '{print $1/4}')
	#filtered_reads=$(cat ref_based_assembly_${name}/${name}.prinseq.fastq | wc -l | awk '{print $1/4}')
	mapped_reads=$(samtools view -@ $threads -c -F 4 ref_based_assembly_${name}/${name}.srt.bam "${line}")
	cov_1x=$(samtools depth -a -r ${line} ref_based_assembly_${name}/${name}.srt.bam | awk '$3 >=1' | wc -l)
	perc_cov_1x=$(awk -v rl=${ref_length} -v rc=${cov_1x} 'BEGIN {print (rc/rl*100)}')
	cov_20x=$(samtools depth -a -r ${line} ref_based_assembly_${name}/${name}.srt.bam | awk '$3 >=20' | wc -l)
	perc_cov_20x=$(awk -v rl=${ref_length} -v rc=${cov_20x} 'BEGIN {print (rc/rl*100)}')
	mean_depth=$(samtools depth -d 0 ref_based_assembly_${name}/${name}.srt.bam | awk '{sum+=$3} END { print sum/NR}')

	echo ${line} ${ref_length} ${reads} ${mapped_reads} ${cov_1x} ${perc_cov_1x} ${cov_20x} ${perc_cov_20x} ${mean_depth} | tr " " "\t" > ref_based_assembly_${name}/alignment_stats/${name}-${line}-stats.temp
done

cat ref_based_assembly_${name}/alignment_stats/*.temp > ref_based_assembly_${name}/alignment_stats/${name}-stats.txt
sed -i '1i Reference\tSize\tReads\tMapped_reads\tCoverage_1X\tPercent_cov_1X\tCoverage_20X\tPercent_cov_20X\tMean_depth' ref_based_assembly_${name}/alignment_stats/${name}-stats.txt


grep ">" ${ref} | cut -c2- | while read -r line;
do
	# get depth per position on clipped file; also include references that were not used
	samtools depth -d 0 -a -r ${line} ref_based_assembly_${name}/${name}.srt.bam > ref_based_assembly_${name}/alignment_stats/${name}.${line}.per-base-dep.txt
done


# clean-up; remove temp files
rm ref_based_assembly_${name}/alignment_stats/*.temp
rm ref_based_assembly_${name}/alignment_stats/${name}_alignment_summary_stats.txt
#rm ref_based_assembly_${name}/${name}.srt.bam*
rm ref_based_assembly_${name}/*.fastq
rm ref_based_assembly_${name}/*.mmi
rm ref_based_assembly_${name}/*.fai


echo "";
echo -e "\033[31m Analysis finished!\033[0m"
echo "";
