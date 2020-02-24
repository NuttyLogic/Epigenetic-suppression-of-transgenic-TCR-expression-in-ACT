#!/bin/bash

# methylation calling SGE shell
while getopts f:d:o:i: option
do
 case "${option}"
 in
 f) files_to_process=${OPTARG};;
 d) file_directory=${OPTARG};;
 o) output_folder=${OPTARG};;
 i) index=${OPTARG};;

 esac
done

mkdir -p $output_folder

cd ${output_folder}

. /u/local/Modules/default/init/modules.sh
module load python/3.7.2
module load samtools

sample=$(cat $files_to_process | head -${SGE_TASK_ID} | tail -1 )

echo samtools fixmate -p -m ${file_directory}${sample}.bam ${file_directory}${sample}.fixmates.bam
echo samtools sort -@ 8 -o ${file_directory}${sample}.sorted.bam ${file_directory}${sample}.fixmates.bam
echo rm ${file_directory}${sample}.fixmates.bam
echo samtools markdup ${file_directory}${sample}.sorted.bam ${file_directory}${sample}.dup.bam
echo rm ${file_directory}${sample}.sorted.bam
echo python3 -m BSBolt CallMethylation -min-qual 25 -t 8 -min 5 -0 ${file_directory}${sample} -I ${file_directory}${sample}.dup.bam -DB ${index} 


samtools fixmate -p -m ${file_directory}${sample}.bam ${file_directory}${sample}.fixmates.bam
samtools sort -@ 4 -o ${file_directory}${sample}.sorted.bam ${file_directory}${sample}.fixmates.bam
rm ${file_directory}${sample}.fixmates.bam
samtools markdup ${file_directory}${sample}.sorted.bam ${file_directory}${sample}.dup.bam
rm ${file_directory}${sample}.sorted.bam
samtools index ${file_directory}${sample}.dup.bam
python3 -m BSBolt CallMethylation -t 8 -min-qual 25 -min 5 -O ${file_directory}${sample} -I ${file_directory}${sample}.dup.bam -DB ${index} > ${file_directory}${sample}.meth.log
