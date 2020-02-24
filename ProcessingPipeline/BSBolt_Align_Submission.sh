#!/bin/bash

# BSBolt SGE alignment shell
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

echo python3 -m  BSBolt Align -D -BT2-p 8 -DB ${index} -F1 ${file_directory}${sample}_1.fastq -F2 ${file_directory}${sample}_2.fastq -O ${output_folder}${sample} -BT2-local -BT2-score-min L,160,0 -discord


python3 -m BSBolt Align -D -BT2-p 8 -DB ${index} -F1 ${file_directory}${sample}_1.fastq -F2 ${file_directory}${sample}_2.fastq -O ${output_folder}${sample} -BT2-local -discord  -BT2-score-min L,160,0 > ${output_folder}${sample}.log
