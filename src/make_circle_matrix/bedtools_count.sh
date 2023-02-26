#!/usr/bin/bash
bam=''
name=''

# extract options and their arguments into variables
while getopts "b:n:" options; do
    case ${options} in 
        b) bam=${OPTARG};;
        n) name=${OPTARG};;
        *) exit 1;;
    esac
done
mkdir new_count

bedtools intersect \
           -a for_count.bed -b ${bam} \
           -c \
           -wa  > new_count/${name}.bed
Rscript read.r new_count/${name}.bed 
cut -d "," -f4 new_count/${name}.bed  > new_count/${name}_new.bed
sed -i "1s/^/${name}\n/" new_count/${name}_new.bed

