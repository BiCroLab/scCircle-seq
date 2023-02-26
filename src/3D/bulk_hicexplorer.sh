#!/usr/bin/bash
GENOME="/share/home/chenjx/index/bwa_index/hg38/hg38.fa"
fastq1=""
fastq2=""
OUTPUT=""
nThreads=8
PREFIX=""
digest='/share/home/chenjx/index/bwa_index/hg38/MboI_site_positions.bed'

# extract options and their arguments into variables
while getopts "f:r:o:t:p:" options; do
    case ${options} in 
        f) fastq1=${OPTARG};;
        r) fastq2=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        t) nThreads=${OPTARG};;
        p) PREFIX=${OPTARG};;
        *) exit 1;;
    esac
done

[[ -z ${OUTPUT} ]] && { echo "OUTPUT: NONE!"; exit 1; }
[ -d ${OUTPUT} ] || { mkdir -p ${OUTPUT}; }
#HiCExplorer
bwa mem -A1 -B4  -E50 -L0 -t ${nThreads} ${GENOME} \
    ${fastq1} 2>> ${OUTPUT}/${PREFIX}_R1.log | samtools view -Shb - > ${OUTPUT}/${PREFIX}_R1.bam
bwa mem -A1 -B4  -E50 -L0 -t ${nThreads} ${GENOME} \
    ${fastq2} 2>> ${OUTPUT}/${PREFIX}_R2.log | samtools view -Shb - > ${OUTPUT}/${PREFIX}_R2.bam
mkdir ${OUTPUT}/${PREFIX}_10k
hicBuildMatrix --samFiles ${OUTPUT}/${PREFIX}_R1.bam ${OUTPUT}/${PREFIX}_R2.bam \
                 --binSize 10000 \
                 --restrictionSequence GATC \
                 --danglingSequence GATC \
                 --restrictionCutFile ${digest} \
                 --threads ${nThreads} \
                 --inputBufferSize 100000 \
                 --outBam ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hic.bam \
                 -o ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hic_10k_matrix.h5 \
                 --QCfolder ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hicQC
hicCorrectMatrix correct -m ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hic_10k_matrix.h5 --filterThreshold -1.5 5 -o ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hic_10k_corrected_matrix.h5
mkdir ${OUTPUT}/${PREFIX}_20k
hicMergeMatrixBins -m ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hic_10k_matrix.h5 -o ${OUTPUT}/${PREFIX}_20k/${PREFIX}_hic_20k_matrix.h5 -nb 2
hicCorrectMatrix correct -m ${OUTPUT}/${PREFIX}_20k/${PREFIX}_hic_20k_matrix.h5 --filterThreshold -1.5 5 -o ${OUTPUT}/${PREFIX}_20k/${PREFIX}_hic_20k_corrected_matrix.h5
mkdir ${OUTPUT}/${PREFIX}_100k
hicMergeMatrixBins -m ${OUTPUT}/${PREFIX}_10k/${PREFIX}_hic_10k_matrix.h5 -o ${OUTPUT}/${PREFIX}_100k/${PREFIX}_hic_100k_matrix.h5 -nb 10
hicCorrectMatrix correct -m ${OUTPUT}/${PREFIX}_100k/${PREFIX}_hic_100k_matrix.h5 --filterThreshold -1.5 5 -o ${OUTPUT}/${PREFIX}_100k/${PREFIX}_hic_100k_corrected_matrix.h5
