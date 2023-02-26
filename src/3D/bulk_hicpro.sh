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
mkdir -p ${OUTPUT}/${PREFIX}_HiCPro
r1_format=${${fastq1#*_}%%.*}
r2_format=${${fastq2#*_}%%.*}
data_path=${fastq1%%/*}
#HiCExplorer
cat << EOF > ${OUTPUT}/config-hicpro.txt
# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = ${nThreads}
SORT_RAM = 768M
LOGFILE = hicpro.log

JOB_NAME = ${PREFIX}_HiCPro
JOB_MEM = 20G
JOB_WALLTIME = 1
JOB_QUEUE = 1
JOB_MAIL = 1

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _${r1_format}
PAIR2_EXT = _${r2_format}

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10

BOWTIE2_IDX_PATH = /share/home/chenjx/index/bowtie2_index/GRCh38/hg38
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = hg38
GENOME_SIZE = /share/home/chenjx/index/bwa_index/hg38/hg38.chrom.sizes

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Capture Hi-C analysis
#######################################################################

CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /share/home/chenjx/index/bwa_index/hg38/MboI_site_positions.bed
LIGATION_SITE = 
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 0
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 20000 40000 150000 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1

EOF

cd ${OUTPUT}
HiC-Pro -i $data_path -o ${PREFIX}_HiCPro -c config-hicpro.txt -p
sbatch HiCPro_step1.sh
# sbatch HiCPro_step2.sh
