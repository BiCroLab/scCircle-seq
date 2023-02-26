#!/usr/bin/bash
genome='GRCh38'
fastq1=""
fastq2=""
OUTPUT=""
nThreads=8
PREFIX=""
Deconvolve=''
threshold=0
trimmomatic='/home/jinxin/anaconda3/bin/trimmomatic'
picard='/home/jinxin/tools/picard.jar'
linear_spike_in='/home/jinxin/projects/ecDNA_Project/scripts/src/spike_in/linear_spike_in/lamda_dna.fa'
circular_spike_in='/home/jinxin/projects/ecDNA_Project/scripts/src/spike_in/circular_spike_in/RBD_plasmid.fa'
myc_ecDNA='/home/jinxin/projects/ecDNA_Project/scripts/src/colo320dm/COLO320DM_Version2_ecDNA1_cycles_cycle1.bed'
Circlebed2AAbed='/home/jinxin/projects/ecDNA_Project/scripts/src/Circlebed2AAbed.R'
circleEnrichFilter='/home/jinxin/projects/ecDNA_Project/scripts/src/call_circle_region.sh'
pgltools='/home/jinxin/tools/pgltools/sh/pgltools' # v2.2.0; https://github.com/billgreenwald/pgltools
bedtools='/usr/local/bin/bedtools'
# extract options and their arguments into variables
while getopts "f:r:o:t:p:e:D:g:" options; do
    case ${options} in 
        f) fastq1=${OPTARG};;
        r) fastq2=${OPTARG};;
        o) OUTPUT=${OPTARG};;
        t) nThreads=${OPTARG};;
        p) PREFIX=${OPTARG};;
        D) Deconvolve=${OPTARG};;
        e) threshold=${OPTARG};;
        g) genome=${OPTARG};;
        *) exit 1;;
    esac
done

if [ ${genome} == 'hg19' ];then
   GENOME='/media/bs2-pro/jinxin/reference/hg19/hg19.fa.gz'
   PREFIX=hg19_${PREFIX}
   OUTPUT=hg19_${OUTPUT}
   ref='hg19'
elif [ ${genome} == 'GRCh38' ];then
   GENOME='/media/bs2-pro/jinxin/reference/hg38/ncbi-genomes-2022-10-13/GCF_000001405.40_GRCh38.p14_genomic.fna.gz'
   ref='GRCh38'
   PREFIX=hg38_${PREFIX}
   OUTPUT=hg38_${OUTPUT}
fi


[[ -z ${OUTPUT} ]] && { echo "OUTPUT: NONE!"; exit 1; }
[ -d ${OUTPUT} ] || { mkdir -p ${OUTPUT}; }
mkdir -p ${OUTPUT}/${PREFIX}_trimed_fastq

# Trim the nextera adaptor in case it case false-positive edge reads
# ${trimmomatic} PE -threads ${nThreads} -phred33 ${fastq1} ${fastq2} \
# ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R1.fq.gz ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R1.unmapped.fq.gz ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R2.fq.gz ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R2.unmapped.fq.gz \
# ILLUMINACLIP:/home/jinxin/tools/Trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:51

# Read Mapping 
#conda activate AA
# bwa mem -q -t ${nThreads} ${GENOME} ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R1.fq.gz ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R2.fq.gz | gzip > ${OUTPUT}/${PREFIX}.sam.gz
bwa mem -q -t ${nThreads} ${linear_spike_in} ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R1.fq.gz ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R2.fq.gz | gzip > ${OUTPUT}/${PREFIX}_linear_spike_in.sam.gz
bwa mem -q -t ${nThreads} ${circular_spike_in} ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R1.fq.gz ${OUTPUT}/${PREFIX}_trimed_fastq/${PREFIX}_R2.fq.gz | gzip > ${OUTPUT}/${PREFIX}_circular_spike_in.sam.gz
samtools sort -@ ${nThreads} -o ${OUTPUT}/${PREFIX}.Sort.bam ${OUTPUT}/${PREFIX}.sam.gz
samtools index -@ ${nThreads} ${OUTPUT}/${PREFIX}.Sort.bam ${OUTPUT}/${PREFIX}.Sort.bam.bai
java -jar ${picard} MarkDuplicates \
I=${OUTPUT}/${PREFIX}.Sort.bam \
O=${OUTPUT}/${PREFIX}.dedup.bam \
M=${OUTPUT}/${PREFIX}.dedup.txt
samtools index -@ ${nThreads} ${OUTPUT}/${PREFIX}.dedup.bam ${OUTPUT}/${PREFIX}.dedup.bam.bai

# Spike-in QC
samtools sort -@ ${nThreads} -o ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam ${OUTPUT}/${PREFIX}_linear_spike_in.sam.gz
samtools index -@ ${nThreads} ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam.bai
samtools flagstat -@ ${nThreads} ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam > ${OUTPUT}/${PREFIX}_linear_spike_in.log
samtools depth -@ ${nThreads} -o ${OUTPUT}/${PREFIX}_linear_spike_in.depth  ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam 

samtools sort -@ ${nThreads} -o ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam ${OUTPUT}/${PREFIX}_circular_spike_in.sam.gz
samtools index -@ ${nThreads} ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam.bai
samtools flagstat -@ ${nThreads} ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam > ${OUTPUT}/${PREFIX}_circular_spike_in.log
samtools depth -@ ${nThreads} -o ${OUTPUT}/${PREFIX}_circular_spike_in.depth  ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam 

rm ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam 
rm ${OUTPUT}/${PREFIX}_circular_spike_in.sam.gz
rm ${OUTPUT}/${PREFIX}_circular_spike_in.Sort.bam.bai
rm ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam 
rm ${OUTPUT}/${PREFIX}_linear_spike_in.sam.gz
rm ${OUTPUT}/${PREFIX}_linear_spike_in.Sort.bam.bai

# Basic QC
   #genome coverage
samtools flagstat -@ ${nThreads} ${OUTPUT}/${PREFIX}.dedup.bam > ${OUTPUT}/${PREFIX}_map.log
bedtools genomecov -bg -ibam ${OUTPUT}/${PREFIX}.dedup.bam > ${OUTPUT}/cov.log
a=$(cat ${OUTPUT}/cov.log | awk -F '\t' '{print($3-$2)}' | awk '{sum +=$1};END {print sum}')
genome_cov_percent=$(echo "scale=6; ($a/3088458224)*100" | bc)
cat << EOF > ${OUTPUT}/${PREFIX}_genome_coverage_percent.log 
genome_coverage_percent(%): ${genome_cov_percent}
EOF
rm ${OUTPUT}/cov.log
   #mitochondria percentage
MT=$(echo "scale=6; 100 *($(samtools view -c ${OUTPUT}/${PREFIX}.dedup.bam chrM) / $(samtools view -c ${OUTPUT}/${PREFIX}.dedup.bam))" | bc -l)
cat << EOF > ${OUTPUT}/${PREFIX}_MT.log 
MT(%): ${MT}
EOF

# Call circle region&circle region QC
   #call circle region based on coverage
cd ${PWD}/${OUTPUT}/
/media/bs2-pro/jinxin/ecDNA_Project/scripts/src/call_circle_region.sh ${PWD}/${PREFIX}.dedup.bam ${PWD}/${PREFIX}_circle_filter_threshold_${threshold} ${threshold}
cd ..
   #circle read enrichment ratio
b=$(cat ${PWD}/${OUTPUT}/${PREFIX}_circle_filter_threshold_${threshold}/Regions_${PREFIX}.dedup.enriched.merged.b2bRefined.Counts.bed | awk -F '\t' '{print($5)}' | awk '{sum +=$1};END {print sum}')
c=$(cat ${OUTPUT}/${PREFIX}_map.log | awk "NR==1{print}" | cut -d ' ' -f 1) 
circle_read_percent=$(echo "scale=6;($b/$c)*100" | bc)
cat << EOF > ${OUTPUT}/${PREFIX}_circle_read_percent.log 
circle_read_percent(%): ${circle_read_percent}
EOF
  #extract the circle-region reads
samtools view -hb -L ${PWD}/${OUTPUT}/${PREFIX}_circle_filter_threshold_${threshold}/Regions_${PREFIX}.dedup.enriched.merged.b2bRefined.Counts.bed \
${PWD}/${OUTPUT}/${PREFIX}.dedup.bam  > ${PWD}/${OUTPUT}/${PREFIX}_target.bam
frac=$( samtools idxstats ${PWD}/${OUTPUT}/${PREFIX}_target.bam | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=1500000/total; if (frac > 1) {print 1} else {print frac}}' )
if [ $frac < 1];then
    samtools view -h -s $frac ${PWD}/${OUTPUT}/${PREFIX}_target.bam > ${PWD}/${OUTPUT}/${PREFIX}_target_subsample.bam
else
    samtools view -h ${PWD}/${OUTPUT}/${PREFIX}_target.bam > ${PWD}/${OUTPUT}/${PREFIX}_target_subsample.bam
fi

if [${Deconvolve} == 'TRUE' ];then
# Deconvolve the circle by Circle-Map
   #pre-define the file
source activate cnvkit
bam=${OUTPUT}/${PREFIX}.dedup.bam
qname_bam=${bam/%.bam/_qname.bam}
circular_read_candidates=${bam/%.bam/_circular_read_candidates.bam}
circular_read_candidates_sort=${bam/%.bam/_circular_read_candidates.Sort.bam}
circular_read_candidates_bw=${bam/%.bam/_circular_read_candidates.bw}
my_unknown_circle=${bam/%.bam/_CircleMap.bed}
chimericBedpe=${bam/%.bam/.complete_Chimeric.bedpe}
edgethreshOut=${OUTPUT}/${PREFIX}_circle_filter_threshold_${threshold}/*.b2bRefined.Counts.ThreshFinal.bed
   #Circle-map pipeline
samtools sort -@ ${nThreads} -n -o $qname_bam $bam
Circle-Map ReadExtractor -i $qname_bam -o $circular_read_candidates
samtools sort -@ ${nThreads} -o $circular_read_candidates_sort $circular_read_candidates
   #generate bigwig file for circle-supporting reads
samtools index -@ ${nThreads} $circular_read_candidates_sort
bamCoverage --extendReads 0 --minMappingQuality 20 --ignoreDuplicates --binSize 25 --numberOfProcessors ${nThreads} --scaleFactor 1 --normalizeUsing None -b $circular_read_candidates_sort -o $circular_read_candidates_bw
Circle-Map Realign -t ${nThreads} -i $circular_read_candidates_sort -qbam $qname_bam -sbam $bam -fasta ${GENOME} -o $my_unknown_circle
   #generate bedpe containging all discordant reads and split reads
$bedtools intersect -u -a $circular_read_candidates_sort -b $edgethreshOut | $samtools view - | awk -v OFS='\t' 'match($0, /SA:Z:[^[:space:]]+/) { SAfield=substr($0,RSTART+5,RLENGTH-5); split(SAfield,chim,/\;/); for (i=1; i <= length(chim); i++) { split(chim[i], coordinfo,/,/); if( coordinfo[5] >= 20) { print $3,$4,$4+10, coordinfo[1], coordinfo[2], coordinfo[2]+10} } }' | $pgltools formatbedpe | uniq | awk '{print $0 "\t."}' | $pgltools merge -noH -stdInA -d 500 -c 7 -o count | awk '{ if($7 >= 2) print $0 }' > $chimericBedpe
fi
