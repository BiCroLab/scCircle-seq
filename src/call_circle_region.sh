#!/bin/bash

## The code below is to find enriched regions in Circle-seq data from a paired end Illumina run
## (Long reads, both PacBio and Nanopore, are processed differenty, as are PE reads from WGS data.)

## The only inputs are
## 1) position-sorted bam file from an alignment using bwa mem
## 2) output directory (can be just the name or the full path+name)


# With 10 threads as above, most deeply sequenced Circle-seq bam files take ~30 minutes to run everything below.


###### Software dependencies.
### Here, including paths to minimize confusion with alt installs within the CER.
samtools=/usr/bin/samtools
java=/usr/bin/java # version 1.8 for compatibility with picard
picard=/home/jinxin/tools/picard.jar # using 2.18 on our cluster
homer=/home/jinxin/tools/homer
bedtools=/usr/local/bin/bedtools
pgltools=/home/jinxin/tools/pgltools/sh/pgltools # v2.2.0; https://github.com/billgreenwald/pgltools

# The above are called directly using the name, but deeptools must be in path
# Functions used: bamCoverage, computeMatrix, plotProfile, plotHeatmap
deeptools=~/anaconda3/bin/deeptools

###### Params set for neuroblastoma Circle-seq data, but can be adjusted as needed
nthreads=40 # for processes which can be multithreaded
qfilt=20 # MAPQ filter for reads
enrichFDR=0.001 # FDR for initial segment enrichment step
mergedist=3000 # bp for initial enriched region merging; based on tests in NB data
circEdgeThresh=$3 # threshold for the number of circle-supporting (split and/or outward facing) reads at edge of each putative circle 


echo "=========================================================="
echo "Starting on : $(date)"
echo "Current directory : $(pwd)"
echo "=========================================================="


###### Input params: bam file + folder

inbam=$1
inbamnodir=${inbam##*/} # rm path for renaming
inbambase=${inbamnodir/%.bam/}    ###/% / is for cutting some character out from the string

outdir=$2

mkdir -p $outdir
cd $outdir;


###### Find enriched regions
### First, simply find blocks, no split/outward facing read requirements

# Tag Dir
tagdir=TagDir_${inbambase}

mkdir $tagdir
$homer/bin/makeTagDirectory $tagdir -format sam -mapq $qfilt -rmsoft $inbam

# Enriched regions file name:
peakfile=Regions_${inbambase}.txt

$homer/bin/findPeaks $tagdir -o $peakfile -style histone -fdr $enrichFDR -regionRes 20

# Output enriched regions in bed format
# Mostly used for QC with normalized (to 1e7) tag count included
peakbed=${peakfile/%.txt/.orig.bed}
grep -v "^#" $peakfile | cut -f2-4,6 > $peakbed

# Output merged regions of enrichment
# Edges may not be perfect due to where enrichment block falls, but works in most cases and can be corrected below
mergedbed=${peakfile/%.txt/.enriched.merged.orig.bed}
grep -v "^#" $peakfile | cut -f2-4 | grep -v "random\|chrUn" | sort -k1,1 -k2,2n | $bedtools merge -d $mergedist -i stdin > $mergedbed


# Output merged reads
# (sometimes more accurate delineation of circle junctions, but can fall into trap of cascading reads outside of enriched region)
# First get bam to bed
bam2bed=${inbamnodir/%.bam/.bam2bed.bed}
$samtools view --threads $nthreads -bq $qfilt $inbam | $bedtools bamtobed -i stdin -splitD | cut -f1-3 | grep -v "random\|chrUn" | sort -k1,1 -k2,2n | bedtools merge -d 250 -i stdin > $bam2bed

# Then get bam2bed overlap with enriched blocks, but only keep segments with >5 reads
# (this is for edge fine-tuning, not circle calling)
enrichedbam2bed=${peakfile/%.txt/.enrichedbam2bed.merged.bed}
$bedtools intersect -u -a $bam2bed -b $mergedbed | $bedtools multicov -q $qfilt -bams $inbam -bed stdin | awk '{ if($4 >= 5) print }'  > $enrichedbam2bed



### Fine-tune circle edges
# Note: this is often required due to enrichment probabilities tapering off at edges when calling enlarged blocks.
# Will then use outward-facing and split reads below to threshold

# Use bedtools closest so need only get small window on left and right sides
# 'Left' edge first
# bam2bed
bam2bedLeft=${enrichedbam2bed/%.bed/.tmpLeft.bed}
cut -f1,2 $enrichedbam2bed | awk -v OFS='\t' '{ print $1,$2,$2+10}' > $bam2bedLeft

# enriched left
mergedbedLeft=${mergedbed/%.bed/.tmpLeft.bed}
cut -f1,2 $mergedbed | awk -v OFS='\t' '{ print $1,$2,$2+10}' > $mergedbedLeft

# refine left
mergedbedLeftb2b=${mergedbedLeft/%.bed/.B2B.bed}
$bedtools closest -d -t first -a $mergedbedLeft -b $bam2bedLeft | awk -v OFS='\t' '{ print $1,$5 }' > $mergedbedLeftb2b

# 'Right' edge next
# bam2bed
bam2bedRight=${enrichedbam2bed/%.bed/.tmpRight.bed}
cut -f1,3 $enrichedbam2bed | awk -v OFS='\t' '{ print $1,$2-10,$2}' > $bam2bedRight

# enriched right
mergedbedRight=${mergedbed/%.bed/.tmpRight.bed}
cut -f1,3 $mergedbed | awk -v OFS='\t' '{ print $1,$2-10,$2}' > $mergedbedRight

# refine right
mergedbedRightb2b=${mergedbedRight/%.bed/.B2B.bed}
$bedtools closest -d -t last -a $mergedbedRight -b $bam2bedRight | awk -v OFS='\t' '{ print $1,$6 }' > $mergedbedRightb2b

# combine to get enriched circle region with modified edges:
mergedbedRefine=${mergedbed/%.orig.bed/.b2bRefined.bed}
paste $mergedbedLeftb2b $mergedbedRightb2b | cut -f1,2,4 | awk -v OFS='\t' '{dist=$3-$2; if (dist < 0) { print $1,$3,$2 } else { print $0 }}' | sort -k1,1 -k2,2n | bedtools merge -d 1000 -i stdin > $mergedbedRefine




###### Filter bam for potential circle-supporting reads
### First, outward facing

# Note about filter below: samtools flag -F is an OR statement (unlike -f, which is AND)
# Instead break down components within awk, using bitwise operator ('and') and direct reading/testing of sam flag
# Skip if following true: 1 = paired, 16 = read on rev, 32 = mate on rev
# This should get rid of R1R2 and R2R1 read pairs.

# Temp bams for storage of each class
# F1 of R2F1; flag 97, neg
v97f1=${inbamnodir/%.bam/.F1neg.bam}
# R2 of R2F1; flag 145, pos
v145r2=${inbamnodir/%.bam/.R2pos.bam}
# R1 of R1F2; flag 81, pos
v81r1=${inbamnodir/%.bam/.R1pos.bam}
# F2 of R1F2; flag 161, neg
v161f2=${inbamnodir/%.bam/.F2neg.bam}

#R2,F1; here read 2 with NEG TLEN
# F1=v97f1
$samtools view --threads $nthreads -h -f 97 -q $qfilt $inbam | awk '{if ($1 ~ /^@/) { print $_ } else if (and($2,1) && and($2,16) && and($2,32)) { next } else if ($9 < 0) {print $_} }' | $samtools view --threads $nthreads -bS -o $v97f1 -

# next v145r2
$samtools view --threads $nthreads -h -f 145 -q $qfilt $inbam | awk '{if ($1 ~ /^@/) { print $_ } else if (and($2,1) && and($2,16) && and($2,32)) { next } else if ($9 > 0) {print $_} }' | $samtools view --threads $nthreads -bS -o $v145r2 -

# v81r1
$samtools view --threads $nthreads -h -f 81 -q $qfilt $inbam | awk '{if ($1 ~ /^@/) { print $_ } else if (and($2,1) && and($2,16) && and($2,32)) { next } else if ($9 > 0) {print $_} }' | $samtools view --threads $nthreads -bS -o $v81r1 -

# v161f2
$samtools view --threads $nthreads -h -f 161 -q $qfilt $inbam | awk '{if ($1 ~ /^@/) { print $_ } else if (and($2,1) && and($2,16) && and($2,32)) { next } else if ($9 < 0) {print $_} }' | $samtools view --threads $nthreads -bS -o $v161f2 -


### Second, add split reads
# note sam flag distinction:
# secondary ('not primary') reads = the same part of sequence aligns to multiple locations
# supplementary (chimeric/split reads) = where (mostly) non-overlapping parts of a sequence align to multiple locations
# supp sam flag value = 2048, with any combination thereof

vSplit=${inbamnodir/%.bam/.split.bam}

$samtools view --threads $nthreads -h -q $qfilt $inbam | awk '{ if(/^@/) print; else if (/SA:Z:/) print }' | $samtools view --threads $nthreads -bS -o $vSplit -




###### Merge all circle-supporting reads into a single bam file
# output name
vmerge=${inbamnodir/%.bam/.SplitOutface.merged.bam}

# infile bams to be merged
$java -Xms500m -Xmx6g -jar $picard MergeSamFiles AS=true USE_THREADING=true OUTPUT=$vmerge INPUT=$v97f1 INPUT=$v145r2 INPUT=$v81r1 INPUT=$v161f2 INPUT=$vSplit

# Index:
$samtools index $vmerge


# Cleanup 1:
# Note: if merged bam created successfully, should be able to delete interim bams:
# -s = file exists and not empty
if [ -s "$vmerge" ]; then
    # alignment files
    rm $v97f1
    rm $v145r2
    rm $v81r1
    rm $v161f2
    rm $vSplit
    # tag dir (stepwise for safer removal, avoid rf flag
    rm $tagdir/*
    rmdir $tagdir
    # region files and tagdir
    rm $peakbed $bam2bed $bam2bedLeft $mergedbedLeft $mergedbedLeftb2b $bam2bedRight $mergedbedRight $mergedbedRightb2b
fi



###### Get circle junction-supporting reads

# Circle 'edge' coordinates
# here use +/-N bp for reads around putative junctions
enrichededges=${mergedbedRefine/%.bed/.edges.bed}
edgewin=300 # +/-N bp around junctions
awk -v OFS='\t' "{if(\$2 - $edgewin < 0) start1=0; else start1=\$2 - $edgewin; print \$1,start1,\$2+$edgewin; print \$1,\$3 - $edgewin,\$3 + $edgewin"} $mergedbedRefine > $enrichededges

bamedge=${vmerge/%.bam/.Edges.bam}
$bedtools intersect -u -a $vmerge -b $enrichededges > $bamedge

$samtools index $bamedge





###### Threshold based on junction-supporting split and outward facing reads
# intersect entire enriched region with circle-supporting edge bam from above
# also count all reads over this segment, to examine ratio downstream
# So edge counts = column 4, all counts = column 5

# Note multicov discrepancy when compared with featureCounts:
# multicov won't collapse counts by read name, which may be
# desirable here since multiple split and outward-facing reads
# at different circle junctions can lend support in different ways.
# If not desired, can create SAF from bed and run featureCounts in its place. 
edgecount=${mergedbedRefine/%.bed/.Counts.bed}
$bedtools multicov -q $qfilt -bams $bamedge $inbam -bed $mergedbedRefine > $edgecount

edgethreshOut=${edgecount/%.bed/.ThreshFinal.bed}
# uses $circEdgeThresh from above (usually 2)
awk "{ if(\$4 >= $circEdgeThresh) print }" $edgecount > $edgethreshOut






###### Create bedpe file of all chimeric reads within circle-enriched windows
# Note: reads must pass qfilt 20 threshold on both ends
# Requires at least 2 separate chimeric reads within 500 bp to score junction
# This utlizes the pgltools package for paired genomic loci
# Note: split reads are marked at both primary read and officially documented supplementary read
# => so to avoid counting twice, collapse all duplicates with uniq after bedpe formatting.

chimericBedpe=${mergedbedRefine/%.bed/.Chimeric.bedpe}
chimericBedpe_chr8=${chimericBedpe/%.Chimeric.bedpe/.Chimeric_chr8.bedpe}
$bedtools intersect -u -a $vmerge -b $edgethreshOut | $samtools view - | awk -v OFS='\t' 'match($0, /SA:Z:[^[:space:]]+/) { SAfield=substr($0,RSTART+5,RLENGTH-5); split(SAfield,chim,/\;/); for (i=1; i <= length(chim); i++) { split(chim[i], coordinfo,/,/); if( coordinfo[5] >= 20) { print $3,$4,$4+10, coordinfo[1], coordinfo[2], coordinfo[2]+10} } }' | $pgltools formatbedpe | uniq | awk '{print $0 "\t."}' | $pgltools merge -noH -stdInA -d 500 -c 7 -o count | awk '{ if($7 >= 2) print $0 }' > $chimericBedpe
cat ${chimericBedpe} | grep "chr8" > ${chimericBedpe_chr8}





###### Visualization of enriched regions and circle-supporting reads (calls to deepTools)

### Create bigwigs of all reads as well as circle-junction supporting reads
# all reads
bwall=${inbamnodir/%.bam/.bw}
bamCoverage --extendReads 0 --minMappingQuality $qfilt --ignoreDuplicates --binSize 25 --numberOfProcessors $nthreads --scaleFactor 1 --normalizeUsing None -b $inbam -o $bwall

# circle edge reads
bwedge=${bamedge/%.bam/.bw}
bamCoverage --extendReads 0 --minMappingQuality $qfilt --ignoreDuplicates --binSize 25 --numberOfProcessors $nthreads --scaleFactor 1 --normalizeUsing None -b $bamedge -o $bwedge

# all outward facing and split reads
# NOTE: while not always at edges, these reads are critical for matching internal circle rearrangements
bwcirc=${vmerge/%.bam/.bw}

bamCoverage --extendReads 0 --minMappingQuality $qfilt --ignoreDuplicates --binSize 25 --numberOfProcessors $nthreads --scaleFactor 1 --normalizeUsing None -b $vmerge -o $bwcirc


### Meta plots and heat maps of above

## First, chrM plots:
# get chrM coords for this genome assembly, based on bam header:
chrMcoord="chrM.bed"
$samtools view -H $inbam | grep "SN:chrM" | awk -v OFS='\t' '{ split($3,a,":"); print "chrM",1,a[2]}' > $chrMcoord

# meta plots
# create data matrix (one site, so here a simple vector)
datamatChrM="DataMat_${inbamnodir}.chrM.tab.gz"
computeMatrix scale-regions -R chrM.bed -S $bwall $bwedge --beforeRegionStartLength 0 --regionBodyLength 3000 --afterRegionStartLength 0 --missingDataAsZero -out $datamatChrM

# plot chrM: all and circle-supporting reads
outChrMsep="Fig_${inbamnodir}.chrM.meta.sep.pdf"
outChrMgroup="Fig_${inbamnodir}.chrM.meta.group.pdf"

plotProfile -m $datamatChrM --yMax 5000 500 --colors dodgerblue red --samplesLabel AllReads CircEdges --startLabel Start --endLabel End --regionsLabel '' --yAxisLabel readCount --plotTitle chrM -out $outChrMsep
plotProfile -m $datamatChrM --perGroup --colors dodgerblue red --samplesLabel AllReads CircEdges --startLabel Start --endLabel End --regionsLabel '' --yAxisLabel readCount --plotTitle chrM -out $outChrMgroup


## Next: plots of all enriched regions passing filter
# create data matrix
datamatAllPF="DataMat_${inbamnodir}.Enriched.tab.gz"
computeMatrix scale-regions -R $edgethreshOut -S $bwall $bwedge --beforeRegionStartLength 500 --regionBodyLength 2000 --afterRegionStartLength 500 --missingDataAsZero -out $datamatAllPF

# meta plots
outCrcsep="Fig_${inbamnodir}.Circles.meta.sep.pdf"
plotProfile -m $datamatAllPF --yMax 30 10 --colors dodgerblue firebrick --samplesLabel AllReads CircEdges --startLabel Start --endLabel End --regionsLabel '' --yAxisLabel readCount --plotTitle AllCircles -out $outCrcsep
outCrcgroup="Fig_${inbamnodir}.Circles.meta.group.pdf"
plotProfile -m $datamatAllPF --perGroup --colors dodgerblue firebrick --samplesLabel AllReads CircEdges --startLabel Start --endLabel End --regionsLabel '' --yAxisLabel readCount --plotTitle AllCircles -out $outCrcgroup

# heat maps
outCrcHeat="Fig_${inbamnodir}.Circles.Heatmap.pdf"
plotHeatmap -m $datamatAllPF  --missingDataColor 'white'  --colorList 'white,dodgerblue' 'white,firebrick' --zMax 20 5 --heatmapWidth 6 --samplesLabel AllReads CircEdges --startLabel Start --endLabel End --regionsLabel '' --yAxisLabel readCount -out $outCrcHeat





####### Clean up 2
if [ -s "$edgecount" ]; then
    rm $enrichededges $mergedbedRefine $chrMcoord
fi

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="



