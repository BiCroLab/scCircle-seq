library(cotools)
library(data.table)

merge <- read.table('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_stability/293T/293T_merge.soted.bed',sep='\t',header=FALSE) # input the files containing all the circle regions called from single cells
pseudobulk <- read.table('/dshare/xielab/analysisdata/ThreeDGene/chenjx/ecDNA_project/analysis/circle_probability/293T/293T_merge_sort/Regions_293T_merge.sort.enriched.merged.b2bRefined.Counts.ThreshFinal.bed',sep='\t',header=FALSE) # input the files containing circle regions called from merged pseudobulk sample
metadata <- read.csv('/gshare/xielab/chenjx/ecDNA_Project/ecDNA_metadata.csv',sep=',',header=TRUE,row.names=1) # input metadata table for filtering
cell_select <- rownames(metadata[which(metadata$Circle_Read_Percent... > 50 & metadata$Mapping_Ratio... > 80 & metadata$Cell_type == '293T' & metadata$Drug == 'No'),])
merge <- merge[which(merge$V6 %in% cell_select),]
total <- length(unique(merge$V6))

into_grange <- function(genomic_region_data_frame){
    grange <- makeGRangesFromDataFrame(genomic_region_data_frame,
                            keep.extra.columns=FALSE,
                            ignore.strand=TRUE,
                            seqnames.field='V1',
                            start.field='V2',
                            end.field='V3')
    grange
}

find_overlap <- function(single_region,multi_region){
    groupA <- data.table(
    chr = single_region$V1,
    start = single_region$V2,
    end = single_region$V3)
    setkey(groupA, chr, start, end)
    minoverlap = floor(0.1*(groupA$end - groupA$start))
    single_region_grange <- into_grange(single_region)

    groupB <- data.table(
        chr = multi_region$V1,
        start = multi_region$V2,
        end = multi_region$V3)
    setkey(groupB, chr, start, end)
    
    over <- foverlaps(groupA, groupB, which=TRUE)
    if (is.na(over$yid)) {
        print('No hit found')
        overlap <- NA
    } else {
        overlap <- multi_region[over$yid,]
        name <- unique(overlap$V6)
        for (i in c(1:length(name))){
            overlap_of_specific_cell <- overlap[which(overlap$V6 == name[i]),]
            grange_overlap_of_specific_cell <- into_grange(overlap_of_specific_cell)
            #calculate overlapping length
            hits <- findOverlaps(single_region_grange, grange_overlap_of_specific_cell)
            p <- Pairs(single_region_grange, grange_overlap_of_specific_cell, hits = hits)
            intersect <- pintersect(p)
            overlap_length <- sum(intersect@ranges@width)
            if (overlap_length < minoverlap){
                overlap <- overlap[which(overlap$V6 != name[i]),]
            }
        }
        if (nrow(overlap) == 0){
            overlap <- NA
        }
    }
    overlap
}

calculate_frequency <- function(overlap){
    if (is.na(overlap)){
        frequency = 0
    }else{
        frequency = length(unique(overlap$V6))/total 
    }
    frequency 
}

calculate_jaccard_similarity <- function(single_region,overlap){
    if (is.na(overlap)){
        jaccard_index = 0
    }else{ 
        single_region_grange <- into_grange(single_region)
        name <- unique(overlap$V6)
        if (length(name) > 1) {
            jaccard_index = 0
            for (i in (1: length(name))){
                overlap_of_specific_cell <- overlap[which(overlap$V6 == name[i]),]
                grange_overlap_of_specific_cell <- into_grange(overlap_of_specific_cell)
                #calculate overlapping length
                jaccard_index_current_grange <- genomicCorr.jaccard(single_region_grange, grange_overlap_of_specific_cell, restrict = NULL)
                jaccard_index = jaccard_index + jaccard_index_current_grange
                # print(jaccard_index_current_grange)
            }
        jaccard_index = jaccard_index/length(name)
        }else {
            overlap_grange <- into_grange(overlap)
            jaccard_index <- genomicCorr.jaccard(overlap_grange, single_region_grange, restrict = NULL)/total
        }
    }
    jaccard_index
}

main_function <- function(pseudobulk,merge){
    pseudobulk$frequency <- 0
    pseudobulk$jaccard_index <- 0
    for (i in c(1:nrow(pseudobulk))){
        print(i)
        overlap <- find_overlap(pseudobulk[i,],merge)
        pseudobulk$frequency[i] <- calculate_frequency(overlap)
        pseudobulk$jaccard_index[i] <- calculate_jaccard_similarity(pseudobulk[i,],overlap)
    }
    pseudobulk
}

pseudobulk_new <- main_function(pseudobulk,merge)
pseudobulk_new <- pseudobulk_new[order(pseudobulk_new$frequency),]
pseudobulk_new <- pseudobulk_new[order(pseudobulk_new$jaccard_index),]
