#!/usr/bin/bash
Rscript peak2gtf.r human_GRCh38_windows.bed
sed -i '1d' name.matrix
sed -i '1d' for_count.bed
