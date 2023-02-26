#!/usr/bin/bash
for i in {1..4};do
for j in {1..2};do
sbatch <<EOT
#!/usr/bin/bash
#SBATCH -J Cp${i}-${j}
#SBATCH -c 20
#SBATCH --mem 20G
#SBATCH --partition=compute_new
#SBATCH -o logs/Cp${i}-${j}.log
./bedtools_count.sh -b /gshare/xielab/chenjx/data/ecDNA/6th_k562_colobulk/Cp${i}-${j}/Cp${i}-${j}.dedup.bam -n Cp${i}-${j}
EOT
done
done

