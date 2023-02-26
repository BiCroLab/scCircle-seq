#!/bin/bash
mkdir logs
for i in {1..12}; do
sbatch <<EOT
#!/bin/bash
#SBATCH -J P${i}e
#SBATCH -c 10
#SBATCH --mem 20G
#SBATCH -p compute_new
#SBATCH -o logs/P${i}e.log
/gshare/xielab/chenjx/ecDNA_Project/scripts/circle_map.sh -f /gshare/xielab/chenjx/ecDNA_Project/data/2rd_batch/raw_data/ecDNA/P${i}e*_1.fq.gz -r /gshare/xielab/chenjx/ecDNA_Project/data/2rd_batch/raw_data/ecDNA/P${i}e*_2.fq.gz -o P${i}e -t 10 -p P${i}e -A pc3 -g hg19 -e 0 -R FALSE
EOT
done

for i in {1..24}; do
sbatch <<EOT
#!/bin/bash
#SBATCH -J H${i}e
#SBATCH -c 10
#SBATCH --mem 20G
#SBATCH -p compute_new
#SBATCH -o logs/H${i}e.log
/gshare/xielab/chenjx/ecDNA_Project/scripts/circle_map.sh -f /gshare/xielab/chenjx/ecDNA_Project/data/2rd_batch/raw_data/ecDNA/H${i}e*_1.fq.gz -r /gshare/xielab/chenjx/ecDNA_Project/data/2rd_batch/raw_data/ecDNA/H${i}e*_2.fq.gz -o H${i}e -t 10 -p H${i}e -A hela -g hg19 -e 0 -R FALSE
EOT
done

for i in {1..11}; do
sbatch <<EOT
#!/bin/bash
#SBATCH -J C${i}e
#SBATCH -c 10
#SBATCH --mem 20G
#SBATCH -p compute_new
#SBATCH -o logs/C${i}e.log
/gshare/xielab/chenjx/ecDNA_Project/scripts/circle_map.sh -f /gshare/xielab/chenjx/ecDNA_Project/data/2rd_batch/raw_data/ecDNA/C${i}e*_1.fq.gz -r /gshare/xielab/chenjx/ecDNA_Project/data/2rd_batch/raw_data/ecDNA/C${i}e*_2.fq.gz -o C${i}e -t 10 -p C${i}e -A colo -g hg19 -e 0 -R FALSE
EOT
done


