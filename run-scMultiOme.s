#!/bin/bash
#SBATCH --out=scMultiOme_%j.out
#SBATCH --error=scMultiOme_%j.err



source configFile.txt

module load r/4.1.2

cellranger_out=$(head -n${SLURM_ARRAY_TASK_ID} ${inputs} | tail -n1)

echo ${cellranger_out}

out_dir=$(head -n${SLURM_ARRAY_TASK_ID} ${outs} | tail -n1)
mkdir -p ${out_dir}

echo ${out_dir}

Rscript scATAC-WNN.R ${cellranger_out} ${out_dir} ${species} 
