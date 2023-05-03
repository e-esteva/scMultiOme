#!/bin/bash
#SBATCH --time=0-1
#SBATCH --mem=10GB
#SBATCH --out=scMultiOme-meta_%j.out
#SBATCH --error=scMultiOme-meta_%j.err


source configFile.txt

# run multiOme preProcessing
multiOme function {
	local job_id=($(sbatch --export=configFile=$configFile --mem=${module1_mem} --time=${module1_time}  --array=1-${batch_size} ${module1_Path}))
	echo ${job_id[3]}
}


module_1_job_id=$(multiOme)
