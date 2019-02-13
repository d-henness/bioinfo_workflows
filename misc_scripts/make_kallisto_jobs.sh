boiler_plate(){
job_string="#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00

echo \"start: \$(date)\"
set -e
module load kallisto
kallisto_index=~/hg38_data/cDNA/Homo_sapiens.GRCh38.cdna.all.kallisto.index
fq1=$1
fq2=$2
dir=$3
kallisto quant -b 500 --index=\$kallisto_index --output-dir=kallisto/\$dir --threads=4 \$fq1 \$fq2 &> kallisto/\$dir/\$dir\_kallisto.log
kallisto h5dump -o=kallisto/\$dir/h5dump kallisto/\$dir/abundance.h5
python3 ~/aux/get_mean_exp.py kallisto/\$dir/h5dump > kallisto/\$dir/\$dir\_mean_exp.tsv
echo \"end: \$(date)\""
echo "$job_string"
}

dir=$1
dir=$(echo "$dir" | sed 's/\/$//')
if [[ ( "$dir" == '' ) || !(-d $dir ) ]]; then
	echo "Directory not found"
	exit
fi
fq1=$(ls $dir/*1.fq.gz)
fq2=$(ls $dir/*2.fq.gz)
if [[ ( "$fq1" == '' ) || ( "$fq2" == '' ) ]]; then
	echo "Matching fq.gz files not found"
	exit
fi

if [[ !(-d "kallisto/$dir") ]]; then
	mkdir -p "kallisto/$dir"
fi

boiler_plate $fq1 $fq2 $dir > job_$dir.sh
