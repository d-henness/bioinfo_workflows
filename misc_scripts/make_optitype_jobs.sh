boiler_plate(){
job_string="#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=1:00:00
set -e
module load miniconda3
source activate optitype

fq1=\"$1\"
fq2=\"$2\"
dir=\"$3\"
echo \"Running Optitype on \$fq1\"

gzip -dc \$fq1 | parallel -j 16 --pipe -L4 -N10000000 \"cat > \$fq1\_{#}.fastq; razers3 -i 95 -m 1 -dr 0 -o \$fq1\_{#}.bam \$CONDA_PREFIX/share/optitype-1.3.2-1/data/hla_reference_dna.fasta \$fq1\_{#}.fastq; rm \$fq1\_{#}.fastq\"
gzip -dc \$fq2 | parallel -j 16 --pipe -L4 -N10000000 \"cat > \$fq2\_{#}.fastq; razers3 -i 95 -m 1 -dr 0 -o \$fq2\_{#}.bam \$CONDA_PREFIX/share/optitype-1.3.2-1/data/hla_reference_dna.fasta \$fq2\_{#}.fastq; rm \$fq2\_{#}.fastq\"
ls \$fq1*bam > \$fq1\_concat.txt
ls \$fq2*bam > \$fq2\_concat.txt
samtools cat -b \$fq1\_concat.txt -o \$fq1\_concat.bam
samtools cat -b \$fq2\_concat.txt -o \$fq2\_concat.bam
samtools bam2fq \$fq1\_concat.bam > \$fq1\_filtered.fastq
samtools bam2fq \$fq2\_concat.bam > \$fq2\_filtered.fastq
rm \$fq1*bam
rm \$fq2*bam
python \$CONDA_PREFIX/bin/OptiTypePipeline.py -i \$fq1\_filtered.fastq \$fq2\_filtered.fastq --dna -v -o \$dir -p \$(echo \$dir | sed 's/\///') -c \$CONDA_PREFIX/share/optitype-1.3.2-1/config.ini
sed -rn '/[A-C]\*/{s/(\t([A-C])\*)/\tHLA-\2/g;p}' \$dir/*tsv > \$dir/\$dir.hlatype"
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

boiler_plate $fq1 $fq2 $dir > job_$dir.sh
