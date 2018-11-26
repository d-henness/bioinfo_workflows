# Bioinformatics workflows
currenty supported workflows:
* GATK preprocessing
* MuTect2
* ichorCNA
* TitanCNA
* Strelka
## Example command to run MuTect2 calculations
```
snakemake -s path/to/mutect2.snakefile  \
          --configfile path/to/configfile.yaml  \
          --cores 16  \
          --cluster 'sbatch --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb} --exclude={params.exclude_list}' \
          -p
```
