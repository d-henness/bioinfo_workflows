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
          --use-conda \
          --cluster 'sbatch --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb} --exclude={params.exclude_list}' \
          -p
```
* -s: gives the snakefile for the pipeline we want
* --configfile: the path to the job specific yaml config file, see below for an example and formatting
* --cores: gives the maximum number of cores that any particular job can use, 16 is currently the max that any job uses so just set this to 16
* --use-conda: allows snakemake to use conda enviroments to install the needed programs
* --cluster: needed for running the calculations on a cluster, here I am assuming that the cluster is using the slurm scheduler
* (optional) -p: print the shell command that runs a particular job, this is useful for debugging but is not actully needed
* (optional) -n: print out all steps that will be run but won't actually run them, this is useful to make sure that the jobs to be run match what you expect

## Example configfile.yaml file
snakemake needs a config file to tell it where job specific data is located and also how the data is related.
Here is an example yaml file for a mutect run. There are three different libraries from one patient to be run: 
one normal sample (normal_0) and two tumor samples (tumor_0 and tumor_1)
```
libraries:
  normal_0_0: ['/path/to/normal_0_1.fg.gz', '/path/to/normal_0_2.fg.gz']
  tumor_0_0: ['/path/to/tumor_0_0_1.fg.gz', '/path/to/tumor_0_0_2.fg.gz']
  tumor_0_1: ['/path/to/tumor_0_1_1.fg.gz', '/path/to/tumor_0_1_2.fg.gz']
  tumor_1_0: ['/path/to/tumor_1_0_1.fg.gz', '/path/to/tumor_1_0_2.fg.gz']
merge_libs:
  normal_0: [normal_0_0]
  tumor_0:  [tumor_0_0, tumor_0_1]
  tumor_1:  [tumor_1_0]
pairs:
  tumor_0:  normal_0
  tumor_1:  normal_0
```
