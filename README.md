# TODO
* automatically download other git repos when needed

# Bioinformatics workflows
currently supported workflows:
* GATK preprocessing
* MuTect2
* ichorCNA
* TitanCNA
* Strelka
* PhyloWGS
* VarScan
* MuPeXI (kind of)

~~Note that only MuTect2 calculations can currently be run in a simple way. Many workflows will need to be manually editted inorder to work on a different machine/user. This will be addressed in the future.~~ Most workflows have been completely hooked together, but not all have good time and memory estimates.

## Example command to run MuTect2 calculations
```
snakemake -s path/to/mutect2.snakefile  \
          --configfile path/to/configfile.yaml  \
          --jobs xx  \
          --use-conda \
          --cluster 'sbatch --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb} --time={resources.time_min}' \
          -p \
          --restart-times 3 \
          --cluster-status /path/to/scripts_dir/slurm_cluster_status.py
```
* -s: gives the snakefile for the pipeline we want.
* --configfile: the path to the job specific yaml config file, see below for an example and formatting.
* --jobs: the maximum number of jobs that snakemake will try to submit at once. Check with sys admins to figure out this limit.
* --use-conda: allows snakemake to use conda enviroments to install the needed programs.
* --cluster: needed for running the calculations on a cluster, here I am assuming that the cluster is using the slurm scheduler.
* --cluster-status: needed so that snakemake can detect job timeouts. This will only work for the slurm scheduler.
* (optional) -p: print the shell command that runs a particular job, this is useful for debugging but is not actully needed.
* (optional) -n: print out all steps that will be run but won't actually run them, this is useful to make sure that the jobs to be run match what you expect.

One liner
```
snakemake -s path/to/mutect2.snakefile    --configfile path/to/configfile.yaml    --jobs xx    --use-conda   --cluster 'sbatch --cpus-per-task={threads} --mem-per-cpu={resources.mem_mb} --time={resources.time_min}'   -p   --restart-times 3   --cluster-status /path/to/scripts_dir/slurm_cluster_status.py
```

## Example configfile.yaml file
snakemake needs a config file to tell it where job specific data is located and also how the data is related.
Here is an example yaml file for a mutect run. There are three different libraries from one patient to be run:
one normal sample (normal_0) and two tumor samples (tumor_0 and tumor_1). The data for tumor_0 was broken up into
two different Illumina runs. Note that the input files must be in the fq.gz format.

```
libraries:
  normal_0_0: ['/path/to/normal_0_1.fq.gz', '/path/to/normal_0_2.fq.gz']
  tumor_0_0: ['/path/to/tumor_0_0_1.fq.gz', '/path/to/tumor_0_0_2.fq.gz']
  tumor_0_1: ['/path/to/tumor_0_1_1.fq.gz', '/path/to/tumor_0_1_2.fq.gz']
  tumor_1_0: ['/path/to/tumor_1_1.fq.gz', '/path/to/tumor_1_2.fq.gz']
merge_libs:
  normal_0: [normal_0_0]
  tumor_0:  [tumor_0_0, tumor_0_1]
  tumor_1:  [tumor_1_0]
pairs:
  tumor_0:  normal_0
  tumor_1:  normal_0
```

The libraries and merge_libs sections can be genreated automatically using the make_config_2.py script (in scripts_dir).
It takes one or more arguments which are a strings that matche something that all libraries we want to work on contain.
For instance we want to do runs on the fq.gz files in the directories called RG4-2, RG4-3, RG4_PBMC, RG5-1, RG5-2, and RG5_PBMC.
We would run the following command in a directory that contains all of these other directories:

```
python3 /path/to/bioinfo_workflows/scripts_dir/make_config_2.py RG > config.yaml
```

This creates a file called config.yaml which would contain something like the following.

```
libraries:
  RG4-2_0: ['path/to/RG4-2/RG4_2_USPD16082184-2_H3TGLCCXY_L3_1.fq.gz', 'path/to/RG4-2/RG4_2_USPD16082184-2_H3TGLCCXY_L3_2.fq.gz']
  RG4-3_0: ['path/to/RG4-3/RG4_3_USPD16082184-3_H3TGLCCXY_L3_1.fq.gz', 'path/to/RG4-3/RG4_3_USPD16082184-3_H3TGLCCXY_L3_2.fq.gz']
  RG4_PBMC_0: ['path/to/RG4_PBMC/RG4_PBMC_USPD16083698-A61_HGJJFBBXX_L8_1.fq.gz', 'path/to/RG4_PBMC/RG4_PBMC_USPD16083698-A61_HGJJFBBXX_L8_2.fq.gz']
  RG5-1_0: ['path/to/RG5-1/RG5_1_USPD16082185-4_H3TGLCCXY_L4_1.fq.gz', 'path/to/RG5-1/RG5_1_USPD16082185-4_H3TGLCCXY_L4_2.fq.gz']
  RG5-2_0: ['path/to/RG5-2/RG5_2_USPD16082185-5_H3TGLCCXY_L4_1.fq.gz', 'path/to/RG5-2/RG5_2_USPD16082185-5_H3TGLCCXY_L4_2.fq.gz']
  RG5_PBMC_0: ['path/to/RG5_PBMC/RG5_PBMC_USPD16087877_HMGG5BBXX_L5_1.fq.gz', 'path/to/RG5_PBMC/RG5_PBMC_USPD16087877_HMGG5BBXX_L5_2.fq.gz']
merge_libs:
  RG4-2: [RG4-2_0]
  RG4-3: [RG4-3_0]
  RG4_PBMC: [RG4_PBMC_0]
  RG5-1: [RG5-1_0]
  RG5-2: [RG5-2_0]
  RG5_PBMC: [RG5_PBMC_0]
```

The pairs section needs to be made manually.
