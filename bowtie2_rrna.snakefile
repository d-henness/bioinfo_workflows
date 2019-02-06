rna_libs = config["libraries"]

rule all:
  input:
    expand("bowtie2_rrna/{rna_lib}/{rna_lib}.sam", rna_lib = rna_libs)

rule bowtie2:
  input:
    fq1 = lambda wildcards: rna_libs[wildcards.rna_lib][0],
    fq2 = lambda wildcards: rna_libs[wildcards.rna_lib][1]
  output:
    sam_file = "bowtie2_rrna/{rna_lib}/{rna_lib}.sam"
  conda:
    "/usr/local/bioinfo_workflows/envs_dir/bowtie2.yaml"
  params:
    regions = "/data/shared/hg38/ftp.ensembl.org/pub/release-95/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.rRNA",
  resources:
    mem_mb = 4000
  threads: 1
  log:
    "bowtie2_rrna/{rna_lib}/{rna_lib}.log"
  shell:
    """
      bowtie2 \
        --no-unal \
        -x {params.regions} \
        -1 {input.fq1} \
        -2 {input.fq2} \
        -S {output.sam_file} &> {log}
    """
