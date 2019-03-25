configfile: "{}/ref.yaml".format(workflow.basedir)

rule run_varscan:
  input:
    expand("varscan_runs/{tumor}/varscan/{tumor}_snp.vcf", zip, sample = config["merge_libs"], tumor = config["pairs"])

rule make_pileup:
  input:
    bam_in = "runs/{sample}/ApplyBQSR/recal.bam",
  output:
    pileup_out = "varscan_runs/{sample}/mpileup/{sample}.pileup",
  conda:
    "envs_dir/varscan.yaml"
  log:
    "varscan_runs/log/{sample}_pileup.log",
  benchmark:
    "varscan_runs/benchmark/{sample}_pileup.benchmark",
  params:
    ref_fasta = config["ref_fasta"],
  threads: 1
  resources:
    mem_mb = 4000,
  shell:
    """
      samtools mpileup -q 1 -f {params.ref_fasta} {input.bam_in} > {output.pileup_out} 2> {log}
    """

rule varscan:
  input:
    tumor = "varscan_runs/{tumor}/mpileup/{tumor}.pileup",
    normal = lambda wildcards: "varscan_runs/" + config["pairs"][wildcards.tumor] + "/mpileup/" + config["pairs"][wildcards.tumor] + ".pileup",
  output:
    snp = "varscan_runs/{tumor}/varscan/{tumor}_snp.vcf",
    indel = "varscan_runs/{tumor}/varscan/{tumor}_indel.vcf",
  conda:
    "envs_dir/varscan.yaml"
  log:
    "varscan_runs/log/{tumor}_varscan.log",
  benchmark:
    "varscan_runs/benchmark/{tumor}_varscan.benchmark",
  threads: 1
  resources:
    mem_mb = 4000,
  shell:
    """
      varscan somatic \
        {input.normal} \
        {input.tumor} \
        --output-vcf 1 \
        --output-snp {output.snp} \
        --output-indel {output.indel} 2> {log}
    """
