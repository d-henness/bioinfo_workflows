include: "fastp_rna.snakefile"
configfile: "{}/ref.yaml".format(workflow.basedir)


rule mixcr_all:
    input:
      expand("mixcr/{rna_lib}/mixcr/{rna_lib}_clones.txt", rna_lib = config["rna_libraries"])

alignment_threads = 4
rule mixcr_align:
  input:
    fq1 = rules.fastp_paired.output.fq1_out,
    fq2 = rules.fastp_paired.output.fq2_out,
  output:
    alignments = "mixcr/{rna_lib}/mixcr/alignments.vdjca",
  conda: "envs_dir/mixcr.yaml",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 7 * 1024,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  threads: 1 + alignment_threads
  benchmark: "mixcr/benchmark/{rna_lib}_mixcr_align.benchmark"
  log: "mixcr/log/{rna_lib}_mixcr_align.log",
  shell:
    """
      mixcr align -p rna-seq -s hsa -t {alignment_threads} -r {log} -OallowPartialAlignments=true {input.fq1} {input.fq2} {output.alignments}
    """

rule mixcr_assemble:
  input:
    alignments = rules.mixcr_align.output.alignments
  output:
    rescued_1 = "mixcr/{rna_lib}/mixcr/rescued_1.vdjca",
    rescued_2 = "mixcr/{rna_lib}/mixcr/rescued_2.vdjca",
    clones = "mixcr/{rna_lib}/mixcr/clones.clna",
  threads: 4
  log: "mixcr/log/{rna_lib}_mixcr_assemble.log",
  benchmark: "mixcr/benchmark/{rna_lib}_mixcr_assemble.benchmark"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 1 * 1024,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  conda: "envs_dir/mixcr.yaml"
  shell:
    """
      mixcr assemblePartial {input.alignments} {output.rescued_1}
      mixcr assemblePartial {output.rescued_1} {output.rescued_2}
      mixcr assemble -t {threads} -r {log} {output.rescued_2} {output.clones}
    """

rule mixcr_exportclones:
  input:
    clones = rules.mixcr_assemble.output.clones,
  output:
    tra_summary = "mixcr/{rna_lib}/mixcr/{rna_lib}_tra_clones.txt",
    trb_summary = "mixcr/{rna_lib}/mixcr/{rna_lib}_trb_clones.txt",
    trg_summary = "mixcr/{rna_lib}/mixcr/{rna_lib}_trg_clones.txt",
    all_summary = "mixcr/{rna_lib}/mixcr/{rna_lib}_clones.txt",
  conda: "envs_dir/mixcr.yaml"
  log: "mixcr/log/{rna_lib}_mixcr_exportclones.log",
  benchmark: "mixcr/benchmark/{rna_lib}_mixcr_exportclones.benchmark"
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 1 * 1024,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,  # time in minutes
  shell:
    """
      mixcr exportClones -c TRA {input.clones} {output.tra_summary}
      mixcr exportClones -c TRB {input.clones} {output.trb_summary}
      mixcr exportClones -c TRG {input.clones} {output.trg_summary}
      mixcr exportClones {input.clones} {output.all_summary}
    """
