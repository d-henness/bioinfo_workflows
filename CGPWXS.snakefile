configfile: "{}/ref.yaml".format(workflow.basedir)

include: "pre_pro_af_merge_alt_bed.snakefile"

tumors = config["pairs"]
normals = [config["pairs"][tumor] for tumor in tumors]

rule CGPWXS_all:
  input:
    expand("CGPWXS/{tumor}/results/results/run.params", tumor = tumors),

rule make_ref_dir:
  output:
    fai = "CGPWXS/ref/Homo_sapiens_assembly38.fasta.fai",
    fasta = "CGPWXS/ref/Homo_sapiens_assembly38.fasta"
  params:
    ref_fai = config["ref_fasta_index"],
    ref_fasta = config["ref_fasta"],
    ref_dir = "CGPWXS/ref",
  log:
    "CGPWXS/logs/make_ref_dir.log",
  benchmark:
    "CGPWXS/benchmarks/make_ref_dir.benchmark"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
  threads: 1
  shell:
    """
      if [[ !(-d {params.ref_dir}) ]]; then
        mkdir -p {params.ref_dir}
      fi

      cp {params.ref_fai} {output.fai}
      cp {params.ref_fasta} {output.fasta}
    """

rule make_bas:
  input:
    fai = rules.make_ref_dir.output.fai,
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
  output:
    tumor_bas = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bas",
  params:
    ref_fai = config["ref_fasta_index"],
    ref_dir = "CGPWXS/ref",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
  log:
    "CGPWXS/logs/{tumor}_make_bas.log",
  conda:
    "envs_dir/CGPWXS.yaml"
  benchmark:
    "CGPWXS/benchmarks/{tumor}_make_bas.benchmark"
  shell:
    """
      export CGPWXS_VER=3.1.6

      docker run -i --rm \
      --read-only --tmpfs /tmp \
      --env HOME=/var/spool/results \
      -v $PWD/{params.ref_dir}:/var/spool/ref:ro \
      -v $PWD/GATK_runs/{wildcards.tumor}/ApplyBQSR:/var/spool/results:rw \
      quay.io/wtsicgp/dockstore-cgpwxs:$CGPWXS_VER \
      bam_stats \
        -i /var/spool/results/$(basename {input.tumor_bam}) \
        -o /var/spool/results/$(basename {output.tumor_bas}) \
        -r /var/spool/ref/$(basename {params.ref_fai}) \
        -@ 1

    """

rule caveman:
  input:
    fai = rules.make_ref_dir.output.fai,
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
    tumor_bas = rules.make_bas.output.tumor_bas,
    normal_bam = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam",
    normal_bai = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bai",
    normal_bas = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bas",
  output:
    "CGPWXS/{tumor}/results/results/run.params"
  params:
    ref_fai = config["ref_fasta_index"],
    ref_dir = "CGPWXS/ref",
    out_dir = "$PWD/CGPWXS/{tumor}/results",
  threads: 4
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 8000,
  log:
    "CGPWXS/logs/{tumor}_caveman.log",
  conda:
    "envs_dir/CGPWXS.yaml"
  benchmark:
    "CGPWXS/benchmarks/{tumor}_caveman.benchmark"
  shell:
    """
      if [[ !(-d {params.out_dir}/data) ]]; then
        mkdir -p {params.out_dir}/data
      fi

      if [[ !(-d {params.out_dir}/results) ]]; then
        mkdir -p {params.out_dir}/results
      fi

      cp {input.tumor_bam} {params.out_dir}/data
      cp {input.tumor_bai} {params.out_dir}/data
      cp {input.tumor_bas} {params.out_dir}/data
      cp {input.normal_bam} {params.out_dir}/data
      cp {input.normal_bai} {params.out_dir}/data
      cp {input.normal_bas} {params.out_dir}/data


      export CGPWXS_VER=3.1.6

      docker run -i --rm \
        --read-only --tmpfs /tmp \
        --env HOME=/var/spool/results \
        -v $PWD/{params.ref_dir}:/var/spool/ref:ro \
        -v {params.out_dir}/data:/var/spool/data:ro \
        -v {params.out_dir}/results:/var/spool/results:rw \
        quay.io/wtsicgp/dockstore-cgpwxs:$CGPWXS_VER \
        ds-cgpwxs.pl \
        -sp "HOMO SAPIENS" \
        -r /var/spool/ref/core_ref_GRCh38_hla_decoy_ebv.tar.gz \
        -a /var/spool/ref/VAGrENT_ref_GRCh38_hla_decoy_ebv_ensembl_91.tar.gz \
        -si /var/spool/ref/SNV_INDEL_ref_GRCh38_hla_decoy_ebv-fragment.tar.gz \
        -e 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV' \
        -t /var/spool/data/$(basename {input.tumor_bam}) \
        -tidx /var/spool/data/$(basename {input.tumor_bai}) \
        -n /var/spool/data/$(basename {input.normal_bam}) \
        -nidx /var/spool/data/$(basename {input.normal_bai}) \
        -c 4 \
        -o /var/spool/results &> {log}

    """
