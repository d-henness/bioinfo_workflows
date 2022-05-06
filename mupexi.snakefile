configfile: "{}/ref.yaml".format(workflow.basedir)

include: "./mutect2_alt_bed.snakefile"
include: "kallisto.snakefile"

def make_normal_set(tumors):
  normals = []
  for tumor in tumors:
    if tumors[tumor] not in normals:
      normals.append(tumors[tumor])
  return normals

def kallisto_runs(tumor):
  try:
    return f"kallisto/{config['rna_pairs'][tumor]}/kallisto/{config['rna_pairs'][tumor]}_mean_exp.tsv"
  except:
    return f"mupexi_runs/{config['pairs'][tumor]}/optitype/{config['pairs'][tumor]}_result.tsv"


#noramls = make_normal_set(config['pairs'])
rule run_mupexi:
  input:
#   expand("mupexi_runs/{lib}/optitype/{lib}_result.tsv", lib = noramls),
    expand("mupexi_runs/{tumor}/mupexi/{tumor}_merged.mupexi", tumor = config['pairs']),

rule razers3_1:
  input:
    fq1_in = lambda wildcards: config["libraries"][config["merge_libs"][wildcards.lib][0]][0],
  output:
    fq1_out = temp("mupexi_runs/{lib}/razers3/{lib}_1_filtered.fq"),
    bam1_out = temp("mupexi_runs/{lib}/razers3/{lib}_1_filtered.bam"),
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    razers3_log_1 = "mupexi_runs/log/{lib}_razers3_1.log",
  benchmark: "mupexi_runs/benchmark/{lib}_razers3_1.benchmark"
  threads: 2
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (16 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  shell:
    """
      gzip -dc {input.fq1_in} | parallel -j {threads} --pipe -L4 -N1000000 "cat > mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_1_{{#}}.fastq; razers3 -i 95 -m 1 -dr 0 -o mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_1_{{#}}.bam $CONDA_PREFIX/bin/data/hla_reference_dna.fasta mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_1_{{#}}.fastq; rm mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_1_{{#}}.fastq"
      ls mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_1_*.bam > mupexi_runs/{wildcards.lib}/razers3/all_bam.txt
      samtools cat -b mupexi_runs/{wildcards.lib}/razers3/all_bam.txt > {output.bam1_out}
      samtools bam2fq {output.bam1_out} > {output.fq1_out}
    """

rule razers3_2:
  input:
    fq2_in = lambda wildcards: config["libraries"][config["merge_libs"][wildcards.lib][0]][1],
  output:
    fq2_out = temp("mupexi_runs/{lib}/razers3/{lib}_2_filtered.fq"),
    bam2_out = temp("mupexi_runs/{lib}/razers3/{lib}_2_filtered.bam"),
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    razers3_log_2 = "mupexi_runs/log/{lib}_razers3_2.log",
  benchmark: "mupexi_runs/benchmark/{lib}_razers3_2.benchmark"
  threads: 2
  resources:
    mem_mb = lambda wildcards, attempt: attempt * (16 * 1024),
    time_min = lambda wildcards, attempt: attempt * (24 * 60),	# time in minutes
  shell:
    """
      gzip -dc {input.fq2_in} | parallel -j {threads} --pipe -L4 -N1000000 "cat > mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_2_{{#}}.fastq; razers3 -i 95 -m 1 -dr 0 -o mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_2_{{#}}.bam $CONDA_PREFIX/bin/data/hla_reference_dna.fasta mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_2_{{#}}.fastq; rm mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_2_{{#}}.fastq"
      ls mupexi_runs/{wildcards.lib}/razers3/{wildcards.lib}_2_*.bam > mupexi_runs/{wildcards.lib}/razers3/all_bam.txt
      samtools cat -b mupexi_runs/{wildcards.lib}/razers3/all_bam.txt > {output.bam2_out}
      samtools bam2fq {output.bam2_out} > {output.fq2_out}
    """

rule optitype:
  input:
    fq1_in = rules.razers3_1.output.fq1_out,
    fq2_in = rules.razers3_2.output.fq2_out,
  output:
    hlas = "mupexi_runs/{lib}/optitype/{lib}_result.tsv",
  conda:
    "./envs_dir/optitype.yaml"
  log:
    optitype_log = "mupexi_runs/log/{lib}_optitype.log",
  benchmark: "mupexi_runs/benchmark/{lib}_optitype.benchmark"
  params:
    out_dir = "mupexi_runs/{lib}/optitype",
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 16000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  shell:
    """
      if [[ !(-d mupexi_runs/{wildcards.lib}/optitype) ]]; then
        mkdir -p mupexi_runs/{wildcards.lib}/optitype
      fi
      cp $CONDA_PREFIX/bin/config.ini mupexi_runs/{wildcards.lib}/optitype/config.ini
      sed -i 's/razers3.*/razers3=razers3/' mupexi_runs/{wildcards.lib}/optitype/config.ini
      sed -i 's/threads.*/threads={threads}/' mupexi_runs/{wildcards.lib}/optitype/config.ini

      OptiTypePipeline.py -i {input.fq1_in} {input.fq2_in} --dna -v -o {params.out_dir} -p {wildcards.lib} -c mupexi_runs/{wildcards.lib}/optitype/config.ini 2> {log.optitype_log}
    """

rule split_vcf_files:
  input:
    vcf = rules.FilterByOrientationBias.output.vcf,
  output:
    chr1 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr1.vcf",
    chr2 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr2.vcf",
    chr3 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr3.vcf",
    chr4 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr4.vcf",
    chr5 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr5.vcf",
    chr6 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr6.vcf",
    chr7 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr7.vcf",
    chr8 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr8.vcf",
    chr9 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr9.vcf",
    chr10 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr10.vcf",
    chr11 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr11.vcf",
    chr12 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr12.vcf",
    chr13 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr13.vcf",
    chr14 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr14.vcf",
    chr15 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr15.vcf",
    chr16 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr16.vcf",
    chr17 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr17.vcf",
    chr18 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr18.vcf",
    chr19 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr19.vcf",
    chr20 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr20.vcf",
    chr21 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr21.vcf",
    chr22 = "mupexi_runs/{tumor}/mupexi/{tumor}_chr22.vcf",
    chrX = "mupexi_runs/{tumor}/mupexi/{tumor}_chrX.vcf",
    chrY = "mupexi_runs/{tumor}/mupexi/{tumor}_chrY.vcf",
  log:
    razers3_log_2 = "mupexi_runs/log/{tumor}_split_vcf_files.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_split_vcf_files.benchmark"
  params:
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
    out_dir = "mupexi_runs/{tumor}/mupexi",
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 1000,
    time_min = lambda wildcards, attempt: attempt * 10,	# time in minutes
  shell:
    """
      python3 {params.bioinfo_workflows_path}/scripts_dir/split_up_vcf_by_chr.py {input.vcf} {params.out_dir} {wildcards.tumor}
    """

rule mupexi_chr1:
  input:
    vcf = rules.split_vcf_files.output.chr1,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr1/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr1.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr1.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr1/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 72 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 24,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr2:
  input:
    vcf = rules.split_vcf_files.output.chr2,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr2/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr2.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr2.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr2/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 16,
  shell:
    """
      echo "--------------here---------------------"
      which python
      /usr/bin/env python --version
      echo "--------------here---------------------"
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr3:
  input:
    vcf = rules.split_vcf_files.output.chr3,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr3/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr3.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr3.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr3/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 48 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 16,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr4:
  input:
    vcf = rules.split_vcf_files.output.chr4,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr4/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr4.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr4.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr4/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 12 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr5:
  input:
    vcf = rules.split_vcf_files.output.chr5,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr5/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr5.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr5.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr5/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 12,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr6:
  input:
    vcf = rules.split_vcf_files.output.chr6,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr6/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr6.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr6.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr6/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 48 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 12,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr7:
  input:
    vcf = rules.split_vcf_files.output.chr7,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr7/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr7.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr7.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr7/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 16,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr8:
  input:
    vcf = rules.split_vcf_files.output.chr8,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr8/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr8.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr8.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr8/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr9:
  input:
    vcf = rules.split_vcf_files.output.chr9,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr9/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr9.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr9.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr9/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 12 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr10:
  input:
    vcf = rules.split_vcf_files.output.chr10,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr10/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr10.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr10.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr10/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 12,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr11:
  input:
    vcf = rules.split_vcf_files.output.chr11,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr11/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr11.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr11.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr11/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 16,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr12:
  input:
    vcf = rules.split_vcf_files.output.chr12,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr12/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr12.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr12.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr12/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 12,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr13:
  input:
    vcf = rules.split_vcf_files.output.chr13,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr13/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr13.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr13.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr13/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr14:
  input:
    vcf = rules.split_vcf_files.output.chr14,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr14/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr14.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr14.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr14/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr15:
  input:
    vcf = rules.split_vcf_files.output.chr15,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr15/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr15.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr15.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr15/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr16:
  input:
    vcf = rules.split_vcf_files.output.chr16,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr16/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr16.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr16.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr16/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 16,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr17:
  input:
    vcf = rules.split_vcf_files.output.chr17,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr17/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr17.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr17.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr17/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 20,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr18:
  input:
    vcf = rules.split_vcf_files.output.chr18,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr18/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr18.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr18.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr18/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 20,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr19:
  input:
    vcf = rules.split_vcf_files.output.chr19,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr19/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr19.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr19.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr19/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 36 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 20,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr20:
  input:
    vcf = rules.split_vcf_files.output.chr20,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr20/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr20.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr20.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr20/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 24 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr21:
  input:
    vcf = rules.split_vcf_files.output.chr21,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr21/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr21.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr21.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr21/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 12 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chr22:
  input:
    vcf = rules.split_vcf_files.output.chr22,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chr22/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chr22.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chr22.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chr22/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 12 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chrX:
  input:
    vcf = rules.split_vcf_files.output.chrX,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chrX/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chrX.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chrX.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chrX/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 12 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule mupexi_chrY:
  input:
    vcf = rules.split_vcf_files.output.chrY,
    hlas = lambda wildcards: f"mupexi_runs/{config['pairs'][wildcards.tumor]}/optitype/{config['pairs'][wildcards.tumor]}_result.tsv",
    rna = lambda wildcards: kallisto_runs(wildcards.tumor),
  output:
    mupexi_out = "mupexi_runs/{tumor}/mupexi/chrY/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_mupexi_chrY.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi_chrY.benchmark"
  params:
    mupexi_path = f"{config['bioinfo_workflows_path']}/MuPeXI/MuPeXI.py",
    mupexi_config = f"{config['bioinfo_workflows_path']}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}/mupexi/chrY/",
  threads: 1
  resources:
    time_min = lambda wildcards, attempt: attempt * 12 * 60,
    mem_mb = lambda wildcards, attempt: attempt * 1024 * 8,
  shell:
    """
      hla_string=$(sed -En '/[A-C]\*/{{s/(\\t([A-C])\*)/\\tHLA-\\2/g;p}}' {input.hlas} | sed -E 's/\s/\,/g; s/^[^,]*,//; s/,[^,]*,[^,]*$//')

      if [[ "$hla_string" == '' ]]; then
        echo Problem with hla file
        exit 1
      fi

      rna_var=''
      if [[ "{input.rna}" == "kallisto"* ]]; then
        rna_var="-e {input.rna} "
      fi

      {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} $rna_var -t -d {params.out_dir} -p {wildcards.tumor} -l 8-11 -f &> {log.mupexi_log}
    """

rule merge_mupexi:
  input:
    chr1 = rules.mupexi_chr1.output.mupexi_out,
    chr2 = rules.mupexi_chr2.output.mupexi_out,
    chr3 = rules.mupexi_chr3.output.mupexi_out,
    chr4 = rules.mupexi_chr4.output.mupexi_out,
    chr5 = rules.mupexi_chr5.output.mupexi_out,
    chr6 = rules.mupexi_chr6.output.mupexi_out,
    chr7 = rules.mupexi_chr7.output.mupexi_out,
    chr8 = rules.mupexi_chr8.output.mupexi_out,
    chr9 = rules.mupexi_chr9.output.mupexi_out,
    chr10 = rules.mupexi_chr10.output.mupexi_out,
    chr11 = rules.mupexi_chr11.output.mupexi_out,
    chr12 = rules.mupexi_chr12.output.mupexi_out,
    chr13 = rules.mupexi_chr13.output.mupexi_out,
    chr14 = rules.mupexi_chr14.output.mupexi_out,
    chr15 = rules.mupexi_chr15.output.mupexi_out,
    chr16 = rules.mupexi_chr16.output.mupexi_out,
    chr17 = rules.mupexi_chr17.output.mupexi_out,
    chr18 = rules.mupexi_chr18.output.mupexi_out,
    chr19 = rules.mupexi_chr19.output.mupexi_out,
    chr20 = rules.mupexi_chr20.output.mupexi_out,
    chr21 = rules.mupexi_chr21.output.mupexi_out,
    chr22 = rules.mupexi_chr22.output.mupexi_out,
    chrX = rules.mupexi_chrX.output.mupexi_out,
    chrY = rules.mupexi_chrY.output.mupexi_out,
  output:
    merged_mupexi = "mupexi_runs/{tumor}/mupexi/{tumor}_merged.mupexi"
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}_merge_mupexi.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_merge_mupexi.benchmark"
  params:
    dirnm = "mupexi_runs/{tumor}/mupexi/",
    workflow_path = config['bioinfo_workflows_path'],
  threads: 1
  resources:
    mem_mb = 4000,
  shell:
    """
      python3 {params.workflow_path}/scripts_dir/merge_mupexi_files.py {params.dirnm} {wildcards.tumor} &> {log.mupexi_log}
    """
