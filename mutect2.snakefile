configfile: "{}/ref.yaml".format(workflow.basedir)

include: "pre_pro_af_merge.snakefile"

tumors = config["pairs"]
normals = [config["pairs"][tumor] for tumor in tumors]

rule mutect_all:
  input:
#   expand("GATK_runs/{tumor}/M2_{normal}/out.vcf", zip, tumor = tumors, normal = normals)
#   expand("GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_detail_metrics", tumor = tumors)
#   expand("GATK_runs/{tumor}/CalculateContamination_{normal}/seg_tab.table", tumor = tumors, normal = normals)
#   expand("GATK_runs/{tumor}/MergeVCFs_{normal}/out.vcf", tumor = tumors, normal = normals)
#   expand("GATK_runs/{tumor}/Filter_{normal}/out.vcf", tumor = tumors, normal = normals)
    expand("GATK_runs/{tumor}/FilterByOrientationBias/{tumor}.vcf", tumor = tumors)
#   expand("GATK_runs/{tumor}/FuncotateMaf_{normal}/out.vcf.maf.annotated", zip, tumor = tumors, normal = normals) There is an issue with funcotator, uncomment this if a fix is found


rule CollectF1R2Counts:
  input:
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
  output:
    alt_table = temp("GATK_runs/{tumor}/CollectF1R2Counts/alt_tab.tsv"),
    ref_hist = temp("GATK_runs/{tumor}/CollectF1R2Counts/alt_tab.metrics"),
    ref_hists = temp("GATK_runs/{tumor}/CollectF1R2Counts/alt-depth1.metrics"),
    tumor_name = temp("GATK_runs/{tumor}/CollectF1R2Counts/tumor_name.txt"),
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/CollectF1R2Counts/out.log",
  conda:
    "envs_dir/pre_proc.yaml"
  benchmark:
    "benchmarks/{tumor}.CollectF1R2Counts.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m GetSampleName -R {params.ref_fasta} -I {input.tumor_bam} -O {output.tumor_name} -encode
      tumor_name=$(head -n 1 {output.tumor_name})

			gatk --java-options -Xmx4000m CollectF1R2Counts \
      -I {input.tumor_bam} -R {params.ref_fasta} \
      -L {params.interval} \
      -alt-table "{output.alt_table}" \
      -ref-hist "{output.ref_hist}" \
			-alt-hist "{output.ref_hists}" &> {log}
    """

rule LearnReadOrientationModel:
  conda:
    "envs_dir/pre_proc.yaml"
  input:
    alt_table = rules.CollectF1R2Counts.output.alt_table,
    ref_hist = rules.CollectF1R2Counts.output.ref_hist,
    ref_hists = rules.CollectF1R2Counts.output.ref_hists,
  output:
    temp("GATK_runs/{tumor}/LearnReadOrientationModel/art_tab.tsv"),
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/LearnReadOrientationModel/log.log",
  benchmark:
    "benchmarks/{tumor}.LearnReadOrientationModel.benchmark.txt"
  shell:
    """
			gatk --java-options -Xmx4000m LearnReadOrientationModel \
			-alt-table {input.alt_table} \
			-ref-hist {input.ref_hist} \
			-alt-hist {input.ref_hists} \
			-O "{output}" &> {log}
    """

rule M2:
  conda:
    "envs_dir/pre_proc.yaml"
  input:
    normal_bam = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam",
    normal_bai = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bai",
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
    art_tab = rules.LearnReadOrientationModel.output
  output:
    vcf = temp("GATK_runs/{tumor}/M2/out.vcf"),
    vcf_index = temp("GATK_runs/{tumor}/M2/out.vcf.idx"),
    bam_out = temp("GATK_runs/{tumor}/M2/out.bam"),
    tumor_name = temp("GATK_runs/{tumor}/M2/tumor_name.txt"),
    normal_name = temp("GATK_runs/{tumor}/M2/normal_name.txt"),
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    exclude_list = 'neuron,biolx95'
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60 * 30,	# time in minutes
  threads: 4 # default for mutect asks for 4 threads
  log:
    "GATK_runs/{tumor}/M2/out.log"
  benchmark:
    "benchmarks/{tumor}.M2.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m GetSampleName -R {params.ref_fasta} -I {input.tumor_bam} -O {output.tumor_name} -encode
      tumor_command_line="-I {input.tumor_bam} -tumor `cat {output.tumor_name}`"

      gatk --java-options -Xmx4000m GetSampleName -R {params.ref_fasta} -I {input.normal_bam} -O {output.normal_name} -encode 
      normal_command_line="-I {input.normal_bam} -normal `cat {output.normal_name}`"

      gatk --java-options -Xmx4000m Mutect2 \
      -R {params.ref_fasta} \
      ${{tumor_command_line}} \
      ${{normal_command_line}} \
      --germline-resource {params.gnomad} \
      -L {params.interval} \
      -O "{output.vcf}" \
      --bam-output {output.bam_out} \
      --orientation-bias-artifact-priors {input.art_tab} &> {log}
    """

rule CollectSequencingArtifactMetrics:
  input:
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
  output:
    bait_det_met = "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.bait_bias_detail_metrics",
    bait_sum_met = "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.bait_bias_summary_metrics",
    err_sum_met = "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.error_summary_metrics",
    ada_det_met = "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_detail_metrics",
    ada_sum_met = "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_summary_metrics",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.log"
  conda:
    "envs_dir/pre_proc.yaml"
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}.CollectSequencingArtifactMetrics.benchmark.txt"
  shell:
      """
        gatk --java-options -Xmx4000m CollectSequencingArtifactMetrics \
        -I {input.tumor_bam} -O GATK_runs/{wildcards.tumor}/CollectSequencingArtifactMetrics/gatk -R {params.ref_fasta} -VALIDATION_STRINGENCY LENIENT &> {log}
      """

rule CalculateContamination:
  input:
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
    normal_bam = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam",
    normal_bai = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bai",
  output:
    normal_pil_tab = temp("GATK_runs/{tumor}/CalculateContamination/normal_pileups.table"),
    pil_tab = temp("GATK_runs/{tumor}/CalculateContamination/pileups.table"),
    con_tab = temp("GATK_runs/{tumor}/CalculateContamination/con_tab.table"),
    seg_tab = temp("GATK_runs/{tumor}/CalculateContamination/seg_tab.table"),
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/CalculateContamination/out.log"
  conda:
    "envs_dir/pre_proc.yaml"
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}.CalculateContamination.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m GetPileupSummaries -I {input.normal_bam} --interval-set-rule INTERSECTION -L {params.interval} \
        -V {params.variants_con} -L {params.variants_con} -O {output.normal_pil_tab} &> {log}
      NORMAL_CMD="-matched {output.normal_pil_tab}"

      gatk --java-options -Xmx4000m GetPileupSummaries -R {params.ref_fasta} -I {input.tumor_bam} --interval-set-rule INTERSECTION -L {params.interval} \
        -V {params.variants_con} -L {params.variants_con} -O {output.pil_tab} &>> {log}
      gatk CalculateContamination -I {output.pil_tab} -O {output.con_tab} --tumor-segmentation {output.seg_tab} ${{NORMAL_CMD}} &>> {log}
    """

rule MergeVCFs:
  input:
    vcfs = rules.M2.output.vcf,
    vcf_index = rules.M2.output.vcf_index,
  output:
    merged_vcf = temp("GATK_runs/{tumor}/MergeVCFs/out.vcf"),
    merged_vcf_index = temp("GATK_runs/{tumor}/MergeVCFs/out.vcf.idx"),
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/MergeVCFs/out.log",
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}.MergeVCFs.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m MergeVcfs -I {input.vcfs} -O {output.merged_vcf} &> {log}
    """

rule Filter:
  input:
    vcf = rules.MergeVCFs.output.merged_vcf,
    vcf_index = rules.MergeVCFs.output.merged_vcf_index,
    con_tab = rules.CalculateContamination.output.con_tab,
    maf_seg = rules.CalculateContamination.output.seg_tab,
  output:
    filt_vcf = temp("GATK_runs/{tumor}/Filter/out.vcf"),
    filt_vcf_index = temp("GATK_runs/{tumor}/Filter/out.vcf.idx"),
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/Filter/out.log"
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}.Filter.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m FilterMutectCalls -V {input.vcf} \
        -O {output.filt_vcf} \
      	--contamination-table {input.con_tab} \
      	--tumor-segmentation {input.maf_seg} &> {log} \
    """

rule FilterByOrientationBias:
  input:
    vcf = rules.Filter.output.filt_vcf,
    vcf_index = rules.Filter.output.filt_vcf_index,
    ada_det_met = rules.CollectSequencingArtifactMetrics.output.ada_det_met
  output:
    vcf = "GATK_runs/{tumor}/FilterByOrientationBias/{tumor}.vcf",
    vcf_index = "GATK_runs/{tumor}/FilterByOrientationBias/{tumor}.vcf.idx",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/FilterByOrientationBias/out.log"
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}.FilterByOrientationBias.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m FilterByOrientationBias \
        -V {input.vcf} \
        -P {input.ada_det_met} \
        -O {output.vcf} &> {log}
    """
# removed final_artifact_modes

rule FuncotateMaf:
  input:
    vcf = rules.FilterByOrientationBias.output.vcf,
    vcf_index = rules.FilterByOrientationBias.output.vcf_index,
  output:
    final_out_name = "GATK_runs/{tumor}/FuncotateMaf_{normal}/out.vcf.maf.annotated",
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 5000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/FuncotateMaf_{normal}/out.log",
  conda:
    "envs_dir/pre_proc.yaml"
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    data_sources = config["data_sources"],
    ref_ver = "hg38",
    output_format = "MAF",
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}_{normal}.FuncotateMaf.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m Funcotator \
        --data-sources-path {params.data_sources} \
        --ref-version {params.ref_ver} \
        --output-file-format {params.output_format} \
        -R {params.ref_fasta} \
        -V {input.vcf} \
        -O {output.final_out_name} \
        -L {params.interval} \
        --annotation-default normal_barcode:{config['pairs'][wildcards.tumor]} \
        --annotation-default tumor_barcode:{wildcards.tumor} \
    """
# took out
#         ${"--transcript-selection-mode " + transcript_selection_mode} \
#         ${"--transcript-list " + transcript_selection_list} \
#        --annotation-default Center:${default="Unknown" sequencing_center} \
#        --annotation-default source:${default="Unknown" sequence_source} \
#         ${filter_funcotations_args} \
