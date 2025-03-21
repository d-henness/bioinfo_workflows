configfile: "{}/ref.yaml".format(workflow.basedir)

include: "pre_pro_af_merge_alt_bed.snakefile"

tumors = config["pairs"]
#normals = [config["pairs"][tumor] for tumor in tumors]

rule mutect_all:
  input:
    expand("easy_transfer/{tumor}_{f_score_thresh}_VEP_parsed.tsv", tumor = tumors, f_score_thresh = [0.8]),
    expand("GATK_runs/{tumor}_{f_score_thresh}/VEP/{tumor}_{f_score_thresh}.vcf", tumor = tumors, f_score_thresh = [0.8])


#rule CollectF1R2Counts:
#  input:
#    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
#    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
#  output:
#    alt_table = temp("GATK_runs/{tumor}/CollectF1R2Counts/alt_tab.tsv"),
#    ref_hist = temp("GATK_runs/{tumor}/CollectF1R2Counts/alt_tab.metrics"),
#    ref_hists = temp("GATK_runs/{tumor}/CollectF1R2Counts/alt-depth1.metrics"),
#    tumor_name = temp("GATK_runs/{tumor}/CollectF1R2Counts/tumor_name.txt"),
#  params:
#    ref_fasta = config["ref_fasta"],
#    ref_fai = config["ref_fasta_index"],
#    ref_dict = config["ref_dict"],
#    gnomad = config["gnomad"],
#    interval = config["alt_bed"],
#    exclude_list = ''
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 4000,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
#  log:
#    "GATK_runs/{tumor}/CollectF1R2Counts/out.log",
#  conda:
#    "envs_dir/pre_proc.yaml"
#  benchmark:
#    "benchmarks/{tumor}.CollectF1R2Counts.benchmark.txt"
#  shell:
#    """
#      gatk --java-options -Xmx4000m GetSampleName -R {params.ref_fasta} -I {input.tumor_bam} -O {output.tumor_name} -encode
#      tumor_name=$(head -n 1 {output.tumor_name})
#
#			gatk --java-options -Xmx4000m CollectF1R2Counts \
#      -I {input.tumor_bam} -R {params.ref_fasta} \
#      -L {params.interval} \
#      -alt-table "{output.alt_table}" \
#      -ref-hist "{output.ref_hist}" \
#			-alt-hist "{output.ref_hists}" &> {log}
#    """
#

rule M2:
  conda:
    "envs_dir/pre_proc.yaml"
  input:
#    normal_bam = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam",
#    normal_bai = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bai",
#    normal_bam = lambda wildcards: "".join(["GATK_runs/",config['pairs'][wildcards.tumor],"/ApplyBQSR/",config['pairs'][wildcards.tumor],"_recal.bam"]),
#    normal_bai = lambda wildcards: "".join(["GATK_runs/",config['pairs'][wildcards.tumor],"/ApplyBQSR/",config['pairs'][wildcards.tumor],"_recal.bam.bai"]),
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
  output:
    vcf = temp("GATK_runs/{tumor}/M2/out.vcf"),
    vcf_index = temp("GATK_runs/{tumor}/M2/out.vcf.idx"),
    bam_out = temp("GATK_runs/{tumor}/M2/out.bam"),
    tumor_name = temp("GATK_runs/{tumor}/M2/tumor_name.txt"),
#    normal_name = temp("GATK_runs/{tumor}/M2/normal_name.txt"),
    f1r2_tar_gz = temp("GATK_runs/{tumor}/M2/f1r2.tar.gz"),
#  input:
#    normal_bam = "GATK_runs/{normal}/ApplyBQSR/{normal}_recal.bam",
#    normal_bia = "GATK_runs/{normal}/ApplyBQSR/{normal}_recal.bam.bai",
#    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
#    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
#    art_tab = rules.LearnReadOrientationModel.output
#  output:
#    vcf = temp("GATK_runs/{tumor}/M2_{normal}/out.vcf"),
#    vcf_index = temp("GATK_runs/{tumor}/M2_{normal}/out.vcf.idx"),
#    bam_out = temp("GATK_runs/{tumor}/M2_{normal}/out.bam"),
#    tumor_name = temp("GATK_runs/{tumor}/M2_{normal}/tumor_name.txt"),
#    normal_name = temp("GATK_runs/{tumor}/M2_{normal}/normal_name.txt"),
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    ref_pon = "/home/dylan/ref_data/hg38/1000g_pon.hg38.vcf.gz",
    interval = config["alt_bed"],
    exclude_list = 'neuron,biolx95'
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
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

      gatk --java-options -Xmx4000m Mutect2 \
      -R {params.ref_fasta} \
      ${{tumor_command_line}} \
      --germline-resource {params.gnomad} \
      -L {params.interval} \
      -O "{output.vcf}" \
      --bam-output {output.bam_out} \
      -pon {params.ref_pon} \
      --f1r2-tar-gz {output.f1r2_tar_gz} &> {log}
    """

rule LearnReadOrientationModel:
  conda:
    "envs_dir/pre_proc.yaml"
  input:
    f1r2_tar_gz = rules.M2.output.f1r2_tar_gz,
  output:
    art_tab = temp("GATK_runs/{tumor}/LearnReadOrientationModel/art_tab.tsv.tar.gz"),
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    exclude_list = ''
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}/LearnReadOrientationModel/log.log",
  benchmark:
    "benchmarks/{tumor}.LearnReadOrientationModel.benchmark.txt"
  shell:
    """
			gatk --java-options -Xmx4000m LearnReadOrientationModel \
			-I "{input.f1r2_tar_gz}" \
			-O "{output.art_tab}" &> {log}
    """

#rule CollectSequencingArtifactMetrics:
#  input:
#    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
#    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
#  output:
#    bait_det_met = temp("GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.bait_bias_detail_metrics"),
#    bait_sum_met = temp("GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.bait_bias_summary_metrics"),
#    err_sum_met = temp("GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.error_summary_metrics"),
#    ada_det_met = temp("GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_detail_metrics"),
#    ada_sum_met = temp("GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_summary_metrics"),
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 4000,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
#  log:
#    "GATK_runs/{tumor}/CollectSequencingArtifactMetrics/gatk.log"
#  conda:
#    "envs_dir/pre_proc.yaml"
#  params:
#    ref_fasta = config["ref_fasta"],
#    ref_fai = config["ref_fasta_index"],
#    ref_dict = config["ref_dict"],
#    gnomad = config["gnomad"],
#    interval = config["alt_bed"],
#    exclude_list = ''
#  benchmark:
#    "benchmarks/{tumor}.CollectSequencingArtifactMetrics.benchmark.txt"
#  shell:
#      """
#        gatk --java-options -Xmx4000m CollectSequencingArtifactMetrics \
#        -I {input.tumor_bam} -O GATK_runs/{wildcards.tumor}/CollectSequencingArtifactMetrics/gatk -R {params.ref_fasta} -VALIDATION_STRINGENCY LENIENT &> {log}
#      """

rule CalculateContamination:
  input:
    tumor_bam = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam",
    tumor_bai = "GATK_runs/{tumor}/ApplyBQSR/{tumor}_recal.bam.bai",
#    normal_bam = lambda wildcards: "".join(["GATK_runs/",config['pairs'][wildcards.tumor],"/ApplyBQSR/",config['pairs'][wildcards.tumor],"_recal.bam"]),
#    normal_bai = lambda wildcards: "".join(["GATK_runs/",config['pairs'][wildcards.tumor],"/ApplyBQSR/",config['pairs'][wildcards.tumor],"_recal.bam.bai"]),
#    normal_bam = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam",
#    normal_bai = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bai",
#    normal_bam = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam",
#    normal_bai = lambda wildcards: f"GATK_runs/{config['pairs'][wildcards.tumor]}/ApplyBQSR/{config['pairs'][wildcards.tumor]}_recal.bam.bai",
#    normal_bam = "GATK_runs/{normal}/ApplyBQSR/{normal}_recal.bam",
#    normal_bia = "GATK_runs/{normal}/ApplyBQSR/{normal}_recal.bam.bai",
  output:
#    normal_pil_tab = temp("GATK_runs/{tumor}/CalculateContamination/normal_pileups.table"),
    pil_tab = temp("GATK_runs/{tumor}/CalculateContamination/pileups.table"),
    con_tab = temp("GATK_runs/{tumor}/CalculateContamination/con_tab.table"),
    seg_tab = temp("GATK_runs/{tumor}/CalculateContamination/seg_tab.table"),
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
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
    interval = config["alt_bed"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}.CalculateContamination.benchmark.txt"
  shell:
    """

      gatk --java-options -Xmx4000m GetPileupSummaries -R {params.ref_fasta} -I {input.tumor_bam} --interval-set-rule INTERSECTION -L {params.interval} \
        -V {params.variants_con} -L {params.variants_con} -O {output.pil_tab} &>> {log}
      gatk CalculateContamination -I {output.pil_tab} -O {output.con_tab} --tumor-segmentation {output.seg_tab} &>> {log}
    """

#rule MergeVCFs:
#  input:
#    vcfs = rules.M2.output.vcf,
#    vcf_index = rules.M2.output.vcf_index,
#  output:
#    merged_vcf = temp("GATK_runs/{tumor}/MergeVCFs/out.vcf"),
#    merged_vcf_index = temp("GATK_runs/{tumor}/MergeVCFs/out.vcf.idx"),
#  conda:
#    "envs_dir/pre_proc.yaml"
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 4000,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
#  log:
#    "GATK_runs/{tumor}/MergeVCFs/out.log",
#  params:
#    ref_fasta = config["ref_fasta"],
#    ref_fai = config["ref_fasta_index"],
#    ref_dict = config["ref_dict"],
#    gnomad = config["gnomad"],
#    interval = config["alt_bed"],
#    variants_con = config["variants_for_contamination"],
#    variants_con_index = config["variants_for_contamination_index"],
#    exclude_list = ''
#  benchmark:
#    "benchmarks/{tumor}.MergeVCFs.benchmark.txt"
#  shell:
#    """
#      gatk --java-options -Xmx4000m MergeVcfs -I {input.vcfs} -O {output.merged_vcf} &> {log}
#    """

rule Filter:
  input:
    vcf = rules.M2.output.vcf,
    vcf_index = rules.M2.output.vcf_index,
    con_tab = rules.CalculateContamination.output.con_tab,
    maf_seg = rules.CalculateContamination.output.seg_tab,
    art_tab = rules.LearnReadOrientationModel.output.art_tab,
  output:
    filt_vcf = "GATK_runs/{tumor}_{f_score_thresh}/Filter/{tumor}_{f_score_thresh}.vcf",
    filt_vcf_index = "GATK_runs/{tumor}_{f_score_thresh}/Filter/{tumor}_{f_score_thresh}.vcf.idx",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 4000,
    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
  log:
    "GATK_runs/{tumor}_{f_score_thresh}/Filter/out.log"
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["alt_bed"],
    variants_con = config["variants_for_contamination"],
    variants_con_index = config["variants_for_contamination_index"],
    exclude_list = ''
  benchmark:
    "benchmarks/{tumor}_{f_score_thresh}.Filter.benchmark.txt"
  shell:
    """
      gatk --java-options -Xmx4000m FilterMutectCalls -V {input.vcf} \
        -O {output.filt_vcf} \
	-R {params.ref_fasta} \
      	--contamination-table {input.con_tab} \
	--ob-priors {input.art_tab} \
	-f-score-beta {wildcards.f_score_thresh} \
      	--tumor-segmentation {input.maf_seg} &> {log} \
    """

#rule FilterByOrientationBias:
#  input:
#    vcf = rules.Filter.output.filt_vcf,
#    vcf_index = rules.Filter.output.filt_vcf_index,
#    ada_det_met = rules.CollectSequencingArtifactMetrics.output.ada_det_met
#  output:
#    vcf = "GATK_runs/{tumor}/FilterByOrientationBias/{tumor}.vcf",
#    vcf_index = "GATK_runs/{tumor}/FilterByOrientationBias/{tumor}.vcf.idx",
#  conda:
#    "envs_dir/pre_proc.yaml"
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 4000,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
#  log:
#    "GATK_runs/{tumor}/FilterByOrientationBias/out.log"
#  params:
#    ref_fasta = config["ref_fasta"],
#    ref_fai = config["ref_fasta_index"],
#    ref_dict = config["ref_dict"],
#    gnomad = config["gnomad"],
#    interval = config["alt_bed"],
#    variants_con = config["variants_for_contamination"],
#    variants_con_index = config["variants_for_contamination_index"],
#    exclude_list = ''
#  benchmark:
#    "benchmarks/{tumor}.FilterByOrientationBias.benchmark.txt"
#  shell:
#    """
#      gatk --java-options -Xmx4000m FilterByOrientationBias \
#        -V {input.vcf} \
#        -P {input.ada_det_met} \
#        -O {output.vcf} &> {log}
#    """

rule VEP:
  input:
    vcf = rules.Filter.output.filt_vcf,
  output:
    vcf_out = "GATK_runs/{tumor}_{f_score_thresh}/VEP/{tumor}_{f_score_thresh}.vcf",
    vcf_out_zip = "vep/{tumor}_{f_score_thresh}/VEP/{tumor}_{f_score_thresh}_VEP.vcf.gz",
    summary = "GATK_runs/{tumor}_{f_score_thresh}/VEP/{tumor}_{f_score_thresh}.vcf_summary.html",
    parsed_output = "GATK_runs/{tumor}_{f_score_thresh}/VEP/{tumor}_{f_score_thresh}_VEP_parsed.tsv",
    easy_transfer = "easy_transfer/{tumor}_{f_score_thresh}_VEP_parsed.tsv",
  conda: "envs_dir/vep.yaml"
  log: "GATK_runs/log/{tumor}_{f_score_thresh}_VEP.log"
  benchmark: "GATK_runs/benchmark/{tumor}_{f_score_thresh}_snp_VEP.benchmark"
  params:
    ref_fasta = config["ref_fasta"],
    vep_cache = config["vep_cache"],
    vep_plugins = config["vep_plug_dir"],
    dbNSFP_config = config["dbNSFP_config"],
    condel_config = config["condel_config"],
    loftool_config = config["loftool_config"],
    bioinfo_workflows_path = config["bioinfo_workflows_path"],
  threads: 1
  resources:
    mem_mb = lambda wildcards, attempt: attempt * 3000,
  shell:
    """
      vep \
        --input_file {input.vcf} \
        --output_file {output.vcf_out} \
        --format vcf \
        --vcf \
        --symbol \
        --terms SO \
        --tsl \
        --hgvs \
        --hgvsg \
        --fasta {params.ref_fasta} \
        --offline --cache --dir_cache {params.vep_cache} \
        --dir_plugins {params.vep_plugins} \
        --plugin Downstream \
        --plugin dbNSFP,{params.dbNSFP_config},ALL \
        --plugin Condel,{params.condel_config},b,2 \
        --plugin LoFtool,{params.loftool_config} \
        --plugin Blosum62 \
        --pick \
        --sift b \
        --polyphen b \
        --transcript_version &> {log}
      bgzip -c {output.vcf_out} > {output.vcf_out_zip}
      bcftools index -f {output.vcf_out_zip}

      python3 {params.bioinfo_workflows_path}/scripts_dir/vep_vcf_parser.py \
        -f SYMBOL Gene Consequence SIFT PolyPhen Condel LoFtool BLOSUM62 Protein_position Amino_acids Codons \
        -v {output.vcf_out} -p > {output.parsed_output}

      if [[ !(-d easy_transfer) ]]; then
        mkdir easy_transfer
      fi

      cp {output.vcf_out} {output.parsed_output} easy_transfer
    """
# removed final_artifact_modes

#rule FuncotateMaf:
#  input:
#    vcf = rules.FilterByOrientationBias.output.vcf,
#    vcf_index = rules.FilterByOrientationBias.output.vcf_index,
#  output:
#    final_out_name = "GATK_runs/{tumor}/FuncotateMaf_{normal}/out.vcf.maf.annotated",
#  resources:
#    mem_mb = lambda wildcards, attempt: attempt * 5000,
#    time_min = lambda wildcards, attempt: attempt * 24 * 60,	# time in minutes
#  log:
#    "GATK_runs/{tumor}/FuncotateMaf_{normal}/out.log",
#  conda:
#    "envs_dir/pre_proc.yaml"
#  params:
#    ref_fasta = config["ref_fasta"],
#    ref_fai = config["ref_fasta_index"],
#    ref_dict = config["ref_dict"],
#    gnomad = config["gnomad"],
#    interval = config["alt_bed"],
#    variants_con = config["variants_for_contamination"],
#    variants_con_index = config["variants_for_contamination_index"],
#    data_sources = config["data_sources"],
#    ref_ver = "hg38",
#    output_format = "MAF",
#    exclude_list = ''
#  benchmark:
#    "benchmarks/{tumor}_{normal}.FuncotateMaf.benchmark.txt"
#  shell:
#    """
#      gatk --java-options -Xmx4000m Funcotator \
#        --data-sources-path {params.data_sources} \
#        --ref-version {params.ref_ver} \
#        --output-file-format {params.output_format} \
#        -R {params.ref_fasta} \
#        -V {input.vcf} \
#        -O {output.final_out_name} \
#        -L {params.interval} \
#        --annotation-default normal_barcode:{wildcards.normal} \
#        --annotation-default tumor_barcode:{wildcards.tumor} \
#    """
## took out
##         ${"--transcript-selection-mode " + transcript_selection_mode} \
##         ${"--transcript-list " + transcript_selection_list} \
##        --annotation-default Center:${default="Unknown" sequencing_center} \
##        --annotation-default source:${default="Unknown" sequence_source} \
##         ${filter_funcotations_args} \
