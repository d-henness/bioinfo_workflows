configfile: "{}/ref.yaml".format(workflow.basedir)

include: "pre_pro_af_merge.snakefile"

tumors = config["pairs"]
normals = [config["pairs"][tumor] for tumor in tumors]

rule all:
  input:
#   expand("runs/{tumor}/M2_{normal}/out.vcf", zip, tumor = tumors, normal = normals)
#   expand("runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_detail_metrics", tumor = tumors)
#   expand("runs/{tumor}/CalculateContamination_{normal}/seg_tab.table", tumor = tumors, normal = normals)
#   expand("runs/{tumor}/MergeVCFs_{normal}/out.vcf", tumor = tumors, normal = normals)
#   expand("runs/{tumor}/Filter_{normal}/out.vcf", tumor = tumors, normal = normals)
    expand("runs/{tumor}/FilterByOrientationBias_{normal}/out.vcf", zip, tumor = tumors, normal = normals)
#   expand("runs/{tumor}/FuncotateMaf_{normal}/out.vcf.maf.annotated", zip, tumor = tumors, normal = normals)


rule CollectF1R2Counts:
  input:
    tumor_bam = "runs/{tumor}/ApplyBQSR/recal.bam",
    tumor_bai = "runs/{tumor}/ApplyBQSR/recal.bai",
  output:
    alt_table = "runs/{tumor}/CollectF1R2Counts/alt_tab.tsv",
    ref_hist = "runs/{tumor}/CollectF1R2Counts/alt_tab.metrics",
    ref_hists = "runs/{tumor}/CollectF1R2Counts/alt-depth1.metrics",
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    exclude_list = ''
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/CollectF1R2Counts/out.log",
  conda:
    "envs_dir/pre_proc.yaml"
  benchmark:
    "benchmarks/{tumor}.CollectF1R2Counts.benchmark.txt"
  shell:
    """
      gatk GetSampleName -R {params.ref_fasta} -I {input.tumor_bam} -O tumor_name.txt -encode
      tumor_name=$(head -n 1 tumor_name.txt)

			gatk CollectF1R2Counts \
      -I {input.tumor_bam} -R {params.ref_fasta} \
      -L {params.interval} \
      -alt-table "{output.alt_table}" \
      -ref-hist "{output.ref_hist}" \
			-alt-hist "{output.ref_hists}"
    """

rule LearnReadOrientationModel:
  conda:
    "envs_dir/pre_proc.yaml"
  input:
    alt_table = rules.CollectF1R2Counts.output.alt_table,
    ref_hist = rules.CollectF1R2Counts.output.ref_hist,
    ref_hists = rules.CollectF1R2Counts.output.ref_hists,
  output:
    "runs/{tumor}/LearnReadOrientationModel/art_tab.tsv",
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    exclude_list = ''
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/LearnReadOrientationModel/log.log",
  benchmark:
    "benchmarks/{tumor}.LearnReadOrientationModel.benchmark.txt"
  shell:
    """
			gatk LearnReadOrientationModel \
			-alt-table {input.alt_table} \
			-ref-hist {input.ref_hist} \
			-alt-hist {input.ref_hists} \
			-O "{output}"
    """

rule M2:
  conda:
    "envs_dir/pre_proc.yaml"
  input:
    normal_bam = "runs/{normal}/ApplyBQSR/recal.bam",
    normal_bia = "runs/{normal}/ApplyBQSR/recal.bai",
    tumor_bam = "runs/{tumor}/ApplyBQSR/recal.bam",
    tumor_bai = "runs/{tumor}/ApplyBQSR/recal.bai",
    art_tab = rules.LearnReadOrientationModel.output
  output:
    vcf = "runs/{tumor}/M2_{normal}/out.vcf",
    vcf_index = "runs/{tumor}/M2_{normal}/out.vcf.idx",
    bam_out = "runs/{tumor}/M2_{normal}/out.bam",
  params:
    ref_fasta = config["ref_fasta"],
    ref_fai = config["ref_fasta_index"],
    ref_dict = config["ref_dict"],
    gnomad = config["gnomad"],
    interval = config["exom_padded_primary"],
    exclude_list = 'neuron,biolx95'
  resources:
    mem_mb = 5000,
  threads: 4 # default for mutect asks for 4 threads
  log:
    "runs/{tumor}/M2_{normal}/out.log"
  benchmark:
    "benchmarks/{tumor}_{normal}.M2.benchmark.txt"
  shell:
    """
      gatk GetSampleName -R {params.ref_fasta} -I {input.tumor_bam} -O tumor_name.txt -encode
      tumor_command_line="-I {input.tumor_bam} -tumor `cat tumor_name.txt`"

      gatk GetSampleName -R {params.ref_fasta} -I {input.normal_bam} -O normal_name.txt -encode 
      normal_command_line="-I {input.normal_bam} -normal `cat normal_name.txt`"

      gatk --java-options -Xmx4000m Mutect2 \
      -R {params.ref_fasta} \
      ${{tumor_command_line}} \
      ${{normal_command_line}} \
      --germline-resource {params.gnomad} \
      -L {params.interval} \
      -O "{output.vcf}" \
      --bam-output {output.bam_out} \
      --orientation-bias-artifact-priors {input.art_tab}  \
    """

rule CollectSequencingArtifactMetrics:
  input:
    tumor_bam = "runs/{tumor}/ApplyBQSR/recal.bam",
    tumor_bai = "runs/{tumor}/ApplyBQSR/recal.bai",
  output:
    bait_det_met = "runs/{tumor}/CollectSequencingArtifactMetrics/gatk.bait_bias_detail_metrics",
    bait_sum_met = "runs/{tumor}/CollectSequencingArtifactMetrics/gatk.bait_bias_summary_metrics",
    err_sum_met = "runs/{tumor}/CollectSequencingArtifactMetrics/gatk.error_summary_metrics",
    ada_det_met = "runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_detail_metrics",
    ada_sum_met = "runs/{tumor}/CollectSequencingArtifactMetrics/gatk.pre_adapter_summary_metrics",
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/CollectSequencingArtifactMetrics/gatk.log"
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
        gatk CollectSequencingArtifactMetrics \
        -I {input.tumor_bam} -O runs/{wildcards.tumor}/CollectSequencingArtifactMetrics/gatk -R {params.ref_fasta} -VALIDATION_STRINGENCY LENIENT 
      """

rule CalculateContamination:
  input:
    tumor_bam = "runs/{tumor}/ApplyBQSR/recal.bam",
    tumor_bai = "runs/{tumor}/ApplyBQSR/recal.bai",
    normal_bam = "runs/{normal}/ApplyBQSR/recal.bam",
    normal_bia = "runs/{normal}/ApplyBQSR/recal.bai",
  output:
    normal_pil_tab = "runs/{tumor}/CalculateContamination_{normal}/normal_pileups.table",
    pil_tab = "runs/{tumor}/CalculateContamination_{normal}/pileups.table",
    con_tab = "runs/{tumor}/CalculateContamination_{normal}/con_tab.table",
    seg_tab = "runs/{tumor}/CalculateContamination_{normal}/seg_tab.table",
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/CalculateContamination_{normal}/out.log"
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
    "benchmarks/{tumor}_{normal}.CalculateContamination.benchmark.txt"
  shell:
    """
      gatk GetPileupSummaries -I {input.normal_bam} --interval-set-rule INTERSECTION -L {params.interval} \
        -V {params.variants_con} -L {params.variants_con} -O {output.normal_pil_tab}
      NORMAL_CMD="-matched {output.normal_pil_tab}"

      gatk GetPileupSummaries -R {params.ref_fasta} -I {input.tumor_bam} --interval-set-rule INTERSECTION -L {params.interval} \
        -V {params.variants_con} -L {params.variants_con} -O {output.pil_tab}
      gatk CalculateContamination -I {output.pil_tab} -O {output.con_tab} --tumor-segmentation {output.seg_tab} ${{NORMAL_CMD}}
    """

rule MergeVCFs:
  input:
    vcfs = rules.M2.output.vcf,
    vcf_index = rules.M2.output.vcf_index,
  output:
    merged_vcf = "runs/{tumor}/MergeVCFs_{normal}/out.vcf",
    merged_vcf_index = "runs/{tumor}/MergeVCFs_{normal}/out.vcf.idx",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/MergeVCFs_{normal}/out.log",
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
    "benchmarks/{tumor}_{normal}.MergeVCFs.benchmark.txt"
  shell:
    """
      gatk MergeVcfs -I {input.vcfs} -O {output.merged_vcf}
    """
    

rule Filter:
  input:
    vcf = rules.MergeVCFs.output.merged_vcf,
    vcf_index = rules.MergeVCFs.output.merged_vcf_index,
    con_tab = rules.CalculateContamination.output.con_tab,
    maf_seg = rules.CalculateContamination.output.seg_tab,
  output:
    filt_vcf = "runs/{tumor}/Filter_{normal}/out.vcf",
    filt_vcf_index = "runs/{tumor}/Filter_{normal}/out.vcf.idx",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/Filter_{normal}/out.log"
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
    "benchmarks/{tumor}_{normal}.Filter.benchmark.txt"
  shell:
    """
      gatk FilterMutectCalls -V {input.vcf} \
        -O {output.filt_vcf} \
      	--contamination-table {input.con_tab} \
      	--tumor-segmentation {input.maf_seg} \
    """

rule FilterByOrientationBias:
  input:
    vcf = rules.Filter.output.filt_vcf,
    vcf_index = rules.Filter.output.filt_vcf_index,
    ada_det_met = rules.CollectSequencingArtifactMetrics.output.ada_det_met
  output:
    vcf = "runs/{tumor}/FilterByOrientationBias_{normal}/out.vcf",
    vcf_index = "runs/{tumor}/FilterByOrientationBias_{normal}/out.vcf.idx",
  conda:
    "envs_dir/pre_proc.yaml"
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/FilterByOrientationBias_{normal}/out.log"
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
    "benchmarks/{tumor}_{normal}.FilterByOrientationBias.benchmark.txt"
  shell:
    """
      gatk FilterByOrientationBias \
        -V {input.vcf} \
        -P {input.ada_det_met} \
        -O {output.vcf}
    """
# removed final_artifact_modes

rule FuncotateMaf:
  input:
    vcf = rules.FilterByOrientationBias.output.vcf,
    vcf_index = rules.FilterByOrientationBias.output.vcf_index,
  output:
    final_out_name = "runs/{tumor}/FuncotateMaf_{normal}/out.vcf.maf.annotated",
  resources:
    mem_mb = 5000
  log:
    "runs/{tumor}/FuncotateMaf_{normal}/out.log",
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
      gatk Funcotator \
        --data-sources-path {params.data_sources} \
        --ref-version {params.ref_ver} \
        --output-file-format {params.output_format} \
        -R {params.ref_fasta} \
        -V {input.vcf} \
        -O {output.final_out_name} \
        -L {params.interval} \
        --annotation-default normal_barcode:{wildcards.normal} \
        --annotation-default tumor_barcode:{wildcards.tumor} \
    """
# took out
#         ${"--transcript-selection-mode " + transcript_selection_mode} \
#         ${"--transcript-list " + transcript_selection_list} \
#        --annotation-default Center:${default="Unknown" sequencing_center} \
#        --annotation-default source:${default="Unknown" sequence_source} \
#         ${filter_funcotations_args} \
