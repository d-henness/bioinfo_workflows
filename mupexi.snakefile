configfile: "{}/ref.yaml".format(workflow.basedir)

rule run_mupexi:
  input:
    expand("mupexi_runs/{tumor}/{tumor}.mupexi", tumor = config["vcf"])


rule mupexi:
  input:
    vcf = lambda wildcards: config["vcf"][wildcards.tumor],
    hla = lambda wildcards: config["hla"][wildcards.tumor],
  output:
    mupexi_out = "mupexi_runs/{tumor}/{tumor}.mupexi",
  conda:
    "./envs_dir/mupexi.yaml"
  log:
    mupexi_log = "mupexi_runs/log/{tumor}.log",
  benchmark: "mupexi_runs/benchmark/{tumor}_mupexi.benchmark"
  params:
    mupexi_path = "{workflow.basedir}/MuPeXI/MuPeXI.py",
    mupexi_config = "{workflow.basedir}/MuPeXI/config.ini",
    out_dir = "mupexi_runs/{tumor}",
  threads: 1
  resources:
    mem_mb = 4000,
  shell:
    """
    hla_string=$(sed -E 's/^\S*\s*(([A-Za-z0-9\:\-]*\s*){{6}}).*/\1/; s/\s+([A-Za-z0-9\:\-]+)/\,\1/g' {input.hla})
    if [[ "$hla_string" == '' ]]; then
      echo Problem with hla file
      exit 1
    fi
    {params.mupexi_path} -v {input.vcf} -a $hla_string -c {params.mupexi_config} -t -d {params.out_dir} -l 8-11 -f
    """
