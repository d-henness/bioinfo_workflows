import sys

path = sys.argv[1]
sample_name = sys.argv[2]
purity = sys.argv[3]
print(f"""
# Specifies working directory for analysis. All paths in the rest of the file are relative to this.
working_dir: {path}

# Where the trace (output) from the PyClone MCMC analysis will be written.
trace_dir: trace

# Specifies which density will be used to model read counts. Most people will want pyclone_beta_binomial or pyclone_binomial
density: pyclone_beta_binomial

# Number of iterations of the MCMC chain.
num_iters: 10000

# Specifies parameters in Beta base measure for DP. Most people will want the values below.
base_measure_params:
  alpha: 1
  beta: 1

# Specifies initial values and prior parameters for the prior on the concentration (alpha) parameter in the DP. If the prior node is not set the concentration will not be estimated and the specified value will be used.
concentration:
# Initial value if prior is set, or fixed value otherwise for concentration parameter.
  value: 1.0

# Specifies the parameters in the Gamma prior over the concentration parameter.
  prior:
    shape: 1.0
    rate: 0.001

beta_binomial_precision_params:
  # Starting value
  value: 1000

  # Parameters for Gamma prior distribution
  prior:
    shape: 1.0
    rate: 0.0001

  # Precision of Gamma proposal function for MH step
  proposal:
    precision: 0.01

samples:
  # Unique sample ID
  {sample_name}:
    # Path where YAML formatted mutations file for the sample is placed.
    mutations_file: {sample_name}_parse_input.yaml

    tumour_content:
      # The predicted tumour content for the sample. If you have no estimate set this to 1.0.
      value: {purity}

    # Expected sequencing error rate for sample
    error_rate: 0.001
    """
    )
