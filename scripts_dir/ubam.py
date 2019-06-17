import sys
from xopen import xopen
import subprocess
import shlex

read1 = sys.argv[1]
read2 = sys.argv[2]
library = sys.argv[3]
output_bam = sys.argv[4]

if '/SRR' in read1:                 # check to see if the file stars with SRR
    flow_cell = 'spoof_flow_cell'   # if it does then assume it came from ncbi
    lane = 0                        # and spoof the information we want
else:
    with xopen(read1, 'r') as gatk_file:
        line = gatk_file.readline()
    if '@SRR' in line:
        flow_cell = 'spoof_flow_cell'   # above case might not catch everything, this should get remaining
        lane = 0
    else:
        line_split = line.split(':')
        flow_cell = line_split[2]
        lane = line_split[3]
ubam_command = shlex.split(f'picard FastqToSam FASTQ={read1} FASTQ2={read2} OUTPUT={output_bam} READ_GROUP_NAME="{flow_cell}.{lane}.{library}" SAMPLE_NAME="{library}" LIBRARY_NAME="{library}" PLATFORM_UNIT="{flow_cell}.{lane}.{library}" PLATFORM=ILLUMINA TMP_DIR="{flow_cell}.{lane}.{library}_tmp"')
subprocess.check_call(ubam_command)
