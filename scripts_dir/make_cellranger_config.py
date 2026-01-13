import os
import argparse
import yaml

parser = argparse.ArgumentParser()
parser.add_argument('dirs', nargs = '+')
parser.add_argument('-y', '--yaml_file', default = "config.yaml")
parser.add_argument('-s', '--sample_string', default = None)
args = parser.parse_args()

samples = {}

for thing in os.listdir("."):
    if os.path.isdir(thing):
        if thing in args.dirs:
            if args.sample_string is None:
                samples[thing] = thing 
            else:
                samples[thing] = args.sample_string

config = {'samples': samples}
with open(args.yaml_file, 'w') as data:
    yaml.dump(config, data, default_flow_style = False)

