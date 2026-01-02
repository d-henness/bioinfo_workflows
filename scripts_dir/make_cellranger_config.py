import os
import argparse
import yaml

parser = argparse.ArgumentParser()
parser.add_argument('dir_strings', nargs = '+')
parser.add_argument('-y', '--yaml_file', default = "config.yaml")
args = parser.parse_args()

samples = []

for thing in os.listdir("."):
    if os.path.isdir(thing):
        for dir_string in args.dir_strings:
            if thing.startswith(dir_string):
                samples.append(thing)

config = {'samples': samples}
with open(args.yaml_file, 'w') as data:
    yaml.dump(config, data, default_flow_style = False)

