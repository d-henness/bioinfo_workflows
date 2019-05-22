#!/bin/bash
#TODO add in setup of needed git repos

echo "Enter path to bioinfo_workflows directory, leave blank for default: $PWD"
read -r bio_path
if [[ "$bio_path" == '' ]]; then
  bio_path=$PWD
fi
echo "Enter path to references directory"
read -r data_path

cp ref.yaml-base ref.yaml
sed -i 's,/path/to/bioinfo_workflows,'"$bio_path"',' ref.yaml
sed -i 's,/path/to/refs,'"$data_path"',' ref.yaml
