#!/bin/bash
#TODO add in setup of needed git repos

echo "Enter path to references directory"
read -r data_path

cp ref.yaml-base ref.yaml
sed -i 's,/path/to/refs,'"$data_path"',' ref.yaml
