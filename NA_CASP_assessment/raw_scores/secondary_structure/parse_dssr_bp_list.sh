#!/bin/bash

# ./parse_dssr_bp_list.sh DSSR.out DSSR_bp_list.csv

input_file="$1"
output_file="$2"

# get the base pair text
awk '/List of [0-9]+ base pairs/,/^$/' "$input_file" | \
  sed '1d;s/^[[:space:]]*//' | \
  awk '{print $0}' > $output_file
# format into valid csv
sed -i '1s/^/index,/g' $output_file
sed -i 's/^[[:space:]]*//' $output_file
sed -i 's/[[:space:]]\+/,/g' $output_file
