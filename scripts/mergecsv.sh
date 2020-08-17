#!/usr/bin/env bash

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# DIR=/cluster/home/trinhhe/outputs
DIR=/home/henry/Documents/ETH/BA/outputs

GRAPHS=("barabasi_albert_9" "erdos_renyi_Ntimes5" "watts_strogatz_24_0.3" "planted_partition_128_0.5_1div4n")
# GRAPHS=('watts_strogatz_24_0.3')
MINCUT=('mincut0' 'mincut' 'stoer' 'karger')
# MINCUT=('mincut')
# HEADER="nodes, edges, mincut estimate, packing treshold, \"H size, p, packing time, packing value, packing size\", sample size, trees generator time, cut values of sampled trees, end result vs correct result"     # define header
# MERGED="ba.csv" 

for mincut in "${MINCUT[@]}";
do
    for g in "${GRAPHS[@]}";
    do
        cat ${DIR}/${mincut}/unweighted/${g}/*.csv > ${mincut}_unweighted_${g}.csv
    done
    
done

                    

