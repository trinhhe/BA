#!/usr/bin/env bash

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# GRAPHS=("barabasi_albert_9" "erdos_renyi_Ntimes5" "watts_strogatz_24_0.3" "random_partition_graph_0.25_1div4n")
GRAPHS=("barabasi_albert_9")
# MINCUT=('mincut1' 'mincut2' 'mincut3' 'mincut4')
MINCUT=("mincut3")
HEADER="nodes, edges, mincut estimate, packing treshold, \"H size, p, packing time, packing value, packing size\", sample size, trees generator time, cut values of sampled trees, end result vs correct result"     # define header
# MERGED="ba.csv" 

for mincut in "${MINCUT[@]}";
do
    for g in "${GRAPHS[@]}";
    do
        echo $HEADER > ${i%_*:(-1)}.csv
        for i in `ls -v ../test/out/${mincut}/unweighted/${g}/*.csv`; 
        do
            sed "/$HEADER/d" $i >> ${i%_*}.csv
        done
    done
    
done

                    

#TODO MAKE IT WORK PROPERLY
