#!/usr/bin/env bash

#execute all the graphs

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."

# GRAPHS=("barabasi_albert_9" "erdos_renyi_Ntimes5" "watts_strogatz_24_0.3" "random_partition_graph_0.25_1div4n")
GRAPHS=("random_partition_graph_0.25_1div4n")
for G in "${GRAPHS[@]}";
do
    for file in ${DIR}/test/graphs/unweighted/${G}/*
    do
        OUT=${file##*/}
        # ${DIR}/mincut1 -sf "$file" -n1 > ../test/out/mincut1/unweighted/${G}/${OUT%.wel}.csv
        # ${DIR}/mincut2 -sf "$file" -n1 > ../test/out/mincut2/unweighted/${G}/${OUT%.wel}.csv
        ${DIR}/mincut3 -sf "$file" -n1  > ../test/out/mincut3/unweighted/${G}/${OUT%.wel}.csv
        ${DIR}/mincut4 -sf "$file" -n1  > ../test/out/mincut4/unweighted/${G}/${OUT%.wel}.csv

    done
done

# for G in "${GRAPHS[@]}";
# do
#     for file in ${DIR}/test/graphs/weighted/${G}/*
#     do
#         ${DIR}/mincut1 -sf "$file" -n1 > ${DIR}/test/out/mincut1/weighted/${G}/"$file"
#     done
# done