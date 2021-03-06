#!/usr/bin/env bash

#execute all the graphs

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
DIR=/cluster/home/trinhhe/gapbs

# GRAPHS=("barabasi_albert_9" "erdos_renyi_Ntimes5" "watts_strogatz_24_0.3" "planted_partition_128_0.5_1div4n")
GRAPHS=("watts_strogatz_24_0.3")
for G in "${GRAPHS[@]}";
do
    for file in ${DIR}/../graphs/unweighted/${G}/*
    do
        OUT=${file##*/}
        ${DIR}/mincut0 -sf "$file" -n1 > ${DIR}/../outputs/mincut0/unweighted/${G}/${OUT%.graph}.csv
        ${DIR}/mincut -sf "$file" -n1 > ${DIR}/../outputs/cores/48cores/unweighted/${G}/${OUT%.graph}.csv
        ${DIR}/stoer_wagner -sf "$file" -n1 > ${DIR}/../outputs/stoer/unweighted/${G}/${OUT%.graph}.csv
        ${DIR}/karger_stein "$file" 1 > ${DIR}/../outputs/karger/unweighted/${G}/${OUT%.graph}.csv
    done
    for file in ${DIR}/../graphs/weighted_100_5/${G}/*
    do
        OUT=${file##*/}
        ${DIR}/mincut0 -sf "$file" -n1 > ${DIR}/../outputs/mincut0/weighted_100_5/${G}/${OUT%.graph}.csv
        ${DIR}/mincut -sf "$file" -n1 > ${DIR}/../outputs/cores/48cores/weighted_100_5/${G}/${OUT%.graph}.csv
        ${DIR}/stoer_wagner -sf "$file" -n1 > ${DIR}/../outputs/stoer/weighted_100_5/${G}/${OUT%.graph}.csv
        ${DIR}/karger_stein "$file" 1 > ${DIR}/../outputs/karger/weighted_100_5/${G}/${OUT%.graph}.csv
    done
    for file in ${DIR}/../graphs/weighted_100_30/${G}/*
    do
        OUT=${file##*/}
        ${DIR}/mincut0 -sf "$file" -n1 > ${DIR}/../outputs/mincut0/weighted_100_30/${G}/${OUT%.graph}.csv
        ${DIR}/mincut -sf "$file" -n1 > ${DIR}/../outputs/cores/48cores/weighted_100_30/${G}/${OUT%.graph}.csv
        ${DIR}/stoer_wagner -sf "$file" -n1 > ${DIR}/../outputs/stoer/weighted_100_30/${G}/${OUT%.graph}.csv
        ${DIR}/karger_stein "$file" 1 > ${DIR}/../outputs/karger/weighted_100_30/${G}/${OUT%.graph}.csv
    done
done

# for G in "${GRAPHS[@]}";
# do
#     for file in ${DIR}/test/graphs/weighted/${G}/*
#     do
#         ${DIR}/mincut1 -sf "$file" -n1 > ${DIR}/test/out/mincut1/weighted/${G}/"$file"
#     done
# done