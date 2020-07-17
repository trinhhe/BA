#!/usr/bin/env bash

#execute all the graphs

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."

# GRAPHS=("barabasi_albert_3" "barabasi_albert_9" "erdos_renyi_Ntimes15" "erdos_renyi_Ntimes60" "watts_strogatz_8_0.3" "watts_strogatz_24_0.3")
GRAPHS=("barabasi_albert_9")
for G in "${GRAPHS[@]}";
do
    for file in ${DIR}/test/graphs/weighted/${G}/*
    do
        OUT=${file##*/}
        # ${DIR}/mincut2 -sf "$file" -n1 > ../test/out/mincut2/weighted/${G}/${OUT%.wel}
        ${DIR}/mincut1 -sf "$file" -n1 > ../test/out/mincut1/weighted/${G}/${OUT%.wel}
        # ${DIR}/mincut3 -sf "$file" -n1 > ../test/out/mincut3/weighted/${G}/${OUT%.wel}

        # ${DIR}/mincut1 -sf "$file" -n1 > ../test/out/mincut1/weighted/${G}/${OUT%.wel}
        # ${DIR}/mincut2 -sf "$file" -n1 > ../test/out/mincut2/weighted/${G}/${OUT%.wel}

    done
done

# for G in "${GRAPHS[@]}";
# do
#     for file in ${DIR}/test/graphs/weighted/${G}/*
#     do
#         ${DIR}/mincut1 -sf "$file" -n1 > ${DIR}/test/out/mincut1/weighted/${G}/"$file"
#     done
# done