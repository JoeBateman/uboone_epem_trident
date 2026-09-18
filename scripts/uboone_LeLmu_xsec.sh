#!/bin/bash

gprime=1e-3

mZprime_list=(0.01 ) # 0.05 0.1 0.250 0.5 0.75 1.0

for mZprime in ${mZprime_list[@]};
do
    echo "Running for mZ' = $mZprime GeV, g' = $gprime"
    ./TEG_v2 <<EOF
7
Ar
2
LmuLtau
$gprime
$mZprime
CrossSection
EOF
done