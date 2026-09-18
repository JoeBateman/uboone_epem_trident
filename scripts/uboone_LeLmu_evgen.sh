#!/bin/bash

gprime=0.001
mZprime=0.01 # 0.1 0.250 0.5 0.75 1.0

echo "Running for mZ' = $mZprime GeV, g' = $gprime"
./TEG_v2 <<EOF
7
Ar
2
LmuLtau
$gprime
$mZprime
GenerateEvents
250
outputs/0.001_0.01_epem
2
