#!/bin/bash


# Loop over a list of energy values (in GeV)
beam_energy=62


echo "Getting the SM Cross-Section"
    ./TEG_v2 <<EOF
7 
Ar
2
/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/numi/g4lbne_RHC_ND_${beam_energy}_GeV.root
SM
CrossSection
EOF

mZ_list=( 0.001 0.00158 0.00251 0.00398 0.00631 0.01 0.0158 0.0251 0.0398 0.0631 0.1 0.158 0.251 0.398 0.631 1.0 1.58 2.51 3.98 6.31 10.0 )
for mZ in ${mZ_list[@]};
do

# e+ e- final state, Ar target, BSM interactions
echo "Running for mZ = $mZ GeV"
    ./TEG_v2 <<EOF
7 
Ar
2
/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/data/flux/numi/g4lbne_RHC_ND_62_GeV.root
LmuLtau
8e-4
$mZ
CrossSection
EOF
done