NEvents=$1
outputDir="$DATA/DarkNews/Trident/outputs/root/"
OutName="tfg_hepmc_${NEvents}.root"

lar -c fcls/TextFileGen.fcl -o ${outputDir}$OutName -n $NEvents
lar -c wirecell_g4_uboone.fcl -s ${outputDir}$OutName -o ${outputDir}g4_$OutName -n $NEvents
OutName=g4_$OutName
lar -c wirecell_detsim_uboone.fcl -s ${outputDir}$OutName -o ${outputDir}detsim_$OutName -n $NEvents
OutName=detsim_$OutName
lar -c reco_uboone_mcc9_8_driver_stage1.fcl -s ${outputDir}$OutName -o ${outputDir}reco1_$OutName -n $NEvents
OutName=reco1_$OutName
lar -c reco_uboone_mcc9_8_driver_stage2_fullMC.fcl -s ${outputDir}$OutName -o ${outputDir}reco2_$OutName -n $NEvents

rm *.db *.log *.root *.pndr # Clean up intermediate files