# 7 -> nu_mu -> nu_mu e+ e-
# Ar target
# 2 -> use uboone flux file
# Use standard model
# Generate events
# 50 -> number of events to generate
# Output filename: numu_epem_50
# 2 -> Hepevt format

./TEG_v2 <<EOF
7
Ar
3
SM
GenerateEvents
5000
outputs/numu_numi_5000
2
EOF

./TEG_v2 <<EOF
1
Ar
3
SM
GenerateEvents
5000
outputs/nue_numi_5000
2
EOF