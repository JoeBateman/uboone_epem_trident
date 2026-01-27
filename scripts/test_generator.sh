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
2
SM
GenerateEvents
5000
outputs/numu_epem_5000
2
EOF