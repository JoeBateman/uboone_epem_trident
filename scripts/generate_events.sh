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
4
SM
GenerateEvents
20000
outputs/epem_isotropic_20k
2
EOF
