# Selecting a process:
# [1] nu_e -> nu_e e+ e-  [7] nu_mu -> nu_mu e+ e- 
# [2] nu_e -> nu_e mu+ mu- [8] nu_mu -> nu_mu mu+ mu- 
# [3] nu_e -> nu_mu mu+ e- [9] nu_mu -> nu_e e+ mu- 

# Choosing target material:
# [Ar] Argon [Fe] Iron [proton] free proton [neutron] free neutron

# Choosing flux to use:
# [1] fixed energy [2] uboone (BNB) flux [3] uboone (NuMI) flux [4] BNB isotropic [5] import a flux file

# Stanadard Model or BSM 
# 4F - set vector and axial coupling modifications, 
# LmuLtau - gauged LaLb Z' model, where a/b flavours dicated by the choice of process 
# SM - Standard model trident production
# [4F], [SM], [LmuLtau]

# [GenerateEvents], [CrossSection]

# Optional (if generating events):
# N events to generate
# Output filename (without file extension)
# [1] original "teg" output [2] Hepevt format [3] HepMC3 format

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
