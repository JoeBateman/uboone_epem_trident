# Selecting a process:
# [1] nu_e -> nu_e e+ e-  [7] nu_mu -> nu_mu e+ e- 
# [2] nu_e -> nu_e mu+ mu- [8] nu_mu -> nu_mu mu+ mu- 
# [3] nu_e -> nu_mu mu+ e- [9] nu_mu -> nu_e e+ mu- 

# Choosing target material:
# [Ar] Argon [Fe] Iron [proton] free proton [neutron] free neutron

# Choosing flux to use:
# [1] fixed energy [2] uboone (BNB) flux [3] uboone (NuMI FHC) flux [4] uboone (NuMI RHC) flux [5] BNB isotropic [6] import a flux file

# Stanadard Model or BSM 
# 4F - set vector and axial coupling modifications, 
# LmuLtau - gauged LaLb Z' model, where a/b flavours dicated by the choice of process - set gprime and mZp before choosing CrossSection or GenerateEvents
# SM - Standard model trident production
# [4F], [SM], [LmuLtau]

# [GenerateEvents], [CrossSection]

# Optional (if generating events):
# N events to generate
# Output filename (without file extension)
# [1] original "teg" output [2] Hepevt format [3] HepMC3 format

# SM example
./TEG_v2 <<EOF
7
Ar
2
SM
GenerateEvents
1000
outputs/sm_trident_bnb
2
EOF

# # LmuLtau example
# ./TEG_v2 <<EOF
# 7
# Ar
# 2
# LmuLtau
# 0.001
# 0.1
# GenerateEvents
# 1000
# outputs/0.001_0.1_epem
# 2
# EOF