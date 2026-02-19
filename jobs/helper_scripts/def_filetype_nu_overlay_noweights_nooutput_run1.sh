#! /bin/bash
#------------------------------------------------------------------
#
# Purpose: This script is intended to update the configuration
#          of the fhicl using the initilization source hook. It
#          makes a wrapper fcl file that overrides certain fcl
#          parameters used by module WireCellAnaTree and
#          WireCellEventWeightTree, specifically:
#        
#          physics.analyzers.wcpselection.FileType: "data_extbnb_run3_G1"
#          physics.analyzers.wcpweights.FileType: "data_extbnb_run3_G1"
#
# Created: Wenqiang Gu, 14-Jan-2021
#
#------------------------------------------------------------------

# Make sure batch environment variables needed by this script are defined.

if [ x$FCL = x ]; then
  echo "Variable FCL not defined."
  exit 1
fi

# Make sure fcl file $FCL exists.

if [ ! -f $FCL ]; then
  echo "Fcl file $FCL does not exist."
  exit 1
fi

# Rename the existing fcl file $FCL to something else.

mv $FCL mix_wrapper.fcl

# Generate wrapper.

cat <<EOF > $FCL
#include "mix_wrapper.fcl"

physics.analyzers.wcpselection.FileType: "prodgenie_bnb_nu_overlay_run1"
physics.analyzers.wcpselection.SaveWeights: false
physics.analyzers.wcpweights.SaveWeights: false
physics.analyzers.wcpweights.SaveFullWeights: false
physics.analyzers.wcpweights.SaveGenieWeights: false
#physics.analyzers.wcpweights.FileType: "prodgenie_bnb_nu_overlay_run1"

physics.end_paths: [ stream1, ana ]

services.TFileService.fileName: 	"reco_stage_2_hist.root"
microboone_tfile_metadata.JSONFileName: ["reco_stage_2_hist.root.json"]

physics.stream1:[]

EOF

# Make sure IFDH service is configured in fcl file.

if ! lar --debug-config=/dev/stdout -c $FCL | grep -q IFDH:; then
  cat <<EOF >> $FCL
services.IFDH:
{
}

EOF
fi
