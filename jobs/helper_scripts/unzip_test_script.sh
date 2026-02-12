## This script is intended for use when grid job to add overlays to the simulated event in HEPevt format
## It needs to be passed to grid node, and run on grid node before the start of any job.
## What it does stepwise: 
## 	1. Unzip the tarball provided (variable TARname)
##      2. Read the HEPevt file list  (variable FileListName)
##      3. Based on Process number assigned to given node, read the path of HEPevt file from the file list
##      4. Replace the placeholder in template fhicl by the HEPevt file name and save it as the fhicl that will be run
## 	   to turn HEPevt into artroot and add overlays
## 
## Configurable parameters:
##   1. TARname     2. FileListName
##   3. TemplateFhicl    4. FinalFhicl
##
## writen by Guanqun, in Oct 2021
#!/bin/bash

## configurable parameter
TARname="numu_epem_50evts.tar"
FileListName="numu_epem/numu_epem_5000_splited_file_50evts.list"   # name of the list of HEPevt files
TemplateFhicl="Simulation_AddOverlay_batch_TEMPLATE.fcl"
FinalFhicl="Simulation_AddOverlay_batch_LOCAL.fcl"


## unzip the tarball with HEPevt files
## unzip command isn't found on grid node, use tar instead
## which will give you a folder which has list of HEPevt file list, and a subfolder containing all HEPevt files
## by default, assume the folder after unzipping has the same name as the tarball
echo "untaring "$TARname"..."
tar -xvf $TARname

## directory name and the HEPevt file list that will be used
DirName=$(basename -s .tar $TARname)
# FileList=$DirName"/"$FileListName # getting wrong path if TARname doesn't match the highest level directory name within the archive
FileList=$FileListName


## use $PROCESS number of grid node to tell it which hepevt file to use
#declare -i LINECNT=1
LINECNT=$((PROCESS+1))
hepevtfile=$(sed "${LINECNT}q;d" $FileList) #getting the file name to run on from ext-unbiased
echo "Line: "$LINECNT
echo "HEPevt file used is: "$hepevtfile


## replace the PlaceHolder in template fcl by the hepevt file
sed -e "s:TEMPTEMPTEMP:$hepevtfile:g" $TemplateFhicl > $FinalFhicl 