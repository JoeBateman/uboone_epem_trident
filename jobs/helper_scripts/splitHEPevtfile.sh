## The script to split HEPevt file into smaller files with configurable number of events per file
## and then generating file list (both in full path, and in relative path with smaller batches) under the same output directoty
## written by Guanqun

#!/bin/bash

CurrentDir=$PWD

## --------- CONFIGURABLE -------------------------

## number of event you'd like smaller HEPevt file to contain
declare -i BatchSize=50
particlesPerEvent=6  ## number of particles per event in HEPevt file

## the directory where large HEPevt text file sits
LargeFilePath="/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/outputs"
FileName="numu_epem_5000.txt"

## the directory where HEPevt text file with smaller size will sit
OutputFilePath="/exp/uboone/app/users/jbateman/workdir/DarkNews/Trident/jobs/HEPevtfile_list/numu_epem"

## --------- CONFIGURABLE -------------------------

# N initial/final state particles, plus one line for each event header
LNCountPerFile=$(((particlesPerEvent+1)*BatchSize)) 
SmallFileSubDir="splited_file_"$BatchSize"evts"
mkdir -p $OutputFilePath"/"$SmallFileSubDir

## ======= loop over each large HEPevt file, and split it into smaller files and save to output directory
cd $LargeFilePath
for file in $FileName; do
    echo $file
    split -l $LNCountPerFile $file $OutputFilePath"/"$SmallFileSubDir"/"$file 
done


## generate a filelist with full path --> used for simulating events with no overlay
ls $OutputFilePath"/"$SmallFileSubDir"/"* > $OutputFilePath"/"gntp_HEPevt_all_splited_file_$BatchSize"evts".list 


## ======= generate filelist in smaller batches and zip the whole folder ============
##  --> used for simulating events with overlays
## each filelist contain relative path to smaller HEPevt files splited from a certain larger HEPevt file

cd $OutputFilePath
OutputDirName=$(basename $OutputFilePath)  # name of the directory, used for zip command
cd ../
for file in $LargeFilePath"/"$FileName;do

   ## remove the path infront and suffix
   base_name=$(basename -s .txt $file)
   ls $OutputDirName"/"$SmallFileSubDir"/"$base_name* > $OutputDirName"/"$base_name"_splited_file_"$BatchSize"evts.list"
done

## unzip command is not available on grid node
## use tar instead
#zip -r $OutputDirName".zip" $OutputDirName
tar -cvf $OutputDirName"_"$BatchSize"evts.tar" $OutputDirName/*$BatchSize"evts"*

## =======  end of generating filelist in smaller batches and zipping the whole folder ================


## go back to current dir
cd $CurrentDir