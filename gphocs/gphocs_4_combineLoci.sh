#!/bin/bash -l

FILE_ID=$1

SOURCEFOLDER=analyses/gphocs/input_prep
TARGETFOLDER=analyses/gphocs/input

date
echo "Script: gphocs_combineLoci.sh"
echo "File ID: $FILE_ID"
echo "Output file: $TARGETFOLDER/$FILE_ID.gphocsInput.txt"

## Remove files with no sequences:
echo "Files with no sequences (will be removed):"; find $SOURCEFOLDER -maxdepth 1 -name "*$FILE_ID*input" -size -100c
find $SOURCEFOLDER -maxdepth 1 -name "*input" -size -100c -print0 | xargs -0 rm

## Count and print number of loci:
NLOCI=$(grep "$FILE_ID" $SOURCEFOLDER/$FILE_ID* | wc -l) 
echo "Number of loci: $NLOCI"
printf "${NLOCI}\n\n" > $TARGETFOLDER/$FILE_ID.gphocsInput.txt

## Add loci to file:
cat $SOURCEFOLDER/$FILE_ID*input >> $TARGETFOLDER/$FILE_ID.gphocsInput.txt

echo "Done with script."
date