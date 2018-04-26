#!/bin/bash -l

## Command-line args
FILE_ID=$1 
MEM=$2 # Memory

## Other scripts and files:
PGDS=/proj/cmarlab/users/jelmer/software/PGDSpider_2.1.1.0/PGDSpider2-cli.jar
SPID=scripts/conversion/fasta2nexus.spid

## Other variables:
INPUT=seqdata/fasta/$FILE_ID.fasta
OUTPUT=seqdata/nexus/$FILE_ID.nexus

## Report:
date
echo "Script: vcf2geno.sh"
echo "File ID: $FILE_ID"
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo "SPID file: $SPID"

## Convert:
java -Xmx${MEM}G -Xms${MEM}G -jar $PGDS -inputfile $INPUT -inputformat FASTA -outputfile $OUTPUT.tmp -outputformat NEXUS -spid $SPID

## Remove taxon sets block:
cat $OUTPUT.tmp | head -n -1 | grep -v "BEGIN SETS" | grep -v "TaxSet" | grep -v "TaxPartition" | grep -v pop_1 > $OUTPUT
rm $OUTPUT.tmp

rm *log

echo "Done with script."
date