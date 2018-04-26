IND=$1
shift
BAM_FOLDER_IN=$2
shift
BAM_FOLDER_OUT=$3

TO_SORT='true'
TO_DEDUP='true'
TO_INDEX='true'
CHECK_CVG='true'

while getopts 'SDIC' flag; do
  case "${flag}" in
    S) TO_SORT='false' ;;
    D) TO_DEDUP='false' ;;
    I) TO_INDEX='false' ;;
    C) CHECK_CVG='false' ;;
    *) error "Unexpected option ${flag}" ;;
  esac
done

PICARD=/nas02/apps/picard-2.2.4/picard-tools-2.2.4/picard.jar

echo "Sorting, marking duplicates and indexing bam file for: $IND"
echo "To sort: $TO_SORT"
echo "To dedup: $TO_DEDUP"
echo "To index: $TO_INDEX"
echo "Check coverage: $CHECK_CVG"

if [ ! -d "tmp" ]; then mkdir tmp; fi

if [ $TO_SORT == true ]
  then
  echo "Step 1: Sorting bam file..." # 43 minutes
  java -Xmx2g -Djava.io.tmpdir=tmp -jar $PICARD SortSam SORT_ORDER=coordinate INPUT=$BAM_FOLDER_IN/$IND.sam OUTPUT=$BAM_FOLDER_OUT/$IND.sort.bam TMP_DIR=tmp
fi

if [ $TO_DEDUP == true ]
  then
  echo "Step 2: Marking duplicates with Picard..." # 38 minutes, needs 8GB memory
  java -Xmx14g -Djava.io.tmpdir=tmp -jar $PICARD MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=$BAM_FOLDER_OUT/$IND.sort.bam OUTPUT=$BAM_FOLDER_OUT/$IND.dedup.bam TMP_DIR=tmp METRICS_FILE=$BAM_FOLDER_OUT/dedupMetrics.$IND.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
fi

if [ $TO_INDEX == true ]
  then
  echo "Step 3: Indexing deduplicated bam file..." # 5 minutes, needs at least 8GB memory
  java -Xmx14g -Djava.io.tmpdir=tmp -jar $PICARD BuildBamIndex INPUT=$BAM_FOLDER_OUT/$IND.dedup.bam TMP_DIR=tmp
fi

echo "Removing sorted bam file (keeping only sorted+dedupped)..."
rm $BAM_FOLDER_OUT/$IND.sort.bam

if [ $CHECK_CVG == true ]
  then
  echo "Checking coverage..."
  DEPTH=`samtools depth $BASEPATH/processed/$IND.dedup.bam | awk '{sum+=$3} END {print sum/NR}'`
  echo "Coverage for $IND is: $DEPTH"
  echo "$IND $DEPTH" > depth.$IND.txt
fi

echo "Done with script."