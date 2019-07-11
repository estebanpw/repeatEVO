#!/bin/bash


if [ "$#" -ne 3 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <repeatScout repeats> <repeatID> <dna sequence>"
   echo ""
   exit -1
fi

REPEATS=$1
ID=$2
DNA=$3

# Converts multifasta line to single-line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $REPEATS > $REPEATS.oneline


IDextract=`expr $2 + 1`
IDextract=`expr $IDextract \* 2`
IDextract=`expr $IDextract - 1`

# Extracts the repetition corresponding to the id
awk -v id="$IDextract" 'BEGIN{n=0;} { if(n==id) print $0 ; if(n==(id+1)) print $0 ; n=n+1;}' $REPEATS.oneline > repeat-$ID.fasta

# Use blast to get alignments with gaps
cp $DNA .
dbDNA=$(basename $DNA)
makeblastdb -dbtype nucl -in $dbDNA

# Run blast
blastn -query repeat-$ID.fasta -db $dbDNA -outfmt 6 -out all-repeats-$ID.blast

# Extract DNA from repetitions
/home/estebanpw/software/seqExtractor/seqExtractorNatural all-repeats-$ID.blast $dbDNA > all-repeats-from-$ID.fasta

rm $REPEATS.oneline

