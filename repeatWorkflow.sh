#!/bin/bash


if [ "$#" -ne 4 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <repeatScout repeats> <repeatID> <dna sequence> <border>"
   echo ""
   exit -1
fi

REPEATS=$1
ID=$2
DNA=$3
BORDER=$4

# Converts multifasta line to single-line fasta
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $REPEATS > $REPEATS.oneline


IDextract=`expr $2 + 1`
IDextract=`expr $IDextract \* 2`
IDextract=`expr $IDextract - 1`

# Extracts the repetition corresponding to the id
awk -v id="$IDextract" 'BEGIN{n=0;} { if(n==id) print $0 ; if(n==(id+1)) print $0 ; n=n+1;}' $REPEATS.oneline > repeat-$ID.fasta

# Get length of repeat
length=$(tail -n +2 repeat-$ID.fasta | wc -c)
echo "Length of repeat is $length"

# Use blast to get alignments with gaps
cp $DNA .
dbDNA=$(basename $DNA)
makeblastdb -dbtype nucl -in $dbDNA

# Run blast
blastn -query repeat-$ID.fasta -db $dbDNA -outfmt 6 -out all-repeats-$ID.blast

awk -v OFS="\t" -v border="$BORDER" '{ if($9 < $10) { print $1,$2,$3,$4,$5,$6,$7,$8,$9-border,$10+border,$11,$12} else { print $1,$2,$3,$4,$5,$6,$7,$8,$9+border,$10-border,$11,$12  }  }' all-repeats-$ID.blast > all-repeats-$ID-border.blast


# Extract DNA from repetitions (minimum a 70% of the original rep)
percentageLength=$(echo "" |  awk -v num="$length" '{ printf("%d", num*0.7) }')
echo "Looking for repeats of size (minimum) $percentageLength"
/home/estebanpw/software/seqExtractor/seqExtractorNatural all-repeats-$ID-border.blast $dbDNA $percentageLength > all-repeats-from-$ID.fasta

# Align the repetitions
clustalw -INFILE=all-repeats-from-$ID.fasta

# Activate virtual env because python is stoopid
source /home/estebanpw/software/repeatEVO/plot/plotenv/bin/activate

# Run conservation analysis
python /home/estebanpw/software/repeatEVO/plot/clustal2image.py -i all-repeats-from-$ID.aln

# Deactivate unwanted programs (I dont want virtualenv)
deactivate

# Remove unwanted
rm $REPEATS.oneline

# Info

echo "[[[[[ DEAR USER ]]]]] Remember to run your gene analysis using blastx and swissprot. Just like this:"
echo "blastx -query repeat.fasta -db /home/estebanpw/data/swissprot/swissprot.fasta -out genesearch.blast"

