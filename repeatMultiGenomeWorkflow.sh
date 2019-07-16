#!/bin/bash


if [ "$#" -lt 4 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <repeatScout repeats> <repeatID> <dna sequences FOLDER> <border> [%minLen]"
   echo ""
   exit -1
fi

REPEATS=$1
ID=$2
DNAFOLDER=$3
BORDER=$4
PERCENTAGE=0.7

if [ "$#" -eq "5" ]
then

	PERCENTAGE=$5

fi

LASTFOLDER=$(basename $DNAFOLDER)
mkdir dna-$LASTFOLDER

echo "Using percentage $PERCENTAGE"

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

for i in $DNAFOLDER/* ; do

	# Copy sequence
	cp $i dna-$LASTFOLDER/

	# Get chromo name
	IFS='.', read -a splits <<< "$i"
	dimension=${#splits[@]}
	IDENTIFIER=${splits[$dimension-2]}

	# Use blast to get alignments with gaps
	dbDNA=$(basename $i)
	makeblastdb -dbtype nucl -in dna-$LASTFOLDER/$dbDNA

	# Run blast
	blastn -query repeat-$ID.fasta -db dna-$LASTFOLDER/$dbDNA -outfmt 6 -out all-repeats-$ID-$IDENTIFIER.blast

	# Calculate min length
	percentageLength=$(echo "" |  awk -v num="$length" -v perc="$PERCENTAGE" '{ printf("%d", num*perc) }')
	echo "Looking for repeats of size (minimum) $percentageLength"

	cat all-repeats-$ID-$IDENTIFIER.blast | awk -v OFS="\t" -v minlen="$percentageLength" '{ if($4 > minlen) { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}  }' > all-repeats-$ID-$IDENTIFIER-filtered.blast

	rm all-repeats-$ID-$IDENTIFIER.blast
	mv all-repeats-$ID-$IDENTIFIER-filtered.blast all-repeats-$ID-$IDENTIFIER.blast

	# Extract with border (but only if they had minimum length)
	awk -v OFS="\t" -v border="$BORDER" '{ if($9 < $10) { print $1,$2,$3,$4,$5,$6,$7,$8,$9-border,$10+border,$11,$12} else { print $1,$2,$3,$4,$5,$6,$7,$8,$9+border,$10-border,$11,$12  }  }' all-repeats-$ID-$IDENTIFIER.blast > all-repeats-$ID-$IDENTIFIER-border.blast

	# Extract DNA from repetitions (minimum a 70% of the original rep (or given param))
	/home/estebanpw/software/seqExtractor/seqExtractorNatural all-repeats-$ID-$IDENTIFIER-border.blast dna-$LASTFOLDER/$dbDNA $percentageLength | awk -v id="$IDENTIFIER" '{ if(substr($0,1,1) == ">") {printf("%s_%s\n",$0,id); } else {print $0}  }' >  all-repeats-from-$ID-$IDENTIFIER.fasta

	cat all-repeats-from-$ID-$IDENTIFIER.fasta >> concatenated-repeats-from-$ID.fasta

done


# Align the repetitions
clustalw -INFILE=concatenated-repeats-from-$ID.fasta



# Activate virtual env because python is stoopid
source /home/estebanpw/software/repeatEVO/plot/plotenv/bin/activate

# Run conservation analysis
python /home/estebanpw/software/repeatEVO/plot/clustal2image.py -i concatenated-all-repeats-from-$ID.aln -o .

# Deactivate unwanted programs (I dont want virtualenv)
deactivate

# Remove unwanted
rm $REPEATS.oneline

# Info

echo "[[[[[ DEAR USER ]]]]] Remember to run your gene analysis using blastx and swissprot. Just like this:"
echo "blastx -query repeat.fasta -db /home/estebanpw/data/swissprot/swissprot.fasta -out genesearch.blast"

