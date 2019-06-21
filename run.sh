#!/bin/bash

if [ $# != 5 ]; then
	echo "***ERROR*** Use: $0 input.fasta output.fasta gecko-len gecko-sim borderSize"
	exit -1
fi

INPUT=$1
FREQS=$1.freqs
REPS=$2
GECKOLEN=$3
GECKOSIM=$4
BORDERSIZE=$5

# Paths go here
REPEATSCOUT=$SOFT/repeatScout/RepeatScout-1/
GECKO=$SOFT/gecko/bin
REPEATEVO=$SOFT/repeatEVO/

# Names and extensions

dirNameY=$(readlink -f $1 | xargs dirname)
seqYName=$(basename "$1")
extensionY="${seqYName##*.}"
seqYName="${seqYName%.*}"

# Run repeatScout to generate k-mers and their frequency
$REPEATSCOUT/build_lmer_table -l 20 -sequence $INPUT -freq $FREQS

# Run repeatScout to find the repeating elements 
$REPEATSCOUT/RepeatScout -sequence $INPUT -output $REPS -freq $FREQS -l 20

# Names and extensions part 2
dirNameX=$(readlink -f $REPS | xargs dirname)
seqXName=$(basename "$REPS")
extensionX="${seqXName##*.}"
seqXName="${seqXName%.*}"

# Run gecko to get all the repetitions (repeatScout reports only one)
$GECKO/workflow.sh $REPS $INPUT $GECKOLEN $GECKOSIM 32 1

# Get borders
$REPEATEVO/getRepeatsAndBorders.sh results/$seqXName-$seqYName.csv $REPS $INPUT $REPS.all $BORDERSIZE




