#!/bin/bash

if [ "$#" -lt 5 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 fragsFILE.csv fastaX fastaY outputRepeats borderSize"
   echo ""
   exit -1
fi

FRAGS=$1
FASTAX=$2
FASTAY=$3
ALIGN=$4
BORDER=$5




BINDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# generate indices

$BINDIR/indexmaker $FASTAX $FASTAX.idx
$BINDIR/indexmaker $FASTAY $FASTAY.idx

$BINDIR/reverseComplement $FASTAY $FASTAY.rev
$BINDIR/indexmaker $FASTAY.rev $FASTAY.rev.idx

# extract frags

sed 's/,/ /g' $FRAGS > $FRAGS.fix

$BINDIR/csvExtractRepeatsAndBorders $FRAGS.fix $FASTAX $FASTAY $FASTAY.rev $FASTAX.idx $FASTAY.idx $FASTAY.rev.idx $ALIGN $BORDER


rm $FRAGS.fix
rm $FASTAX.idx
rm $FASTAY.idx
rm $FASTAY.rev
rm $FASTAY.rev.idx

