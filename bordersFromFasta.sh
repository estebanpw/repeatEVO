#!/bin/bash

if [ "$#" -ne 2 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 <Repeats fasta> <border size>"
   echo ""
   exit -1
fi

FASTA=$1
BORDER=$2

# Extract repetitions
# Left border
tail -n +2 $FASTA | awk -v border="$BORDER" 'BEGIN{n=0;} { if((n%2)==1){  print substr($0, 0, border+1); } else { print $0  } ; n=n+1;}' > $FASTA.leftborders

# Right border
tail -n +2 $FASTA | awk -v border="$BORDER" 'BEGIN{n=0;} { if((n%2)==1){  print substr($0, length($0)-border);   } else { print $0  } ; n=n+1;}' > $FASTA.rightborders

NAMELEFT=$(basename $FASTA.leftborders)
NAMERIGHT=$(basename $FASTA.rightborders)



cp $FASTA.leftborders ./$NAMELEFT
cp $FASTA.rightborders ./$NAMERIGHT

SECONDNAMELEFT=${NAMELEFT%".leftborders"}
SECONDNAMERIGHT=${NAMERIGHT%".rightborders"}



mv $NAMELEFT  leftborders_$SECONDNAMELEFT
mv $NAMERIGHT rightborders_$SECONDNAMERIGHT

NAMELEFT=leftborders_$SECONDNAMELEFT
NAMERIGHT=rightborders_$SECONDNAMERIGHT

# Align the repetitions
clustalw -INFILE=$NAMELEFT

clustalw -INFILE=$NAMERIGHT

# Activate virtual env because python is stoopid
source /home/estebanpw/software/repeatEVO/plot/plotenv/bin/activate

WORDTOREMOVE=.fasta
SUPERNAMELEFT=$(printf '%s\n' "${NAMELEFT//$WORDTOREMOVE/}")
SUPERNAMERIGHT=$(printf '%s\n' "${NAMERIGHT//$WORDTOREMOVE/}")

# Run conservation analysis
python /home/estebanpw/software/repeatEVO/plot/clustal2image.py -i $SUPERNAMELEFT.aln
python /home/estebanpw/software/repeatEVO/plot/clustal2image.py -i $SUPERNAMERIGHT.aln

# Deactivate unwanted programs (I dont want virtualenv)
deactivate

