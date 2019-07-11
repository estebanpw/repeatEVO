#!/bin/bash

INPUTCSV=$1
THRESHOLD=$2

# First, sort by x start and x end (y not necessary since its only one sequence)
tail -n +18 $INPUTCSV | sort -g -t , -k2 -k3 > $INPUTCSV.sorted

# Second, join the blocks

BLOCKID=0
XP=-1
YP=-1
XCURR=0
YCURR=0

rm $INPUTCSV.blocks
rm $INPUTCSV.synteny

while IFS= read -r line
do
	XCURR=$(echo "$line" | awk -F , '{print $2}')
	YCURR=$(echo "$line" | awk -F , '{print $3}')

	if [ "$XP" -eq "-1" ]
	then
		XP=$XCURR
		YP=$YCURR
		echo "$line, $BLOCKID" >> $INPUTCSV.blocks
		continue
	fi

	DIFFX=`expr $XCURR - $XP`
	DIFFY=`expr $YCURR - $YP`

	if [ "$DIFFX" -gt "0" ] && [ "$DIFFY" -gt "0" ]
	then
		if [ "$DIFFX" -le "$THRESHOLD" ] && [ "$DIFFY" -le "$THRESHOLD" ]
		then
			
			echo "$line,$BLOCKID" >> $INPUTCSV.blocks
			continue

		fi
	fi

	BLOCKID=`expr $BLOCKID + 1`
	echo "$line, $BLOCKID" >> $INPUTCSV.blocks

	XP=$XCURR
	YP=$YCURR

done < "$INPUTCSV.sorted"


PREVLINE="NULL"
CURRLINE="NULL"
INITLINE="NULL"
PREVBLOCK=-1
CURRBLOCK=-1
NUM=0

while IFS= read -r line
do

	# Get new lines
	PREVLINE=$CURRLINE

	if [ "$INITLINE" == "NULL" ]
	then
		INITLINE=$PREVLINE
		ENDINGX=$(echo "$line" | awk -F, '{print $4}')
		ENDINGY=$(echo "$line" | awk -F, '{print $5}')
	fi

	CURRLINE=$line

	# Get block ids
	PREVBLOCK=$CURRBLOCK
	CURRBLOCK=$(echo "$line" | awk -F , '{print $15}')

	#echo "$NUM blocks: ($PREVBLOCK, $CURRBLOCK) -> $line"
	#read -p "[Next]..."


	if [ "$PREVBLOCK" -eq "$CURRBLOCK" ]
	then

		# Just grab the ending X and Y coords
		ENDINGX=$(echo "$line" | awk -F, '{print $4}')
		ENDINGY=$(echo "$line" | awk -F, '{print $5}')

		NUM=`expr $NUM + 1`

	else

		if [ "$PREVBLOCK" -ne "-1" ]
		then

			if [ "$NUM" -gt "0" ]
			then	
				# Write it out
				BLOCK=$(echo "$INITLINE" | awk -F "," -v OFS=',' -v a="$ENDINGX" -v b="$ENDINGY" '{print $1,$2,$3,a,b,$6,$7,a-$2,$9,$10,$11,$12,$13,$14,$15}' )
			
				echo "$BLOCK" >> $INPUTCSV.synteny

			fi

			NUM=0
			INITLINE="NULL"

		fi
	fi
	
	

done < "$INPUTCSV.blocks"

#rm $INPUTCSV.sorted

