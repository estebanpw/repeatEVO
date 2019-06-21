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

while IFS= read -r line
do
	XCURR=$(echo "$line" | awk -F , '{print $2}')
	YCURR=$(echo "$line" | awk -F , '{print $3}')

	if [ "$XP" -eq "-1" ]
	then
		XP=$XCURR
		YP=$YCURR
		echo "$line, $BLOCKID"
		continue
	fi

	DIFFX=`expr $XCURR - $XP`
	DIFFY=`expr $YCURR - $YP`

	if [ "$DIFFX" -gt "0" ] && [ "$DIFFY" -gt "0" ]
	then
		if [ "$DIFFX" -le "$THRESHOLD" ] && [ "$DIFFY" -le "$THRESHOLD" ]
		then
			
			echo "$line, $BLOCKID"
			continue

		fi
	fi

	BLOCKID=`expr $BLOCKID + 1`
	echo "$line, $BLOCKID"

	XP=$XCURR
	YP=$YCURR

done < "$INPUTCSV.sorted"
