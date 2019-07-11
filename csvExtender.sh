#!/bin/bash


if [ "$#" -ne 2 ]; then
   echo " ==== ERROR ... you called this script inappropriately."
   echo ""
   echo "   usage:  $0 right-border-extension csv"
   echo ""
   exit -1
fi




BORDER=$1
CSV=$2

awk -F , -v border="$BORDER" 'BEGIN{OFS=",";} {if($1 == substr($0,0,4)) { print $1,$2,$3,$4,$5+border,$6,$7,$8,$9,$10,$11,$12,$13,$14 } else { print $0}   }' $CSV
