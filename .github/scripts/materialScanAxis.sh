#!/bin/bash

# This test scan the material along the X (or Z) axis using materialScan from DD4hep
# then extract the second-to-last line, which contains a summary of the scan, 
# then checks if X0 is larger or smaller than X0 for air

compactFile=$1
axis=$2
IsMaterialExpectedToBeAir=$( awk '{print tolower($0)}' <<< $3)


#Check if file exists...
[ ! -f "$compactFile" ] && { echo "Compact file $compactFile not found."; exit 2; }


#Check if materialScan command from DD4hep if found
if ! command -v materialScan &> /dev/null
then
    echo "materialScan could not be found"
    exit 3
fi



#Get the total X0 integrated along the axis. X0 is provided as 12th parameter of the second-to-last line of 
x0=""
if [ "x" = $axis ] || [ "X" = $axis ] ; then
        x0=$(materialScan $compactFile 0 0 0 1000 0 0 |  awk -F' ' '{prevlast = last; last = $0} END {if (NR >= 2) print  prevlast}' | awk -F' ' '{ print $12 }') 

elif  [ "z" = $axis ] || [ "Z" = $axis ] ; then
        x0=$(materialScan $compactFile 0 0 0 0 0 1000 |  awk -F' ' '{prevlast = last; last = $0} END {if (NR >= 2) print  prevlast}' | head -n 1 | awk -F' ' '{ print $12 }') 

else 
        echo "Axis " $axis " not supported"
        exit 22
fi

x0=$(expr "$x0" )
echo "x0 is " $x0 "cm"


xair=0.032756
if [ "air" = $IsMaterialExpectedToBeAir ] ; then
        if (( $(echo "$x0 > $xair" |bc -l) )); then
                echo "X0 larger than or equal to air X0 (0.032756 cm)"
                exit 22
        else
                exit 0
        fi
else 
        if (( $(echo "$x0 <= $xair" |bc -l) )); then
                echo "X0 smaller than or equal to air X0 (0.032756 cm)"
                exit 22
        else
                exit 0
        fi
fi