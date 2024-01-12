#!/bin/bash

#Use like: 
#./filelist_script.sh "../../output/HICAR_snow_1000m/icar_out_201*" file_list_HICARsnow250m.txt


pattern=$1
out_file=$2
temp_file="temp1223334332.txt"
echo $pattern
readlink -f $pattern > $temp_file

while read line; do
	echo "\"$line\""
        echo "\"$line\"" >> $out_file
done <$temp_file
rm $temp_file


