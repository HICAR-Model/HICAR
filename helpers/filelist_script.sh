#!/bin/bash

pattern=$1
out_file=$2
temp_file="temp1223334332.txt"
readlink -f $pattern > $temp_file

while read line; do
	echo "\"$line\""
        echo "\"$line\"" >> $out_file
done <$temp_file
rm $temp_file


