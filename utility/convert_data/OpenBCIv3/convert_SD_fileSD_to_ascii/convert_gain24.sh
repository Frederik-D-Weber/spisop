#!/bin/bash
#cd /cinserv/workspace/OpenBCI/convertSD
find -H -name '*.TXT' > temp

#awk 'gsub(/ /, "\\ "); {print "\x27" $0 "\x27"}' temp > temp2
 
FILES=$(cat temp)  
for FILE in $FILES; do
	TO_FILE=$FILE".new.csv"
	#echo $FILE
	gawk -f "convert_gain24.awk" $FILE > $TO_FILE
	#awk -f convert.awk $FILE > $TO_FILE
	#cp $TO_FILE $FILE
	#rm $TO_FILE
done

rm temp
#rm temp2
